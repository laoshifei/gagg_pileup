"""Continuous-charge supervised estimators with explicit estimation risks.

The network is a small NumPy MLP. Its inputs are observation-derived single-fit
innovations; labels are accepted only by training losses. Count scores are
calibrated separately, and are not interpreted as likelihoods or multiplied
by a second event-rate prior.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize
from scipy.special import expit, logsumexp

from .alpha_pair_joint_likelihood import refine
from .alpha_pair_model_sampler import GaussianFactor

RISKS = ('absolute', 'logarithmic', 'relative')
TIME_GRID = np.arange(-384., 384.01, 2.)
INPUT_SIZE = 1021 + len(TIME_GRID) + 2
OUTPUT_SIZE = 3 + 2*len(TIME_GRID)


def charge_risk(predicted, target, kind):
    """Elementwise risk and derivative with respect to predicted Q."""
    p = np.asarray(predicted, float)
    q = np.asarray(target, float)
    if np.any(p <= 0.) or np.any(q <= 0.):
        raise ValueError('Positive charges required')
    if kind == 'absolute':
        return ((p-q)/5000.)**2, 2.*(p-q)/5000.**2
    if kind == 'logarithmic':
        return np.log(p/q)**2, 2.*np.log(p/q)/p
    if kind == 'relative':
        return (p/q-1.)**2, 2.*(p/q-1.)/q
    raise ValueError(kind)


def single_features(response, noise, voltage):
    """Full voltage -> a finite K=1 pilot and full/aligned innovations.

    The pilot is an inferred hypothesis, not knowledge that the frame is a
    single event. Even empty and failed-pilot frames are retained.
    """
    y = np.asarray(voltage, float)
    if y.shape != (1024,) or not np.isfinite(y).all():
        raise ValueError('A finite complete 1024-sample voltage is required')
    observed = y-response.baseline-noise.background_mV
    reference = response.components([256.], [3000.])[0]
    delay = response.sample_grid[1+np.argmax(np.diff(gaussian_filter1d(reference, 1.))) ]-256.
    peak = int(np.argmax(np.diff(gaussian_filter1d(observed, 1.))))+1
    anchor = float(np.clip(response.sample_grid[peak]-delay, 0., response.sample_grid[-1]))
    charge = float(np.clip(observed.sum()*2./response.base.area_per_pe, 50., 50000.))
    pilot = refine(response, noise, y, dict(times_ns=[anchor], charge_pe=[charge]), radius=24., maxiter=35)
    anchor = float(pilot['times_ns'][0]);charge = float(pilot['charge_pe'][0])
    curves = response.components([anchor], [charge])
    z = noise.transform(observed-curves[0])
    factor = GaussianFactor.from_noise(noise, curves, [charge])
    whitened = np.arcsinh(factor.whiten(z)/3.)
    aligned = np.interp(anchor+TIME_GRID, response.sample_grid[3:], whitened, left=0., right=0.)
    feature = np.r_[whitened, aligned, np.log(charge/3000.), anchor/2048.].astype(np.float32)
    assert feature.shape == (INPUT_SIZE,) and np.isfinite(feature).all()
    return feature, pilot


def unpack_output(output, pilot_charge):
    z = np.asarray(output, float)
    total = np.asarray(pilot_charge, float)*np.exp(.1*np.tanh(z[:, 1]))
    fraction = expit(np.clip(z[:, 2], -30., 30.))
    q = total[:, None]*np.column_stack((fraction, 1.-fraction))
    return q, z[:, 3:].reshape(len(z), 2, len(TIME_GRID))


def objective(output, pilot_charge, true_charge, true_times, pilot_time, count, kind):
    """BCE + 3*charge risk + .4*time cross entropy; analytic output gradient."""
    z = np.asarray(output, float);batch = len(z)
    pair = np.asarray(count) == 2;number = max(int(pair.sum()), 1)
    predicted, logits = unpack_output(z, pilot_charge)
    truth = np.maximum(np.asarray(true_charge, float), 1.)
    values, derivative = charge_risk(predicted, truth, kind)
    values = values*pair[:, None];derivative = derivative*pair[:, None]/(2.*number)
    qrisk = float(values.sum()/(2.*number))
    label = pair.astype(float)
    count_loss = float(np.mean(np.logaddexp(0., z[:, 0])-label*z[:, 0]))
    relative = np.asarray(true_times)-np.asarray(pilot_time)[:, None]
    density = np.exp(-.5*((TIME_GRID[None, None, :]-relative[:, :, None])/2.)**2)
    # Inactive rows have no time loss. No truth-dependent inference realignment.
    density /= np.maximum(density.sum(-1, keepdims=True), 1e-300)
    log_prob = logits-logsumexp(logits, axis=-1, keepdims=True)
    tloss = float(-(density*log_prob*pair[:, None, None]).sum()/(2.*number))
    gradient = np.zeros_like(z)
    gradient[:, 0] = (expit(z[:, 0])-label)/batch
    gradient[:, 1] = 3.*np.sum(derivative*predicted, axis=1)*.1*(1.-np.tanh(z[:, 1])**2)
    fraction = expit(np.clip(z[:, 2], -30., 30.));total = predicted.sum(1)
    gradient[:, 2] = 3.*(derivative[:, 0]-derivative[:, 1])*total*fraction*(1.-fraction)*(np.abs(z[:, 2])<30.)
    gradient[:, 3:] = (.4*(np.exp(log_prob)-density)*pair[:, None, None]/(2.*number)).reshape(batch, -1)
    return count_loss+3.*qrisk+.4*tloss, gradient, dict(count=count_loss, charge=qrisk, time=tloss)


class RiskNet:
    """Two SiLU hidden layers; all preprocessing statistics are training-only."""
    def __init__(self, mean, scale, seed, hidden=(128, 96)):
        self.mean=np.asarray(mean, np.float32);self.scale=np.asarray(scale, np.float32)
        rng=np.random.default_rng(seed);widths=(len(mean), *hidden, OUTPUT_SIZE)
        self.parameters={}
        for layer,(before,after) in enumerate(zip(widths[:-1], widths[1:])):
            self.parameters[f'W{layer}']=(rng.normal(0., np.sqrt(2./before), (before, after))*(.1 if layer==2 else 1.)).astype(np.float32)
            self.parameters[f'b{layer}']=np.zeros(after, np.float32)

    def forward(self, features, *, cache=False):
        x=np.clip((np.asarray(features, np.float32)-self.mean)/self.scale, -8., 8.)
        saved=[x]
        for layer in range(2):
            z=x@self.parameters[f'W{layer}']+self.parameters[f'b{layer}']
            probability=expit(z);x=z*probability;saved.extend((z, probability, x))
        output=x@self.parameters['W2']+self.parameters['b2']
        return (output, saved) if cache else output

    def backward(self, gradient, saved):
        g=np.asarray(gradient, np.float32);grads={}
        grads['W2']=saved[6].T@g;grads['b2']=g.sum(0);g=g@self.parameters['W2'].T
        for layer in (1,0):
            z,probability=saved[1+3*layer:3+3*layer]
            g=g*(probability+z*probability*(1.-probability))
            grads[f'W{layer}']=saved[3*layer].T@g;grads[f'b{layer}']=g.sum(0)
            if layer:g=g@self.parameters[f'W{layer}'].T
        return grads

    def save(self, path, metadata=None):
        np.savez_compressed(path, mean=self.mean, scale=self.scale, **self.parameters,
                            metadata=json.dumps(metadata or {}, ensure_ascii=False))

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as values:
            model=cls(values['mean'], values['scale'], 0,
                      hidden=(values['W0'].shape[1], values['W1'].shape[1]))
            model.parameters={name:values[name].copy() for name in model.parameters}
        return model


class Adam:
    def __init__(self, parameters, rate=.0007):
        self.rate=rate;self.step=0
        self.first={k:np.zeros_like(v) for k,v in parameters.items()}
        self.second={k:np.zeros_like(v) for k,v in parameters.items()}

    def update(self, parameters, gradients):
        self.step+=1
        norm=np.sqrt(sum(float(np.sum(g.astype(float)**2)) for g in gradients.values()))
        factor=min(1.,5./max(norm,1e-20))
        for name,g in gradients.items():
            g=g*factor
            self.first[name]=.9*self.first[name]+.1*g
            self.second[name]=.999*self.second[name]+.001*g*g
            first=self.first[name]/(1.-.9**self.step)
            second=self.second[name]/(1.-.999**self.step)
            parameters[name]-=self.rate*first/(np.sqrt(second)+1e-8)
        return float(norm)


def decode(output, pilot_charge, pilot_time):
    """Two candidate slots; count selection remains a separate blind decision."""
    q, logits=unpack_output(output, pilot_charge)
    q=np.clip(q, 50., 50000.)
    offsets=np.empty((len(q),2))
    for i in range(len(q)):
        for j in range(2):
            center=int(np.argmax(logits[i,j]));ids=np.arange(max(0,center-2),min(len(TIME_GRID),center+3))
            weights=np.exp(logits[i,j,ids]-np.max(logits[i,j,ids]))
            offsets[i,j]=float(TIME_GRID[ids]@weights/weights.sum())
    t=np.clip(np.asarray(pilot_time)[:,None]+offsets,0.,2046.25)
    fits=[]
    for i in range(len(q)):
        order=np.argsort(t[i])
        fits.append(dict(times_ns=t[i,order].tolist(),charge_pe=q[i,order].tolist(),
                         network_score=float(output[i,0]),converged=True,truth_inputs=[]))
    return fits


def calibrated_gate(candidates, single_charges, *, edges=(1500.,5000.,12500.), alpha=.025):
    groups=[[] for _ in range(len(edges)+1)]
    for candidate, charge in zip(candidates,single_charges):
        score=float(candidate['network_score']) if min(candidate['charge_pe'])>=500. else -1e12
        groups[int(np.searchsorted(edges,charge,side='right'))].append(score)
    if any(not row for row in groups):raise ValueError('Empty calibration stratum')
    return dict(edges_pe=list(edges),alpha=alpha,scores=[sorted(row) for row in groups],
                score='Learned two-event logit; invalid sub-500pe pair -> -1e12; not a likelihood')


def select_candidate(baseline, candidate, gate):
    if len(baseline)<2 or baseline[1]['charge_pe'][0]<500. or baseline[1]['cost']>=baseline[0]['cost']:
        return dict(baseline[0],birth_accepted=False,null_rank_p=None)
    charge=sum(baseline[1]['charge_pe']);group=int(np.searchsorted(gate['edges_pe'],charge,side='right'))
    eligible=min(candidate['charge_pe'])>=500.
    score=float(candidate['network_score']) if eligible else -1e12
    reference=np.asarray(gate['scores'][group]);rank=(1+np.count_nonzero(reference>=score))/(1+len(reference))
    accepted=bool(eligible and rank<=gate['alpha'])
    return dict(candidate if accepted else baseline[1],birth_accepted=accepted,null_rank_p=float(rank),calibration_bin=group)


def regularized_refine(response,noise,voltage,initial,reference,sigma,*,regularized=True,maxiter=70):
    """Data-dependent learned-Q penalty, not an independent Bayesian prior."""
    y=np.asarray(voltage,float)-response.baseline-noise.background_mV
    order=np.argsort(initial['times_ns']);times=np.asarray(initial['times_ns'])[order]
    charges=np.asarray(initial['charge_pe'])[order]
    target=np.log(np.clip(reference['charge_pe'],50.,50000.));sigma=np.asarray(sigma,float)
    lower=np.r_[np.maximum(0.,times-24.)/10.,np.full(2,np.log(50.))]
    upper=np.r_[np.minimum(response.sample_grid[-1],times+24.)/10.,np.full(2,np.log(50000.))]
    start=np.clip(np.r_[times/10.,np.log(charges)],lower+1e-8,upper-1e-8)
    zy=noise.transform(y)
    def penalty(z):
        ix=np.argsort(z[:2])
        return float(np.sum(((z[2:][ix]-target)/sigma)**2)) if regularized else 0.
    def cost(z):
        q=np.exp(z[2:]);curves=noise.transform(response.components(z[:2]*10.,q))
        return noise.score(zy-curves.sum(0),curves,q,transformed=True)+penalty(z)
    fit=minimize(cost,start,method='L-BFGS-B',bounds=list(zip(lower,upper)),
                 options=dict(maxiter=maxiter,ftol=1e-9,gtol=1e-4,maxls=30,eps=1e-5))
    chosen=fit.x if fit.fun<=cost(start) else start
    ix=np.argsort(chosen[:2]);t=chosen[:2][ix]*10.;q=np.exp(chosen[2:][ix]);curves=response.components(t,q)
    result=noise.score(y-curves.sum(0),curves,q,details=True)
    return dict(result,times_ns=t.tolist(),charge_pe=q.tolist(),penalty=penalty(chosen),
                penalized_cost=result['cost']+penalty(chosen),converged=bool(fit.success),
                evaluations=int(fit.nfev),message=str(fit.message),network_score=reference['network_score'])
