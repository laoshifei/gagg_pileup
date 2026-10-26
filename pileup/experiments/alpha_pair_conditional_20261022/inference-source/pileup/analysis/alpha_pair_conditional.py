"""Calibrated charge/neighbor-dependent means and observation-fixed GLS.

Truth is accepted only by ``calibrate``. Inference takes voltage and blind
initial estimates. The neighbor state is an inferred response surrogate, not
the simulator's cell state. The scintillation lifetimes remain fixed in the
data generator; this empirical mean does not fit an event's lifetimes.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares

from .alpha_correlated_fit import InnovationsMetric, ResidualCalibration
from .alpha_pair_fitting import PairFitter


def features(grid, charge, other_charge=0., other_delay=0.):
    """Mean-shape regression features, with a causal later-neighbor envelope."""
    s = np.asarray(grid, float)
    x = charge / 10000.
    other = other_charge / 10000.
    past = other if other_delay < 0. else 0.
    p95 = past * np.exp(-abs(other_delay) / 95.)
    p351 = past * np.exp(-abs(other_delay) / 351.)
    u = np.maximum(s - other_delay, 0.)
    future = other if other_delay >= 0. else 0.
    # Rise on the readout scale and decay on the fixed scintillation scales.
    # No later neighbor can change a component before that neighbor arrives.
    f95 = future * (-np.expm1(-u / 8.)) * np.exp(-u / 95.)
    f351 = future * (-np.expm1(-u / 8.)) * np.exp(-u / 351.)
    one = np.ones_like(s)
    return np.column_stack((one, x*one, x*x*one, p95*one,
                            p351*one, x*p351*one, f95, f351, x*f351))


@dataclass
class ConditionalResponse:
    settings: dict
    grid: np.ndarray
    coefficients: np.ndarray
    residual: dict

    def __post_init__(self):
        self.base = PairFitter(self.settings)
        self.baseline = self.base.baseline
        self.sample_grid = np.arange(1024)*2.+.25

    def unit(self, relative_time, charge, other_charge=0., other_delay=0.):
        s = np.asarray(relative_time, float)
        correction = np.sum(features(self.grid, charge, other_charge, other_delay)
                            * self.coefficients, axis=1)
        base = np.interp(s, self.base.grid, self.base.kernel, left=0., right=0.)
        correction = np.interp(s, self.grid, correction, left=0., right=0.)
        # Fixed smooth joins; never extrapolate an empirical correction.
        onset = np.clip(s/8., 0., 1.)
        taper = onset*onset*(3.-2.*onset) * np.clip((self.grid[-1]-s)/120., 0., 1.)
        return np.maximum(0., base + correction*taper)

    def components(self, times, charge):
        t = np.asarray(times, float); q = np.asarray(charge, float)
        if t.ndim != 1 or t.shape != q.shape or len(t)>2:
            raise ValueError('Expected zero, one or two aligned times/charges')
        curves = []
        for j in range(len(t)):
            other = 1-j if len(t)==2 else None
            curves.append(q[j]*self.unit(self.sample_grid-t[j], q[j],
                q[other] if other is not None else 0.,
                t[other]-t[j] if other is not None else 0.))
        return np.asarray(curves).reshape(len(t), len(self.sample_grid))

    def metric(self, voltage, correlated=True):
        if not correlated:
            return None
        # Fixed observed pilot: all K/Q/t candidates in a frame use exactly C.
        pilot = gaussian_filter1d(np.maximum(np.asarray(voltage)-self.baseline, 0.), 3.)
        return InnovationsMetric(pilot, ResidualCalibration(**self.residual))

    def refine(self, voltage, initial, *, metric=None, radius=16.):
        y = np.asarray(voltage, float)-self.baseline
        times = np.asarray(initial['times_ns'], float)
        charge = np.asarray(initial['charge_pe'], float)
        if y.shape != (1024,) or not np.isfinite(y).all():
            raise ValueError('A finite 1024-sample voltage is required')
        apply = (lambda x:x) if metric is None else metric.apply
        if not len(times):
            r = apply(y)
            return dict(times_ns=[], charge_pe=[], cost=float(r@r), converged=True)
        order = np.argsort(times); times=times[order]; charge=charge[order]
        # Optimize scaled times and log Q, no physical truth or true count.
        n=len(times)
        lower=np.r_[np.maximum(0.,times-radius)/10., np.full(n,np.log(50.))]
        upper=np.r_[np.minimum(self.sample_grid[-1],times+radius)/10., np.full(n,np.log(50000.))]
        z0=np.clip(np.r_[times/10.,np.log(np.maximum(charge,50.))],lower+1e-9,upper-1e-9)
        def residual(z):
            model=self.components(z[:n]*10.,np.exp(z[n:])).sum(0)
            return apply(model-y)
        result=least_squares(residual,z0,bounds=(lower,upper),max_nfev=60,
                             ftol=1e-7,xtol=1e-7,gtol=1e-5,diff_step=1e-4)
        z=result.x if np.sum(result.fun**2)<=np.sum(residual(z0)**2) else z0
        t=z[:n]*10.;q=np.exp(z[n:]);order=np.argsort(t)
        return dict(times_ns=t[order].tolist(),charge_pe=q[order].tolist(),
                    cost=float(np.sum(residual(z)**2)),converged=bool(result.success),
                    evaluations=int(result.nfev))

    def save(self, path):
        import json
        np.savez_compressed(path, grid=self.grid, coefficients=self.coefficients,
                            settings=json.dumps(self.settings), residual=json.dumps(self.residual))

    @classmethod
    def load(cls, path):
        import json
        with np.load(path, allow_pickle=False) as f:
            return cls(json.loads(str(f['settings'])), f['grid'], f['coefficients'],
                       json.loads(str(f['residual'])))


def calibrate(settings, components, times, charges, voltage, *, ridge=.03):
    """Development-only response regression; independent parents for validation.

    Components are actual avalanche ancestry waveforms. Nonlinear mean and
    empirical AR(1) residual covariance are separate approximations. The
    latter is frozen from an observed pilot during inference, so its logdet
    is constant across proposals and cancels from model comparisons.
    """
    base=PairFitter(settings);sample_grid=np.arange(1024)*2.+.25
    grid=np.arange(-16.,1720.01,2.)
    x=[];target=[];valid=[]
    for curves,t,q in zip(components,times,charges):
        ids=np.flatnonzero(np.asarray(q)>=500.)
        for j in ids:
            other=next((k for k in ids if k!=j),None)
            x.append(features(grid,q[j],q[other] if other is not None else 0.,
                              t[other]-t[j] if other is not None else 0.))
            aligned=np.interp(grid+t[j],sample_grid,curves[j],left=0.,right=0.)/q[j]
            target.append(aligned-np.interp(grid,base.grid,base.kernel,left=0.,right=0.))
            valid.append((grid+t[j]>=sample_grid[0])&(grid+t[j]<=sample_grid[-1]))
    x=np.asarray(x);target=np.asarray(target);valid=np.asarray(valid)
    # Late pulses do not supply their unobserved tail as a zero-valued target.
    weight=valid.astype(float)/np.maximum(valid.sum(0),1.)[None]
    gram=np.einsum('nsi,nsj,ns->sij',x,x,weight)
    rhs=np.einsum('nsi,ns,ns->si',x,target,weight)
    regularizer=np.eye(x.shape[-1])*ridge;regularizer[0,0]*=.01
    coef=np.linalg.solve(gram+regularizer[None],rhs[...,None])[...,0]
    coef=gaussian_filter1d(coef,2.,axis=0)  # 4 ns smoothness fixed in advance.
    response=ConditionalResponse(settings,grid,coef,dict(electronic_sigma_mV=.1,
        dark_variance_mV2=.01,shot_coefficient_mV=1.,rho=.9))
    residual=[];pilot=[]
    for y,t,q in zip(voltage,times,charges):
        ids=np.flatnonzero(np.asarray(q)>=500.)
        residual.append(np.asarray(y)-base.baseline-response.components(np.asarray(t)[ids],np.asarray(q)[ids]).sum(0))
        pilot.append(gaussian_filter1d(np.maximum(np.asarray(y)-base.baseline,0.),3.))
    residual=np.asarray(residual);pilot=np.asarray(pilot)
    # Robust moment estimates, excluding the baseline from shot calibration.
    mask=pilot>5.
    coefficient=float(np.median((residual[mask]**2)/pilot[mask])/.4549364231195727)
    coefficient=max(.001,coefficient)
    normalized=residual/np.sqrt(.02+coefficient*pilot)
    adjacent=mask[:,1:]&mask[:,:-1]
    a=normalized[:,:-1][adjacent];b=normalized[:,1:][adjacent]
    rho=float(np.clip(a@b/np.sqrt((a@a)*(b@b)),0.,.995))
    response.residual.update(shot_coefficient_mV=coefficient,rho=rho)
    return response,dict(component_rows=len(x),frames=len(times),ridge=ridge,
        smooth_sigma_ns=4.,mean_features=9,residual_calibration=response.residual,
        component_relative_rmse=float(np.sqrt(np.mean((np.einsum('nsi,si->ns',x,coef)-target)**2))),
        truth_used='Development calibration only: component ancestry, Q, first-active time',
        covariance='Empirical heteroscedastic AR(1), frozen from observed pilot per frame')


def select_count(candidates, penalty):
    norm=max(candidates[0]['cost'],1.)
    eligible=[c for c in candidates if all(q>=500. for q in c['charge_pe'])]
    return min(eligible,key=lambda c:c['cost']/norm+penalty*len(c['times_ns']))


def selective_refit(raw, refined, *, gap_threshold_ns, max_log_charge_change):
    """Frozen observable-only gate; never changes a network's event count."""
    if len(raw['times_ns'])!=len(refined['times_ns']):
        raise ValueError('A same-count proposal is required')
    if not len(raw['times_ns']):
        return dict(raw, accepted_refit=False)
    q=np.asarray(raw['charge_pe']);newq=np.asarray(refined['charge_pe'])
    gap=float(np.diff(np.sort(raw['times_ns']))[0]) if len(q)==2 else np.inf
    change=float(np.max(np.abs(np.log(np.maximum(newq,1.)/np.maximum(q,1.)))))
    accept=gap>=gap_threshold_ns and change<=max_log_charge_change and bool(refined.get('converged',True))
    return dict(refined if accept else raw,accepted_refit=bool(accept))
