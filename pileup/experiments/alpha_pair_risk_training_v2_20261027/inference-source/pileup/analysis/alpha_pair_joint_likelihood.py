"""Candidate-dependent Gaussian surrogate for blind alpha-pair inference.

The fixed electronics transform is not a scintillation lifetime fit. Its first
three outputs are marginalized. Covariance includes short-range shot noise,
cross-component noise, and a low-rank charge constraint. All candidate scores
include log(det C); they are not the previous normalized GLS costs.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import cholesky_banded, cho_solve_banded
from scipy.optimize import minimize
from scipy.signal import lfilter

from .alpha_pair_conditional import ConditionalResponse, features


class ConservingResponse(ConditionalResponse):
    """Normalize every conditional component to realized charge times gain."""
    def __post_init__(self):
        super().__post_init__()
        self._area_grid=np.arange(-16.,6000.001,1.)
        self._base_area=float(np.trapezoid(np.interp(self._area_grid,self.base.grid,self.base.kernel),self._area_grid))
        self._area_base=np.interp(self.grid,self.base.grid,self.base.kernel,left=0.,right=0.)
        onset=np.clip(self.grid/8.,0.,1.)
        self._taper=onset*onset*(3.-2.*onset)*np.clip((self.grid[-1]-self.grid)/120.,0.,1.)

    def unit(self,relative_time,charge,other_charge=0.,other_delay=0.):
        correction=np.sum(features(self.grid,charge,other_charge,other_delay)*self.coefficients,axis=1)
        changed=np.maximum(0.,self._area_base+correction*self._taper)-self._area_base
        area=self._base_area+float(np.trapezoid(changed,self.grid))
        if area<=0.:raise ValueError('Nonpositive conditional template area')
        s=np.asarray(relative_time,float);base=np.interp(s,self.base.grid,self.base.kernel,left=0.,right=0.)
        onset=np.clip(s/8.,0.,1.)
        taper=onset*onset*(3.-2.*onset)*np.clip((self.grid[-1]-s)/120.,0.,1.)
        values=np.maximum(0.,base+np.interp(s,self.grid,correction,left=0.,right=0.)*taper)
        return values*self.base.area_per_pe/area


def band_rowsum(band):
    """Row sums of a symmetric matrix in scipy lower-banded storage."""
    out=band[0].copy();n=len(out)
    for k in range(1,len(band)):
        out[k:]+=band[k,:n-k];out[:n-k]+=band[k,:n-k]
    return out


@dataclass
class JointNoise:
    coefficients: np.ndarray
    noise_lags: np.ndarray
    parameters: np.ndarray
    lag1: float
    cross: float
    kappa: float
    background_mV: float = 0.

    def transform(self,x):
        # Explicit short FIR also accepts the empty component axis at K=0.
        x=np.asarray(x,float);order=len(self.coefficients)-1;n=x.shape[-1]
        return sum(c*x[...,order-k:n-k] for k,c in enumerate(self.coefficients))

    def covariance(self,curves,q,*,transformed=False):
        """Banded base minus V V^T, with a guaranteed PSD joint construction."""
        m=np.asarray(curves,float);q=np.asarray(q,float)
        z=m if transformed else self.transform(m)
        n=z.shape[-1];k=len(q)
        base=np.zeros((len(self.noise_lags),n))
        for lag,value in enumerate(self.noise_lags):base[lag,:n-lag]=value
        if k==0:return base,np.empty((n,0))
        if k>2:raise ValueError('Local K <= 2 required')
        # Features depend only on the proposed components and charges.
        a,b,c=self.parameters
        activity=np.maximum(z,0.)
        neighbor=activity[::-1] if k==2 else np.zeros_like(activity)
        scale=np.exp(np.clip(a+b*np.log(np.maximum(q[:,None],50.)/3000.)+c*neighbor/10.,-8.,4.))
        v=scale*activity;root=np.sqrt(v)
        bands=[];row_sums=[]
        for j in range(k):
            band=np.zeros_like(base);band[0]=v[j]
            band[1,:n-1]=self.lag1*root[j,1:]*root[j,:-1]
            bands.append(band);row_sums.append(band_rowsum(band))
            base+=band
        S=np.diag([x.sum() for x in row_sums]);u=np.column_stack(row_sums)
        if k==2 and self.cross:
            # C12 = cross * diag(sqrt(v1)) R diag(sqrt(v2)).
            diag=self.cross*root[0]*root[1]
            lower=self.cross*self.lag1*root[0,1:]*root[1,:-1]
            upper=self.cross*self.lag1*root[0,:-1]*root[1,1:]
            rows=diag.copy();cols=diag.copy()
            rows[1:]+=lower;rows[:-1]+=upper
            cols[:-1]+=lower;cols[1:]+=upper
            base[0]+=2*diag;base[1,:n-1]+=lower+upper
            u[:,0]+=cols;u[:,1]+=rows
            S[0,1]=S[1,0]=rows.sum()
        if not self.kappa:return base,np.empty((n,0))
        # S is the covariance of component totals before conditioning.
        # Cbase - U inv(S) U.T is PSD; 0 <= kappa < 1 keeps it positive.
        S+=np.eye(k)*1e-12
        chol=np.linalg.cholesky(S)
        V=np.linalg.solve(chol,u.T).T*np.sqrt(self.kappa)
        return base,V

    def score(self,residual,curves,q,*,transformed=False,details=False):
        r=np.asarray(residual,float) if transformed else self.transform(residual)
        band,V=self.covariance(curves,q,transformed=transformed)
        factor=cholesky_banded(band,lower=True,check_finite=False)
        solved=cho_solve_banded((factor,True),np.column_stack((r,V)),check_finite=False)
        quadratic=float(r@solved[:,0]);logdet=float(2*np.log(factor[0]).sum())
        if V.shape[1]:
            H=np.eye(V.shape[1])-V.T@solved[:,1:]
            sign,det=np.linalg.slogdet(H)
            if sign!=1:raise np.linalg.LinAlgError('Invalid charge-conditioning Schur complement')
            rhs=V.T@solved[:,0]
            quadratic+=float(rhs@np.linalg.solve(H,rhs));logdet+=float(det)
        return dict(cost=quadratic+logdet,quadratic=quadratic,logdet=logdet) if details else quadratic+logdet

    def save(self,path):
        np.savez_compressed(path,coefficients=self.coefficients,noise_lags=self.noise_lags,
            parameters=self.parameters,lag1=self.lag1,cross=self.cross,kappa=self.kappa,background_mV=self.background_mV)

    @classmethod
    def load(cls,path):
        with np.load(path) as f:
            return cls(f['coefficients'],f['noise_lags'],f['parameters'],float(f['lag1']),float(f['cross']),float(f['kappa']),float(f['background_mV']))


def refine(response,noise,voltage,initial,*,radius=24.,maxiter=70):
    y=np.asarray(voltage,float)-response.baseline-noise.background_mV;n=len(initial['times_ns'])
    if n==0:
        return dict(times_ns=[],charge_pe=[],converged=True,evaluations=1,
            **noise.score(y,np.empty((0,len(y))),[],details=True))
    t=np.asarray(initial['times_ns']);q=np.asarray(initial['charge_pe']);order=np.argsort(t);t=t[order];q=q[order]
    lower=np.r_[np.maximum(0.,t-radius)/10.,np.full(n,np.log(50.))]
    upper=np.r_[np.minimum(response.sample_grid[-1],t+radius)/10.,np.full(n,np.log(50000.))]
    z0=np.clip(np.r_[t/10.,np.log(np.maximum(q,50.))],lower+1e-8,upper-1e-8)
    zy=noise.transform(y)
    def objective(z):
        charges=np.exp(z[n:]);curves=response.components(z[:n]*10.,charges)
        transformed=noise.transform(curves)
        return noise.score(zy-transformed.sum(0),transformed,charges,transformed=True)
    fit=minimize(objective,z0,method='L-BFGS-B',bounds=list(zip(lower,upper)),
        options=dict(maxiter=maxiter,ftol=1e-9,gtol=1e-4,maxls=30,eps=1e-5))
    z=fit.x if fit.fun<=objective(z0) else z0
    t=z[:n]*10.;q=np.exp(z[n:]);order=np.argsort(t);t=t[order];q=q[order]
    curves=response.components(t,q)
    return dict(times_ns=t.tolist(),charge_pe=q.tolist(),converged=bool(fit.success),
        evaluations=int(fit.nfev),message=str(fit.message),**noise.score(y-curves.sum(0),curves,q,details=True))


def waveform_candidates(response,noise,voltage):
    """Voltage-only starts; new NLL ranks and refines a declared finite bank."""
    from .alpha_pair_residual import PrecisionBlend, residual_starts, proposal_from_times
    y=np.asarray(voltage,float)
    raw=response.base.candidates(y)
    zero=refine(response,noise,y,dict(times_ns=[],charge_pe=[]))
    if len(raw)==1:return [zero]
    one=refine(response,noise,y,raw[1])
    two=refine(response,noise,y,raw[2])
    # The legacy weighted tangent scan proposes times only; the NLL selects them.
    metric=PrecisionBlend(response,y,0.,rho=.9)
    births=residual_starts(response,y,one,metric,maximum=4)
    time=float(one['times_ns'][0]);charge=float(one['charge_pe'][0]);split=[]
    observed=y-response.baseline-noise.background_mV
    for gap in (4.,8.,12.,20.,40.,80.,160.,320.):
        for direction in (-1,1):
            times=[time,time+direction*gap]
            if min(times)<0. or max(times)>response.sample_grid[-1]:continue
            group=[]
            for share in (.15,.5,.85):
                initial=proposal_from_times(times,[charge*share,charge*(1.-share)])
                curves=response.components(initial['times_ns'],initial['charge_pe'])
                group.append((noise.score(observed-curves.sum(0),curves,initial['charge_pe']),initial))
            split.append(min(group,key=lambda z:z[0]))
    starts=births+[r[1] for r in sorted(split,key=lambda r:r[0])[:4]]
    trial=[refine(response,noise,y,s,maxiter=20) for s in starts]
    best=sorted([two]+trial,key=lambda r:r['cost'])[:2]
    trial += [refine(response,noise,y,s) for s in best]
    return [zero,one,min([two]+trial,key=lambda r:r['cost'])]


def evidence(candidates):
    if len(candidates)<3 or any(q<500. for q in candidates[2]['charge_pe']):return 0.
    return max(0.,candidates[1]['cost']-candidates[2]['cost'])


def calibrate_gate(single_candidates,*,edges=(1500.,5000.,12500.),alpha=.025):
    scores=[[] for _ in range(len(edges)+1)]
    for candidates in single_candidates:
        q=sum(candidates[1]['charge_pe']) if len(candidates)>1 else 0.
        scores[int(np.searchsorted(edges,q,side='right'))].append(evidence(candidates))
    if any(not values for values in scores):raise ValueError('Every fitted-Q bin needs native singles')
    return dict(edges_pe=list(edges),alpha=alpha,scores=[sorted(x) for x in scores],
        score='max(NLL1-NLL2,0); no division by a possibly negative NLL',
        presence='Q>=500pe and NLL1<NLL0; no true K input')


def select(candidates,gate):
    zero=candidates[0]
    if len(candidates)<2 or candidates[1]['charge_pe'][0]<500. or candidates[1]['cost']>=zero['cost']:
        return dict(zero,birth_accepted=False,null_rank_p=None)
    one=candidates[1];bin_id=int(np.searchsorted(gate['edges_pe'],sum(one['charge_pe']),side='right'))
    score=evidence(candidates);values=np.asarray(gate['scores'][bin_id])
    rank=(1+np.count_nonzero(values>=score))/(1+len(values))
    accepted=bool(len(candidates)==3 and score>0. and rank<=gate['alpha'])
    return dict(candidates[2] if accepted else one,birth_accepted=accepted,null_rank_p=float(rank),
        evidence_score=score,calibration_bin=bin_id)
