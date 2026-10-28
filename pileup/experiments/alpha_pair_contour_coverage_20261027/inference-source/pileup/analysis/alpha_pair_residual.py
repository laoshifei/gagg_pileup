"""Regularized precision and residual proposals for blind local alpha fits.

The response calibration is reused unchanged. A fixed observation defines one
positive precision matrix per frame, shared by every count/time/charge trial.
No truth, source energy, rate or actual microcell state is an inference input.
"""
from __future__ import annotations

import numpy as np
from scipy.signal import find_peaks

from .alpha_correlated_fit import InnovationsMetric, ResidualCalibration


class PrecisionBlend:
    """Convex combination of normalized GLS and ordinary LS precisions.

    This regularizes an empirical precision; it is not an exact likelihood.
    Both normalizations use the fixed observed baseline-subtracted voltage.
    The zero-signal candidate has cost 1e6 for nonempty inputs. This common
    numerical scale keeps the existing optimizer's absolute tolerances useful;
    count selection divides by that zero-model cost again.
    """
    def __init__(self, response, voltage, fraction, *, rho=None):
        if not 0. <= fraction <= 1.:
            raise ValueError('LS fraction must be in [0,1]')
        self.fraction=float(fraction)
        self.correlated=response.metric(voltage)
        if rho is not None:
            params=dict(response.residual,rho=float(rho))
            self.correlated=InnovationsMetric(self.correlated.pilot,ResidualCalibration(**params))
        observed=np.asarray(voltage,float)-response.baseline
        weighted=self.correlated.apply(observed)
        self.gls_scale=np.sqrt(1e6*(1.-fraction)/max(float(weighted@weighted),1.))
        self.ols_scale=np.sqrt(1e6*fraction/max(float(observed@observed),1.))

    def apply(self, values):
        x=np.asarray(values,float)
        if self.fraction==0.:
            return self.gls_scale*self.correlated.apply(x)
        if self.fraction==1.:
            return self.ols_scale*x
        return np.concatenate((self.gls_scale*self.correlated.apply(x),self.ols_scale*x),axis=0)


def proposal_from_times(times,charges):
    order=np.argsort(times)
    return dict(times_ns=np.asarray(times,float)[order].tolist(),
                charge_pe=np.asarray(charges,float)[order].tolist())


def residual_starts(response, voltage, one, metric, *, maximum=4):
    """Screen possible births after projecting out the fitted pulse tangent.

    A fixed linear template only proposes locations; every accepted candidate
    is jointly refitted with the nonlinear conditional response afterwards.
    The existing pulse's amplitude/time directions are nuisance parameters.
    """
    if len(one['times_ns'])!=1:
        raise ValueError('One fitted event is required for a residual scan')
    t=float(one['times_ns'][0]);q=float(one['charge_pe'][0])
    grid=np.arange(0.,response.sample_grid[-1]-32.,4.)
    model=response.components([t],[q]).sum(0)
    residual=np.asarray(voltage,float)-response.baseline-model
    slopes=np.column_stack((
        (response.components([t],[q*1.001]).sum(0)-model)/.001,
        (response.components([t+.05],[q]).sum(0)-model)/.05))
    tangent=metric.apply(slopes)
    u,s,_=np.linalg.svd(tangent,full_matrices=False)
    u=u[:,s>max(float(s[0]),1.)*1e-10]
    bank=metric.apply(response.base.columns(response.sample_grid,grid))
    wr=metric.apply(residual)
    # Partial regression removes an existing event's local fit directions.
    bank-=u@(u.T@bank);wr-=u@(u.T@wr)
    norm=np.einsum('ij,ij->j',bank,bank)
    dot=bank.T@wr
    gain=np.maximum(dot,0.)**2/np.maximum(norm,1e-30)
    gain[np.abs(grid-t)<2.]=0.
    peaks,_=find_peaks(gain,distance=2)
    choices=sorted(peaks,key=lambda j:(-float(gain[j]),int(j)))[:maximum]
    starts=[]
    for j in choices:
        if gain[j]<=0.:continue
        added=float(np.clip(dot[j]/max(norm[j],1e-30),500.,25000.))
        starts.append(proposal_from_times([t,float(grid[j])],[q,added]))
    return starts


def choose(candidates, penalty):
    norm=max(candidates[0]['cost'],1.)
    eligible=[r for r in candidates if all(q>=500. for q in r['charge_pe'])]
    return min(eligible,key=lambda r:(r['cost']/norm+penalty*len(r['times_ns']),len(r['times_ns'])))


def fit_candidates(response, voltage, seeds, *, fraction, rho=None, maximum_births=4, metric=None):
    """Refit existing candidates and search births/splits, up to local K=2.

    The neural probabilities are not used as a count prior. Pure waveform
    candidates and neural-assisted candidates share the same final criterion.
    """
    y=np.asarray(voltage,float)
    if metric is None:
        metric=PrecisionBlend(response,y,fraction,rho=rho)
    original=seeds['candidates']
    zero=response.refine(y,dict(times_ns=[],charge_pe=[]),metric=metric)
    if len(original)==1:
        return dict(existing=[zero],expanded=[zero],waveform_only=[zero],
                    hybrids={'1':zero,'2':zero},birth_starts=[],split_starts=[],fraction=fraction,rho=rho)
    one=response.refine(y,original[1],metric=metric,radius=24.)
    two=response.refine(y,original[2],metric=metric,radius=24.)
    hybrids={str(j):response.refine(y,seeds[f'ml_{j}'],metric=metric,radius=24.) for j in (1,2)}
    existing=[zero,one,two]
    for f in hybrids.values():
        k=len(f['times_ns'])
        if f['cost']<existing[k]['cost']:existing[k]=f
    births=residual_starts(response,y,one,metric,maximum=maximum_births)
    t=float(one['times_ns'][0]);q=float(one['charge_pe'][0])
    split=[]
    # A small discrete search also covers unresolved peaks with no residual
    # maximum after nuisance projection. No true separation is supplied.
    for gap in (4.,8.,12.,20.,40.,80.,160.,320.):
        for direction in (-1,1):
            times=[t,t+direction*gap]
            if min(times)<0. or max(times)>response.sample_grid[-1]:continue
            group=[]
            for share in (.15,.5,.85):
                initial=proposal_from_times(times,[q*share,q*(1.-share)])
                r=metric.apply(response.components(initial['times_ns'],initial['charge_pe']).sum(0)-(y-response.baseline))
                group.append((float(r@r),initial))
            split.append(min(group,key=lambda z:z[0]))
    split=[r[1] for r in sorted(split,key=lambda z:z[0])[:4]]
    # Screen difficult starts with a bounded local solve, then polish only the
    # best two. This is a declared search budget, not a global optimum claim.
    fitted=[response.refine(y,s,metric=metric,radius=24.,max_nfev=20) for s in births+split]
    screened=sorted(fitted,key=lambda r:r['cost'])[:2]
    fitted += [response.refine(y,s,metric=metric,radius=24.) for s in screened]
    best=min([two]+fitted,key=lambda r:r['cost'])
    waveform=[zero,one,best]
    expanded=[dict(r) for r in existing]
    if best['cost']<expanded[2]['cost']:expanded[2]=best
    return dict(existing=existing,expanded=expanded,waveform_only=waveform,hybrids=hybrids,
                birth_starts=births,split_starts=split,fraction=fraction,rho=rho)
