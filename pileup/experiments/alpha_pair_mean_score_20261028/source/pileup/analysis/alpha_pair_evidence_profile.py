"""Observation-only count tradeoffs and finite, conditional K=2 Q2 profiles.

Q2 denotes the chronologically later slot, not the true weaker particle.
The profile is a numerical compatibility diagnostic, not a confidence interval.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from .alpha_pair_physical_gate import decide


def rank_grid(gates):
    values={0.,.025,.05,.1,.25,.5,1.}
    for gate in gates.values():
        for group in gate['scores']:
            values.update(j/(len(group)+1) for j in range(1,len(group)+2))
    return sorted(values)


def decision_endpoints(row,method,gate):
    low=decide(row,method,dict(gate,alpha=0.))
    high=decide(row,method,dict(gate,alpha=1.))
    assert len(low['times_ns'])<2
    return dict(low=low,high=high,eligible=len(high['times_ns'])==2,
                rank=high.get('null_rank_p'),evidence=high.get('evidence_score'))


def choose_endpoint(endpoints,alpha):
    return 'high' if endpoints['eligible'] and endpoints['rank']<=alpha else 'low'


def encode(times,q1,bounds):
    low,high=bounds;t=np.clip(np.sort(times),low,high)
    fraction=(t[1]-t[0])/max(high-t[0],1e-12)
    return np.array([t[0]/10.,100.*fraction,np.log(np.clip(q1,50.,50000.))])


def decode(z,q2,bounds):
    high=bounds[1];t1=float(z[0]*10.);t2=t1+(high-t1)*z[1]/100.
    return [t1,float(t2)],[float(np.exp(z[2])),float(q2)]


def fit_at_charge(response,noise,voltage,q2,seeds,bounds,*,maxiter=50,starts=2):
    if not 50.-1e-8<=q2<=50000.+1e-8:raise ValueError('Q2 outside declared fit domain')
    y=np.asarray(voltage,float)-response.baseline-noise.background_mV;zy=noise.transform(y)
    def objective(z):
        t,q=decode(z,q2,bounds);m=noise.transform(response.components(t,q))
        return float(noise.score(zy-m.sum(0),m,q,transformed=True))
    initial={}
    for seed in seeds:
        z=encode(seed['times_ns'],seed['charge_pe'][0],bounds)
        initial.setdefault(tuple(z.tolist()),z)
    ranked=sorted(((objective(z),z) for z in initial.values()),key=lambda pair:pair[0])[:starts]
    if not ranked:raise ValueError('At least one observation-derived start required')
    results=[];evaluations=len(initial)
    for initial_cost,z0 in ranked:
        fit=minimize(objective,z0,method='L-BFGS-B',bounds=[(bounds[0]/10.,bounds[1]/10.),(0.,100.),(np.log(50.),np.log(50000.))],
            options=dict(maxiter=maxiter,ftol=1e-9,gtol=1e-4,maxls=30,eps=1e-5))
        evaluations+=int(fit.nfev)
        use_fit=np.isfinite(fit.fun) and fit.fun<=initial_cost
        z=fit.x if use_fit else z0;t,q=decode(z,q2,bounds)
        results.append(dict(times_ns=t,charge_pe=q,cost=float(objective(z)),converged=bool(fit.success),
            retained_initial=not use_fit,message=str(fit.message),evaluations=int(fit.nfev),initial_cost=initial_cost))
    best=min(results,key=lambda r:r['cost'])
    return dict(best,search_evaluations=evaluations,trial_converged=[r['converged'] for r in results])


def seeds_for(one,base,learned,q2,bounds,warm=None):
    seeds=[]
    for candidate in (base,learned,warm):
        if candidate is None:continue
        seeds.append(candidate)
        seeds.append(dict(times_ns=candidate['times_ns'],charge_pe=[max(50.,sum(candidate['charge_pe'])-q2),q2]))
    t=one['times_ns'][0];q1=max(50.,sum(one['charge_pe'])-q2)
    for gap in (0.,4.,8.,16.,32.,64.,128.,256.,512.):
        seeds.append(dict(times_ns=[t,min(bounds[1],t+gap)],charge_pe=[q1,q2]))
    return seeds


def compatibility_segments(q,cost,reference,level):
    """Piecewise-linear interpolation in log Q; keep disconnected components."""
    q=np.asarray(q,float);d=np.asarray(cost,float)-reference
    if len(q)<2 or np.any(np.diff(q)<=0):raise ValueError('Sorted unique Q grid required')
    x=np.log(q);pieces=[]
    for i in range(len(q)-1):
        a,b=d[i]<=level,d[i+1]<=level
        if a and b:pieces.append([float(q[i]),float(q[i+1])])
        elif a or b:
            u=float((level-d[i])/(d[i+1]-d[i]));cross=float(np.exp(x[i]+u*(x[i+1]-x[i])))
            pieces.append([float(q[i]),cross] if a else [cross,float(q[i+1])])
    merged=[]
    for a,b in pieces:
        if merged and np.isclose(merged[-1][1],a,rtol=1e-10,atol=1e-10):merged[-1][1]=b
        else:merged.append([a,b])
    return dict(segments=merged,lower=merged[0][0] if merged else None,upper=merged[-1][1] if merged else None,
        lower_boundary=bool(d[0]<=level),upper_boundary=bool(d[-1]<=level),
        log_span=float(sum(np.log(b/a) for a,b in merged)),level=float(level))


def profile(response,noise,voltage,row,*,grid_count=33,maxiter=50,starts=2):
    base=row['sets']['baseline'][2];learned=row['sets']['mlp'][2];one=row['sets']['baseline'][1]
    trigger=one['times_ns'][0]
    bounds=[max(0.,trigger-96.),min(float(response.sample_grid[-1]),trigger+512.)]
    qs=np.unique(np.clip(np.r_[np.geomspace(50.,50000.,grid_count),500.,base['charge_pe'][1],
        np.asarray([.5,.8,1.,1.2,2.])*learned['charge_pe'][1]],50.,50000.))
    fits={};history=[];evaluations=0;nonconverged=0
    for direction in (1,-1):
        warm=None
        for q2 in qs[::direction]:
            seeds=seeds_for(one,base,learned,q2,bounds,warm)
            if float(q2) in fits:seeds.append(fits[float(q2)])
            fit=fit_at_charge(response,noise,voltage,float(q2),seeds,bounds,maxiter=maxiter,starts=starts)
            evaluations+=fit['search_evaluations'];nonconverged+=sum(not v for v in fit['trial_converged'])
            if float(q2) not in fits or fit['cost']<fits[float(q2)]['cost']:fits[float(q2)]=fit
            warm=fits[float(q2)]
        history.append([fits[float(q)]['cost'] for q in qs])
    best=min(fits.values(),key=lambda r:r['cost'])
    # Include an explicit fixed-nuisance slice in the finite envelope. This
    # guarantees that a failed optimizer cannot make a profile worse than it.
    y=np.asarray(voltage,float)-response.baseline-noise.background_mV;slices=[]
    for q2 in qs:
        t=best['times_ns'];q=[best['charge_pe'][0],float(q2)];m=response.components(t,q)
        cost=float(noise.score(y-m.sum(0),m,q));slices.append(cost)
        if cost<fits[float(q2)]['cost']:
            fits[float(q2)]=dict(times_ns=list(t),charge_pe=q,cost=cost,converged=False,retained_slice=True)
    values=[fits[float(q)] for q in qs];reference=min(r['cost'] for r in values)
    assert np.max(np.array([r['cost'] for r in values])-slices)<=1e-8
    return dict(bounds_ns=bounds,q_grid_pe=qs.tolist(),fits=values,reference_cost=reference,
        slice_cost=slices,slice_nuisance=dict(times_ns=best['times_ns'],Q1_pe=best['charge_pe'][0]),
        forward_cost=history[0],reverse_cost=history[1],evaluations=evaluations,nonconverged_trials=nonconverged,
        cached_outside_local_domain=bool(any(t<bounds[0] or t>bounds[1] for t in learned['times_ns'])),
        interpretation='Forced K2, ordered later Q2, trigger-local time domain, finite multistart profile; not a confidence interval')
