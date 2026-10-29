"""Observation-only polishing AFTER the frozen blind event-count decision.

This is an optimization ablation, not a new count selector or a proof of global
optimality. The original fit is always retained as a feasible competitor.
"""
from __future__ import annotations

import copy
import numpy as np
from scipy.optimize import minimize

from .alpha_pair_blind_uncertainty import (
    box_encode, box_decode, domain_for, in_domain, waveform_objective,
)

DEFAULTS = dict(starts=6, maxiter=160, powell_maxiter=100, powell_maxfev=2400,
                ftol=1e-12, gtol=1e-6, eps=1e-5)


def seed_bank(one, candidate, domain):
    """Starts use fitted parameters only; no true gap, charge or event labels."""
    t = np.asarray(candidate['times_ns'], float)
    q = np.asarray(candidate['charge_pe'], float)
    seeds = [dict(times_ns=t.tolist(), charge_pe=q.tolist())]
    for coordinate in (0, 1):
        for step in (-32., -16., -8., -4., -2., -1., 1., 2., 4., 8., 16., 32.):
            tt=t.copy(); tt[coordinate]+=step; order=np.argsort(tt)
            seeds.append(dict(times_ns=np.clip(tt[order], *domain['time']).tolist(),
                              charge_pe=q[order].tolist()))
    anchor=float(one['times_ns'][0]); total=float(sum(one['charge_pe']))
    for gap in (4., 8., 16., 32., 64., 128., 256., 512.):
        for direction in (-1., 1.):
            tt=np.clip(sorted([anchor, anchor+direction*gap]), *domain['time']).tolist()
            for share in (.02, .05, .15, .5, .85, .95, .98):
                seeds.append(dict(times_ns=tt, charge_pe=np.clip(
                    [share*total, (1.-share)*total], *domain['charge']).tolist()))
    return seeds


def polish_objective(objective, one, candidate, domain, settings=None):
    cfg=dict(DEFAULTS if settings is None else settings)
    if not in_domain(candidate['times_ns'],candidate['charge_pe'],domain):
        raise ValueError('Original fit outside observation-derived domain')
    limits=[domain['time'][:], domain['time'][:],
            np.log(domain['charge']).tolist(), np.log(domain['charge']).tolist()]
    bounds=[(limits[0][0]/10.,limits[0][1]/10.),(0.,100.),limits[2],limits[3]]
    calls=0; failures=[]; trials=[]
    def decode(z):
        t,q=box_decode(z,limits)
        return np.clip(t,*domain['time']).tolist(),np.clip(q,*domain['charge']).tolist()
    def fun(z):
        nonlocal calls
        calls+=1
        t,q=decode(z)
        try:
            value=float(objective(t,q))
            return value if np.isfinite(value) else np.inf
        except (ArithmeticError,ValueError,np.linalg.LinAlgError):
            return np.inf
    def point(z, **extra):
        t,q=decode(z); value=fun(z)
        return dict(times_ns=t,charge_pe=q,cost=value,**extra)
    z0=box_encode(candidate['times_ns'],candidate['charge_pe'],limits)
    best=point(z0,converged=bool(candidate.get('converged',False)),origin='original')
    if not np.isfinite(best['cost']): raise ValueError('Nonfinite reference cost')
    original_cost=best['cost']
    def run(z, method, name):
        nonlocal best
        initial=point(z,converged=False,origin=name+'_seed')
        if not np.isfinite(initial['cost']):
            failures.append(dict(stage=name,reason='nonfinite seed')); return
        if initial['cost']<best['cost']: best=initial
        options=(dict(maxiter=cfg['maxiter'],ftol=cfg['ftol'],gtol=cfg['gtol'],
                      eps=cfg['eps'],maxls=40) if method=='L-BFGS-B' else
                 dict(maxiter=cfg['powell_maxiter'],maxfev=cfg['powell_maxfev'],
                      ftol=1e-10,xtol=1e-5))
        try:
            fit=minimize(fun,z,method=method,bounds=bounds,options=options)
            valid=np.isfinite(fit.x).all() and all(a-1e-8<=v<=b+1e-8
                    for v,(a,b) in zip(fit.x,bounds))
            proposed=point(fit.x,converged=bool(fit.success),origin=name) if valid else None
            if proposed and np.isfinite(proposed['cost']) and proposed['cost']<best['cost']:
                best=proposed
            trials.append(dict(stage=name,initial_cost=initial['cost'],
                final_cost=proposed['cost'] if proposed and np.isfinite(proposed['cost']) else None,
                converged=bool(fit.success),evaluations=int(fit.nfev),message=str(fit.message)))
        except (ArithmeticError,ValueError,np.linalg.LinAlgError) as exc:
            failures.append(dict(stage=name,reason=type(exc).__name__+': '+str(exc)))
    run(z0,'L-BFGS-B','restart')
    restart=copy.deepcopy(best)
    unique={}
    for seed in seed_bank(one,candidate,domain)+[restart]:
        z=box_encode(seed['times_ns'],seed['charge_pe'],limits)
        unique.setdefault(tuple(z),z)
    ranked=sorted(((fun(z),z) for z in unique.values()),key=lambda x:x[0])
    chosen=[]
    for value,z in ranked:
        if not np.isfinite(value): continue
        t,q=box_decode(z,limits); position=np.r_[np.asarray(t)/4.,np.log(q)/.15]
        if not chosen or min(np.linalg.norm(position-p) for p in chosen)>=1.:
            chosen.append(position); run(z,'L-BFGS-B',f'multistart_{len(chosen)}')
        if len(chosen)==cfg['starts']: break
    multistart=copy.deepcopy(best)
    z=box_encode(best['times_ns'],best['charge_pe'],limits)
    run(z,'Powell','powell')
    polished=copy.deepcopy(best)
    stages=dict(restart=restart,multistart=multistart,polished=polished)
    costs=[original_cost]+[s['cost'] for s in stages.values()]
    assert all(b<=a+1e-8 for a,b in zip(costs,costs[1:]))
    assert all(in_domain(s['times_ns'],s['charge_pe'],domain) for s in stages.values())
    return dict(stages=stages,trials=trials,failures=failures,objective_calls=calls,
                original_cost=original_cost,domain=domain,truth_inputs=[],
                global_optimum_claim=False)


def waveform_polish(response,noise,voltage,candidates,selected,settings=None):
    """Keep the original unknown-K decision, including any missed/extra events."""
    if len(selected['times_ns'])!=2:
        return dict(stages={n:copy.deepcopy(selected) for n in ('restart','multistart','polished')},
                    searched=False,truth_inputs=[],count_unchanged=True)
    objective=waveform_objective(response,noise,voltage)
    domain=domain_for(response,candidates[1],selected)
    result=polish_objective(objective,candidates[1],selected,domain,settings)
    y=np.asarray(voltage)-response.baseline-noise.background_mV
    for stage in result['stages'].values():
        curves=response.components(stage['times_ns'],stage['charge_pe'])
        detail=noise.score(y-curves.sum(0),curves,stage['charge_pe'],details=True)
        assert abs(detail['cost']-stage['cost'])<1e-7*max(1.,abs(detail['cost']))
        stage.update(detail)
    result.update(searched=True,count_unchanged=True)
    return result
