"""Observation-only charge-allocation coordinates and scalar-profile starts.

The likelihood, pulse physics and individual 50..50000 pe charge domain are
unchanged. No prior, true count, labels or true component enters this module.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from .alpha_pair_joint_likelihood import refine

QMIN, QMAX = 50., 50000.
GAPS = (4., 8., 12., 20., 40., 80., 160., 320.)
SHARES = (.01, .03, .05, .10, .15, .5, .85, .90, .95, .97, .99)
VARIANTS = ('component_local', 'allocation_local', 'component_profile', 'allocation_profile')


def forbid_labels(event,args):
    if event=='open':
        path=str(args[0])
        if '.truth.' in path or path.endswith(('.cases.json','.root','native-records.json')) or '/truth/' in path or '/native-inputs/' in path:
            raise PermissionError('Inference attempted to read labels: '+path)


def charge_interval(total):
    """Exactly retain the old individual-charge box at fixed total Q."""
    return max(QMIN, total-QMAX), min(QMAX, total-QMIN)


def encode_charges(charges):
    q=np.asarray(charges,float)
    if q.shape!=(2,) or not np.isfinite(q).all() or np.any(q<QMIN-1e-8) or np.any(q>QMAX+1e-8):
        raise ValueError('Two charges in the original physical fit domain required')
    total=float(q.sum());low,high=charge_interval(total)
    fraction=float((q[0]-low)/(high-low)) if high>low else .5
    return np.array([100.*np.log(total),10.*np.clip(fraction,0.,1.)])


def decode_charges(coordinates):
    x=np.asarray(coordinates,float)
    if x.shape!=(2,) or not np.isfinite(x).all():
        raise ValueError('Two finite allocation coordinates required')
    total=float(np.clip(np.exp(x[0]/100.),2*QMIN,2*QMAX))
    low,high=charge_interval(total)
    first=low+(high-low)*np.clip(x[1]/10.,0.,1.)
    return np.array([first,total-first])


def refine_allocation(response,noise,voltage,initial,*,radius=24.,maxiter=70):
    """Same objective/time bounds/iteration limits as the component control."""
    if len(initial['times_ns'])!=2:
        raise ValueError('This refinement requires a two-component candidate')
    y=np.asarray(voltage,float)-response.baseline-noise.background_mV
    order=np.argsort(initial['times_ns'])
    t=np.asarray(initial['times_ns'],float)[order]
    q=np.clip(np.asarray(initial['charge_pe'],float)[order],QMIN,QMAX)
    lower=np.r_[np.maximum(0.,t-radius)/10.,100.*np.log(2*QMIN),0.]
    upper=np.r_[np.minimum(response.sample_grid[-1],t+radius)/10.,100.*np.log(2*QMAX),10.]
    z0=np.clip(np.r_[t/10.,encode_charges(q)],lower+1e-8,upper-1e-8)
    transformed_y=noise.transform(y)

    def objective(z):
        charges=decode_charges(z[2:])
        curves=response.components(z[:2]*10.,charges)
        transformed=noise.transform(curves)
        return noise.score(transformed_y-transformed.sum(0),transformed,charges,transformed=True)

    fit=minimize(objective,z0,method='L-BFGS-B',bounds=list(zip(lower,upper)),
        options=dict(maxiter=maxiter,ftol=1e-9,gtol=1e-4,maxls=30,eps=1e-5))
    z=fit.x if fit.fun<=objective(z0) else z0
    t=z[:2]*10.;q=decode_charges(z[2:]);order=np.argsort(t);t=t[order];q=q[order]
    curves=response.components(t,q)
    return dict(times_ns=t.tolist(),charge_pe=q.tolist(),converged=bool(fit.success),
        evaluations=int(fit.nfev),message=str(fit.message),
        **noise.score(y-curves.sum(0),curves,q,details=True))


def profile_starts(response,noise,voltage,one):
    """At each observed-anchor time pair and share, numerically profile total Q.

    Bounded scalar minimization is a finite numerical search, not a proof of
    a global profile. Keep the unprofiled point too. No true timing is used.
    """
    y=np.asarray(voltage,float)-response.baseline-noise.background_mV
    zy=noise.transform(y)
    anchor=float(one['times_ns'][0]);total=float(sum(one['charge_pe']))
    profiles=[];groups={};evaluations=0
    for direction in (-1,1):
        for gap in GAPS:
            times=np.sort([anchor,anchor+direction*gap])
            if times[0]<0. or times[1]>response.sample_grid[-1]:continue
            for share in SHARES:
                # share is the charge fraction of the chronological first pulse.
                lower=QMIN/min(share,1.-share)
                upper=QMAX/max(share,1.-share)
                initial=np.log(np.clip(total,lower,upper))

                def objective(log_total):
                    q=np.exp(log_total)*np.array([share,1.-share])
                    curves=noise.transform(response.components(times,q))
                    return noise.score(zy-curves.sum(0),curves,q,transformed=True)

                initial_cost=float(objective(initial))
                fit=minimize_scalar(objective,method='bounded',bounds=(np.log(lower),np.log(upper)),
                                    options=dict(maxiter=24,xatol=1e-4))
                best,cost=(float(fit.x),float(fit.fun)) if fit.fun<initial_cost else (initial,initial_cost)
                q=np.exp(best)*np.array([share,1.-share])
                proposal=dict(times_ns=times.tolist(),charge_pe=q.tolist())
                row=dict(direction=direction,gap_ns=gap,share_first=share,total_charge_pe=float(q.sum()),
                    cost=cost,unprofiled_cost=initial_cost,converged=bool(fit.success),evaluations=int(fit.nfev)+1,
                    **proposal)
                profiles.append(row);evaluations+=row['evaluations']
                key=(direction,gap)
                if key not in groups or cost<groups[key]['cost']:groups[key]=row
    starts=[]
    for direction in (-1,1):
        ranked=sorted([r for k,r in groups.items() if k[0]==direction],key=lambda r:r['cost'])
        starts += [dict(times_ns=r['times_ns'],charge_pe=r['charge_pe']) for r in ranked[:4]]
    return dict(starts=starts,profiles=profiles,evaluations=evaluations)


def extend_allocation(response,noise,voltage,baseline):
    """Return both coordinates, each with local-only and common-profile search.

    All variants retain the baseline K=2 candidate and exactly the K=0/1
    models. They share profile starts and iteration limits; actual objective
    evaluations and incremental runtime are reported rather than assumed equal.
    """
    import time
    if len(baseline)<3:
        return dict(variants={name:dict(candidates=baseline,evaluations=0,seconds=0.,selected_source='baseline') for name in VARIANTS},
                    profile=dict(starts=[],profiles=[],evaluations=0,seconds=0.))
    started=time.monotonic();bank=profile_starts(response,noise,voltage,baseline[1]);bank['seconds']=time.monotonic()-started
    variants={}
    for coordinate,optimizer in [('component',refine),('allocation',refine_allocation)]:
        started=time.monotonic()
        local=optimizer(response,noise,voltage,baseline[2],radius=24.,maxiter=70)
        best=min([baseline[2],local],key=lambda r:r['cost'])
        variants[coordinate+'_local']=dict(candidates=[*baseline[:2],best],
            evaluations=local['evaluations'],seconds=time.monotonic()-started,
            selected_source='baseline' if best is baseline[2] else 'local',trials=[local])
        started=time.monotonic()
        screened=[optimizer(response,noise,voltage,s,radius=24.,maxiter=20) for s in bank['starts']]
        finalists=sorted([baseline[2],local]+screened,key=lambda r:r['cost'])[:2]
        polished=[optimizer(response,noise,voltage,s,radius=24.,maxiter=70) for s in finalists]
        candidates=[baseline[2],local]+screened+polished
        best=min(candidates,key=lambda r:r['cost'])
        selected=candidates.index(best)
        variants[coordinate+'_profile']=dict(candidates=[*baseline[:2],best],
            evaluations=sum(r['evaluations'] for r in [local]+screened+polished),
            seconds=variants[coordinate+'_local']['seconds']+time.monotonic()-started,
            selected_source=('baseline' if selected==0 else 'local' if selected==1 else 'screen' if selected<2+len(screened) else 'polish'),
            trials=[local]+screened+polished,screened=len(screened),polished=len(polished))
    return dict(variants=variants,profile=bank)
