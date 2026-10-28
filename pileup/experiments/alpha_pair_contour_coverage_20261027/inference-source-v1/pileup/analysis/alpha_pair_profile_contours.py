"""Refine finite charge-profile crossings with an explicit bracket check.

Only evaluated fits determine crossings and stopping. Truth is not an input.
This checks a sampled local contour, not global optimization or coverage.
"""
from __future__ import annotations

import copy
import time
import numpy as np

from .alpha_pair_profile_stability import refine_stage, compare_regions

LEVELS=(1.,4.,9.)


def crossing_cells(curve, levels=LEVELS):
    q=np.asarray(curve['q_grid_pe'],float)
    d=np.array([f['cost'] for f in curve['fits']])-curve['reference_cost']
    out=[]
    for level in levels:
        for i in np.flatnonzero((d[:-1]<=level)!=(d[1:]<=level)):
            out.append(dict(level=float(level),index=int(i),lower=float(q[i]),upper=float(q[i+1]),
                log_width=float(np.log(q[i+1]/q[i]))))
    return out


def local_nodes(curve,levels=LEVELS):
    q=np.asarray(curve['q_grid_pe'],float);ids=set()
    for r in crossing_cells(curve,levels):ids.update((r['index'],r['index']+1))
    best=int(np.argmin([f['cost'] for f in curve['fits']]))
    ids.update(range(max(0,best-1),min(len(q),best+2)))
    return q[sorted(ids)]


def refinement_targets(curve,levels=LEVELS,bracket_fraction=.01):
    q=np.asarray(curve['q_grid_pe'],float);tol=np.log1p(bracket_fraction)
    ids={r['index'] for r in crossing_cells(curve,levels) if r['log_width']>tol}
    best=int(np.argmin([f['cost'] for f in curve['fits']]))
    for i in (best-1,best):
        if 0<=i<len(q)-1 and np.log(q[i+1]/q[i])>tol:ids.add(i)
    return np.array([np.sqrt(q[i]*q[i+1]) for i in sorted(ids)])


def quality(curve,levels=LEVELS,bracket_fraction=.01):
    crossings=crossing_cells(curve,levels);tol=np.log1p(bracket_fraction)
    return dict(crossings=crossings,maximum_bracket_fraction=float(np.expm1(max((r['log_width'] for r in crossings),default=0.))),
        bracket_pass=all(r['log_width']<=tol+1e-12 for r in crossings),
        unresolved_cells=sum(r['log_width']>tol+1e-12 for r in crossings))


def converge_contours(response,noise,voltage,row,previous,settings,progress=None,refiner=refine_stage):
    start=time.monotonic();levels=settings['levels'];current=copy.deepcopy(previous)
    stages={'initial':current};history=[]
    def step(label,targets,refit,options):
        nonlocal current
        old=current
        current=refiner(response,noise,voltage,row,old,targets,refit_existing=refit,label=label,**options)
        before={float(q):f['cost'] for q,f in zip(old['q_grid_pe'],old['fits'])}
        after={float(q):f['cost'] for q,f in zip(current['q_grid_pe'],current['fits'])}
        assert all(after[q]<=cost+1e-8 for q,cost in before.items())
        new=[q for q in after if q not in before]
        deviations=[dict(Q2_pe=q,interpolated_cost=float(np.interp(np.log(q),np.log(old['q_grid_pe']),[f['cost'] for f in old['fits']])),
                         evaluated_cost=after[q]) for q in new]
        changes=[compare_regions(old,current,l,settings['reference_tolerance'],settings['endpoint_fraction']) for l in levels]
        record=dict(stage=label,points=len(after),added=len(new),old_node_max_improvement=max(before[q]-after[q] for q in before),
            reference_gain=old['reference_cost']-current['reference_cost'],comparisons=changes,
            quality=quality(current,levels,settings['bracket_fraction']),new_point_comparison=deviations,
            seconds=current.get('seconds'),trials=current.get('trials',0),nonconverged_trials=current.get('nonconverged_trials',0))
        stages[label]=current;history.append(record)
        if progress:progress(record)
        return record
    step('local_search',local_nodes(current,levels),True,settings['local'])
    stable=False;verified=0
    for iteration in range(1,settings['max_rounds']+1):
        targets=refinement_targets(current,levels,settings['bracket_fraction'])
        if len(targets):step(f'midpoint_{iteration:02d}',targets,False,settings['grid'])
        if quality(current,levels,settings['bracket_fraction'])['bracket_pass'] and not len(refinement_targets(current,levels,settings['bracket_fraction'])):
            record=step(f'verify_{iteration:02d}',local_nodes(current,levels),True,settings['verify']);verified+=1
            stable=record['quality']['bracket_pass'] and all(c['stable'] for c in record['comparisons'])
            if stable:break
    return dict(stages=stages,history=history,final_stage=history[-1]['stage'],stable=bool(stable),verification_passes=verified,
        final_quality=quality(current,levels,settings['bracket_fraction']),seconds=time.monotonic()-start,
        interpretation='Finite known-K2 crossing-bracket and stronger local-search check, not a global or statistical guarantee')
