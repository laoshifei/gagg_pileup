"""Grouped, conditional stability checks for the frozen mean-shape recipe.

Only the incremental correction is cross-fitted. The old response, spatial
prior and nominal influence operator remain frozen and have prior exposure.
"""
from __future__ import annotations

import hashlib
import numpy as np

from .alpha_mean_shape_correction import calibrate_profiles


def mother_folds(identities, groups, seed, folds=4):
    identities=list(identities);groups=np.asarray(groups)
    if len(set(identities))!=len(identities):raise ValueError('Repeated optical mother identity')
    output=np.full(len(identities),-1,int)
    for group in sorted(set(groups)):
        indices=np.flatnonzero(groups==group).tolist()
        if len(indices)%folds:raise ValueError('Groups must divide evenly among folds')
        indices.sort(key=lambda i:hashlib.sha256((str(seed)+':'+identities[i]).encode()).hexdigest())
        output[indices]=np.arange(len(indices))%folds
    return output


def bootstrap_mothers(groups, folds, seed):
    """Resample whole mothers within fixed group/fold strata.

    Every copy of a mother remains in the same fold, including its eight
    digitization repeats. Thus bootstrap duplicates cannot leak into training
    for their own held-out evaluation.
    """
    groups,folds=np.asarray(groups),np.asarray(folds)
    rng=np.random.default_rng(seed);indices=[]
    for group in sorted(set(groups)):
        for fold in sorted(set(folds)):
            pool=np.flatnonzero((groups==group)&(folds==fold))
            if not len(pool):raise ValueError('Empty group/fold stratum')
            indices.extend(rng.choice(pool,len(pool),replace=True))
    return np.asarray(indices,int)


def many_profiles(correction, charge, apply_scale=True):
    charge=np.asarray(charge,float)
    x=np.clip(np.log(charge),*correction.log_charge_knots[[0,-1]])
    j=np.clip(np.searchsorted(correction.log_charge_knots,x,side='right')-1,0,len(correction.log_charge_knots)-2)
    fraction=(x-correction.log_charge_knots[j])/(correction.log_charge_knots[j+1]-correction.log_charge_knots[j])
    values=charge[:,None]*((1-fraction[:,None])*correction.profiles[j]+fraction[:,None]*correction.profiles[j+1])
    return values*(correction.scale if apply_scale else 1.)


def interpolate_many(profiles, grid, target):
    """Vectorized np.interp with zero endpoints, preserving the old recipe."""
    target=np.asarray(target);j=np.searchsorted(grid,target,side='right')-1
    inside=(target>=grid[0])&(target<=grid[-1]);j=np.clip(j,0,len(grid)-2)
    f=(target-grid[j])/(grid[j+1]-grid[j])
    return ((1-f)*profiles[:,j]+f*profiles[:,j+1])*inside


def many_voltage(correction, charge, time, onset=128.):
    return interpolate_many(many_profiles(correction,charge),correction.grid,np.asarray(time)-onset)


def cached_positivity_scale(correction, cache):
    curves=many_profiles(correction,cache['q0'],apply_scale=False)
    bound=1.
    for j,phase in enumerate(cache['phases_ns']):
        delta=interpolate_many(curves,correction.grid,cache['time_ns']-128.-phase)
        negative=delta<0
        if np.any(negative):bound=min(bound,float(np.min(cache['old_mV'][j][negative]/(-delta[negative]))))
    return 1. if bound>=1 else max(0.,.99*bound)


def fit_update(data, train, cache=None):
    """Labels and residuals indexed only at training mothers; fixed recipe."""
    train=np.asarray(train,int)
    if set(data['groups'][train])!={0,1,2,3}:raise ValueError('All four training energy groups are required')
    correction,_=calibrate_profiles(data['time_ns'],data['means_mV'][train],data['old_predictions_mV'][train],
                                    data['q0'][train],data['groups'][train])
    if cache is not None:correction.scale=cached_positivity_scale(correction,cache)
    return correction


def residual_statistics(old, new, groups, operator):
    groups=np.asarray(groups);high=groups==3
    before=100*old[high]@operator.T;after=100*new[high]@operator.T
    projected=[]
    for j,gap in enumerate((8.,40.,200.)):
        a,b=before[:,j],after[:,j]
        projected.append(dict(gap_ns=gap,old_mean=float(a.mean()),corrected_mean=float(b.mean()),
            old_mse=float(np.mean(a*a)),corrected_mse=float(np.mean(b*b)),
            mse_ratio=float(np.mean(b*b)/np.mean(a*a))))
    waveform=[]
    for g in range(4):
        a,b=old[groups==g],new[groups==g]
        waveform.append(dict(group=g,old_mean_mse=float(np.mean(a.mean(0)**2)),
            corrected_mean_mse=float(np.mean(b.mean(0)**2)),
            mean_mse_ratio=float(np.mean(b.mean(0)**2)/np.mean(a.mean(0)**2)),
            old_mother_mse=float(np.mean(a*a)),corrected_mother_mse=float(np.mean(b*b)),
            mother_mse_ratio=float(np.mean(b*b)/np.mean(a*a))))
    return dict(projection=projected,waveform=waveform),before,after


def crossfit(data, cache, folds, sample=None, keep_arrays=False):
    sample=np.arange(len(folds)) if sample is None else np.asarray(sample,int)
    chosen_folds=np.asarray(folds)[sample]
    old=data['means_mV'][sample]-data['old_predictions_mV'][sample]
    new=np.empty_like(old);models=[];profiles=[];knots=[]
    for fold in range(4):
        train=sample[chosen_folds!=fold];testpos=np.flatnonzero(chosen_folds==fold);test=sample[testpos]
        if set(train)&set(test):raise ValueError('Optical mother leaks between train and evaluation')
        correction=fit_update(data,train,cache)
        delta=many_voltage(correction,data['q0'][test],data['time_ns'])
        new[testpos]=old[testpos]-delta
        nominal=many_voltage(correction,[cache['nominal_high_q0']],data['time_ns'])[0]
        models.append(dict(fold=fold,training_rows=len(train),test_rows=len(test),
            unique_training_mothers=len(set(train)),unique_test_mothers=len(set(test)),
            scale=correction.scale,positive_scale_gate=bool(correction.scale>=.5),
            nominal_shift_100logE=(100*nominal@data['operator'].T).tolist()))
        if keep_arrays:profiles.append(correction.profiles);knots.append(correction.log_charge_knots)
    summary,before,after=residual_statistics(old,new,data['groups'][sample],data['operator'])
    summary['models']=models
    if keep_arrays:
        arrays=dict(old_residual_mV=old,corrected_residual_mV=new,old_projection_100logE=before,
            corrected_projection_100logE=after,profiles=np.array(profiles),log_charge_knots=np.array(knots))
        return summary,arrays
    return summary


def projected_variance(repeats, predicted, operator):
    """Finite-repeat correction; retain negative between-mother estimates."""
    residual=repeats-predicted[:,None,:]
    y=100*np.einsum('mrt,jt->mrj',residual,operator)
    count,reps,_=y.shape;mean=y.mean(1)
    within=y.var(1,ddof=1).mean(0);observed=mean.var(0,ddof=1);between=observed-within/reps
    return [dict(gap_ns=gap,mothers=count,repeats=reps,mean=float(mean[:,j].mean()),
        within_variance=float(within[j]),between_variance_unclipped=float(between[j]),
        observed_mean_variance=float(observed[j]),mean_standard_error=float(np.sqrt(observed[j]/count)),
        between_fraction_mean8=float(between[j]/observed[j]) if observed[j]>0 else None,
        planning_mothers_halfwidth10_repeat8=int(np.ceil(1.96**2*observed[j]/100)),
        planning_mothers_halfwidth10_infinite_repeats=int(np.ceil(1.96**2*max(between[j],0)/100)))
        for j,gap in enumerate((8.,40.,200.))]
