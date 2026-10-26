"""Regularized empirical pulse covariance with an observation-only pilot.

Single-event calibration labels align development residuals. Inference uses
only voltage and fitted pilot Q/t. The resulting positive covariance remains
fixed across every K/Q/t candidate in one frame. This is a project covariance
approximation, not a microscopic SiPM likelihood or a copied TES algorithm.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class EmpiricalNoise:
    grid: np.ndarray
    edges: np.ndarray
    variance: np.ndarray
    modes: np.ndarray
    mean_residual: np.ndarray
    counts: np.ndarray
    noise_floor: float = .01

    def save(self, path):
        np.savez_compressed(path, grid=self.grid, edges=self.edges, variance=self.variance,
                            modes=self.modes, mean_residual=self.mean_residual, counts=self.counts,
                            noise_floor=self.noise_floor)

    @classmethod
    def load(cls, path):
        with np.load(path) as f:
            return cls(f['grid'], f['edges'], f['variance'], f['modes'],
                       f['mean_residual'], f['counts'], float(f['noise_floor']))


def calibrate(grid, residuals, true_charge, *, rank=32, edges=(1500.,5000.,12500.), noise_floor=.01):
    """Train covariance on isolated optical parents, distinct from gate data."""
    residuals=np.asarray(residuals,float);charge=np.asarray(true_charge,float)
    if residuals.shape!=(len(charge),len(grid)) or not np.isfinite(residuals).all():
        raise ValueError('Expected finite aligned single-event residuals')
    bins=np.searchsorted(edges,charge,side='right')
    variance=[];modes=[];means=[];counts=[]
    for index in range(len(edges)+1):
        x=residuals[bins==index]
        if len(x)<=rank+1:raise ValueError('Too few independent calibration parents for covariance rank')
        mean=x.mean(0);centered=(x-mean)/np.sqrt(len(x)-1.)
        _,s,vh=np.linalg.svd(centered,full_matrices=False)
        u=(vh[:rank].T*s[:rank])
        variance.append(np.maximum(np.sum(centered**2,axis=0),noise_floor))
        modes.append(u);means.append(mean);counts.append(len(x))
    return EmpiricalNoise(np.asarray(grid,float),np.asarray(edges,float),np.asarray(variance),
                          np.asarray(modes),np.asarray(means),np.asarray(counts),float(noise_floor))


class LowRankPrecision:
    """Exact whitening for C=diag(d)+UU^T, up to a fixed scalar factor."""
    def __init__(self, diagonal, modes, *, pilot=None):
        self.diagonal=np.asarray(diagonal,float)
        self.modes=np.asarray(modes,float)
        if self.diagonal.ndim!=1 or self.modes.ndim!=2 or self.modes.shape[0]!=len(self.diagonal):
            raise ValueError('Incompatible covariance dimensions')
        if not np.isfinite(self.diagonal).all() or not np.isfinite(self.modes).all() or np.any(self.diagonal<=0.):
            raise ValueError('The regularized covariance must be finite and positive')
        self.root_inverse=1./np.sqrt(self.diagonal)
        a=self.root_inverse[:,None]*self.modes
        if a.shape[1]:
            self.left,s,_=np.linalg.svd(a,full_matrices=False)
            self.shrink=1.-1./np.sqrt(1.+s*s)
        else:
            self.left=np.empty((len(self.diagonal),0));self.shrink=np.empty(0)
        self.scale=1.
        if pilot is not None:
            w=self.apply(np.asarray(pilot,float))
            self.scale=np.sqrt(1e6/max(float(w@w),1.))

    def apply(self, values):
        x=np.asarray(values,float)
        w=self.root_inverse*x if x.ndim==1 else self.root_inverse[:,None]*x
        projected=self.left.T@w
        correction=self.shrink*projected if x.ndim==1 else self.shrink[:,None]*projected
        return self.scale*(w-self.left@correction)


def pilot(response, voltage):
    """Estimate a single smooth pulse from voltage, with no labels."""
    raw=response.base.candidates(voltage)
    if len(raw)>1:
        fit=response.refine(voltage,raw[1],radius=24.)
        return raw,fit
    return raw,dict(times_ns=[],charge_pe=[])


def metric_from_pilot(calibration, response, fit, voltage, *, rank=16, fraction=.9):
    if not 0.<=fraction<=1. or not 0<=rank<=calibration.modes.shape[-1]:
        raise ValueError('Invalid covariance rank or fraction')
    if not fit['times_ns']:
        return LowRankPrecision(np.full(1024,calibration.noise_floor),np.empty((1024,0)),
                                pilot=np.asarray(voltage)-response.baseline)
    time=float(fit['times_ns'][0]);charge=float(sum(fit['charge_pe']))
    index=int(np.searchsorted(calibration.edges,charge,side='right'))
    relative=response.sample_grid-time
    grid=calibration.grid
    # Calibration reads the original 20us records, so the complete 2.048us
    # inference record lies within this relative-time grid for normal pilots.
    # Outside the calibrated interval, keep a positive decaying diagonal tail;
    # never extrapolate oscillatory covariance modes or unobserved zeros.
    variance=np.interp(relative,grid,calibration.variance[index],
                       left=calibration.noise_floor,right=calibration.variance[index,-1])
    late=relative>grid[-1]
    variance[late]=calibration.noise_floor+(variance[late]-calibration.noise_floor)*np.exp(
        -(relative[late]-grid[-1])/351.)
    modes=np.column_stack([np.interp(relative,grid,calibration.modes[index,:,j],left=0.,right=0.)
                           for j in range(rank)]) if rank else np.empty((len(relative),0))
    modes*=np.sqrt(fraction)
    diagonal=np.maximum(variance-np.sum(modes*modes,axis=1),calibration.noise_floor)
    smooth=response.components(fit['times_ns'],fit['charge_pe']).sum(0)
    result=LowRankPrecision(diagonal,modes,pilot=smooth)
    result.calibration_bin=index;result.pilot_time_ns=time;result.pilot_charge_pe=charge
    return result


def waveform_candidates(response, calibration, voltage, *, rank=16, fraction=.9):
    from .alpha_pair_residual import fit_candidates
    raw,one=pilot(response,voltage)
    raw=[{k:v for k,v in r.items() if k!='model_mV'} for r in raw]
    start=raw[min(1,len(raw)-1)]
    metric=metric_from_pilot(calibration,response,one,voltage,rank=rank,fraction=fraction)
    seeds=dict(candidates=raw,ml_1=start,ml_2=start)
    return fit_candidates(response,voltage,seeds,fraction=0.,metric=metric)['waveform_only']
