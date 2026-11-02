"""Area-preserving conditional isolated-mean correction, not event nuisance fits.

Calibration uses deposited-energy labels on training mothers. Inference uses
only candidate E/t, a frozen q0(E) map and the learned correction profiles.
The inherited shared-minus-isolated recovery term is deliberately unchanged.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.ndimage import gaussian_filter1d

from .alpha_pair_spatial_recovery import SpatialRecovery


def temporal_window(grid):
    s=np.asarray(grid,float)
    # A2ns guard makes the first interpolation interval causal at fractional
    # ADC phases without clipping after interpolation (which would alter area).
    onset=np.clip((s-2.)/8.,0.,1.)
    end=np.clip((1600.-s)/160.,0.,1.)
    return onset*onset*(3-2*onset)*end*end*(3-2*end)


@dataclass
class MeanShapeCorrection:
    grid: np.ndarray
    log_charge_knots: np.ndarray
    profiles: np.ndarray
    scale: float = 1.

    def __post_init__(self):
        self.grid=np.asarray(self.grid,float)
        self.log_charge_knots=np.asarray(self.log_charge_knots,float)
        self.profiles=np.asarray(self.profiles,float)
        if self.profiles.shape != (len(self.log_charge_knots),len(self.grid)):
            raise ValueError('Charge knots and correction grid disagree')
        if np.any(np.diff(self.log_charge_knots)<=0) or not np.allclose(np.diff(self.grid),2.):
            raise ValueError('Increasing charge knots and2ns grid required')
        if not np.isfinite(self.profiles).all() or not 0<=self.scale<=1:
            raise ValueError('Finite profiles and scale in[0,1] required')
        if np.max(abs(self.profiles.sum(1)))>1e-12:
            raise ValueError('Correction must have zero complete-grid area')
        if np.any(self.profiles[:,(self.grid<=.25)|(self.grid>=1600)]!=0):
            raise ValueError('Correction needs a causal interpolation guard and1600ns support')

    def profile(self, charge):
        q=float(charge)
        if not np.isfinite(q) or q<=0:raise ValueError('Positive frozen q0(E) required')
        x=float(np.clip(np.log(q),self.log_charge_knots[0],self.log_charge_knots[-1]))
        j=int(np.clip(np.searchsorted(self.log_charge_knots,x,side='right')-1,0,len(self.log_charge_knots)-2))
        fraction=(x-self.log_charge_knots[j])/(self.log_charge_knots[j+1]-self.log_charge_knots[j])
        return self.scale*q*((1-fraction)*self.profiles[j]+fraction*self.profiles[j+1])

    def voltage(self, grid, onset, charge):
        return np.interp(np.asarray(grid)-onset,self.grid,self.profile(charge),left=0.,right=0.)

    def save(self,path):
        np.savez_compressed(path,grid=self.grid,log_charge_knots=self.log_charge_knots,
                            profiles=self.profiles,scale=self.scale)

    @classmethod
    def load(cls,path):
        with np.load(path) as f:return cls(f['grid'],f['log_charge_knots'],f['profiles'],float(f['scale']))


def calibrate_profiles(sample_grid, means, predicted, q0, groups):
    """Equal mother weighting; fixed4ns smoothing and zero-area projection.

    Residuals are normalized by q0(true E), never by realized pulse charge.
    The positive old-template compensation shape removes the area component;
    this stage intentionally does not recalibrate the scalar energy map.
    """
    grid=np.arange(-16.,1602.,2.)+.25
    source_grid=np.asarray(sample_grid)-128.
    means,predicted=np.asarray(means),np.asarray(predicted)
    q0=np.asarray(q0,float);groups=np.asarray(groups)
    if means.shape!=predicted.shape or means.shape!=(len(q0),len(source_grid)):
        raise ValueError('Training arrays must share mother and sample axes')
    profiles,knots,records=[],[],[]
    for group in sorted(set(groups)):
        index=np.flatnonzero(groups==group)
        # Interpolation is used on average calibration curves, not as a noise
        # likelihood transform. All requested relative samples are observed.
        residual=np.mean((means[index]-predicted[index])/q0[index,None],axis=0)
        reference=np.mean(predicted[index]/q0[index,None],axis=0)
        residual=np.interp(grid,source_grid,residual,left=0.,right=0.)
        reference=np.maximum(np.interp(grid,source_grid,reference,left=0.,right=0.),0.)
        window=temporal_window(grid)
        correction=gaussian_filter1d(residual,2.,mode='constant')*window
        # Suppress unstable extrapolation where the old mean is nearly zero.
        correction*=reference/(reference+.01*reference.max())
        removed_area=float(2*correction.sum())
        weight=reference*window
        if weight.sum()<=0:raise ValueError('No positive calibration signal')
        correction-=correction.sum()*weight/weight.sum()
        correction-=correction.sum()*weight/weight.sum()
        profiles.append(correction);knots.append(float(np.log(q0[index]).mean()))
        records.append(dict(group=int(group),mothers=len(index),geometric_mean_q0=float(np.exp(knots[-1])),
            removed_area_per_q0=removed_area,
            final_area_per_q0=float(2*correction.sum())))
    order=np.argsort(knots)
    return MeanShapeCorrection(grid,np.array(knots)[order],np.array(profiles)[order]),records


class MeanShapeRecovery(SpatialRecovery):
    def __init__(self,anchor,configuration,library,charge_map,profile,correction):
        super().__init__(anchor,configuration,library,charge_map,profile)
        self.correction=correction

    def isolated(self,times,energies):
        original,charge=super().isolated(times,energies)
        q0=self.anchor.q_isolated(energies)
        changed=original+np.array([self.correction.voltage(self.readout.adc_centers,t,q)
                                  for t,q in zip(times,q0)])
        # Do not clip: clipping would silently destroy area preservation.
        if changed.min(initial=0.) < -1e-8:
            raise ValueError('Frozen mean correction produced negative voltage')
        return changed,charge


def positivity_scale(model,correction,energies,phases):
    """Training-model numerical constraint only; no evaluation outcomes used."""
    bound=1.
    for e in energies:
        q=float(model.anchor.q_isolated([e])[0])
        for phase in phases:
            old=model.isolated([128.+phase],[e])[0][0]
            delta=correction.voltage(model.readout.adc_centers,128.+phase,q)
            negative=delta<0
            if np.any(negative):bound=min(bound,float(np.min(old[negative]/(-delta[negative]))))
    return 1. if bound>=1. else max(0.,.99*bound)
