"""Causal shared-cell correction anchored to calibrated isolated voltages.

The age model supplies G(light1+light2)-G(light1)-G(light2), not an absolute
replacement of a calibrated voltage by a regularized deconvolution estimate.
This anchored difference is a project hypothesis, not Rosado's analytic model.
True count/times are explicit inputs in this experiment; energy remains unknown.
"""
from dataclasses import replace
from types import SimpleNamespace
import ctypes

import numpy as np
from scipy.optimize import least_squares, nnls

from .sipm_readout import ReadoutConfig, synthesize_readout
from .alpha_pair_recovery_forward import CompiledAgeRecovery

METHODS=('isolated_sum','causal_uniform','causal_calibrated','legacy_joint')


class DifferentiableAge(CompiledAgeRecovery):
    def __init__(self,sensor,library,dt_ns=1.):
        super().__init__(sensor,library,dt_ns)
        self.grad_function=self.library.gagg_age_gradient
        pointer=np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS')
        self.grad_function.argtypes=[pointer,pointer,ctypes.c_int,ctypes.c_int,ctypes.c_double,
                                    pointer,pointer,pointer,ctypes.c_int,ctypes.c_double,pointer,pointer]
        self.grad_function.restype=ctypes.c_int

    def with_gradient(self,light,derivatives):
        x=np.ascontiguousarray(light,float);d=np.ascontiguousarray(derivatives,float)
        if d.ndim!=2 or d.shape[1:]!=x.shape or not np.isfinite(d).all():raise ValueError('Invalid input derivatives')
        y=np.empty_like(x);dy=np.empty_like(d)
        status=self.grad_function(x,d.ravel(),len(x),len(d),self.cells,self.recovery,self.pde,self.trigger,
                                  len(self.recovery),self.xt,y,dy.ravel())
        if status:raise ValueError(f'Age sensitivity failed: {status}')
        return y,dy


def light_bin_probabilities(times, dt_ns, duration_ns, fast_fraction=.3):
    """Integrate a causal 95/351 ns two-decay density over absolute time bins.

    No negative-time density, no periodic wrap, no interpolation of noisy
    reconstructed timelines. Fraction is a frozen project optical-shape input.
    Optical transit and spatial illumination remain approximation limits.
    """
    t=np.asarray(times,float)
    if t.ndim!=1 or np.any(t<0) or not np.isfinite(t).all() or not 0<=fast_fraction<=1:
        raise ValueError('Finite nonnegative times and a fraction in [0,1] required')
    u=np.maximum(np.arange(int(round(duration_ns/dt_ns))+1)*dt_ns-t[:,None],0.)
    cumulative=fast_fraction*(-np.expm1(-u/95.))+(1-fast_fraction)*(-np.expm1(-u/351.))
    return np.maximum(np.diff(cumulative,axis=-1),0.)


class CausalReadout:
    """Actual aperture with coarse released charges placed at RIGHT bin edges.

    This introduces a bounded causal delay (<dt), never early information from
    a later part of the bin. Grid convergence is required. The original 0.25 ns
    electronic integration, 0.5 ns aperture and 2 ns ADC period are unchanged.
    """
    def __init__(self,configuration,dt_ns=1.,duration_ns=2048.):
        original=ReadoutConfig(**configuration['readout'])
        if dt_ns not in (.5,1.,2.):raise ValueError('Supported causal grids: .5,1,2 ns')
        self.dt=float(dt_ns);self.duration=float(duration_ns)
        self.n=int(round(duration_ns/dt_ns));self.stride=int(round(original.sample_period_ns/dt_ns))
        self.adc_centers=np.arange(int(duration_ns/original.sample_period_ns))*original.sample_period_ns+original.aperture_ns/2
        self.aperture=original.aperture_ns
        config=replace(original,sample_period_ns=dt_ns)
        full=configuration['sensor']['gain_electrons']*1.602176634e-19
        records=np.array([(0.,full)],dtype=[('time_ns','f8'),('charge_C','f8')])
        avalanches=SimpleNamespace(records=records,window_ns=duration_ns,full_charge_C=full)
        self.kernel=synthesize_readout(avalanches,config,seed=0).arrays['signal_mV']

    def __call__(self, released):
        x=np.asarray(released,float)
        if x.shape!=(self.n,) or not np.isfinite(x).all():raise ValueError('Released-charge grid mismatch')
        # Direct convolution retains exact zero support; no FFT roundoff before
        # the first excitation. The final bin fires outside this observation.
        shifted=np.r_[0.,x[:-1]]
        return np.convolve(shifted,self.kernel)[:self.n][::self.stride]


def make_light_charge_map(age,anchor,*,knots=81,duration_ns=6144.,fast_fraction=.3):
    """Map full-PDE-equivalent primary light to released charge within the model.

    Inference scales it to the frozen isolated mean Q(E_dep). No observed or
    true per-event charge is a fitting regressor; Birks is not applied again.
    """
    light=np.geomspace(250.,40000.,knots)
    shape=light_bin_probabilities([0.],age.dt,duration_ns,fast_fraction)[0]
    charge=np.array([age.forward(a*shape).sum() for a in light])
    if not np.all(np.diff(charge)>0):raise ValueError('Nonmonotonic isolated response map')
    e=np.linspace(*anchor.energy_bounds,81)
    target=np.array([anchor.components([0.],[x])[1][0] for x in e])
    if min(target)<charge[0] or max(target)>charge[-1]:raise ValueError('Light map does not cover inherited E bounds')
    return dict(light_pe=light,charge_pe=charge,dt_ns=age.dt,fast_fraction=fast_fraction,duration_ns=duration_ns)


class CausalRecovery:
    def __init__(self,anchor,configuration,library,charge_map,strength=1.):
        self.anchor=anchor;self.mapping=charge_map
        self.dt=float(charge_map['dt_ns']);self.fraction=float(charge_map['fast_fraction'])
        self.age=DifferentiableAge(configuration['sensor'],library,self.dt)
        self.readout=CausalReadout(configuration,self.dt)
        self.background=anchor.baseline+anchor.background
        self.energy_bounds=anchor.energy_bounds
        self.strength=float(strength)

    def isolated(self,times,energies):
        curves=[];charges=[]
        for t,e in zip(times,energies):
            v,q=self.anchor.components([t],[e]);v=v[0].copy()
            # Remove the inherited sub-grid interpolation leakage only for
            # ADC apertures lying WHOLLY before the event (not just centers).
            v[self.readout.adc_centers+self.readout.aperture/2<=t]=0.
            curves.append(v);charges.append(float(q[0]))
        return np.asarray(curves).reshape(-1,len(self.readout.adc_centers)),np.asarray(charges)

    def prepare(self,times):
        return light_bin_probabilities(times,self.dt,self.readout.duration,self.fraction)

    def components(self,times,energies,*,interaction=True,probabilities=None,strength=None):
        t=np.asarray(times,float);e=np.asarray(energies,float)
        if t.shape!=e.shape or len(t) not in (1,2):raise ValueError('Expected one/two aligned events')
        anchored,q=self.isolated(t,e)
        if not interaction or len(t)==1:return anchored
        p=self.prepare(t) if probabilities is None else probabilities
        a=np.interp(q,self.mapping['charge_pe'],self.mapping['light_pe'])
        light=a[:,None]*p;total=light.sum(0)
        joint=self.age.forward(total)
        isolated=np.array([self.age.forward(x) for x in light])
        shares=np.divide(light,total[None],out=np.zeros_like(light),where=total[None]>0)
        corrections=np.array([self.readout(joint*shares[j]-isolated[j]) for j in range(len(t))])
        scale=self.strength if strength is None else strength
        return anchored+scale*corrections

    def mean_and_gradient(self,times,log_energy,probabilities,*,strength):
        e=np.exp(log_energy);t=np.asarray(times,float)
        anchor,q=self.isolated(t,e);dv=[];dq=[]
        for j in range(len(e)):
            up=min(self.energy_bounds[1],e[j]*np.exp(1e-4));down=max(self.energy_bounds[0],e[j]*np.exp(-1e-4))
            vp,qp=self.isolated([t[j]],[up]);vm,qm=self.isolated([t[j]],[down]);den=np.log(up/down)
            dv.append((vp[0]-vm[0])/den);dq.append((qp[0]-qm[0])/den)
        dv=np.array(dv)
        if len(e)==1 or strength==0:return anchor.sum(0),dv
        knots=self.mapping['charge_pe'];light_knots=self.mapping['light_pe']
        indices=np.clip(np.searchsorted(knots,q,side='right')-1,0,len(knots)-2)
        slopes=np.diff(light_knots)[indices]/np.diff(knots)[indices]
        amplitudes=np.interp(q,knots,light_knots)
        lights=amplitudes[:,None]*probabilities;dlights=(slopes*np.array(dq))[:,None]*probabilities
        joint,djoint=self.age.with_gradient(lights.sum(0),dlights)
        singles=[];dsingles=[]
        for light,dlight in zip(lights,dlights):
            single,dsingle=self.age.with_gradient(light,dlight[None]);singles.append(single);dsingles.append(dsingle[0])
        delta=self.readout(joint-np.sum(singles,0))
        derivatives=np.array([self.readout(djoint[j]-dsingles[j]) for j in range(len(e))])
        return anchor.sum(0)+strength*delta,dv+strength*derivatives

    def initial_energy(self,voltage,times):
        reference=np.full(len(times),2.5)
        curves,q=self.isolated(times,reference)
        estimate=nnls((curves/q[:,None]).T,np.asarray(voltage)-self.background,maxiter=100)[0]
        return self.anchor.initial_energy(estimate)

    def fit(self,voltage,times,initial,metric,*,method):
        if method not in METHODS:raise ValueError(method)
        t=np.asarray(times,float);p=self.prepare(t)
        target=metric.apply(self.anchor.transform(np.asarray(voltage)-self.background))
        cached_x=None;cached_result=None
        def evaluate(log_energy):
            nonlocal cached_x,cached_result
            if cached_x is not None and np.array_equal(log_energy,cached_x):return cached_result
            if method=='legacy_joint':
                energy=np.exp(log_energy);curve=self.anchor.components(t,energy)[0].sum(0);grad=[]
                for j in range(len(t)):
                    up,down=energy.copy(),energy.copy()
                    up[j]=min(self.energy_bounds[1],energy[j]*np.exp(1e-4));down[j]=max(self.energy_bounds[0],energy[j]*np.exp(-1e-4))
                    grad.append((self.anchor.components(t,up)[0].sum(0)-self.anchor.components(t,down)[0].sum(0))/np.log(up[j]/down[j]))
            else:
                strength={'isolated_sum':0.,'causal_uniform':1.,'causal_calibrated':self.strength}[method]
                curve,grad=self.mean_and_gradient(t,log_energy,p,strength=strength)
            value=metric.apply(self.anchor.transform(curve))-target
            jac=np.column_stack([metric.apply(self.anchor.transform(g)) for g in grad])
            cached_x=log_energy.copy();cached_result=(value,jac);return cached_result
        lo,hi=np.log(self.energy_bounds);initial=np.asarray(initial,float);starts=[initial]
        if len(t)==2:
            total=float(self.anchor.q_isolated(initial).sum())
            starts.extend(self.anchor.initial_energy(total*np.array([share,1-share])) for share in (.03,.15,.5,.85,.97))
        fits=[least_squares(lambda x:evaluate(x)[0],np.clip(np.log(start),lo+1e-8,hi-1e-8),jac=lambda x:evaluate(x)[1],bounds=(lo,hi),max_nfev=100,
                            ftol=1e-9,xtol=1e-9,gtol=1e-5,diff_step=2e-5) for start in starts]
        best=min(fits,key=lambda f:float(f.fun@f.fun))
        return dict(energy_MeV=np.exp(best.x).tolist(),times_ns=t.tolist(),cost=float(best.fun@best.fun),
                    converged=bool(best.success),boundary=bool(np.any(best.x-lo<1e-5) or np.any(hi-best.x<1e-5)),
                    starts=len(starts),evaluations=sum(f.nfev for f in fits),metric_energy_MeV=initial.tolist(),method=method)
