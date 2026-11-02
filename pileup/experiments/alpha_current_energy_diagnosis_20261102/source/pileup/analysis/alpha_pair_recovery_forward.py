"""Refold frozen deconvolution-domain alpha means and fit native voltage.

This diagnostic changes fitting space/weights while retaining the age-recovery
physics and light profiles. An optional native-voltage residual mean is fitted
on isolated development events. True count and times remain explicit inputs.
"""
from __future__ import annotations

import ctypes
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.ndimage import gaussian_filter1d

from .alpha_pair_deconvolution import AgeRecovery

METHODS = ('native_ols_centered', 'native_gls_centered', 'native_gls_calibrated')


class CompiledAgeRecovery(AgeRecovery):
    """Exactly the same discrete expectation, accelerated in a small C++ loop."""
    def __init__(self, sensor, library, dt_ns=2.):
        super().__init__(sensor, dt_ns)
        self.library = ctypes.CDLL(str(Path(library).resolve()))
        self.function = self.library.gagg_age_forward
        pointer = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
        self.function.argtypes = [pointer, ctypes.c_int, ctypes.c_double, pointer, pointer,
                                  pointer, ctypes.c_int, ctypes.c_double, pointer]
        self.function.restype = ctypes.c_int
        for name in ('recovery', 'pde', 'trigger'):
            setattr(self, name, np.ascontiguousarray(getattr(self, name)))

    def forward(self, light):
        x = np.ascontiguousarray(light, dtype=np.float64)
        if x.ndim != 1 or not np.isfinite(x).all() or np.any(x < 0):
            raise ValueError('Expected finite nonnegative equivalent light')
        if not len(x): return x.copy()
        result = np.empty_like(x)
        status = self.function(x, len(x), self.cells, self.recovery, self.pde, self.trigger,
                               len(self.recovery), self.xt, result)
        if status: raise ValueError(f'Compiled age response failed ({status})')
        return result


class NativeForward:
    def __init__(self, profiles, convolution, recovery, background_mV, correction=None):
        self.profiles, self.op, self.recovery = profiles, convolution, recovery
        self.background = float(background_mV)
        self.empty_signal = self.op.forward(self.recovery.forward(profiles.background))
        self.correction = correction

    def correction_signal(self, times, energies):
        if self.correction is None: raise ValueError('No independently calibrated voltage correction')
        grid = np.arange(len(self.profiles.background))*2.+.25
        result = np.zeros(len(grid)); knots = self.correction['energy_knots']
        for t,e in zip(times,energies):
            j = int(np.clip(np.searchsorted(knots,e)-1,0,len(knots)-2))
            f = np.clip((e-knots[j])/(knots[j+1]-knots[j]),0.,1.)
            shape = (1-f)*self.correction['profiles'][j]+f*self.correction['profiles'][j+1]
            q = np.interp(e,self.profiles.energy_grid,self.profiles.charge_grid)
            result += q*np.interp(grid-t,self.correction['relative_ns'],shape,left=0.,right=0.)
        return result

    def signal(self, times, energies, *, centered=True, corrected=False):
        light = self.profiles.components(times, energies).sum(0)+self.profiles.background
        voltage = self.op.forward(self.recovery.forward(light))
        voltage = voltage-self.empty_signal if centered else voltage
        if corrected:
            voltage = np.maximum(voltage+self.correction_signal(times,energies),0.)
        return voltage

    def fit(self, voltage, times, initial_energy, *, method, metric=None, transform=None):
        if method not in METHODS: raise ValueError(method)
        times, initial = np.asarray(times, float), np.asarray(initial_energy, float)
        observed = np.asarray(voltage, float)-self.background
        if len(times) not in (1, 2) or initial.shape != times.shape or observed.shape != self.profiles.background.shape:
            raise ValueError('Expected complete voltage and aligned known K/t with voltage-derived initial energies')
        if method.startswith('native_gls'):
            if metric is None or transform is None: raise ValueError('A fixed covariance metric is required')
            whiten = lambda value: metric.apply(transform(value))
        else:
            whiten = lambda value: value
        target = whiten(observed)
        lo, hi = np.log(self.profiles.bounds)
        starts = [initial]
        if len(times) == 2:
            total = float(np.interp(initial, self.profiles.energy_grid, self.profiles.charge_grid).sum())
            for share in (.03, .15, .5, .85, .97):
                starts.append(np.interp(total*np.array([share, 1-share]), self.profiles.charge_grid, self.profiles.energy_grid))
        def residual(log_energy):
            return whiten(self.signal(times, np.exp(log_energy), corrected=method=='native_gls_calibrated'))-target
        fits = [least_squares(residual, np.clip(np.log(start),lo+1e-8,hi-1e-8),bounds=(lo,hi),max_nfev=100,
                             ftol=1e-8,xtol=1e-8,gtol=1e-5,diff_step=2e-5) for start in starts]
        best = min(fits,key=lambda f:float(f.fun@f.fun))
        return dict(times_ns=times.tolist(),energy_MeV=np.exp(best.x).tolist(),cost=float(best.fun@best.fun),
                    converged=bool(best.success),boundary=bool(np.any(best.x-lo<1e-5) or np.any(hi-best.x<1e-5)),
                    starts=len(starts),evaluations=sum(f.nfev for f in fits),method=method,
                    covariance_energy_MeV=initial.tolist() if method.startswith('native_gls') else None)


def calibrate_voltage_residual(response, records, voltages):
    """Four isolated E-dependent residual means, no pair/evaluation labels."""
    grid = np.arange(-64.,1720.01,2.); native = np.arange(len(response.profiles.background))*2.+.25
    knots, profiles = [], []
    for condition in ('weak','low','mid','high'):
        numerator, denominator = np.zeros(len(grid)),np.zeros(len(grid)); energies=[]
        for record,voltage in zip(records,voltages):
            if record['condition'] != condition: continue
            e,t = record['energy_MeV'][0],record['times_ns'][0]
            q = np.interp(e,response.profiles.energy_grid,response.profiles.charge_grid)
            predicted = response.signal([t],[e])
            residual = (np.asarray(voltage)-response.background-predicted)/q
            valid = (grid+t>=native[0]) & (grid+t<=native[-1])
            numerator += np.interp(grid+t,native,residual,left=0.,right=0.)*valid
            denominator += valid; energies.append(e)
        shape = gaussian_filter1d(numerator/np.maximum(denominator,1),2.)
        shape *= np.clip((grid[-1]-grid)/80.,0.,1.)
        profiles.append(shape);knots.append(np.mean(energies))
    return dict(relative_ns=grid,energy_knots=np.asarray(knots),profiles=np.asarray(profiles))
