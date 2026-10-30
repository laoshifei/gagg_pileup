"""Energy-indexed causal response surrogate for local alpha pairs.

The exponential neighbor-history attenuation is a project hypothesis calibrated
on development charges, not Rosado's analytic response formula. Isolated means
and the scalar energy calibration stay frozen. No truth is read by inference.
"""
from __future__ import annotations

from dataclasses import dataclass
import json

import numpy as np
from scipy.linalg import cholesky_banded, cho_solve_banded
from scipy.optimize import least_squares, minimize
from scipy.signal import lfilter


def shape_features(grid, charge, other_charge=0., delay=0.):
    s = np.asarray(grid, float)
    x, other = charge/10000., other_charge/10000.
    past = other if delay < 0 else 0.
    future = other if delay >= 0 else 0.
    u = np.maximum(s-delay, 0.)
    p95, p351 = past*np.exp(-abs(delay)/95.), past*np.exp(-abs(delay)/351.)
    f95 = future*(-np.expm1(-u/8.))*np.exp(-u/95.)
    f351 = future*(-np.expm1(-u/8.))*np.exp(-u/351.)
    one = np.ones_like(s)
    return np.column_stack((one, x*one, x*x*one, p95*one, p351*one,
                            x*p351*one, f95, f351, x*f351))


@dataclass
class RecoveryEnergyResponse:
    arrays: dict
    energy_model: dict
    parameters: dict

    def __post_init__(self):
        for name, value in self.arrays.items():
            setattr(self, name, np.asarray(value))
        self.dt = 2.
        self.observed_samples = 1024
        self.full_grid = np.arange(3072)*self.dt+.25
        self.sample_grid = self.full_grid[:self.observed_samples]
        self.area_per_pe = float(self.gain)
        self.baseline = float(self.baseline_mV)
        self.background = float(self.background_mV)
        self.recovery_ns = float(self.parameters['recovery_ns'])
        self.cells = int(self.parameters['cells'])
        self.strength = float(self.parameters['strength'])
        self._energy_grid = np.asarray(self.energy_model['energy_MeV'])
        self._charge_grid = np.asarray(self.energy_model['charge_pe'])
        self.energy_bounds = (float(self.energy_model['energy_MeV'][0]),
                              float(self.energy_model['energy_MeV'][-1]))
        self._base_on_grid = np.interp(self.correction_grid, self.base_grid, self.base_kernel)
        onset = np.clip(self.correction_grid/8., 0., 1.)
        self._taper = onset*onset*(3.-2.*onset)*np.clip((self.correction_grid[-1]-self.correction_grid)/120., 0., 1.)
        self._innovation_per_pe = self.area_per_pe/self.dt*float(np.sum(self.ar))
        assert self._innovation_per_pe > 0 and self.recovery_ns > 0 and self.cells > 0

    @classmethod
    def load(cls, folder, strength=None):
        from pathlib import Path
        folder = Path(folder)
        with np.load(folder/'template.npz') as f:
            arrays = {k: f[k].copy() for k in f.files}
        model = json.loads((folder/'energy-model.json').read_text())
        parameter_file = folder/'response-parameters.json'
        if not parameter_file.exists():
            if strength is None:
                raise FileNotFoundError('Calibrate and freeze the response before inference')
            parameter_file = folder/'initial-response-parameters.json'
        params = json.loads(parameter_file.read_text())
        if strength is not None:
            params['strength'] = float(strength)
        return cls(arrays, model, params)

    def q_isolated(self, energy):
        energy = np.asarray(energy, float)
        if not np.isfinite(energy).all() or np.any(energy < self.energy_bounds[0]-1e-12) or np.any(energy > self.energy_bounds[1]+1e-12):
            raise ValueError('Proposed energy is outside the frozen calibration domain')
        return np.interp(energy, self._energy_grid, self._charge_grid)

    def unit(self, relative_time, q, neighbor=0., delay=0.):
        """Numerically identical to the frozen charge-conserving template."""
        correction = np.sum(shape_features(self.correction_grid, q, neighbor, delay)*self.mean_coefficients, axis=1)
        changed = np.maximum(0., self._base_on_grid+correction*self._taper)-self._base_on_grid
        area = float(self.base_area)+float(np.trapezoid(changed, self.correction_grid))
        if area <= 0:
            raise ValueError('Nonpositive isolated mean area')
        s = np.asarray(relative_time)
        base = np.interp(s, self.base_grid, self.base_kernel, left=0., right=0.)
        onset = np.clip(s/8., 0., 1.)
        taper = onset*onset*(3.-2.*onset)*np.clip((self.correction_grid[-1]-s)/120., 0., 1.)
        return np.maximum(0., base+np.interp(s, self.correction_grid, correction, left=0., right=0.)*taper)*self.area_per_pe/area

    def legacy_components(self, times, charges):
        t, q = np.asarray(times, float), np.asarray(charges, float)
        return np.asarray([q[j]*self.unit(self.sample_grid-t[j], q[j], q[1-j] if len(q)==2 else 0.,
                          t[1-j]-t[j] if len(q)==2 else 0.) for j in range(len(q))]).reshape(len(q), 1024)

    def drivers(self, times, energies):
        """Fixed-electronics innovations and a causal neighbor-state proxy."""
        times, energies = np.asarray(times, float), np.asarray(energies, float)
        if times.ndim != 1 or times.shape != energies.shape or len(times) > 2:
            raise ValueError('Expected aligned zero/one/two time and energy arrays')
        if not np.isfinite(times).all() or np.any(times < 0) or np.any(times > self.sample_grid[-1]):
            raise ValueError('Proposed time is outside the observed window')
        q = self.q_isolated(energies)
        isolated = np.asarray([qj*self.unit(self.full_grid-tj, qj) for tj, qj in zip(times, q)]).reshape(len(q), len(self.full_grid))
        if not len(q):
            return isolated, isolated, q
        z = lfilter(self.ar, [1.], isolated, axis=-1)
        # Positive innovations provide an equivalent-charge activity proxy.
        # This is not a reconstruction of actual microscopic cell occupancy.
        activity = np.maximum(z, 0.)/self._innovation_per_pe
        history = lfilter([1.], [1., -np.exp(-self.dt/self.recovery_ns)], activity, axis=-1)/self.cells
        neighbor = history[::-1] if len(q)==2 else np.zeros_like(history)
        return z, neighbor, q

    def from_drivers(self, z, neighbor, strength=None):
        k = self.strength if strength is None else float(strength)
        if not np.isfinite(k) or k < 0:
            raise ValueError('Recovery strength must be nonnegative')
        attenuation = np.exp(-k*neighbor)
        shared_z = z*attenuation
        curves = lfilter([1.], self.ar, shared_z, axis=-1) if len(z) else z
        charges = shared_z.sum(-1)/self._innovation_per_pe
        return curves, charges

    def components(self, times, energies, strength=None):
        z, neighbor, _ = self.drivers(times, energies)
        curves, charges = self.from_drivers(z, neighbor, strength)
        return curves[:, :self.observed_samples], charges

    def transform(self, x):
        x = np.asarray(x, float)
        order, n = len(self.ar)-1, x.shape[-1]
        return sum(c*x[..., order-k:n-k] for k, c in enumerate(self.ar))

    def fixed_metric(self, initial):
        """Original banded noise model frozen at blind initial Q/t, kappa=0.

        Charge-total conditioning is inappropriate for an energy-conditioned
        mean. This is a fixed weighting metric, not a calibrated new likelihood.
        """
        q = np.asarray(initial['charge_pe'], float)
        curves = self.legacy_components(initial['times_ns'], q)
        z = self.transform(curves)
        n = z.shape[-1]
        band = np.zeros((len(self.noise_lags), n))
        for lag, value in enumerate(self.noise_lags):
            band[lag, :n-lag] = value
        if len(q):
            activity = np.maximum(z, 0.)
            neighbor = activity[::-1] if len(q)==2 else np.zeros_like(activity)
            a, b, c = self.noise_parameters
            v = activity*np.exp(np.clip(a+b*np.log(np.maximum(q[:, None], 50.)/3000.)+c*neighbor/10., -8., 4.))
            root = np.sqrt(v)
            band[0] += v.sum(0)
            band[1, :-1] += float(self.lag1)*(root[:, 1:]*root[:, :-1]).sum(0)
            if len(q)==2:
                band[0] += 2*float(self.cross)*root[0]*root[1]
                band[1, :-1] += float(self.cross)*float(self.lag1)*(root[0, 1:]*root[1, :-1]+root[1, 1:]*root[0, :-1])
        return cholesky_banded(band, lower=True, check_finite=False)

    def initial_energy(self, charges):
        # Clamping is only an optimizer initialization; boundary hits are reported.
        q = np.asarray(charges, float)
        lo, hi = self.energy_bounds
        e = np.interp(q, self._charge_grid, self._energy_grid)
        return np.clip(e, lo*(1+1e-6), hi*(1-1e-6))


def coupled_charge_inverse(response, times, charges):
    """Charge-only comparator. No change to waveform fit or event count."""
    t, q = np.asarray(times, float), np.asarray(charges, float)
    if not len(t):
        return dict(times_ns=[], charge_pe=[], energy_MeV=[], converged=True, evaluations=0, boundary=False)
    lo, hi = np.log(response.energy_bounds)
    initial = np.log(response.initial_energy(q))
    def residual(x):
        _, predicted = response.components(t, np.exp(x))
        return (predicted-q)/np.maximum(q, 50.)
    fit = least_squares(residual, initial, bounds=(lo, hi), max_nfev=60,
                        ftol=1e-9, xtol=1e-9, gtol=1e-7, diff_step=1e-4)
    e = np.exp(fit.x)
    return dict(times_ns=t.tolist(), charge_pe=q.tolist(), energy_MeV=e.tolist(),
        converged=bool(fit.success), evaluations=int(fit.nfev), boundary=bool(np.any((fit.x-lo)<1e-4) or np.any((hi-fit.x)<1e-4)),
        inversion_relative_residual=float(np.max(np.abs(residual(fit.x)))))


def refine_energy(response, voltage, initial, strength, *, radius=24., maxiter=60):
    """Refine E/t at the event count already selected by a blind baseline."""
    y = np.asarray(voltage, float)
    if y.shape != (1024,) or not np.isfinite(y).all():
        raise ValueError('A finite 1024-sample analog voltage is required')
    times = np.asarray(initial['times_ns'], float)
    q = np.asarray(initial['charge_pe'], float)
    n = len(times)
    if not n:
        return dict(times_ns=[], charge_pe=[], energy_MeV=[], converged=True, evaluations=0, boundary=False)
    order = np.argsort(times); times, q = times[order], q[order]
    factor = response.fixed_metric(dict(times_ns=times, charge_pe=q))
    zy = response.transform(y-response.baseline-response.background)
    loe, hie = np.log(response.energy_bounds)
    low = np.r_[np.maximum(times-radius, 0.)/10., np.full(n, loe)]
    high = np.r_[np.minimum(times+radius, response.sample_grid[-1])/10., np.full(n, hie)]
    e0 = (coupled_charge_inverse(response, times, q)['energy_MeV'] if strength else response.initial_energy(q))
    x0 = np.clip(np.r_[times/10., np.log(e0)], low+1e-8, high-1e-8)
    def objective(x):
        curves, _ = response.components(x[:n]*10., np.exp(x[n:]), strength)
        residual = zy-response.transform(curves.sum(0))
        return float(residual@cho_solve_banded((factor, True), residual, check_finite=False))
    initial_cost = objective(x0)
    fit = minimize(objective, x0, method='L-BFGS-B', bounds=list(zip(low, high)),
        options=dict(maxiter=maxiter, ftol=1e-9, gtol=1e-4, maxls=25, eps=1e-5))
    x = fit.x if np.isfinite(fit.fun) and fit.fun <= initial_cost else x0
    t, e = x[:n]*10., np.exp(x[n:])
    _, realized_q = response.components(t, e, strength)
    order = np.argsort(t)
    return dict(times_ns=t[order].tolist(), charge_pe=realized_q[order].tolist(), energy_MeV=e[order].tolist(),
        isolated_charge_pe=response.q_isolated(e[order]).tolist(), converged=bool(fit.success),
        evaluations=int(fit.nfev), cost=objective(x), initial_cost=initial_cost,
        boundary=bool(np.any(x-low<1e-4) or np.any(high-x<1e-4)),
        message=str(fit.message))
