"""Two-stage voltage unfolding and mean-field SiPM recovery inversion.

The multiplicative updates follow the structure of Bretz et al. (2020),
arXiv:2010.14886, Eqs. 5.8 and 5.12. This is a project adaptation: the readout
kernel is ours, the recovering-cell model is an age-distribution approximation,
and the final energy decomposition uses independently calibrated alpha means.
It is neither the authors' circuit Monte Carlo nor a blind event finder.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares, nnls
from scipy.signal import fftconvolve

from .sipm_readout import ReadoutConfig, synthesize_readout


def spe_kernel(configuration, samples=1024):
    """Average over eight internal-grid phases of one 2 ns latent bin.

    Latent bins represent released avalanche charge in full-cell p.e.; they
    are NOT photon counts. Their centers are 1+2n ns, whereas voltage sample
    centers are 0.25+2n ns. Use the physical aperture at the period's start.
    """
    c = ReadoutConfig(**configuration['readout'])
    phases = np.arange(c.internal_dt_ns/2, c.sample_period_ns, c.internal_dt_ns)
    full_charge = configuration['sensor']['gain_electrons']*1.602176634e-19
    records = np.zeros(len(phases), dtype=[('time_ns', 'f8'), ('charge_C', 'f8')])
    records['time_ns'] = phases
    records['charge_C'] = full_charge/len(phases)
    avalanches = SimpleNamespace(records=records, window_ns=samples*c.sample_period_ns,
                                 full_charge_C=full_charge)
    return np.maximum(synthesize_readout(avalanches, c, seed=0).arrays['signal_mV'], 0.)


class FiniteConvolution:
    """Causal forward/adjoint with the SAME finite observation window."""
    def __init__(self, kernel, samples):
        self.h = np.asarray(kernel, float)[:samples]
        self.n = samples
        if self.h.ndim != 1 or np.any(self.h < 0) or not np.isfinite(self.h).all() or self.h[0] <= 0:
            raise ValueError('Expected a finite nonnegative kernel with observable diagonal')
        self.exposure = self.adjoint(np.ones(samples))

    def forward(self, values):
        return fftconvolve(values, self.h)[:self.n]

    def adjoint(self, values):
        return fftconvolve(np.asarray(values)[::-1], self.h)[:self.n][::-1]

    def unfold(self, voltage, background, checkpoints=(128, 512)):
        """Positive Bayesian/RL update with fixed additive voltage background.

        Voltage is Gaussian-noisy, not a Poisson histogram: this is a working
        I-divergence criterion, not a claim of an exact physical likelihood.
        Fixed early stopping supplies regularization. No truth-dependent cut.
        """
        y = np.asarray(voltage, float)
        b = np.broadcast_to(np.asarray(background, float), y.shape)
        if y.shape != (self.n,) or not np.isfinite(y).all() or np.any(y < 0) or np.any(b <= 0):
            raise ValueError('RL needs nonnegative observed voltage and positive background')
        if not checkpoints or min(checkpoints) < 1:
            raise ValueError('Positive iteration counts required')
        x = np.ones(self.n)  # fixed, independent of event labels/amplitude
        result = {}
        for iteration in range(1, max(checkpoints)+1):
            mu = np.maximum(self.forward(x)+b, 1e-15)
            x *= np.maximum(self.adjoint(y/mu), 0.)/self.exposure
            if iteration in checkpoints:
                result[iteration] = x.copy()
        return result


class AgeRecovery:
    """Uniform-illumination mean response with recovering cell ages.

    At most one discharge per cell per time bin; prompt cross talk uses a
    Poisson mean-field closure rather than the digitizer's neighbour graph.
    The tail bucket is fully recovered (12 tau). Dark history and afterpulses
    are omitted; nonzero afterpulse configuration is rejected. Input is
    full-recovery-PDE-equivalent primary light, not incident photon number.
    """
    def __init__(self, sensor, dt_ns=2.):
        if sensor['afterpulse_probability'] != 0 or sensor['crosstalk_delay_ns'] != 0:
            raise ValueError('This approximation supports prompt cross talk and no afterpulses')
        self.cells = sensor['cells_x']*sensor['cells_y']
        self.dt = float(dt_ns)
        tau = float(sensor['recovery_ns'])
        if self.dt <= 0 or tau <= 0 or self.cells <= 0 or not sensor['recovery_enabled']:
            raise ValueError('Positive grid, cells and enabled recovery required')
        # A cell reset in the preceding bin has recovered for one bin by the
        # next opportunity; same-bin re-firing is excluded in this approximation.
        age = self.dt*np.arange(1, int(np.ceil(12*tau/self.dt))+1)
        self.recovery = -np.expm1(-age/tau)
        self.recovery[-1] = 1.
        pde = np.asarray(sensor['pde_probabilities'], float)
        relative = pde/pde[-1]
        if not np.allclose(relative, relative[:, :1], rtol=1e-10, atol=1e-12):
            raise ValueError('Nonseparable wavelength/voltage PDE needs a spectral state')
        voltage = float(sensor['overvoltage_V'])*self.recovery
        self.pde = np.interp(voltage, sensor['pde_overvoltage_V'], relative[:, 0])
        v = np.maximum(voltage-sensor['trigger_threshold_V'], 0.)
        vfull = sensor['overvoltage_V']-sensor['trigger_threshold_V']
        self.trigger = -np.expm1(-v/sensor['trigger_scale_V'])/-np.expm1(-vfull/sensor['trigger_scale_V'])
        self.xt = float(sensor['crosstalk_probability'])

    def forward(self, light):
        x = np.asarray(light, float)
        if x.ndim != 1 or not np.isfinite(x).all() or np.any(x < 0):
            raise ValueError('Expected nonnegative equivalent primary light per bin')
        state = np.zeros(len(self.recovery)); state[-1] = 1.
        output = np.empty_like(x)
        for i, count in enumerate(x):
            hazard = count/self.cells*self.pde
            cross = 0.
            for _ in range(6):
                fire = -np.expm1(-hazard-cross*self.trigger)
                charge = float(state @ (fire*self.recovery))
                cross = self.xt*charge
            fired = state*fire
            output[i] = self.cells*charge
            survived = state-fired
            state[1:] = survived[:-1]
            state[-1] += survived[-1]
            state[0] = fired.sum()
        return output

    def invert(self, released, max_iterations=40, tolerance=.005):
        """Bretz-style whole-timeline multiplicative recovery correction.

        We use relative L1 forward error, avoiding cancellation in integrated
        signed differences. This stopping rule is our choice, not Eq. 5.14.
        Nonconvergence is retained and reported, never silently discarded.
        """
        y = np.asarray(released, float)
        if y.ndim != 1 or not np.isfinite(y).all() or np.any(y < 0):
            raise ValueError('Expected nonnegative released charge')
        x = y.copy()
        scale = max(float(y.sum()), 1e-15)
        for step in range(max_iterations+1):
            predicted = self.forward(x)
            error = float(np.abs(predicted-y).sum()/scale)
            if error <= tolerance or step == max_iterations:
                return x, dict(iterations=step, relative_l1=error,
                               converged=bool(error <= tolerance),
                               correction=float(x.sum()/scale))
            x *= np.divide(y, predicted, out=np.zeros_like(y), where=predicted > 1e-15)
        raise AssertionError('Unreachable')


@dataclass
class EnergyProfiles:
    """Energy-dependent means learned ONLY from isolated development waves.

    Profiles describe a nonlinear reconstruction's output, including its
    regularization blur. They are not true pre-pileup photon-time profiles.
    Constant normalized-profile extrapolation beyond four energy knots is an
    approximation; the inherited q0(E_dep) calibration supplies the scale.
    """
    relative_ns: np.ndarray
    energy_knots: np.ndarray
    profiles: np.ndarray
    energy_grid: np.ndarray
    charge_grid: np.ndarray
    background: np.ndarray
    exposure: np.ndarray
    dt: float = 2.

    @property
    def bounds(self):
        return (float(self.energy_grid[0]), float(self.energy_grid[-1]))

    def components(self, times, energies):
        curves = []
        grid = (np.arange(len(self.background))+.5)*self.dt
        for t, e in zip(times, energies):
            j = int(np.clip(np.searchsorted(self.energy_knots, e)-1, 0, len(self.energy_knots)-2))
            frac = np.clip((e-self.energy_knots[j])/(self.energy_knots[j+1]-self.energy_knots[j]), 0., 1.)
            shape = (1-frac)*self.profiles[j]+frac*self.profiles[j+1]
            q = np.interp(e, self.energy_grid, self.charge_grid)
            curves.append(q*np.interp(grid-t, self.relative_ns, shape, left=0., right=0.))
        return np.asarray(curves)

    def fit(self, transformed, known_times):
        x, t = np.asarray(transformed, float), np.asarray(known_times, float)
        if x.shape != self.background.shape or len(t) not in (1, 2) or not np.isfinite(x).all():
            raise ValueError('Expected one waveform and explicitly known one/two event times')
        # Working diagonal variance. Unfolding induces correlations; this is
        # NOT GLS with a calibrated covariance or a physical photon likelihood.
        weights = np.sqrt(self.exposure/np.maximum(gaussian_filter1d(x, 5.), 1.))
        target = (x-self.background)*weights
        ref = float(self.energy_knots[1]); qref = np.interp(ref, self.energy_grid, self.charge_grid)
        columns = self.components(t, np.full(len(t), ref)).T/qref
        charges = nnls(columns*weights[:, None], target, maxiter=100)[0]
        initial = np.interp(charges, self.charge_grid, self.energy_grid)
        starts = [initial]
        if len(t) == 2:
            total = max(float(charges.sum()), float(self.charge_grid[0])*2)
            for share in (.03, .15, .5, .85, .97):
                starts.append(np.interp(total*np.array([share, 1-share]), self.charge_grid, self.energy_grid))
        lower, upper = np.log(self.bounds)
        def residual(log_energy):
            return self.components(t, np.exp(log_energy)).sum(0)*weights-target
        fits = [least_squares(residual, np.clip(np.log(start), lower+1e-8, upper-1e-8),
                bounds=(lower, upper), max_nfev=120, ftol=1e-9, xtol=1e-9, gtol=1e-6)
                for start in starts]
        fit = min(fits, key=lambda f: float(f.fun@f.fun))
        return dict(times_ns=t.tolist(), energy_MeV=np.exp(fit.x).tolist(), cost=float(fit.fun@fit.fun),
                    converged=bool(fit.success), starts=len(starts), evaluations=sum(f.nfev for f in fits),
                    boundary=bool(np.any(fit.x-lower < 1e-5) or np.any(upper-fit.x < 1e-5)))


def calibrate_profiles(records, transformed, empty, energy_model, exposure):
    """Align only mean profiles; no interpolation of a claimed noise covariance."""
    background = np.mean(empty, axis=0)
    grid = np.arange(-64., 1720.01, 2.)
    native = (np.arange(len(background))+.5)*2.
    knots, profiles = [], []
    for condition in ('weak', 'low', 'mid', 'high'):
        indices = [i for i, r in enumerate(records) if r['condition'] == condition]
        numerator, denominator = np.zeros(len(grid)), np.zeros(len(grid))
        energies = []
        for i in indices:
            r = records[i]; e = r['energy_MeV'][0]; t = r['times_ns'][0]
            q = np.interp(e, energy_model['energy_MeV'], energy_model['charge_pe'])
            valid = (grid+t >= native[0]) & (grid+t <= native[-1])
            numerator += np.interp(grid+t, native, (transformed[i]-background)/q, left=0., right=0.)*valid
            denominator += valid
            energies.append(e)
        shape = np.maximum(gaussian_filter1d(numerator/np.maximum(denominator, 1.), 2.), 0.)
        # Finite learned support, including pre-onset unfolding blur; only the
        # last 80 ns is tapered. This is an explicit tail approximation.
        shape *= np.clip((grid[-1]-grid)/80., 0., 1.)
        knots.append(float(np.mean(energies))); profiles.append(shape)
    return EnergyProfiles(grid, np.array(knots), np.array(profiles),
                          np.asarray(energy_model['energy_MeV']), np.asarray(energy_model['charge_pe']),
                          background, np.asarray(exposure)/max(exposure))
