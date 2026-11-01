"""Fixed-time energy/noise controls, explicitly conditional on known count/time.

No true energy or charge is accepted by the fitting API. Mean replacement and
Gaussian generation are oracle diagnostics, never a blind reconstruction.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import solve_banded

from .alpha_pair_energy_oracle import EnergyObjective
from .alpha_pair_model_sampler import lower_band_product, lift_innovations


def diagonal_factor(factor):
    """Preserve diag(L L.T), removing only off-diagonal covariance."""
    factor = np.asarray(factor, float)
    diagonal = np.zeros(factor.shape[1])
    for lag, row in enumerate(factor):
        diagonal[lag:] += row[:len(diagonal)-lag]**2
    return np.sqrt(diagonal)[None, :]


class FixedMetricObjective(EnergyObjective):
    def __init__(self, response, voltage, factor):
        self.r = response
        self.signal = np.asarray(voltage, float)-response.baseline-response.background
        self.factor = np.asarray(factor, float)
        if self.signal.shape != (1024,) or not np.isfinite(self.signal).all():
            raise ValueError('Expected a finite full 1024-sample waveform')
        if self.factor.shape[1] != len(response.transform(self.signal)):
            raise ValueError('Metric/transform dimension mismatch')
        self.target = self.whiten(response.transform(self.signal))


def fit_known_times(response, voltage, known_times, factor, weight_initial):
    """Same six-start fixed-time search as the preceding energy-oracle study."""
    times = np.asarray(known_times, float)
    if times.ndim != 1 or len(times) not in (1, 2) or np.any(np.diff(times) < 0):
        raise ValueError('Expected ordered one/two explicitly supplied times')
    objective = FixedMetricObjective(response, voltage, factor)
    starts = [objective.linear_energy_start(times)]
    if len(times) == 2:
        total = max(sum(weight_initial['charge_pe']), 100.)
        starts.extend(response.initial_energy(total*np.array([s, 1-s]))
                      for s in (.03, .15, .5, .85, .97))
    fits = [objective.fit(times, energy, free_times=False) for energy in starts]
    best = min(fits, key=lambda row: row['cost'])
    return dict(best, starts=len(starts), total_evaluations=sum(f['evaluations'] for f in fits))


def replace_conditional_mean(observations, reference, model_mean):
    """Independent reference mean; preserve evaluation residual differences.

    The finite reference-mean error is shared by all corrected observations.
    This requires true parameters for model_mean and is NOT a usable correction.
    """
    y, ref, model = map(lambda x: np.asarray(x, float), (observations, reference, model_mean))
    if y.ndim != 2 or ref.ndim != 2 or y.shape[1:] != model.shape or ref.shape[1:] != model.shape:
        raise ValueError('Incompatible observation/reference/model dimensions')
    if len(ref) < 2:
        raise ValueError('At least two independent reference draws required')
    return y-ref.mean(0)+model


def gaussian_waveforms(response, model_mean, factor, seed, count):
    """Exact covariance in the fitted FIR domain, with zero voltage prefix."""
    rng = np.random.default_rng(seed)
    z = lower_band_product(factor, rng.standard_normal((count, factor.shape[1])))
    noise = np.array([lift_innovations(row, response.ar) for row in z])
    np.testing.assert_allclose(response.transform(noise), z, atol=2e-10, rtol=2e-10)
    return np.asarray(model_mean)+noise


def local_energy_projection(response, times, energies, factor, reference, observations):
    """Truth-point linear diagnostics; neither empirical CIs nor a bound."""
    energy = np.asarray(energies, float)
    mean = response.components(times, energy)[0].sum(0)+response.baseline+response.background
    step = 1e-4
    columns = []
    for j in range(len(energy)):
        up, down = energy.copy(), energy.copy()
        up[j] *= np.exp(step); down[j] *= np.exp(-step)
        columns.append(response.transform((response.components(times, up)[0].sum(0)
                        -response.components(times, down)[0].sum(0))/(2*step)))
    def whiten(x):
        x = np.asarray(x)
        return solve_banded((len(factor)-1, 0), factor, x.T, check_finite=False).T
    jac = whiten(np.asarray(columns)).T
    u, singular, vh = np.linalg.svd(jac, full_matrices=False)
    rank = int(np.sum(singular > singular[0]*1e-10))
    inverse = np.linalg.pinv(jac, rcond=1e-10)
    covariance = inverse@inverse.T
    projected = whiten(response.transform(np.asarray(observations)-mean))@inverse.T
    reference_projection = whiten(response.transform(np.asarray(reference)-mean))@inverse.T
    ref_white = whiten(response.transform(np.asarray(reference)-mean))
    eval_white = whiten(response.transform(np.asarray(observations)-np.asarray(reference).mean(0)))
    return dict(rank=rank, condition_number=float(singular[0]/singular[-1]),
        gaussian_local_logE_sd=np.sqrt(np.diag(covariance)).tolist(),
        gaussian_local_correlation=(covariance/np.sqrt(np.outer(np.diag(covariance), np.diag(covariance)))).tolist(),
        reference_mean_projected_logE=reference_projection.mean(0).tolist(),
        reference_mean_projected_se=(reference_projection.std(0, ddof=1)/np.sqrt(len(reference))).tolist(),
        evaluation_projected_logE_sd=projected.std(0, ddof=1).tolist(),
        whitened_reference_variance_mean=float(ref_white.var(0, ddof=1).mean()),
        whitened_evaluation_variance_mean=float(eval_white.var(0, ddof=1).mean()),
        conditional_mean_white_rms=float(np.sqrt(np.mean(ref_white.mean(0)**2))),
        caveat='Known K/t/E local linearization. Fixed optical mother noise omits transport fluctuations. Finite independent reference mean is not the population mean.')
