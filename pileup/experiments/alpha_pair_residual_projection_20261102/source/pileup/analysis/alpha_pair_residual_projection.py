"""Truth-assisted residual attribution, not an energy reconstruction method.

At a declared energy, project voltage residuals on the log-energy tangent.
Contributions add in log energy; exponentiating each term is not additive.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import lsq_linear

TERMS = ('strong_isolated', 'weak_isolated', 'shared_correction', 'dark_electronics')
WEIGHTS = ('full', 'diagonal', 'unit_fir')


class DiagonalMetric:
    """An explicit fixed metric in the inherited FIR domain."""
    def __init__(self, variance):
        self.sigma = np.sqrt(np.asarray(variance, float))
        if self.sigma.ndim != 1 or not np.all(self.sigma > 0):
            raise ValueError('Expected positive finite marginal variances')
        if not np.isfinite(self.sigma).all():
            raise ValueError('Nonfinite variance')

    def apply(self, values):
        values = np.asarray(values, float)
        return values / self.sigma.reshape((-1,) + (1,) * (values.ndim - 1))


def fixed_metrics(metric):
    variance = metric.marginal_variance()
    return dict(full=metric, diagonal=DiagonalMetric(variance),
                unit_fir=DiagonalMetric(np.ones_like(variance)))


def residual_terms(voltage, shared, isolated, anchor, mean, background, weak):
    """All repeated waveforms must use the same mother, times and seed labels.

    Shared and isolated digitization consume randomness differently: these
    four terms are correlated bookkeeping, not independent causal sources.
    The final term includes background expectation mismatch as well as noise.
    """
    voltage, shared, isolated, anchor, mean = map(
        lambda x: np.asarray(x, float), (voltage, shared, isolated, anchor, mean))
    if shared.shape != isolated.shape or shared.shape[1] != 2:
        raise ValueError('Expected aligned repeated two-ancestor signals')
    if voltage.shape != (shared.shape[0], shared.shape[2]):
        raise ValueError('Voltage and ancestor dimensions differ')
    if anchor.shape != shared.shape[1:] or mean.shape != voltage.shape[1:]:
        raise ValueError('Model dimensions differ')
    result = np.stack((isolated[:, 1-weak]-anchor[1-weak],
                       isolated[:, weak]-anchor[weak],
                       shared.sum(1)-isolated.sum(1)-(mean-anchor.sum(0)),
                       voltage-background-shared.sum(1)), axis=1)
    error = float(np.max(abs(result.sum(1)-(voltage-background-mean))))
    if error > 1e-9:
        raise AssertionError(f'Residual closure failed: {error}')
    return result, error


def project(transform, metric, gradient, terms, energy, bounds, weak):
    """Return an exact linear least-squares influence and bounded comparison."""
    operator = metric.apply(transform)
    jacobian = operator @ np.asarray(gradient).T
    u, s, vt = np.linalg.svd(jacobian, full_matrices=False)
    if s[-1] <= 1e-12*s[0]:
        raise ValueError('Log-energy tangent is rank deficient at fixed cutoff')
    inverse = (vt.T/s) @ u.T
    influence = inverse @ operator
    contributions = np.einsum('et,rkt->rke', influence, terms)
    mean_terms = terms.mean(0)
    delta = influence @ mean_terms.sum(0)
    closure = float(np.max(abs(contributions.mean(0).sum(0)-delta)))
    strong_col, weak_col = jacobian[:, 1-weak], jacobian[:, weak]
    q = weak_col-strong_col*(strong_col@weak_col)/(strong_col@strong_col)
    profile_delta = float(q @ (operator@mean_terms.sum(0)) / (q@q))
    profile_error = abs(profile_delta-delta[weak])
    if max(closure, profile_error) > 1e-9:
        raise AssertionError('Projection/weak-profile additivity failed')
    lower, upper = np.log(bounds[0]/energy), np.log(bounds[1]/energy)
    bounded = lsq_linear(jacobian, operator@mean_terms.sum(0),
                         bounds=(lower, upper), tol=1e-12, max_iter=200)
    if not bounded.success:
        raise AssertionError('Bounded linear comparison did not converge')
    rho = float(strong_col@weak_col / np.linalg.norm(strong_col) / np.linalg.norm(weak_col))
    return dict(
        delta_log_energy=delta.tolist(),
        projected_bias_percent=(100*np.expm1(delta)).tolist(),
        bounded_delta_log_energy=bounded.x.tolist(),
        bounded_bias_percent=(100*np.expm1(bounded.x)).tolist(),
        outside_energy_bounds=bool(np.any(delta<lower) or np.any(delta>upper)),
        contributions_log=contributions.tolist(), column_cosine=rho,
        normalized_column_condition=float(np.sqrt((1+abs(rho))/(1-abs(rho)))),
        profiled_weak_norm=float(np.linalg.norm(q)), singular_values=s.tolist(),
        additivity_error=closure, profile_error=float(profile_error)), influence, jacobian


def paired_bootstrap(values, *, seed, draws=10000):
    """Resample mothers and repeats jointly across models/terms/energies.

    Input axes: mother, repeat, ... . Fixed gradients and weighting are not
    refit. Intervals describe this small, historically exposed collection.
    """
    values = np.asarray(values, float)
    if values.ndim < 2 or not np.isfinite(values).all():
        raise ValueError('Expected finite mother/repeat values')
    rng = np.random.default_rng(seed)
    m, r = values.shape[:2]
    mothers = rng.integers(m, size=(draws, m))
    repeats = rng.integers(r, size=(draws, m, r))
    samples = values[mothers[:, :, None], repeats].mean((1, 2))
    return samples, np.percentile(samples, [2.5, 97.5], axis=0)


def symmetric_change(influence_u, influence_s, residual_u, residual_s):
    """Exact order-symmetric split, descriptive rather than causal attribution."""
    response = .5*(influence_u+influence_s) @ (residual_s-residual_u)
    sensitivity = .5*(influence_s-influence_u) @ (residual_s+residual_u)
    total = influence_s@residual_s-influence_u@residual_u
    if np.max(abs(response+sensitivity-total)) > 1e-9:
        raise AssertionError('Symmetric change closure failed')
    return dict(residual_change_log=response.tolist(),
                sensitivity_change_log=sensitivity.tolist(), total_change_log=total.tolist())
