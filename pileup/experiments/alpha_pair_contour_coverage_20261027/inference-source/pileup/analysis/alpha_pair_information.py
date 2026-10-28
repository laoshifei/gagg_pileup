"""Local information diagnostics for the frozen conditional Gaussian model.

These calculations use true parameters and known K=2. They neither select
pulse count nor establish finite-sample or universal resolution limits.
Coordinates are (t1/10 ns, t2/10 ns, log Q1, log Q2).
"""
from __future__ import annotations

import numpy as np

from .alpha_pair_model_sampler import GaussianFactor


def dense_covariance(band, low_rank):
    """Expand scipy's symmetric lower-band storage and its downdate."""
    band = np.asarray(band, float)
    result = np.diag(band[0])
    n = band.shape[1]
    for k in range(1, len(band)):
        values = band[k, :n-k]
        result += np.diag(values, k) + np.diag(values, -k)
    return result - low_rank @ low_rank.T


def gaussian_information(mean_derivatives, covariance_derivatives, whitener):
    """F_ij = mu_i' C^-1 mu_j + tr(C^-1 C_i C^-1 C_j)/2.

whitener maps vectors to W x, W C W.T = I. Derivative axes come first.
The covariance term is required for this parameter-dependent covariance.
"""
    dm = np.asarray(mean_derivatives, float)
    dc = np.asarray(covariance_derivatives, float)
    W = np.asarray(whitener, float)
    if dm.ndim != 2 or dc.shape != (len(dm), dm.shape[1], dm.shape[1]):
        raise ValueError('Derivative shapes are inconsistent')
    if W.shape != (dm.shape[1], dm.shape[1]):
        raise ValueError('Whitening dimensions do not match')
    linear = dm @ W.T
    transformed = np.asarray([W @ matrix @ W.T for matrix in dc])
    mean_part = linear @ linear.T
    covariance_part = .5 * np.einsum('iab,jab->ij', transformed, transformed)
    return dict(total=mean_part + covariance_part, mean=mean_part,
                covariance=covariance_part, linear_scores=linear,
                quadratic_scores=transformed)


def invert_information(information):
    """Invert an identifiable local model without a pseudoinverse or ridge."""
    F = np.asarray(information, float)
    F = (F + F.T) / 2.
    values = np.linalg.eigvalsh(F)
    if not np.isfinite(values).all() or values[0] <= 0. or values[-1]/values[0] > 1e12:
        raise np.linalg.LinAlgError('Singular or poorly conditioned local information')
    covariance = np.linalg.solve(F, np.eye(len(F)))
    return (covariance + covariance.T) / 2.


def coordinates(times, charges):
    t = np.asarray(times, float)
    q = np.asarray(charges, float)
    if t.shape != (2,) or q.shape != (2,) or np.any(q <= 0.) or not np.isfinite(np.r_[t, q]).all():
        raise ValueError('Two finite times and positive charges are required')
    return np.r_[t / 10., np.log(q)]


def model_distribution(response, noise, z):
    z = np.asarray(z, float)
    curves = response.components(z[:2]*10., np.exp(z[2:]))
    mean = noise.transform(curves.sum(0) + response.baseline + noise.background_mV)
    band, V = noise.covariance(curves, np.exp(z[2:]))
    return mean, dense_covariance(band, V), GaussianFactor.from_covariance(band, V)


def local_information(response, noise, times, charges, *, step=1e-3, keep_scores=False):
    """Central derivatives: step in t/10 and log Q; no fit or truth-free claim."""
    if step <= 0.:
        raise ValueError('Positive derivative step required')
    z = coordinates(times, charges)
    mean, covariance, factor = model_distribution(response, noise, z)
    # whiten operates on final-axis vectors, so whiten(I) has rows of W.T.
    W = factor.whiten(np.eye(len(mean))).T
    dm = []
    dc = []
    for i in range(4):
        shift = np.zeros(4)
        shift[i] = step
        mp, cp, _ = model_distribution(response, noise, z + shift)
        mm, cm, _ = model_distribution(response, noise, z - shift)
        dm.append((mp-mm)/(2.*step))
        dc.append((cp-cm)/(2.*step))
    result = gaussian_information(dm, dc, W)
    F = result['total']
    inverse = invert_information(F)
    weak = int(np.argmin(charges))
    weak_indices = [weak, weak+2]
    weak_joint = inverse[np.ix_(weak_indices, weak_indices)]
    known_strong = invert_information(F[np.ix_(weak_indices, weak_indices)])
    fixed_times = invert_information(F[2:, 2:])
    sigma = np.sqrt(np.diag(inverse))
    result.update(inverse=inverse, correlation=inverse/np.outer(sigma, sigma),
                  weak_joint=weak_joint, known_strong=known_strong,
                  fixed_times_charge=fixed_times, weak_index=weak,
                  sigma_time_ns=sigma[:2]*10., sigma_log_charge=sigma[2:],
                  minimum_eigenvalue=float(np.linalg.eigvalsh(F)[0]),
                  condition_number=float(np.linalg.cond(F)), step=float(step),
                  whitening_relative_error=float(np.linalg.norm(W@covariance@W.T-np.eye(len(mean)), ord='fro')/np.sqrt(len(mean))))
    if not keep_scores:
        result.pop('linear_scores')
        result.pop('quadratic_scores')
    return result


def converged_information(response, noise, times, charges):
    """A predeclared derivative-halving diagnostic, keeping unresolved cases."""
    results = [local_information(response, noise, times, charges, step=s) for s in (1e-3, 5e-4)]
    comparisons = []
    for next_step in (2.5e-4, 1.25e-4, None):
        old, new = results[-2:]
        matrix_change = float(np.linalg.norm(new['total']-old['total'])/np.linalg.norm(new['total']))
        scale_change = float(np.max(abs(np.sqrt(np.diag(new['inverse'])/np.diag(old['inverse']))-1.)))
        comparisons.append(dict(step=new['step'], relative_information_change=matrix_change,
                                maximum_sigma_fractional_change=scale_change))
        if matrix_change <= .02 and scale_change <= .02 or next_step is None:
            break
        results.append(local_information(response, noise, times, charges, step=next_step))
    answer = results[-1]
    answer.update(derivative_checks=comparisons, derivative_converged=bool(
        comparisons[-1]['relative_information_change'] <= .02 and comparisons[-1]['maximum_sigma_fractional_change'] <= .02))
    return answer


def standardized_weak_error(fit, times, charges, information, *, fixed_times=False):
    """Truth-conditioned error normalization, not a blind confidence interval."""
    if len(fit['charge_pe']) != 2:
        raise ValueError('Known-K diagnostic requires two fitted components')
    order = np.argsort(fit['times_ns'])
    fitted = coordinates(np.asarray(fit['times_ns'])[order], np.asarray(fit['charge_pe'])[order])
    truth = coordinates(times, charges)
    error = fitted - truth
    weak = information['weak_index']
    if fixed_times:
        covariance = np.asarray(information['fixed_times_charge'])
        variance = float(covariance[weak, weak])
        return dict(log_charge_error=float(error[weak+2]), charge_fractional_error=float(np.expm1(error[weak+2])),
                    standardized_log_charge_error=float(error[weak+2]/np.sqrt(variance)),
                    time_error_ns=float(error[weak]*10.))
    small = error[[weak, weak+2]]
    covariance = np.asarray(information['weak_joint'])
    standardized = np.linalg.solve(np.linalg.cholesky(covariance), small)
    return dict(time_error_ns=float(small[0]*10.), log_charge_error=float(small[1]),
                charge_fractional_error=float(np.expm1(small[1])),
                standardized=standardized.tolist(), mahalanobis_squared=float(standardized@standardized))
