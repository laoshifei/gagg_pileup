"""Candidate-conditioned covariance for known-count/time alpha energy fits.

This reuses the frozen empirical noise coefficients. Evaluating them on the
recovery-energy mean is a project hypothesis, not a newly calibrated SiPM
likelihood. The fitting API accepts no reference deposited energies or charges.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import cholesky_banded, solve_banded
from scipy.optimize import minimize

from .alpha_pair_energy_noise import FixedMetricObjective, fit_known_times

METHODS = ("reference_gls", "all_legacy_once", "energy_once", "energy_fgls", "energy_nll")


def covariance_band(response, curves, charges):
    """Frozen kappa=0 empirical covariance, now on explicit component curves.

    No conditioning on the observed total charge and no tuning to the current
    repeated digitizations. Units are squared mV after the fixed FIR transform.
    """
    curves, q = np.asarray(curves, float), np.asarray(charges, float)
    if curves.shape != (len(q), response.observed_samples) or len(q) not in (0, 1, 2):
        raise ValueError("Expected zero/one/two aligned component curves and charges")
    if not np.isfinite(curves).all() or not np.isfinite(q).all() or np.any(q < 0):
        raise ValueError("Nonfinite or negative candidate response")
    z = response.transform(curves)
    n = z.shape[-1]
    band = np.zeros((len(response.noise_lags), n))
    for lag, value in enumerate(response.noise_lags):
        band[lag, :n-lag] = value
    if len(q):
        activity = np.maximum(z, 0.)
        neighbor = activity[::-1] if len(q) == 2 else np.zeros_like(activity)
        a, b, c = response.noise_parameters
        variance = activity * np.exp(np.clip(
            a + b*np.log(np.maximum(q[:, None], 50.)/3000.) + c*neighbor/10., -8., 4.))
        root = np.sqrt(variance)
        band[0] += variance.sum(0)
        band[1, :-1] += float(response.lag1)*(root[:, 1:]*root[:, :-1]).sum(0)
        if len(q) == 2:
            band[0] += 2*float(response.cross)*root[0]*root[1]
            band[1, :-1] += float(response.cross)*float(response.lag1)*(
                root[0, 1:]*root[1, :-1] + root[1, 1:]*root[0, :-1])
    return band


def candidate_factor(response, times, energies, *, legacy=False):
    curves, charges = response.components(times, energies)
    if legacy:
        return response.fixed_metric(dict(times_ns=times, charge_pe=charges))
    return cholesky_banded(covariance_band(response, curves, charges), lower=True, check_finite=True)


def gaussian_terms(residual, factor):
    """Twice Gaussian NLL up to the common n*log(2*pi) constant."""
    r, factor = np.asarray(residual, float), np.asarray(factor, float)
    w = solve_banded((len(factor)-1, 0), factor, r.T, check_finite=True).T
    quadratic = np.sum(w*w, axis=-1)
    logdet = float(2*np.log(factor[0]).sum())
    return quadratic, logdet


def _starts(response, voltage, times, factor, initial, previous):
    objective = FixedMetricObjective(response, voltage, factor)
    energies = [np.asarray(previous, float), objective.linear_energy_start(times)]
    if len(times) == 2:
        total = max(sum(initial["charge_pe"]), 100.)
        energies.extend(response.initial_energy(total*np.array([s, 1-s]))
                        for s in (.03, .15, .5, .85, .97))
    return energies


def _fixed_update(response, voltage, times, factor, initial, previous, *, multistart):
    objective = FixedMetricObjective(response, voltage, factor)
    starts = (_starts(response, voltage, times, factor, initial, previous)
              if multistart else [np.asarray(previous, float)])
    fits = [objective.fit(times, e, free_times=False) for e in starts]
    best = min(fits, key=lambda r: r["cost"])
    return dict(best, starts=len(starts), total_evaluations=sum(f["evaluations"] for f in fits))


def fit_candidate_nll(response, voltage, times, starts, *, maxiter=120):
    """Mean and covariance use the same proposed energy and fixed known times."""
    t = np.asarray(times, float)
    target = response.transform(np.asarray(voltage)-response.baseline-response.background)
    lo, hi = np.log(response.energy_bounds)
    calls = 0

    def evaluate(x):
        nonlocal calls
        calls += 1
        curves, q = response.components(t, np.exp(x))
        factor = cholesky_banded(covariance_band(response, curves, q), lower=True, check_finite=False)
        quadratic, logdet = gaussian_terms(target-response.transform(curves.sum(0)), factor)
        return float(quadratic), logdet

    def objective(x):
        a, b = evaluate(x)
        return a+b

    candidates = []
    for energy in starts:
        x0 = np.clip(np.log(energy), lo+1e-9, hi-1e-9)
        cost0 = objective(x0)
        fit = minimize(objective, x0, method="L-BFGS-B", bounds=[(lo, hi)]*len(t),
                       options=dict(maxiter=maxiter, ftol=1e-10, gtol=1e-5,
                                    maxls=30, eps=2e-5))
        use_fit = bool(np.isfinite(fit.fun) and fit.fun <= cost0)
        x = fit.x if use_fit else x0
        quadratic, logdet = evaluate(x)
        e = np.exp(x)
        _, q = response.components(t, e)
        candidates.append(dict(times_ns=t.tolist(), energy_MeV=e.tolist(), charge_pe=q.tolist(),
            cost=quadratic+logdet, quadratic=quadratic, logdet=logdet, initial_cost=cost0,
            converged=bool(fit.success and use_fit), evaluations=int(fit.nfev),
            boundary=bool(np.any(x-lo < 1e-5) or np.any(hi-x < 1e-5)), message=str(fit.message)))
    best = min(candidates, key=lambda r: r["cost"])
    return dict(best, starts=len(starts), total_evaluations=calls,
                per_start=[{k: f[k] for k in ("cost", "converged", "evaluations", "energy_MeV")}
                           for f in candidates])


def compare_known_times(response, voltage, known_times, reference_factor, reference_initial):
    """Observe voltage plus explicitly supplied K/t; no deposited-energy input.

    The reference metric is the preceding study's independent-reference-mean
    one-pulse metric. New metrics use only energies inferred from this waveform.
    All methods stay at the same supplied times. This is not blind detection.
    """
    times = np.asarray(known_times, float)
    if times.ndim != 1 or len(times) not in (1, 2) or not np.isfinite(times).all() or np.any(np.diff(times) < 0):
        raise ValueError("Expected ordered one/two supplied times")
    baseline = fit_known_times(response, voltage, times, reference_factor, reference_initial)
    rows = {"reference_gls": baseline}
    e0 = np.asarray(baseline["energy_MeV"])
    for name, legacy in (("all_legacy_once", True), ("energy_once", False)):
        factor = candidate_factor(response, times, e0, legacy=legacy)
        fit = _fixed_update(response, voltage, times, factor, reference_initial, e0, multistart=True)
        rows[name] = dict(fit, covariance_energy_MeV=e0.tolist(), weight_updates=1)

    previous = np.asarray(rows["energy_once"]["energy_MeV"])
    history = [dict(iteration=1, covariance_energy_MeV=e0.tolist(),
                    energy_MeV=previous.tolist(), converged=rows["energy_once"]["converged"],
                    cost=rows["energy_once"]["cost"], evaluations=rows["energy_once"]["total_evaluations"])]
    stable = float(np.max(abs(np.log(previous/e0)))) <= 1e-3
    last = rows["energy_once"]
    for iteration in range(2, 5):
        if stable:
            break
        point = previous.copy()
        factor = candidate_factor(response, times, point)
        last = _fixed_update(response, voltage, times, factor, reference_initial, previous, multistart=False)
        previous = np.asarray(last["energy_MeV"])
        change = float(np.max(abs(np.log(previous/point))))
        history.append(dict(iteration=iteration, covariance_energy_MeV=point.tolist(),
                            energy_MeV=previous.tolist(), converged=last["converged"],
                            cost=last["cost"], evaluations=last["total_evaluations"], log_energy_change=change))
        stable = change <= 1e-3
    rows["energy_fgls"] = dict(last, covariance_energy_MeV=history[-1]["covariance_energy_MeV"],
        weight_updates=len(history), outer_stable=stable, history=history,
        total_evaluations=sum(h["evaluations"] for h in history))
    # No cross-method cost comparison: every GLS metric has a different scale.
    starts = _starts(response, voltage, times, reference_factor, reference_initial, e0)
    starts.append(np.asarray(rows["energy_once"]["energy_MeV"]))
    rows["energy_nll"] = fit_candidate_nll(response, voltage, times, starts)
    return rows
