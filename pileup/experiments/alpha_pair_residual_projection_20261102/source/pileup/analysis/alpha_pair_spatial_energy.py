"""Paired known-K/t energy fits; no truth energy or photon data in inference.

Uniform and spatial means share the same voltage-derived weighting, bounds,
three starts and optimizer budget. The frozen covariance is a working metric,
not a newly calibrated spatial-model likelihood.
"""
from __future__ import annotations

import time
import numpy as np
from scipy.optimize import least_squares, minimize_scalar

from .alpha_pair_energy_calibrated import learned_metric


def common_inputs(model, voltage, times, noise):
    """Build a single common metric and starts using voltage and known t only."""
    y = np.asarray(voltage, float)
    t = np.asarray(times, float)
    if y.shape != (model.anchor.observed_samples,) or not np.isfinite(y).all():
        raise ValueError('Expected the complete finite 2048 ns analog window')
    if t.ndim != 1 or len(t) not in (1, 2) or np.any(np.diff(t) <= 0):
        raise ValueError('Expected one/two ordered declared event times')
    initial = model.initial_energy(y, t)
    metric = learned_metric(model.anchor, t, initial, noise)
    starts = [initial]
    if len(t) == 2:
        total = float(model.anchor.q_isolated(initial).sum())
        starts += [model.anchor.initial_energy(total * np.array([s, 1-s])) for s in (.1, .9)]
    return dict(initial=initial, metric=metric, starts=np.asarray(starts))


def fit_energy(model, voltage, times, common, *, max_nfev=80):
    """Optimize E only, with exactly the same input/metric for both means."""
    started = time.monotonic()
    t = np.asarray(times, float)
    metric = common['metric']
    target = metric.apply(model.anchor.transform(np.asarray(voltage) - model.background))
    lo, hi = np.log(model.energy_bounds)
    if len(t) == 1:
        def scalar(loge):
            curve = model.isolated(t, [np.exp(loge)])[0][0]
            residual = metric.apply(model.anchor.transform(curve)) - target
            return float(residual @ residual)
        fit = minimize_scalar(scalar, bounds=(lo, hi), method='bounded',
                              options={'xatol': 1e-10, 'maxiter': 150})
        return dict(energy_MeV=[float(np.exp(fit.x))], cost=float(fit.fun),
                    converged=bool(fit.success), boundary=bool(min(fit.x-lo, hi-fit.x) < 1e-5),
                    starts=1, evaluations=int(fit.nfev), initial_MeV=common['initial'].tolist(),
                    start_energy_MeV=common['starts'].tolist(), seconds=time.monotonic()-started)
    probabilities = model.prepare(t)
    cached_x, cached = None, None

    def evaluate(loge):
        nonlocal cached_x, cached
        if cached_x is not None and np.array_equal(loge, cached_x):
            return cached
        mean, gradient = model.mean_and_gradient(t, loge, probabilities, strength=1.)
        residual = metric.apply(model.anchor.transform(mean)) - target
        jacobian = np.column_stack([metric.apply(model.anchor.transform(g)) for g in gradient])
        cached_x, cached = loge.copy(), (residual, jacobian)
        return cached

    candidates = []
    for start in common['starts']:
        fit = least_squares(lambda x: evaluate(x)[0],
            np.clip(np.log(start), lo+1e-8, hi-1e-8), jac=lambda x: evaluate(x)[1],
            bounds=(lo, hi), max_nfev=max_nfev, ftol=1e-8, xtol=1e-8, gtol=1e-5)
        candidates.append(dict(energy_MeV=np.exp(fit.x).tolist(), cost=float(fit.fun @ fit.fun),
            converged=bool(fit.success), boundary=bool(np.any(fit.x-lo < 1e-5) or np.any(hi-fit.x < 1e-5)),
            evaluations=int(fit.nfev), gradient_evaluations=int(fit.njev),
            optimality=float(fit.optimality), status=int(fit.status), message=str(fit.message)))
    best_index = min(range(len(candidates)), key=lambda i: candidates[i]['cost'])
    best = dict(candidates[best_index])
    best.update(starts=len(candidates), evaluations=sum(c['evaluations'] for c in candidates),
                selected_start=best_index, candidates=candidates,
                initial_MeV=common['initial'].tolist(), start_energy_MeV=common['starts'].tolist(),
                seconds=time.monotonic()-started)
    return best
