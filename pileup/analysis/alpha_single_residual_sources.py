"""Grouped repeat diagnostics and voltage-only shape prediction helpers.

Mother events, not repeated digitizations, define training/validation units.
Negative finite-sample between-mother variance estimates are retained.
This module does not perform blind pileup reconstruction.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar


def variance_components(values):
    """Balanced random-effects moments; axes mother, repeat, [...]."""
    x = np.asarray(values, float)
    if x.ndim < 2 or min(x.shape[:2]) < 2 or not np.isfinite(x).all():
        raise ValueError('Need finite values for at least two mothers/repeats')
    m, r = x.shape[:2]
    means = x.mean(1)
    within = x.var(1, ddof=1).mean(0)
    mean_variance = means.var(0, ddof=1)
    between = mean_variance - within/r
    empirical = ((x-x.mean((0, 1)))**2).mean((0, 1))
    identity = (m-1)/m*mean_variance + (r-1)/r*within
    return dict(mean=x.mean((0, 1)), within=within, between=between,
                variance_of_repeat_mean=mean_variance,
                grand_mean_sampling_variance=mean_variance/m,
                empirical_variance=empirical,
                closure_error=float(np.max(abs(empirical-identity))))


def component_intervals(values, *, draws=4000, seed=735201):
    """Mother-cluster percentile intervals, fixed collected repetitions."""
    x = np.asarray(values, float)
    if x.ndim != 3:
        raise ValueError('Expected mother, repeat, scalar diagnostic axes')
    m, r = x.shape[:2]
    idx = np.random.default_rng(seed).integers(m, size=(draws, m))
    means, within = x.mean(1), x.var(1, ddof=1)
    sampled_mean = means[idx].mean(1)
    sampled_within = within[idx].mean(1)
    sampled_between = means[idx].var(1, ddof=1)-sampled_within/r
    return {k: np.percentile(v, [2.5, 97.5], axis=0) for k, v in
            [('mean', sampled_mean), ('within', sampled_within), ('between', sampled_between)]}


def opposite_half_targets(values):
    """Each input repeat targets the mean of the *other* half of its mother."""
    x = np.asarray(values, float)
    if x.ndim < 2 or x.shape[1] < 4 or x.shape[1] % 2:
        raise ValueError('Expected an even repeat count >=4')
    h = x.shape[1]//2
    return np.concatenate((np.repeat(x[:, h:].mean(1, keepdims=True), h, axis=1),
                           np.repeat(x[:, :h].mean(1, keepdims=True), h, axis=1)), axis=1)


def prefix_features(model, voltage, horizon_ns, *, onset_ns=128., noise_floor_mV=.1):
    """Known-onset K=1 diagnostic using only complete apertures in a prefix.

    Return log voltage-fitted energy plus four normalized residual block means.
    True E, source energy, photons and ancestor signals are not accepted.
    The full-window result is an isolated-signal information control only.
    """
    y = np.asarray(voltage, float)
    centers = model.readout.adc_centers
    half = model.readout.aperture/2
    mask = (centers-half >= onset_ns) & (centers+half <= onset_ns+horizon_ns)
    indices = np.flatnonzero(mask)
    if y.shape != centers.shape or len(indices) < 4 or not np.isfinite(y[indices]).all():
        raise ValueError('At least four finite fully observed apertures required')
    target = (y-model.background)[indices]
    lo, hi = np.log(model.energy_bounds)

    def curve(loge):
        return model.isolated([onset_ns], [np.exp(loge)])[0][0][indices]

    fit = minimize_scalar(lambda e: float(np.sum((curve(e)-target)**2)),
                          bounds=(lo, hi), method='bounded',
                          options={'xatol': 1e-9, 'maxiter': 150})
    predicted = curve(fit.x)
    residual = target-predicted
    scale = max(float(np.sqrt(np.mean(predicted**2))), noise_floor_mV)
    feature = np.r_[fit.x, [a.mean()/scale for a in np.array_split(residual, 4)]]
    return feature, dict(energy_MeV=float(np.exp(fit.x)), converged=bool(fit.success),
                        boundary=bool(min(fit.x-lo, hi-fit.x)<1e-5),
                        samples=len(indices), aperture_end_ns=float(centers[indices[-1]]+half))


def ridge_fit(x, y, penalty):
    x = np.asarray(x, float).reshape(-1, x.shape[-1])
    y = np.asarray(y, float).reshape(-1)
    center = x.mean(0)
    scale = np.maximum(x.std(0), 1e-12)
    z = (x-center)/scale
    intercept = float(y.mean())
    beta = np.linalg.solve(z.T@z/len(z)+penalty*np.eye(z.shape[1]), z.T@(y-intercept)/len(z))
    return dict(center=center.tolist(), scale=scale.tolist(), beta=beta.tolist(), intercept=intercept,
                penalty=float(penalty))


def ridge_predict(fit, x):
    return ((np.asarray(x)-fit['center'])/fit['scale'])@np.asarray(fit['beta'])+fit['intercept']


def select_ridge(x, y, *, penalties=(.01, .1, 1., 10.)):
    """Leave-one-mother-out selection; repeats always stay with their mother."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.ndim != 3 or y.shape != x.shape[:2] or len(x) < 4:
        raise ValueError('Need grouped features and targets')
    losses = []
    for penalty in penalties:
        loss = []
        for i in range(len(x)):
            keep = np.arange(len(x)) != i
            fit = ridge_fit(x[keep], y[keep], penalty)
            loss.append(float(np.mean((ridge_predict(fit, x[i])-y[i])**2)))
        losses.append(float(np.mean(loss)))
    chosen = int(np.argmin(losses))
    return ridge_fit(x, y, penalties[chosen]), dict(penalties=list(penalties),
        leave_one_mother_out_mse=losses, selected=chosen, mother_count=len(x))


def paired_loss_ratio(candidate, baseline, *, seed=735301, draws=4000):
    """Inputs are per-mother MSE; CI conditions on the frozen training model."""
    a, b = np.asarray(candidate, float), np.asarray(baseline, float)
    if a.shape != b.shape or a.ndim != 1 or len(a) < 2 or np.mean(b) <= 0:
        raise ValueError('Expected aligned nonzero-baseline mother losses')
    idx = np.random.default_rng(seed).integers(len(a), size=(draws, len(a)))
    ratios = a[idx].mean(1)/b[idx].mean(1)
    return dict(ratio=float(a.mean()/b.mean()), ci95=np.percentile(ratios, [2.5, 97.5]).tolist(),
                candidate_rmse=float(np.sqrt(a.mean())), baseline_rmse=float(np.sqrt(b.mean())))
