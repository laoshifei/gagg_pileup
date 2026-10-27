"""Explicitly truth-assisted local profiles; not a blind reconstruction API."""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from .alpha_pair_joint_likelihood import refine


def cost_at(response, noise, voltage, times, charges):
    curves = response.components(times, charges)
    residual = np.asarray(voltage, float) - response.baseline - noise.background_mV - curves.sum(0)
    return noise.score(residual, curves, charges, details=True)


def profile_weak(response, noise, voltage, weak_time, weak_charge, strong_starts,
                 strong_bounds, maxiter=60):
    """Hold the labelled weak Q/t fixed; optimize the other Q/t jointly.

    This is a bounded local nuisance profile, not a certified global optimum.
    The neighbourhood, fixed count and weak label are diagnostic inputs.
    Covariance and log determinant are recomputed for every nuisance proposal.
    """
    zy = noise.transform(np.asarray(voltage, float) - response.baseline - noise.background_mV)
    lower = np.array([strong_bounds[0] / 10., np.log(50.)])
    upper = np.array([strong_bounds[1] / 10., np.log(50000.)])

    def objective(z):
        charges = [np.exp(z[1]), weak_charge]
        curves = response.components([z[0] * 10., weak_time], charges)
        transformed = noise.transform(curves)
        return noise.score(zy - transformed.sum(0), transformed, charges, transformed=True)

    fits = []
    for strong_time, strong_charge in strong_starts:
        initial = np.clip([strong_time / 10., np.log(max(strong_charge, 50.))], lower + 1e-9, upper - 1e-9)
        fit = minimize(objective, initial, method='L-BFGS-B', bounds=list(zip(lower, upper)),
                       options=dict(maxiter=maxiter, ftol=1e-9, gtol=1e-4, maxls=30, eps=1e-5))
        point = fit.x if fit.fun <= objective(initial) else initial
        fits.append(dict(times_ns=[float(point[0] * 10.), float(weak_time)],
                         charge_pe=[float(np.exp(point[1])), float(weak_charge)],
                         cost=float(objective(point)), converged=bool(fit.success),
                         evaluations=int(fit.nfev), message=str(fit.message)))
    best = min(fits, key=lambda r: r['cost'])
    return dict(best, total_evaluations=sum(r['evaluations'] for r in fits),
                diagnostic_truth=['K=2', 'weak role', 'profile region'],
                scope='Bounded local nuisance minimization; no Wilks confidence interpretation')


def truth_diagnosis(response, noise, voltage, times, charges, baseline):
    """Known K, true initialization and fixed-time Q fits, kept separate."""
    t = np.asarray(times, float); q = np.asarray(charges, float)
    initial = dict(times_ns=t.tolist(), charge_pe=q.tolist())
    assisted = refine(response, noise, voltage, initial, radius=48., maxiter=200)
    assisted['diagnostic_truth'] = ['K=2', 'initial_Q_t', 'time_bounds']
    zy = noise.transform(np.asarray(voltage, float) - response.baseline - noise.background_mV)

    def objective(z):
        qq = np.exp(z); curves = noise.transform(response.components(t, qq))
        return noise.score(zy - curves.sum(0), curves, qq, transformed=True)

    total = sum(baseline[1]['charge_pe'])
    starts = [np.array([total * f, total * (1. - f)]) for f in (.01, .03, .10, .5, .9, .97, .99)]
    fixed = []
    for qq in starts:
        z0 = np.log(np.clip(qq, 50.0001, 49999.99))
        fit = minimize(objective, z0, method='L-BFGS-B', bounds=[(np.log(50.), np.log(50000.))] * 2,
                       options=dict(maxiter=100, ftol=1e-9, gtol=1e-4, maxls=30, eps=1e-5))
        z = fit.x if fit.fun <= objective(z0) else z0
        fixed.append(dict(times_ns=t.tolist(), charge_pe=np.exp(z).tolist(), cost=float(objective(z)),
                          converged=bool(fit.success), evaluations=int(fit.nfev),
                          diagnostic_truth=['K=2', 'fixed_true_times']))
    return dict(truth_point=cost_at(response, noise, voltage, t, q),
                truth_start=assisted, fixed_times=min(fixed, key=lambda r: r['cost']),
                baseline_double=baseline[2], truth_used=True,
                warning='A lower cost than the exact truth is expected in noisy data; it does not alone diagnose model bias')


def profile_map(response, noise, voltage, times, charges, baseline, *, dense=False):
    t = np.asarray(times, float); q = np.asarray(charges, float)
    weak = int(np.argmin(q)); strong = 1 - weak
    # Fixed protocol, not chosen after looking at which frames succeeded.
    offsets = np.array([-48., -32., -24., -16., -8., -4., 0., 4., 8., 16., 24., 32., 48.])
    factors = np.array([.25, .4, .6, .8, 1., 1.2, 1.5, 2., 3.])
    if dense:
        offsets = np.sort(np.unique(np.r_[offsets, (offsets[:-1] + offsets[1:]) / 2.]))
        factors = np.sort(np.unique(np.r_[factors, np.sqrt(factors[:-1] * factors[1:])]))
    tt = t[weak] + offsets
    qq = np.maximum(50., q[weak] * factors)
    strong_bounds = (max(0., t[strong] - 48.), min(response.sample_grid[-1], t[strong] + 48.))
    one = baseline[1]
    # The diagnostic permits a true strong start; the second uses only the old fit.
    starts = [(t[strong], q[strong]), (one['times_ns'][0], one['charge_pe'][0])]
    costs = np.zeros((len(qq), len(tt))); found_t = costs.copy(); found_q = costs.copy()
    converged = np.zeros(costs.shape, bool); evaluations = np.zeros(costs.shape, int)
    for j, qw in enumerate(qq):
        for i, tw in enumerate(tt):
            fit = profile_weak(response, noise, voltage, tw, qw, starts, strong_bounds)
            costs[j, i] = fit['cost']; found_t[j, i] = fit['times_ns'][0]; found_q[j, i] = fit['charge_pe'][0]
            converged[j, i] = fit['converged']; evaluations[j, i] = fit['total_evaluations']
    j, i = np.unravel_index(np.argmin(costs), costs.shape)
    initial = dict(times_ns=[float(found_t[j, i]), float(tt[i])], charge_pe=[float(found_q[j, i]), float(qq[j])])
    polished = refine(response, noise, voltage, initial, radius=48., maxiter=200)
    polished['diagnostic_truth'] = ['K=2', 'profile_region', 'weak_role', 'nuisance_start']
    return dict(weak_times_ns=tt.tolist(), weak_charge_pe=qq.tolist(), costs=costs.tolist(),
                strong_times_ns=found_t.tolist(), strong_charge_pe=found_q.tolist(),
                converged=converged.tolist(), evaluations=evaluations.tolist(),
                weak_truth_index=weak, true_t=t.tolist(), true_q=q.tolist(),
                best_grid=initial, polished=polished, dense=bool(dense),
                grid_min_cost=float(costs.min()), truth_used=True,
                scope='Truth-centred local diagnostic profile; not blind performance or a global identifiability bound')
