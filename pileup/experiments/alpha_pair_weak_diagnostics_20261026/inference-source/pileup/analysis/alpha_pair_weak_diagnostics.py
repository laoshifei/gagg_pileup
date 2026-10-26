"""Truth-assisted development diagnostics, never a blind reconstruction API.

The supplied times, charges and ancestry components are simulation labels.
They diagnose a frozen response/search on exposed data and do not constitute
achievable blind performance or a theoretical resolution bound.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from .alpha_pair_residual import PrecisionBlend


def model_pilot_metric(response, voltage):
    """Reproduce the previous observation-only pilot, independently of labels."""
    raw = response.base.candidates(voltage)
    if len(raw) > 1:
        one = response.refine(voltage, raw[1], radius=24.)
        smooth = response.components(one['times_ns'], one['charge_pe']).sum(0) + response.baseline
        return PrecisionBlend(response, smooth, 0., rho=.9)
    return PrecisionBlend(response, voltage, 0., rho=.9)


def fit_cost(response, voltage, fit, metric=None):
    residual = response.components(fit['times_ns'], fit['charge_pe']).sum(0)
    residual = residual - (np.asarray(voltage, float) - response.baseline)
    if metric is not None:
        residual = metric.apply(residual)
    return float(residual @ residual)


def fixed_time_charge_fit(response, voltage, true_times, total_start, *, metric=None,
                          max_nfev=120):
    """Known K and true times; Q starts use an observed total, not true Q."""
    times = np.sort(np.asarray(true_times, float))
    if times.ndim != 1 or len(times) > 2 or not np.isfinite(times).all():
        raise ValueError('Expected zero, one or two finite true times')
    y = np.asarray(voltage, float) - response.baseline
    if not len(times):
        r = y if metric is None else metric.apply(y)
        return dict(times_ns=[], charge_pe=[], cost=float(r @ r), converged=True,
                    evaluations=0, diagnostic_truth=['K', 'times'])
    total = max(float(total_start), 100.)
    starts = [[total]] if len(times) == 1 else [
        [total * share, total * (1.-share)] for share in (.03, .15, .5, .85, .97)]
    def residual(z):
        r = response.components(times, np.exp(z)).sum(0) - y
        return r if metric is None else metric.apply(r)
    fits = []
    for q in starts:
        z0 = np.log(np.clip(q, 50.0001, 49999.99))
        result = least_squares(residual, z0, bounds=(np.log(50.), np.log(50000.)),
                               max_nfev=max_nfev, ftol=1e-9, xtol=1e-9,
                               gtol=1e-7, diff_step=1e-4)
        z = result.x if result.fun @ result.fun <= residual(z0) @ residual(z0) else z0
        r = residual(z)
        fits.append(dict(times_ns=times.tolist(), charge_pe=np.exp(z).tolist(),
                         cost=float(r @ r), converged=bool(result.success),
                         evaluations=int(result.nfev), diagnostic_truth=['K', 'times']))
    return min(fits, key=lambda f: f['cost'])


def oracle_component_refit(response, voltage, true_times, true_charge, true_components,
                           *, metric, radius=24., max_nfev=120):
    """Remove the exact other component, then fit the weaker component.

    Other-event Q/t remain true conditioning variables. Both weak parameters
    start at truth. Dark events and electronics noise remain in the data.
    The returned full pair also contains the fixed true other event; its
    success must therefore never be counted as a blind pair recovery rate.
    """
    t = np.asarray(true_times, float)
    q = np.asarray(true_charge, float)
    curves = np.asarray(true_components, float)
    if t.shape != (2,) or q.shape != (2,) or curves.shape != (2, 1024):
        raise ValueError('Oracle component diagnosis needs exactly two labeled components')
    weak = int(np.argmin(q)); other = 1-weak
    y = np.asarray(voltage, float) - response.baseline - curves[other]
    lower = [max(0., t[weak]-radius)/10., np.log(50.)]
    upper = [min(response.sample_grid[-1], t[weak]+radius)/10., np.log(50000.)]
    z0 = np.clip([t[weak]/10., np.log(q[weak])], np.asarray(lower)+1e-9, np.asarray(upper)-1e-9)
    def residual(z):
        time = z[0]*10.; charge = np.exp(z[1])
        model = charge*response.unit(response.sample_grid-time, charge, q[other], t[other]-time)
        return metric.apply(model-y)
    result = least_squares(residual, z0, bounds=(lower, upper), max_nfev=max_nfev,
                           ftol=1e-9, xtol=1e-9, gtol=1e-7, diff_step=1e-4)
    z = result.x if result.fun @ result.fun <= residual(z0) @ residual(z0) else z0
    found_t = t.copy(); found_q = q.copy()
    found_t[weak] = z[0]*10.; found_q[weak] = np.exp(z[1])
    order = np.argsort(found_t); r = residual(z)
    return dict(times_ns=found_t[order].tolist(), charge_pe=found_q[order].tolist(),
                cost=float(r @ r), converged=bool(result.success), evaluations=int(result.nfev),
                weak_truth_index=weak, diagnostic_truth=['K', 'other_Q_t', 'other_component', 'initial_Q_t'],
                cost_scope='Only the oracle-subtracted weak-component residual; not a full-frame cost')


def diagnose(response, voltage, true_times, true_charge, true_components, previous_candidates):
    """Return fixed-model oracles and cost comparisons on a labeled frame."""
    y = np.asarray(voltage, float)
    ids = np.flatnonzero(np.asarray(true_charge) >= 500.)
    t = np.asarray(true_times)[ids]; q = np.asarray(true_charge)[ids]
    metric = model_pilot_metric(response, y)
    replay_errors = []
    for fit in previous_candidates:
        cost = fit_cost(response, y, fit, metric)
        error = abs(cost-fit['cost'])/max(1., abs(fit['cost']))
        replay_errors.append(float(error))
        if error > 1e-7:
            raise AssertionError(('The frozen metric was not replayed', error))
    initial = dict(times_ns=t.tolist(), charge_pe=q.tolist())
    total = sum(previous_candidates[min(1, len(previous_candidates)-1)]['charge_pe'])
    fixed = fixed_time_charge_fit(response, y, t, total, metric=metric)
    fixed_ols = fixed_time_charge_fit(response, y, t, total, metric=None)
    joint = response.refine(y, initial, metric=metric, radius=24., max_nfev=120)
    joint['diagnostic_truth'] = ['K', 'initial_Q_t']
    joint_ols = response.refine(y, initial, metric=None, radius=24., max_nfev=120)
    joint_ols['diagnostic_truth'] = ['K', 'initial_Q_t']
    truth_cost = fit_cost(response, y, initial, metric)
    assert joint['cost'] <= truth_cost*(1.+1e-10)+1e-7
    kfit = previous_candidates[min(len(ids), len(previous_candidates)-1)]
    result = dict(fixed_true_times=fixed, fixed_true_times_ols=fixed_ols,
                  truth_start_joint=joint, truth_start_joint_ols=joint_ols,
                  known_count_previous=kfit, exact_truth_cost=truth_cost,
                  metric_replay_relative_error=max(replay_errors, default=0.),
                  cost_truth_start_minus_previous=joint['cost']-kfit['cost'])
    if len(ids) == 2:
        result['oracle_subtract_other'] = oracle_component_refit(
            response, y, t, q, np.asarray(true_components)[ids], metric=metric)
    return result
