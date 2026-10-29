"""Observation-only joint-parameter compatibility and ambiguity diagnostics.

The reference is the FROZEN K2 candidate, not a claimed global optimum. An
explicit low-cost exterior point proves ambiguity within the chosen cost set;
failure of a finite search to find one does not prove set containment.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.optimize import minimize

DEFAULTS = dict(starts=2, maxiter=100, verify_starts=4, verify_maxiter=180,
                stability_cost=.05, q_tolerance=.20, t_tolerance_ns=10.)


def domain_for(response, one, candidate):
    anchor = float(one['times_ns'][0])
    times = candidate['times_ns']
    return dict(time=[max(0., min(anchor-96., min(times))),
                      min(float(response.sample_grid[-1]), max(anchor+512., max(times)))],
                charge=[50., 50000.])


def in_domain(times, charges, domain):
    t, q = np.asarray(times, float), np.asarray(charges, float)
    return bool(t.shape == (2,) and q.shape == (2,) and np.isfinite(t).all() and
        np.isfinite(q).all() and domain['time'][0] <= t[0] <= t[1] <= domain['time'][1]
        and np.all(q >= domain['charge'][0]) and np.all(q <= domain['charge'][1]))


def exterior_boxes(candidate, domain, q_tolerance=.2, t_tolerance_ns=10.):
    """Eight closed exterior boxes, including the accuracy boundary.

    The recovery error divides by TRUTH, so a possible true Q compatible with
    a fixed estimate Qhat lies in [Qhat/1.2,Qhat/.8], not Qhat*(1+/-0.2).
    Each box is intersected with t1<=t2 and the declared fitting domain.
    """
    t, q = np.asarray(candidate['times_ns']), np.asarray(candidate['charge_pe'])
    if not in_domain(t, q, domain):
        raise ValueError('Reference candidate outside its declared domain')
    low_t, high_t = domain['time']; low_q, high_q = domain['charge']
    limits = [list(domain['time']), list(domain['time']),
              [math.log(low_q), math.log(high_q)], [math.log(low_q), math.log(high_q)]]
    boundaries = [t[0]-t_tolerance_ns, t[0]+t_tolerance_ns,
                  t[1]-t_tolerance_ns, t[1]+t_tolerance_ns,
                  math.log(q[0]/(1+q_tolerance)), math.log(q[0]/(1-q_tolerance)),
                  math.log(q[1]/(1+q_tolerance)), math.log(q[1]/(1-q_tolerance))]
    boxes = []
    for j, (coordinate, side) in enumerate((i, s) for i in range(4) for s in ('lower', 'upper')):
        box = [x.copy() for x in limits]
        if side == 'lower': box[coordinate][1] = min(box[coordinate][1], boundaries[j])
        else: box[coordinate][0] = max(box[coordinate][0], boundaries[j])
        box[0][1] = min(box[0][1], box[1][1])
        box[1][0] = max(box[1][0], box[0][0])
        possible = all(a <= b for a, b in box)
        boxes.append(dict(name=('t1', 't2', 'Q1', 'Q2')[coordinate]+'_'+side,
                          coordinate=coordinate, side=side, limits=box, possible=possible))
    return boxes


def box_encode(times, charges, limits):
    t = np.sort(np.asarray(times, float)); q = np.asarray(charges, float)
    t1 = float(np.clip(t[0], *limits[0]))
    lower2 = max(t1, limits[1][0]); upper2 = limits[1][1]
    t2 = float(np.clip(t[1], lower2, upper2))
    fraction = (t2-lower2)/max(upper2-lower2, 1e-12)
    return np.array([t1/10., 100.*fraction,
                     np.clip(np.log(q[0]), *limits[2]), np.clip(np.log(q[1]), *limits[3])])


def box_decode(z, limits):
    t1 = float(10.*z[0]); lower2 = max(t1, limits[1][0])
    t2 = lower2+(limits[1][1]-lower2)*float(z[1])/100.
    return [t1, float(t2)], np.exp(np.asarray(z[2:])).tolist()


def observed_seeds(one, candidate, domain):
    seeds = [candidate]
    anchor, total = one['times_ns'][0], sum(one['charge_pe'])
    for direction in (-1., 1.):
        for gap in (8., 40., 200.):
            times = sorted([anchor, anchor+direction*gap])
            times = np.clip(times, *domain['time']).tolist()
            for share in (.03, .1, .5, .9, .97):
                seeds.append(dict(times_ns=times, charge_pe=[max(50., total*share), max(50., total*(1-share))]))
    return seeds


def search_box(objective, box, seeds, starts, maxiter):
    if not box['possible']:
        return dict(name=box['name'], possible=False, best=None, trials=[])
    limits = box['limits']; unique = {}
    for seed in seeds:
        z = box_encode(seed['times_ns'], seed['charge_pe'], limits)
        unique.setdefault(tuple(z), z)
    def fun(z):
        t, q = box_decode(z, limits)
        return float(objective(t, q))
    ranked = sorted(((fun(z), z) for z in unique.values()), key=lambda p: p[0])
    chosen = []
    for value, z in ranked:
        t, q = box_decode(z, limits); position = np.r_[np.array(t)/8., np.log(q)/.25]
        if not chosen or min(np.linalg.norm(position-p[2]) for p in chosen) >= 1.:
            chosen.append((value, z, position))
        if len(chosen) == starts: break
    trial_bounds = [(limits[0][0]/10., limits[0][1]/10.), (0., 100.), limits[2], limits[3]]
    trials = []
    for initial_cost, initial, _ in chosen:
        fit = minimize(fun, initial, method='L-BFGS-B', bounds=trial_bounds,
            options=dict(maxiter=maxiter, ftol=1e-10, gtol=1e-5, maxls=30, eps=1e-5))
        use = bool(np.isfinite(fit.fun) and fit.fun <= initial_cost)
        z = fit.x if use else initial; t, q = box_decode(z, limits)
        trials.append(dict(times_ns=t, charge_pe=q, cost=fun(z), initial_cost=initial_cost,
            converged=bool(use and fit.success), evaluations=int(fit.nfev), message=str(fit.message)))
    return dict(name=box['name'], possible=True, best=min(trials, key=lambda x: x['cost']), trials=trials)


def exterior_search(objective, one, candidate, domain, settings=None):
    cfg = dict(DEFAULTS if settings is None else settings)
    boxes = exterior_boxes(candidate, domain, cfg['q_tolerance'], cfg['t_tolerance_ns'])
    seeds = observed_seeds(one, candidate, domain)
    first = [search_box(objective, b, seeds, cfg['starts'], cfg['maxiter']) for b in boxes]
    extras = [r['best'] for r in first if r['best'] is not None]
    second = [search_box(objective, b, seeds+extras, cfg['verify_starts'], cfg['verify_maxiter']) for b in boxes]
    retained = []
    for a, b in zip(first, second):
        if b['best'] is not None and a['best']['cost'] < b['best']['cost']:
            b = dict(b, best=a['best'], retained_first=True)
        retained.append(b)
    a = min((x['best'] for x in first if x['best']), key=lambda r: r['cost'])
    b = min((x['best'] for x in retained if x['best']), key=lambda r: r['cost'])
    reference = float(objective(candidate['times_ns'], candidate['charge_pe']))
    gain = float(a['cost']-b['cost'])
    assert gain >= -1e-7
    replay = max(abs(objective(r['best']['times_ns'], r['best']['charge_pe'])-r['best']['cost'])
                 / max(1., abs(r['best']['cost'])) for r in retained if r['best'])
    edge = bool(any(not box['possible'] for box in boxes))
    return dict(domain=domain, reference_cost=reference, first=first, verified=retained,
        best_exterior=b, exterior_delta=float(b['cost']-reference), verification_gain=gain,
        stable=bool(gain <= cfg['stability_cost']), cost_replay_error=float(replay),
        all_exterior_minima_converged=all(r['best']['converged'] for r in retained if r['best']),
        tolerance_box_hits_domain=edge, truth_inputs=[],
        meaning='Finite exterior search; low-cost witness proves ambiguity, absent witness is not a proof of containment')


def waveform_objective(response, noise, voltage):
    observed = np.asarray(voltage, float)-response.baseline-noise.background_mV
    transformed = noise.transform(observed)
    def objective(times, charges):
        curves = noise.transform(response.components(times, charges))
        return float(noise.score(transformed-curves.sum(0), curves, charges, transformed=True))
    return objective


def empirical_level(scores, quantile=.95, minimum=20):
    """Prespecified empirical upper quantile, not a coverage theorem.

    Out-of-domain calibration truths must be encoded as +inf and retained.
    """
    values = sorted(float(s) for s in scores)
    if any(np.isnan(values)) or any(np.isneginf(values)) or not 0 < quantile <= 1:
        raise ValueError('Finite or +inf scores and a valid quantile required')
    rank = int(math.ceil(quantile*len(values))) if values else 0
    value = values[rank-1] if len(values) >= minimum else math.inf
    return dict(level=max(0., value) if np.isfinite(value) else None,
                status='finite' if np.isfinite(value) else 'insufficient_or_unbounded',
                samples=len(values), order=rank, empirical_target=quantile,
                minimum_samples=minimum, guarantee='No nominal or uniform coverage claim')


def diagnostic_label(row, level):
    if len(row['selected']['times_ns']) != 2: return 'not_selected_K2'
    if level is None: return 'uncalibrated'
    search = row['exterior']
    if search['exterior_delta'] <= level: return 'alternative_found'
    if (not search['stable'] or not search['all_exterior_minima_converged'] or
        search['tolerance_box_hits_domain'] or not row['selected'].get('converged', True)):
        return 'numerically_unresolved'
    return 'no_alternative_found'
