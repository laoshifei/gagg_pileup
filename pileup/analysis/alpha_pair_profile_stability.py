"""Finite numerical checks of conditional K=2 later-charge profiles.

All starts, grids and refinement decisions use observed fits only. This module
does not load labels or assign confidence levels to cost contours.
"""
from __future__ import annotations

import copy
import time
import numpy as np

from .alpha_pair_evidence_profile import (
    encode, decode, fit_at_charge, seeds_for, compatibility_segments,
)


def midpoint_grid(q):
    q = np.asarray(q, float)
    if len(q) < 2 or np.any(q <= 0) or np.any(np.diff(q) <= 0):
        raise ValueError('A positive, sorted, unique charge grid is required')
    return np.exp((np.log(q[:-1]) + np.log(q[1:])) / 2.)


def contour_midpoints(curve, levels=(1., 4., 9.)):
    """Split every sampled crossing and both cells adjacent to the minimum."""
    q = np.asarray(curve['q_grid_pe'], float)
    d = np.array([f['cost'] for f in curve['fits']]) - curve['reference_cost']
    mask = np.zeros(len(q) - 1, bool)
    for level in levels:
        mask |= (d[:-1] <= level) != (d[1:] <= level)
    best = int(np.argmin(d))
    if best: mask[best-1] = True
    if best < len(q)-1: mask[best] = True
    return midpoint_grid(q)[mask]


def point_cost(response, noise, voltage, fit, q2=None):
    q = list(fit['charge_pe'])
    if q2 is not None: q[1] = float(q2)
    y = np.asarray(voltage, float) - response.baseline - noise.background_mV
    m = response.components(fit['times_ns'], q)
    return float(noise.score(y-m.sum(0), m, q))


def diverse_starts(response, noise, voltage, q2, seeds, bounds, count):
    """Cost-ranked, separated observed starts; distance is a numerical choice."""
    unique = {}
    for seed in seeds:
        z = encode(seed['times_ns'], seed['charge_pe'][0], bounds)
        t, q = decode(z, q2, bounds)
        fit = dict(times_ns=t, charge_pe=q)
        unique.setdefault(tuple(z.tolist()), fit)
    ranked = sorted(unique.values(), key=lambda f: point_cost(response, noise, voltage, f))
    if not ranked: raise ValueError('At least one observed start is required')
    def position(f):
        return np.r_[np.asarray(f['times_ns'])/8., np.log(f['charge_pe'][0])/.25]
    chosen = []
    for f in ranked:
        if not chosen or min(np.linalg.norm(position(f)-position(g)) for g in chosen) >= 1.:
            chosen.append(f)
        if len(chosen) == count: break
    if len(chosen) < count:
        chosen.extend(f for f in ranked if all(f is not g for g in chosen))
    return chosen[:count]


def observed_starts(row, q2, bounds, previous, warm=None):
    one = row['sets']['baseline'][1]
    seeds = seeds_for(one, row['sets']['baseline'][2], row['sets']['mlp'][2], q2, bounds, warm)
    fits = previous['fits']
    seeds += sorted(fits, key=lambda f: abs(np.log(f['charge_pe'][1]/q2)))[:4]
    best = min(fits, key=lambda f: f['cost'])
    seeds.append(best)
    trigger = one['times_ns'][0]
    q1 = float(np.clip(sum(one['charge_pe'])-q2, 50., 50000.))
    for offset in (-8., 0., 8.):
        for gap in (0., 4., 8., 16., 32., 64., 128., 256., 512.):
            seeds.append(dict(times_ns=[trigger+offset, trigger+offset+gap], charge_pe=[q1, q2]))
    for factor in (.75, 1.25):
        seeds.append(dict(times_ns=best['times_ns'], charge_pe=[best['charge_pe'][0]*factor, q2]))
    return seeds


def refine_stage(response, noise, voltage, row, previous, new_q, *, refit_existing,
                 starts, maxiter, label):
    """Retain every old computed point, and add a finite search/slice envelope."""
    began = time.monotonic()
    bounds = previous['bounds_ns']
    fits = {float(q): copy.deepcopy(f) for q, f in zip(previous['q_grid_pe'], previous['fits'])}
    targets = np.unique(np.asarray(new_q, float))
    if np.any(targets < 50.-1e-8) or np.any(targets > 50000.+1e-8):
        raise ValueError('New charges outside the fixed domain')
    if not refit_existing: targets = np.array([q for q in targets if q not in fits])
    evaluations = trials = failed = 0
    for direction in (1, -1):
        warm = None
        for q2 in targets[::direction]:
            seeds = observed_starts(row, q2, bounds, previous, warm)
            if q2 in fits: seeds.append(fits[q2])
            nearby = sorted(fits, key=lambda q: abs(np.log(q/q2)))[:2]
            seeds.extend(fits[q] for q in nearby)
            chosen = diverse_starts(response, noise, voltage, q2, seeds, bounds, starts)
            candidate = fit_at_charge(response, noise, voltage, float(q2), chosen, bounds,
                                      maxiter=maxiter, starts=len(chosen))
            evaluations += candidate['search_evaluations'] + len(seeds)
            trials += len(candidate['trial_converged'])
            failed += sum(not x for x in candidate['trial_converged'])
            if q2 not in fits or candidate['cost'] < fits[q2]['cost']:
                fits[float(q2)] = candidate
            warm = fits[float(q2)]
    # Preserve an actual evaluated nuisance slice over all old and new nodes.
    # This can also improve old nodes when a new charge reveals a better basin.
    best = min(fits.values(), key=lambda f: f['cost'])
    slice_cost = []
    for q2 in sorted(fits):
        candidate = dict(times_ns=list(best['times_ns']), charge_pe=[best['charge_pe'][0], q2])
        cost = point_cost(response, noise, voltage, candidate)
        slice_cost.append(cost)
        if cost < fits[q2]['cost']:
            fits[q2] = dict(candidate, cost=cost, converged=False, retained_slice=True)
    q = sorted(fits)
    return dict(bounds_ns=list(bounds), q_grid_pe=q, fits=[fits[x] for x in q],
                reference_cost=min(f['cost'] for f in fits.values()), slice_cost=slice_cost,
                slice_nuisance=dict(times_ns=best['times_ns'], Q1_pe=best['charge_pe'][0]),
                stage=label, optimized_nodes=len(targets), new_nodes=len(q)-len(previous['q_grid_pe']),
                starts_per_node_per_pass=starts, maxiter=maxiter, passes=2,
                evaluations_upper_estimate=evaluations, trials=trials, nonconverged_trials=failed,
                seconds=time.monotonic()-began)


def compare_regions(old, new, level, reference_tolerance=.05, endpoint_fraction=.05):
    def region(curve):
        return compatibility_segments(curve['q_grid_pe'], [f['cost'] for f in curve['fits']],
                                      curve['reference_cost'], level)
    a, b = region(old), region(new)
    topology = len(a['segments']) == len(b['segments'])
    boundary = (a['lower_boundary'], a['upper_boundary']) == (b['lower_boundary'], b['upper_boundary'])
    endpoint_log_change = None
    if topology and a['segments']:
        endpoint_log_change = float(np.max(np.abs(np.log(np.array(b['segments'])/np.array(a['segments'])))))
    reference_gain = float(old['reference_cost']-new['reference_cost'])
    stable = bool(topology and boundary and endpoint_log_change is not None and
                  endpoint_log_change <= np.log1p(endpoint_fraction) and reference_gain <= reference_tolerance)
    return dict(level=float(level), old=a, new=b, same_topology=topology, same_boundaries=boundary,
                endpoint_log_change=endpoint_log_change, reference_gain=reference_gain, stable=stable)


def profile_stability(response, noise, voltage, row, previous, settings):
    original = copy.deepcopy(previous)
    original['stage'] = 'original'
    search = refine_stage(response, noise, voltage, row, original, original['q_grid_pe'],
                          refit_existing=True, label='search', **settings['search'])
    dense = refine_stage(response, noise, voltage, row, search, midpoint_grid(search['q_grid_pe']),
                         refit_existing=False, label='dense', **settings['grid'])
    final = refine_stage(response, noise, voltage, row, dense, contour_midpoints(dense, settings['levels']),
                         refit_existing=False, label='contour', **settings['grid'])
    stages = dict(original=original, search=search, dense=dense, contour=final)
    comparisons = {}
    for before, after in (('original', 'search'), ('search', 'dense'), ('dense', 'contour')):
        comparisons[before+'_'+after] = [compare_regions(stages[before], stages[after], level,
            settings['reference_tolerance'], settings['endpoint_fraction']) for level in settings['levels']]
    return dict(stages=stages, comparisons=comparisons,
                interpretation='Known K2 numerical audit on exposed data; no truth inputs, no confidence calibration')
