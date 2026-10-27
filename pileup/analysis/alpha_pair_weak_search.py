"""Voltage-only extra starts for the frozen, charge-conserving joint NLL."""
from __future__ import annotations

import numpy as np

from .alpha_pair_joint_likelihood import refine, waveform_candidates

GAPS = (4., 8., 12., 20., 40., 80., 160., 320.)
SHARES = {
    "legacy_extra": (.15, .5, .85),
    "weak_extra": (.01, .03, .05, .10, .15, .5, .85, .90, .95, .97, .99),
}


def split_proposals(time, total_charge, maximum_time, shares):
    """Keep total charge, chronological order, and the optimizer's Q domain."""
    if total_charge < 100.:
        return []
    proposals = []
    for direction in (-1, 1):
        for gap in GAPS:
            times = np.array([time, time + direction * gap])
            if times.min() < 0. or times.max() > maximum_time:
                continue
            order = np.argsort(times)
            seen = set()
            for share in shares:
                first = float(np.clip(total_charge * share, 50., total_charge - 50.))
                charges = np.array([first, total_charge - first])[order]
                key = tuple(charges)
                if key in seen:
                    continue
                seen.add(key)
                proposals.append(dict(direction=direction, gap_ns=gap,
                                      times_ns=times[order].tolist(), charge_pe=charges.tolist()))
    return proposals


def extend_candidates(response, noise, voltage, variant, baseline=None):
    """The optional cache contains only candidates computed from this voltage.

    There is no true-count, time, charge, particle or component argument. Both
    variants screen at most eight starts and polish two, retaining the old
    best double candidate. The one/zero models and physics are unchanged.
    """
    if variant not in SHARES:
        raise ValueError(variant)
    old = waveform_candidates(response, noise, voltage) if baseline is None else baseline
    if len(old) < 3:
        return dict(candidates=old, search=dict(proposals=0, screened=0, polished=0, evaluations=0))
    y = np.asarray(voltage, float) - response.baseline - noise.background_mV
    one = old[1]
    bank = split_proposals(one['times_ns'][0], sum(one['charge_pe']),
                           response.sample_grid[-1], SHARES[variant])
    groups = {}
    for proposal in bank:
        curves = response.components(proposal['times_ns'], proposal['charge_pe'])
        cost = noise.score(y - curves.sum(0), curves, proposal['charge_pe'])
        key = (proposal['direction'], proposal['gap_ns'])
        if key not in groups or cost < groups[key][0]:
            groups[key] = (float(cost), proposal)
    starts = []
    for direction in (-1, 1):
        ranked = sorted((v for k, v in groups.items() if k[0] == direction), key=lambda item: item[0])
        starts.extend(proposal for _, proposal in ranked[:4])
    screened = [refine(response, noise, voltage, p, maxiter=20) for p in starts]
    finalists = sorted([old[2]] + screened, key=lambda r: r['cost'])[:2]
    polished = [refine(response, noise, voltage, p, maxiter=70) for p in finalists]
    best = min([old[2]] + screened + polished, key=lambda r: r['cost'])
    return dict(candidates=[old[0], old[1], best], search=dict(
        proposals=len(bank), screened=len(screened), polished=len(polished),
        evaluations=sum(p['evaluations'] for p in screened + polished),
        before_cost=old[2]['cost'], after_cost=best['cost'],
        screened_starts=starts,
        interpretation='Extra optimization cost only; the cached baseline search is additional'))
