"""Empirical single-event calibration of blind residual pulse births.

The score includes the entire fixed candidate search. Calibration samples
contain isolated events; inference receives only candidate fits and frozen
score distributions. This is an empirical gate, not a chi-square likelihood
test, a posterior probability, or a guarantee under detector-model shift.
"""
from __future__ import annotations

import numpy as np

from .alpha_pair_residual import PrecisionBlend, choose, fit_candidates


def waveform_candidates(response, voltage, *, pilot='observed'):
    """Exactly reuse the previous waveform-only search, without a network."""
    raw = [{k: v for k, v in fit.items() if k != 'model_mV'}
           for fit in response.base.candidates(voltage)]
    # The old search also computes two hybrid ablations. Identical observed
    # single starts make these harmless; neither affects waveform_only.
    initial = raw[min(1, len(raw)-1)]
    seeds = dict(candidates=raw, ml_1=initial, ml_2=initial)
    metric = None
    if pilot == 'model' and len(raw) > 1:
        # Fit one smooth mean from the observation, without truth or true K.
        # Freeze its weights for every subsequent zero/one/two-event trial.
        one = response.refine(voltage, initial, radius=24.)
        smooth = response.components(one['times_ns'], one['charge_pe']).sum(0)+response.baseline
        metric = PrecisionBlend(response, smooth, 0., rho=.9)
    elif pilot != 'observed' and pilot != 'model':
        raise ValueError('Unknown precision pilot')
    return fit_candidates(response, voltage, seeds, fraction=0., rho=.9, metric=metric)['waveform_only']


def evidence(candidates, mode='relative_residual'):
    """Gain of the best eligible double over the fitted single, after search."""
    if len(candidates) < 3:
        return 0.
    zero, one, two = candidates
    if any(q < 500. for q in two['charge_pe']):
        return 0.
    gain = max(0., one['cost'] - two['cost'])
    if mode == 'relative_residual':
        return gain / max(one['cost'], 1e-12)
    if mode == 'relative_signal':
        return gain / max(zero['cost'], 1.)
    raise ValueError('Unknown evidence score: ' + mode)


def bin_index(candidates, edges):
    charge = sum(candidates[1]['charge_pe']) if len(candidates) > 1 else 0.
    return int(np.searchsorted(edges, charge, side='right'))


def calibrate(single_candidates, *, edges, mode, alpha=.025, penalty=.001, require_birth_penalty=True):
    """Use observable single-fit Q bins fixed before score calibration.

    The caller must provide independent isolated-event calibration parents.
    Empty bins are rejected, not filled by optimistic interpolated tails.
    """
    if not 0. < alpha < 1. or any(b <= a for a, b in zip(edges, edges[1:])):
        raise ValueError('Invalid alpha or charge-bin edges')
    scores = [[] for _ in range(len(edges)+1)]
    for candidates in single_candidates:
        scores[bin_index(candidates, edges)].append(evidence(candidates, mode))
    if any(not values for values in scores):
        raise ValueError('Every observable charge bin needs calibration data')
    return dict(edges_pe=list(edges), mode=mode, alpha=float(alpha), penalty=float(penalty),
                scores=[sorted(values) for values in scores], truth_inputs_at_inference=[],
                require_birth_penalty=bool(require_birth_penalty))


def select(candidates, calibration):
    """Retain the old zero/single decision; require calibrated birth evidence.

    Rank uses >= for ties and the plus-one finite-sample correction. A value
    near zero is not a probability that the emitted second event is false.
    """
    base = choose(candidates[:2], calibration['penalty'])
    if calibration.get('require_birth_penalty', True):
        original = choose(candidates, calibration['penalty'])
        eligible = len(original['times_ns']) == 2
    else:
        # The model-pilot study replaces the birth penalty by the calibrated
        # rank. Zero-vs-single detection remains unchanged. This can recover
        # a genuine weak birth whose gain is below the old fixed penalty.
        eligible = (len(candidates) == 3 and len(base['times_ns']) == 1
                    and all(q >= 500. for q in candidates[2]['charge_pe']))
        original = candidates[2] if eligible else base
    if not eligible:
        return dict(original, null_rank_p=None, birth_accepted=False)
    index = bin_index(candidates, calibration['edges_pe'])
    values = np.asarray(calibration['scores'][index])
    observed = evidence(candidates, calibration['mode'])
    p = (1. + int(np.count_nonzero(values >= observed))) / (len(values) + 1.)
    accepted = bool(p <= calibration['alpha'])
    result = original if accepted else base
    return dict(result, null_rank_p=float(p), birth_accepted=accepted,
                evidence_score=float(observed), calibration_bin=index)
