"""Reuse observation-only candidates while separating count decision from search."""
from __future__ import annotations

from copy import deepcopy

import numpy as np

from .alpha_pair_joint_likelihood import calibrate_gate,select,evidence
from .alpha_pair_risk_learning import calibrated_gate,select_candidate

METHODS=('baseline','learned_gate','local_physical','mlp_physical','mlp_fixed_gate')
PRIMARY=METHODS[:-1]


def compact(candidate,*,learned_score=False):
    fields=('times_ns','charge_pe','cost','quadratic','logdet','converged','evaluations','message')
    result={k:deepcopy(candidate[k]) for k in fields if k in candidate}
    if learned_score:result['network_score']=float(candidate['network_score'])
    return result


def assemble(baseline,hybrid):
    """The local control uses trial 0; the hybrid includes both previous starts.

    Do not choose a different minimum merely because it passes Q>=500. That
    would change the original physical decision rule being controlled here.
    """
    base=deepcopy(baseline)
    if len(base)<3:
        return dict(baseline=base,local=deepcopy(base),mlp=deepcopy(base)),compact(hybrid,learned_score=True)
    assert len(base)==3 and len(hybrid['trials'])==2
    local=min((base[2],hybrid['trials'][0]),key=lambda c:c['cost'])
    assert hybrid['penalty']==0. and hybrid['cost']<=local['cost']+1e-7
    sets=dict(baseline=base,local=deepcopy(base[:2])+[compact(local)],mlp=deepcopy(base[:2])+[compact(hybrid)])
    assert all(group[:2]==base[:2] for group in sets.values())
    assert sets['local'][2]['cost']<=base[2]['cost']+1e-7
    return sets,compact(hybrid,learned_score=True)


def gates_from_controls(controls,baseline_gate,previous_learned_gate):
    assert len(controls)==280
    original=calibrate_gate([r['sets']['baseline'] for r in controls])
    assert original==baseline_gate
    learned=calibrated_gate([r['learned'] for r in controls],
        [sum(r['sets']['baseline'][1]['charge_pe']) for r in controls])
    assert learned==previous_learned_gate
    return dict(baseline=original,learned_gate=learned,
        local_physical=calibrate_gate([r['sets']['local'] for r in controls]),
        mlp_physical=calibrate_gate([r['sets']['mlp'] for r in controls]),mlp_fixed_gate=deepcopy(original))


def bank_for(method):
    return 'baseline' if method=='baseline' else 'local' if method=='local_physical' else 'mlp'


def decide(row,method,gate):
    base=row['sets']['baseline'];bank=row['sets'][bank_for(method)]
    fit=select_candidate(base,row['learned'],gate) if method=='learned_gate' else select(bank,gate)
    if len(fit['times_ns'])<2:
        original=base[len(fit['times_ns'])]
        assert fit['times_ns']==original['times_ns'] and fit['charge_pe']==original['charge_pe']
    return fit


def rejection(row,method,fit):
    base=row['sets']['baseline'];bank=row['sets'][bank_for(method)]
    if len(base)<2 or base[1]['charge_pe'][0]<500. or base[1]['cost']>=base[0]['cost']:return 'no_presence'
    if len(bank)<3:return 'no_pair'
    if min(bank[2]['charge_pe'])<500.:return 'subthreshold_Q'
    if method!='learned_gate' and evidence(bank)<=0.:return 'nonpositive_evidence'
    return 'accepted' if len(fit['times_ns'])==2 else 'rank_rejected'


def replay_cost(response,noise,voltage,candidate):
    """Use float64 before offset subtraction, like the frozen physical fitter."""
    y=np.asarray(voltage,dtype=float)-response.baseline-noise.background_mV
    curves=response.components(candidate['times_ns'],candidate['charge_pe'])
    return float(noise.score(y-curves.sum(0),curves,candidate['charge_pe']))
