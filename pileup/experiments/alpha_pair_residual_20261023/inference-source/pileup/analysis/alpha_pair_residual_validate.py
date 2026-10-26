"""Scientific integrity checks for the frozen residual-search experiment."""
from __future__ import annotations

import ast
from datetime import datetime
import io
import json
import unittest

import numpy as np

from . import alpha_pair_residual_study as study
from .alpha_pair_residual_report import assemble


def main():
    frozen = study.check_freeze()
    manifest = study.load(study.EXPERIMENT / 'manifest.json')
    runs = manifest['runs']
    old_runs = study.load(study.PARENT / 'manifest.json')['runs']
    assert len(runs) == len(manifest['planned_runs']) == 240
    assert len({r['seed'] for r in runs}) == len({r['root'] for r in runs}) == 240
    assert all(r['status'] == 'completed' and r['validation']['runinfo_macro_verified']
               and r['validation']['seed_verified'] for r in runs)
    assert not {r['seed'] for r in runs} & {r['seed'] for r in old_runs}
    assert not {r['root'] for r in runs} & {r['root'] for r in old_runs}
    independence = study.load(study.RECORD / 'parent-independence-audit.json')
    assert independence['passed'] and independence['new_seeds'] == 240
    assert not independence['collisions']
    digitizer_audit = study.load(study.RECORD / 'digitizer-independence-audit.json')
    assert digitizer_audit['passed'] and digitizer_audit['new_base_seeds'] == 280
    assert digitizer_audit['collisions'] == []
    starts = sorted(datetime.fromisoformat(r['started']).timestamp() for r in runs)
    assert starts[0] > datetime.fromisoformat(frozen['frozen_at']).timestamp()
    assert min(np.diff(starts)) >= 1.
    for row in runs:
        assert (study.PROJECT / row['root']).stat().st_size == row['root_bytes']
        assert study.digest(study.PROJECT / row['macro']) == row['macro_sha256']
    runtime = study.load(study.RECORD / 'neural-runtime-preflight.json')
    assert runtime['passed'] and runtime['torch'] == '2.14.1+cpu'
    assert sum(r['frames'] for r in runtime['rows']) == 128
    assert all(all(v == 0. for v in r['max_absolute_difference'].values()) for r in runtime['rows'])
    protocol = study.load(study.EXPERIMENT / 'protocol.json')
    for name, path in (
        ('development_input_sha256', study.EXPERIMENT / 'development.json'),
        ('development_predictions_sha256', study.OUTPUT / 'validation.predictions.jsonl'),
        ('neural_runtime_preflight_sha256', study.RECORD / 'neural-runtime-preflight.json'),
    ):
        assert study.digest(path) == protocol[name], name
    for key in ('test-inference.json', 'test-test_seeds-access-audit.json'):
        audit = study.load(study.RECORD / key)
        assert audit['passed'] and audit['label_access_blocked'] and audit['truth_inputs'] == []
    inference = study.load(study.RECORD / 'test-inference.json')
    assert inference['events'] == inference['requested'] == 2320
    assert inference['prediction_sha256'] == study.digest(study.OUTPUT / 'test.predictions.jsonl')
    assert inference['waveform_sha256'] == study.digest(study.OUTPUT / 'test.waveforms.npz')
    with np.load(study.OUTPUT / 'test.waveforms.npz', allow_pickle=False) as f:
        assert f.files == ['voltage_mV']
        y = f['voltage_mV']
    with np.load(study.OUTPUT / 'test.truth.npz', allow_pickle=False) as f:
        q, t = f['charge_pe'], f['times_ns']
    assert y.shape == (2320, 1024) and np.isfinite(y).all()
    assert q.shape == t.shape == (2320, 2) and np.isfinite(q).all() and np.isfinite(t).all()
    methods, cases = assemble('test')
    assert len(cases) == 2320 and len(methods) == 14
    for i, case in enumerate(cases):
        target = {'empty': 0, 'single': 1, 'pair': 2}[case['kind']]
        assert int(np.sum(q[i] >= 500.)) == target
        if target == 2:
            assert abs(float(t[i, 1] - t[i, 0]) - case['gap_ns']) < 1e-8
    for fits in methods.values():
        assert len(fits) == 2320
        for fit in fits:
            tt, qq = np.asarray(fit['times_ns']), np.asarray(fit['charge_pe'])
            assert len(tt) == len(qq) <= 2
            assert np.isfinite(tt).all() and np.isfinite(qq).all()
            assert (qq >= 0).all() and (np.diff(tt) >= 0).all()
    selection = study.load(study.EXPERIMENT / 'residual-selection.json')
    assert selection['chosen']['eligible']
    assert selection['chosen_policy'] == protocol['selected_policy']
    report = study.load(study.OUTPUT / 'test-report.json')
    assert report['event_rows'] == 2320 * 14 and len(report['summary']) == 54 * 14
    for row in report['summary']:
        assert row['frames'] == (40 if row['kind'] == 'pair' else 80)
        assert 0 <= row['success'] <= row['frames']
        assert 0 <= row['strict_recovered'] <= row['target']
        assert 0 <= row['recovered'] <= row['target']
        assert row['budget_pass'] == (row['extra_fraction'] <= .05)
    physical = study.load(study.RECORD / 'physical-validation.json')
    assert physical['frames'] == 2320 and physical['max_closure_error_mV'] < 1e-8

    # A predetermined fresh pair is replayed with the unoptimized scalar code.
    # This check never feeds truth or results back into reconstruction.
    from . import alpha_pair_resolution_data as data, sipm_response
    assert sipm_response.simulate_avalanches.__module__ == 'pileup.analysis.sipm_response'
    checkpoint = study.cache() / 'physical-cases/00000.npz'
    if checkpoint.exists():
        with np.load(checkpoint) as f:
            case = json.loads(str(f['case']))
            expected = [f[k] for k in ('y', 'q', 't', 'components')]
        precision = 'float64 checkpoint'
    else:
        case = cases[0]
        with np.load(study.OUTPUT / 'test.truth.npz') as f:
            components = f['components_mV'][0]
        expected = [y[0], q[0], t[0], components]
        precision = 'persistent float32 voltage/components; float64 Q/time'
    assert case['parents'] == [0, 1] and case['gap_ns'] == 4.
    data._DONORS = []
    for i in (0, 1):
        with np.load(study.OUTPUT / 'photons' / f'{i:04}.npz') as f:
            data._DONORS.append({k: f[k] for k in f.files})
    records = study.load(study.OUTPUT / 'native-records.json')['test_single']
    data._CONFIG = study.load(study.PROJECT / records[0]['metadata'])['config']
    replay = data.joint_case(case)
    for a, b in zip(expected, replay[:4]):
        np.testing.assert_array_equal(a, np.asarray(b, dtype=a.dtype))

    # Calibration and development parents precede the new independent test.
    development = study.load(study.EXPERIMENT / 'development.json')
    assert development['previous_test_now_development'] and development['validation_frames'] == 1320
    calibration_runs = study.load(study.PROJECT / development['calibration_parent'] / 'manifest.json')['runs']
    assert not {r['root'] for r in calibration_runs} & {r['root'] for r in old_runs + runs}
    old_frozen = study.load(study.PARENT / 'freeze.json')
    for name, sha in old_frozen['files'].items():
        assert study.digest(study.PARENT / name) == sha
    for path in (study.EXPERIMENT / 'inference-source/pileup/analysis').glob('*alpha_pair_residual*.py'):
        ast.parse(path.read_text(), filename=str(path))
    suite = unittest.defaultTestLoader.loadTestsFromNames([
        'pileup.analysis.test_alpha_pair_conditional', 'pileup.analysis.test_alpha_pair_residual'])
    log = io.StringIO()
    result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)
    (study.RECORD / 'targeted-tests.txt').write_text(log.getvalue())
    assert result.wasSuccessful(), log.getvalue()
    answer = dict(passed=True, roots=240, frames=2320, methods=14, score_rows=2320*14,
        freeze_unchanged=True, previous_freeze_unchanged=True, production_after_freeze=True,
        parent_seeds_disjoint=True, digitizer_seeds_disjoint_from_old_cases=True,
        previous_manifests_checked=len(independence['manifests']),
        restored_neural_replay_exact_frames=128, runinfo_validated_during_production=True,
        root_sizes_unchanged=True, observation_only_files=True, inference_label_guards_passed=True,
        component_closure_max_mV=physical['max_closure_error_mV'],
        fresh_scalar_digitizer_replay_exact=True, scalar_replay_comparison_precision=precision,
        targeted_tests=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        report_sha256=study.digest(study.OUTPUT / 'test-report.json'),
        scope='Pure alpha, fixed double lifetime, local K<=2; not energy calibration or rate limit')
    study.save(study.RECORD / 'science-validation.json', answer)
    print(json.dumps(answer, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
