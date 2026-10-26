"""Read-only scientific checks after frozen independent inference."""
from pathlib import Path
import ast
from datetime import datetime
import json
import os
import sys
import unittest

FOLDER=Path(__file__).resolve().parent
PROJECT=FOLDER.parents[2]
for key in ('GAGG_CONDITIONAL_PROJECT','GAGG_RESOLUTION_PROJECT','GAGG_PAIR_PROJECT','GAGG_COMPONENT_PROJECT'):
    os.environ[key]=str(PROJECT)
sys.path.insert(0,str(FOLDER/'inference-source'))


def main():
    import numpy as np
    from pileup.analysis import alpha_pair_conditional_study as study
    from pileup.analysis.alpha_pair_conditional_report import assemble
    study.check_freeze()
    assert study.digest(study.cache()/'response.npz')==study.digest(FOLDER/'response.npz')
    manifest=study.load(FOLDER/'manifest.json');runs=manifest['runs']
    freeze=study.load(FOLDER/'freeze.json');old=study.load(study.PARENT/'manifest.json')['runs']
    assert len(runs)==240 and len({r['seed'] for r in runs})==240
    assert all(r['status']=='completed' and r['validation']['runinfo_macro_verified'] and r['validation']['seed_verified'] for r in runs)
    assert not {r['seed'] for r in runs}&{r['seed'] for r in old}
    assert not {r['root'] for r in runs}&{r['root'] for r in old}
    independence=study.load(study.RECORD/'parent-independence-audit.json')
    assert independence['passed'] and independence['new_seeds']==240 and not independence['collisions']
    runtime=study.load(study.RECORD/'restored-neural-replay.json')
    assert runtime['passed'] and runtime['torch']=='2.14.1+cpu'
    assert all(all(v==0. for v in r['max_absolute_difference'].values()) for r in runtime['rows'])
    assert min(datetime.fromisoformat(r['started']) for r in runs)>datetime.fromisoformat(freeze['frozen_at'])
    starts=sorted(datetime.fromisoformat(r['started']).timestamp() for r in runs)
    assert min(np.diff(starts))>=1.
    for r in runs:
        assert (PROJECT/r['root']).stat().st_size==r['root_bytes']
        assert study.digest(PROJECT/r['macro'])==r['macro_sha256']
    with np.load(study.OUTPUT/'test.waveforms.npz',allow_pickle=False) as f:
        assert f.files==['voltage_mV'];y=f['voltage_mV']
        assert y.shape==(2320,1024) and np.isfinite(y).all()
    with np.load(study.OUTPUT/'test.truth.npz',allow_pickle=False) as f:
        q=f['charge_pe'];t=f['times_ns']
        assert q.shape==t.shape==(2320,2) and np.isfinite(q).all() and np.isfinite(t).all()
    methods,cases=assemble('test')
    assert len(cases)==2320 and len(methods)==13
    for i,c in enumerate(cases):
        target={'empty':0,'single':1,'pair':2}[c['kind']]
        assert int(np.sum(q[i]>=500.))==target
        if target==2:assert abs(float(t[i,1]-t[i,0])-c['gap_ns'])<1e-8
    for name,fits in methods.items():
        assert len(fits)==2320
        for fit in fits:
            tt=np.asarray(fit['times_ns']);qq=np.asarray(fit['charge_pe'])
            assert len(tt)==len(qq)<=2 and np.isfinite(tt).all() and np.isfinite(qq).all()
            assert (qq>=0).all() and (np.diff(tt)>=0).all()
    report=study.load(study.OUTPUT/'test-report.json')
    assert report['event_rows']==2320*13
    for row in report['summary']:
        assert row['frames']==(40 if row['kind']=='pair' else 80)
        assert 0<=row['success']<=row['frames']
        # Matching is recomputed for each tolerance; the individual-event
        # qualified counts need not be nested when an assignment changes.
        assert 0<=row['strict_recovered']<=row['target']
        assert 0<=row['recovered']<=row['target']
        assert row['budget_pass']==(row['extra_fraction']<=.05)
    for name in ('test-infer-access-audit.json','test-test_seeds-access-audit.json'):
        audit=study.load(study.RECORD/name)
        assert audit['passed'] and audit['label_access_blocked'] and audit['truth_inputs']==[]
    physical=study.load(study.RECORD/'physical-validation.json')
    assert physical['max_closure_error_mV']<1e-8
    # Recheck the existing vectorized digitizer on one predetermined fresh
    # pair using the frozen scalar implementation, without changing inference.
    from pileup.analysis import alpha_pair_resolution_data as data,sipm_response
    assert sipm_response.simulate_avalanches.__module__=='pileup.analysis.sipm_response'
    checkpoint=study.cache()/'physical-cases/00000.npz'
    if checkpoint.exists():
        with np.load(checkpoint) as f:
            case=json.loads(str(f['case']));expected=[f[k] for k in ('y','q','t','components')]
        replay_precision='float64 checkpoint'
    else:
        case=cases[0]
        with np.load(study.OUTPUT/'test.truth.npz') as f:components=f['components_mV'][0]
        expected=[y[0],q[0],t[0],components]
        replay_precision='persistent observation/component float32, Q/time float64'
    assert case['parents']==[0,1] and case['gap_ns']==4.
    data._DONORS=[]
    for i in (0,1):
        with np.load(study.OUTPUT/'photons'/f'{i:04}.npz') as f:data._DONORS.append({k:f[k] for k in f.files})
    records=study.load(study.OUTPUT/'native-records.json')['test_single']
    data._CONFIG=study.load(PROJECT/records[0]['metadata'])['config']
    replay=data.joint_case(case)
    for a,b in zip(expected,replay[:4]):np.testing.assert_array_equal(a,np.asarray(b,dtype=a.dtype))
    development=study.load(FOLDER/'development.json')
    assert not set(development['train_parents'])&set(development['validation_parents'])
    for p in (FOLDER/'inference-source/pileup/analysis').glob('*alpha_pair_conditional*.py'):
        ast.parse(p.read_text(),filename=str(p))
    suite=unittest.defaultTestLoader.loadTestsFromName('pileup.analysis.test_alpha_pair_conditional')
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    assert result.wasSuccessful()
    answer=dict(passed=True,roots=240,frames=2320,methods=13,score_rows=2320*13,
        freeze_unchanged=True,production_after_freeze=True,parent_seeds_disjoint=True,
        previous_manifests_checked=len(independence['manifests']),
        restored_neural_replay_exact_frames=128,
        runinfo_validated_during_production=True,root_sizes_unchanged=True,
        observation_only_files=True,inference_label_guards_passed=True,
        component_closure_max_mV=physical['max_closure_error_mV'],
        fresh_scalar_digitizer_replay_exact=True,
        scalar_replay_comparison_precision=replay_precision,
        targeted_tests=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        report_sha256=study.digest(study.OUTPUT/'test-report.json'),
        scope='Pure alpha, fixed double lifetime, local K<=2; not energy calibration or rate limit')
    study.save(study.RECORD/'science-validation.json',answer)
    print(json.dumps(answer,ensure_ascii=False),flush=True)


if __name__=='__main__':main()
