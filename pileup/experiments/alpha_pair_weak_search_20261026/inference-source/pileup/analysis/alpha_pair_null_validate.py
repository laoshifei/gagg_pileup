"""Validate provenance, isolation, scoring and fixed-example delivery inputs."""
from __future__ import annotations

from datetime import datetime
import gzip
import json
from pathlib import Path
import re
import unittest

import numpy as np

from . import alpha_pair_null_study as study
from .alpha_pair_null_report import assemble, summarize


def validate():
    study.check_freeze()
    frozen=study.load(study.EXPERIMENT/'freeze.json')
    freeze_time=datetime.fromisoformat(frozen['frozen_at'])
    old_freezes=[]
    for name in ('alpha_pair_conditional_20261022','alpha_pair_residual_20261023'):
        folder=study.PROJECT/'pileup/experiments'/name
        old=study.load(folder/'freeze.json')
        for file, sha in old['files'].items(): assert study.digest(folder/file)==sha,(name,file)
        old_freezes.append(dict(experiment=name,files=len(old['files']),unchanged=True))
    roots=[]; seeds=[]; stage_counts={}
    for stage in ('pilot','pilot-refinement','development','test'):
        manifest=study.load(study.EXPERIMENT/stage/'manifest.json')
        plans={r['label']:r for r in manifest['planned_runs']}
        assert len(manifest['runs'])==len(plans)
        assert len({r['label'] for r in manifest['runs']})==len(plans)
        stage_counts[stage]=len(plans)
        for r in manifest['runs']:
            assert r['status']=='completed' and r['seed']==plans[r['label']]['seed']
            root=study.PROJECT/r['root']; macro=study.PROJECT/r['macro']
            assert root.parent==study.PROJECT/'g4/data'
            assert re.fullmatch(r'\d{8}_\d{2}h\d{2}m\d{2}s\.root',root.name)
            assert study.digest(root)==r['root_sha256'],r['root']
            assert study.digest(macro)==r['macro_sha256'],r['macro']
            meta=study.load(study.OUTPUT/stage/'native-inputs/metadata'/(root.stem+'.json'))
            assert meta['root_sha256']==r['root_sha256'] and meta['seed']==r['seed']
            assert meta['macro_contents']==macro.read_text()
            assert meta['primaries']==1 and meta['duration_ns']==20000.
            assert r['validation']['runinfo_macro_verified'] and r['validation']['seed_verified']
            if stage=='test': assert datetime.fromisoformat(r['started'])>freeze_time
            roots.append(r['root']);seeds.append(r['seed'])
    assert len(roots)==len(set(roots))==len(seeds)==len(set(seeds))==416
    assert study.load(study.RECORD/'test/digitizer-independence.json')['passed']
    stage='test';output=study.OUTPUT/stage
    with np.load(output/'observations.npz') as f:
        assert f.files==['voltage_mV']
        assert f['voltage_mV'].shape==(3200,1024)
        assert np.isfinite(f['voltage_mV']).all()
    for variant,filename in [('observed','inference.json'),('model','model-inference.json')]:
        info=study.load(study.RECORD/stage/filename)
        assert info['frames']==info['requested']==3200 and info['truth_inputs']==[]
        assert info['worker_label_access_blocked'] and info['waveform_sha256']==study.digest(output/'observations.npz')
        predfile='predictions.jsonl' if variant=='observed' else 'predictions.model.jsonl'
        assert info['prediction_sha256']==study.digest(output/predfile)
    for stage in ('development','test'):
        pre=study.load(study.RECORD/stage/'physical-preflight.json')
        physical=study.load(study.RECORD/stage/'physical-validation.json')
        assert pre['passed'] and pre['scalar_digitizer_exact']
        assert physical['max_closure_error_mV']<1e-8 and not physical['photon_scaling']
    result=study.load(output/'report.json')
    assert result['primary']=='model_'+study.load(study.EXPERIMENT/'model-null-selection.json')['chosen']
    with gzip.open(output/'scores.jsonl.gz','rt') as stream:
        scores=[json.loads(line) for line in stream]
    assert len(scores)==16000 and len(result['summary'])==375
    groups={}
    for r in scores:
        assert 0<=r['strict_recovered']<=r['target']<=2
        assert 0<=r['recovered']<=r['matched']<=r['target']
        assert r['missed']+r['matched']==r['target']
        assert r['extra']+r['matched']==r['count']
        assert r['count']<=2
        assert r['known_count']==(r['method']=='known_count_model')
        groups.setdefault((r['method'],r['kind'],r['condition'],r['gap_ns']),[]).append(r)
    for row in result['summary']:
        key=(row['method'],row['kind'],row['condition'],row['gap_ns'])
        assert summarize(groups[key])==row,key
    rows,fits=assemble('test')
    examples=study.load(output/'examples.json')
    spec=study.load(study.EXPERIMENT/'test-plan.json')['examples_predeclared']
    for example,expected in zip(examples,spec):
        assert all(example['case'][key]==value for key,value in expected.items())
        i=example['index']
        assert example['new']==fits[result['primary']][i]
        assert example['old']==fits['observed_penalty'][i]
        assert example['true_q']==rows[i]['charge'].tolist()
        assert example['true_t']==rows[i]['times'].tolist()
    assert len(examples)==len(spec)==2
    images={p.name:study.digest(p) for p in (study.OUTPUT/'figures').glob('*.png')}
    assert set(images)=={s+'.png' for s in ('0004-49-1','0004-50-1','0004-50-2','0004-50-3','0004-51-1')}
    from PIL import Image
    for name in images:
        with Image.open(study.OUTPUT/'figures'/name) as image: image.verify()
    # Check the audit guard itself without depending on test labels.
    from .alpha_pair_conditional_study import forbid_labels
    blocked=[]
    for name in ('a.root','a.truth.npz','a.cases.json','/truth/a.npz','native-records.json'):
        try: forbid_labels('open',(name,'r',0))
        except PermissionError: blocked.append(name)
        else: raise AssertionError(name)
    suite=unittest.defaultTestLoader.loadTestsFromNames([
        'pileup.analysis.test_alpha_pair_null_gate','pileup.analysis.test_alpha_pair_residual'])
    with (study.RECORD/'targeted-tests.log').open('w') as stream:
        tests=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    assert tests.wasSuccessful()
    report=dict(passed=True,stage_counts=stage_counts,root_hashes_verified=len(roots),
        frozen_files=len(frozen['files']),old_freezes=old_freezes,independent_test_after_freeze=True,
        frames=3200,score_rows=len(scores),summary_rows=len(result['summary']),
        scalar_digitizer_checks=['development','test'],label_guard_cases=blocked,
        targeted_tests=tests.testsRun,figures=images,report_sha256=study.digest(output/'report.json'),
        independent_control_requirements_passed=result['all_control_requirements_passed'],
        distinction='Validation of implementation/provenance is separate from performance acceptance')
    study.save(study.RECORD/'science-validation.json',report)
    print('SCIENCE_VALIDATION_PASSED',report,flush=True)


if __name__=='__main__':validate()
