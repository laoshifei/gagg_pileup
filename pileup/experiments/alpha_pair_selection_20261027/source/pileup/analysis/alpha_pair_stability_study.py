"""Frozen numerical stability experiment, separate from statistical coverage."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import importlib.metadata
import multiprocessing
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import unittest
import json
import numpy as np

from .alpha_pair_information_study import load, save, digest, rows
from .alpha_pair_evidence_study import paths
from .alpha_pair_profile_stability import profile_stability, observed_starts, diverse_starts, point_cost
from .alpha_pair_evidence_profile import fit_at_charge

PROJECT=Path(os.environ.get('GAGG_STABILITY_PROJECT',Path(__file__).resolve().parents[2])).resolve()
NAME='alpha_pair_profile_stability_20261027'
EXPERIMENT=PROJECT/'pileup/experiments'/NAME
OUTPUT=PROJECT/'jupyter/artifacts'/NAME
RECORD=PROJECT/'jupyter/run_records'/NAME
PREVIOUS=PROJECT/'pileup/experiments/alpha_pair_evidence_profile_20261027'
INPUT=PROJECT/'jupyter/artifacts/alpha_pair_evidence_profile_20261027'
OLD_RECORD=PROJECT/'jupyter/run_records/alpha_pair_evidence_profile_20261027'
NOTEBOOK=PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
SETTINGS=dict(search=dict(starts=6,maxiter=150),grid=dict(starts=4,maxiter=150),
              levels=[1.,4.,9.],reference_tolerance=.05,endpoint_fraction=.05)


def copy_chunks(source,target):
    with Path(source).open('rb') as r,Path(target).open('xb') as w:
        while data:=r.read(256*1024): assert w.write(data)==len(data)
    assert digest(source)==digest(target)


def prepare(resume=False):
    prior=load(OLD_RECORD/'final-validation.json');assert prior['passed']
    assert digest(NOTEBOOK)==prior['notebook_sha256']
    for name,value in prior['source_sha256'].items():assert digest(PROJECT/name)==value,name
    for name,value in prior['documentation_sha256'].items():assert digest(PROJECT/name)==value,name
    if resume:
        assert all(p.is_dir() for p in (EXPERIMENT,OUTPUT,RECORD))
        assert not (EXPERIMENT/'freeze.json').exists() and not (EXPERIMENT/'inference-source').exists()
    else:
        for folder in (EXPERIMENT,OUTPUT,RECORD):folder.mkdir(parents=True,exist_ok=False)
    for src,name in ((NOTEBOOK,'notebook.before.ipynb'),(PROJECT/'README.md','README.before.md'),
                     (PROJECT/'docs/literature-use.md','literature.before.md')):
        if resume:assert digest(src)==digest(RECORD/name)
        else:copy_chunks(src,RECORD/name)
    all_tasks=load(PREVIOUS/'tasks.json');assert len(all_tasks)==160
    tasks=[all_tasks[i] for i in range(0,160,5)]
    assert len(tasks)==32
    protocol=dict(created=datetime.now(timezone.utc).isoformat(),question='Are numerical conditional-Q2 contours stable to more starts and nested grids?',
        scope='Known K2 pure-alpha diagnostic; fixed95/351ns and original covariance+logdet; ample dynamic range,1024x2ns waveform, no blind refit changes',
        data='Systematic exposed development subset, eight of each original40-frame group; no new independent samples or random draws',
        selection='Take positions0,5,10,15,20,25,30,35 in each of four prespecified groups, before new calculation; no result-based selection',
        groups=['native40ns','native200ns','Gaussian24-64ns','Gaussian128-256ns'],tasks=32,
        stages=['original frozen curve','same grid, six diverse observed starts,150 iterations,forward/reverse',
                'insert geometric midpoint of every cell, four starts,150 iterations,forward/reverse',
                'insert midpoints at every1/4/9 crossing and cells adjacent to minimum, four starts,150 iterations,forward/reverse'],
        settings=SETTINGS,time_domain='Original observed K1 trigger[-96,+512]ns clipped to[0,2046], ordered times; unchanged',
        charge_domain_pe=[50.,50000.],seed_rule='Cost-ranked seeds separated in(t1/8ns,t2/8ns,lnQ1/.25); observed cached/previous fits, trigger offsets-8/0/+8ns and fixed gaps0..512ns; no true values',
        envelope='Retain old points and new stage best-nuisance slice over all nodes; grid stages can improve old points through this slice. Not an equal-compute comparison.',
        numerical_rule='Last dense->contour step stable only if reference gain<=.05, every matched interval endpoint changes by factor<=1.05, same number of components and same boundary flags, separately at1/4/9. All three needed per frame.',
        interpretation='Finite stability screen, not global-optimum proof or confidence-coverage calibration; small systematic subset cannot certify unexamined frames.',
        probe='Only after inference hash frozen, a separate scoring probe fixes true laterQ2 and refits nuisance. No truth times/Q1 or feedback.',
        pilot=[tasks[i] for i in (0,8,16,24)],examples=[tasks[i] for i in (0,8,16,24)],workers_limit=8,blas_threads=1,
        defaults='All blind candidates,gates,7040 original decisions and prior160 profiles remain unchanged.',
        next='If unresolved numerical changes remain, report them before attempting independent coverage calibration; no confidence-level claims from this audit.',
        figures=['0004-103-1','0004-104-1','0004-104-2','0004-105-1','0004-106-1'])
    if resume:
        old=load(EXPERIMENT/'protocol.json');protocol['created']=old['created'];assert protocol==old
        assert tasks==load(EXPERIMENT/'tasks.json')
    else:
        save(EXPERIMENT/'protocol.json',protocol);save(EXPERIMENT/'tasks.json',tasks)
    for name in ('joint-response.npz','joint_conserved_nll-noise.npz'):
        if resume:assert digest(PREVIOUS/name)==digest(EXPERIMENT/name)
        else:copy_chunks(PREVIOUS/name,EXPERIMENT/name)
    provenance=load(OLD_RECORD/'input-provenance.json')
    inputs=dict(provenance['inputs'])
    for name in ('profiles.jsonl','truth-probes.jsonl','report.json','profile-scores.json','decision-records.jsonl'):
        inputs[str((INPUT/name).relative_to(PROJECT))]=digest(INPUT/name)
    save(RECORD/'input-provenance.json',dict(inputs=inputs,prior_validation=prior,
        prior_freezes={name:digest(PREVIOUS/name) for name in ('freeze.json','report-freeze.json','delivery-freeze.json')},
        inherited_provenance_sha256=digest(OLD_RECORD/'input-provenance.json'),mother_groups=provenance['mother_groups'],
        notebook_sha256=digest(NOTEBOOK),new_ROOTs=0,new_training=0,new_random_draws=0))
    try:
        tracked_diff=subprocess.check_output(['git','-c','core.fileMode=false','diff','--stat','--','README.md','docs/literature-use.md'],
                                             cwd=PROJECT,text=True,timeout=60)
    except subprocess.TimeoutExpired:
        tracked_diff='Scoped Git stat exceeded60s on shared storage; exact prior documentation backups and frozen sources record content.'
    save(RECORD/'environment.json',dict(python=sys.version,executable=sys.executable,
        dependencies={k:importlib.metadata.version(k) for k in ('numpy','scipy','matplotlib','pandas','nbformat')},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        tracked_diff=tracked_diff,git_stat_scope=['README.md','docs/literature-use.md'],
        uncommitted='Earlier information/allocation/risk/gate/evidence studies retained; new numerical audit in separate experiment; no commit or push',
        workers_limit=8,blas_threads=1,seed=None,reason_no_seed='Deterministic reanalysis; input seeds inherited through provenance'))
    shutil.copytree(PREVIOUS/'inference-source',EXPERIMENT/'inference-source')
    for name in ('alpha_pair_profile_stability.py','alpha_pair_stability_study.py','test_alpha_pair_profile_stability.py'):
        copy_chunks(PROJECT/'pileup/analysis'/name,EXPERIMENT/'inference-source/pileup/analysis'/name)
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for key in ('STABILITY','EVIDENCE','GATE','RISK','INFORMATION'):os.environ['GAGG_'+key+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_stability_study import main
    main()
''')
    files={str(p.relative_to(EXPERIMENT)):digest(p) for p in (EXPERIMENT/'inference-source').rglob('*.py')}
    for name in ('run.py','tasks.json','protocol.json','joint-response.npz','joint_conserved_nll-noise.npz'):files[name]=digest(EXPERIMENT/name)
    save(EXPERIMENT/'freeze.json',dict(files=files,created=protocol['created']))
    print('STABILITY_PREPARED',len(files),len(tasks),flush=True)


def check_freeze():
    for name,value in load(EXPERIMENT/'freeze.json')['files'].items():assert digest(EXPERIMENT/name)==value,name


def models():
    from .alpha_pair_joint_likelihood import ConservingResponse,JointNoise
    return ConservingResponse.load(EXPERIMENT/'joint-response.npz'),JointNoise.load(EXPERIMENT/'joint_conserved_nll-noise.npz')


def initialize(load_new=False):
    global RESPONSE,NOISE,OBS,CACHE,OLD,NEW,CONFIG
    from .alpha_pair_charge_allocation import forbid_labels
    RESPONSE,NOISE=models();sys.addaudithook(forbid_labels)
    OBS={};CACHE={}
    for domain in ('native','continuous'):
        with np.load(paths(domain)['observation']) as f:OBS[domain]=f['voltage_mV']
        CACHE[domain]={r['index']:r for r in rows(paths(domain)['candidates'])}
    OLD={(r['domain'],r['index']):r for r in rows(INPUT/'profiles.jsonl')}
    NEW={(r['domain'],r['index']):r for r in rows(OUTPUT/'profiles.jsonl')} if load_new else {}
    CONFIG=load(EXPERIMENT/'protocol.json')['settings']


def run_task(task):
    domain,i=task['domain'],task['index'];began=time.monotonic()
    result=profile_stability(RESPONSE,NOISE,OBS[domain][i],CACHE[domain][i],OLD[domain,i],CONFIG)
    error=0.;last=None
    for name,curve in result['stages'].items():
        assert curve['reference_cost']==min(f['cost'] for f in curve['fits'])
        assert np.all(np.diff(curve['q_grid_pe'])>0)
        current=dict(zip(curve['q_grid_pe'],curve['fits']))
        if last is not None:
            for q,f in last.items():assert current[q]['cost']<=f['cost']+1e-9
        for f in curve['fits']:
            actual=point_cost(RESPONSE,NOISE,OBS[domain][i],f)
            error=max(error,abs(actual-f['cost'])/max(1.,abs(actual)))
            assert curve['bounds_ns'][0]-1e-8<=f['times_ns'][0]<=f['times_ns'][1]<=curve['bounds_ns'][1]+1e-8
        last=current
    assert error<1e-9,error
    return dict(result,**task,seconds=time.monotonic()-began,cost_replay_error=error)


def run_probe(task):
    domain,i,q=task['domain'],task['index'],task['true_later_Q_pe'];curve=NEW[domain,i]['stages']['contour']
    seeds=observed_starts(CACHE[domain][i],q,curve['bounds_ns'],curve)
    seeds.extend(min(c['fits'],key=lambda f:f['cost']) for c in NEW[domain,i]['stages'].values())
    chosen=diverse_starts(RESPONSE,NOISE,OBS[domain][i],q,seeds,curve['bounds_ns'],6)
    began=time.monotonic();fit=fit_at_charge(RESPONSE,NOISE,OBS[domain][i],q,chosen,curve['bounds_ns'],maxiter=150,starts=len(chosen))
    return dict(**task,fit=fit,delta_cost=fit['cost']-curve['reference_cost'],seconds=time.monotonic()-began,
        truth_assistance='True laterQ2 only, knownK2 scoring probe; never fed back to profiles or gates')


def tests():
    check_freeze();start=time.monotonic()
    suite=unittest.defaultTestLoader.loadTestsFromName('pileup.analysis.test_alpha_pair_profile_stability')
    result=unittest.TextTestRunner(verbosity=2).run(suite);assert result.wasSuccessful()
    save(RECORD/'numerical-tests.json',dict(passed=True,tests=result.testsRun,seconds=time.monotonic()-start,
        source_sha256=digest(EXPERIMENT/'inference-source/pileup/analysis/test_alpha_pair_profile_stability.py'),
        core_sha256=digest(EXPERIMENT/'inference-source/pileup/analysis/alpha_pair_profile_stability.py')))


def calculate(command,workers):
    check_freeze();assert 1<=workers<=8
    for name,value in load(RECORD/'input-provenance.json')['inputs'].items():assert digest(PROJECT/name)==value,name
    assert load(RECORD/'numerical-tests.json')['passed']
    tasks=load(EXPERIMENT/'protocol.json')['pilot'] if command=='pilot' else load(EXPERIMENT/'tasks.json')
    prior=[];probe=command=='probe'
    if command!='pilot':assert load(RECORD/'pilot-validation.json')['passed']
    if command=='infer':
        prior=rows(OUTPUT/'pilot-profiles.jsonl');done={(r['domain'],r['index']) for r in prior}
        tasks=[t for t in tasks if (t['domain'],t['index']) not in done]
    if probe:
        assert load(RECORD/'inference.json')['passed']
        assert digest(OUTPUT/'profiles.jsonl')==load(RECORD/'inference.json')['sha256']
        truth={}
        for domain in ('native','continuous'):
            with np.load(paths(domain)['truth']) as f:truth[domain]=f['charge_pe']
        tasks=[dict(t,true_later_Q_pe=float(truth[t['domain']][t['index'],1])) for t in tasks]
    target=OUTPUT/('truth-probes.jsonl' if probe else 'pilot-profiles.jsonl' if command=='pilot' else 'profiles.jsonl')
    start=time.monotonic();results=list(prior)
    with target.open('x') as stream:
        for row in prior:stream.write(json.dumps(row,allow_nan=False)+'\n')
        with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn'),initializer=initialize,initargs=(probe,)) as pool:
            futures=[pool.submit(run_probe if probe else run_task,t) for t in tasks]
            for future in as_completed(futures):
                row=future.result();results.append(row);stream.write(json.dumps(row,allow_nan=False)+'\n');stream.flush()
                print('STABILITY_'+command.upper(),len(results),len(tasks)+len(prior),row['domain'],row['index'],round(row['seconds'],2),flush=True)
    assert len(results)==(4 if command=='pilot' else 32)
    name='probe-validation' if probe else 'pilot-validation' if command=='pilot' else 'inference'
    save(RECORD/(name+'.json'),dict(passed=True,frames=len(results),workers=workers,seconds=time.monotonic()-start,
        sha256=digest(target),primary_profiles_sha256=digest(OUTPUT/'profiles.jsonl') if probe else None,
        primary_truth_inputs=[],label_files_blocked_in_workers=True,probe_truth_input='Fixed true laterQ2 only' if probe else None,
        maximum_cost_replay_error=None if probe else max(r['cost_replay_error'] for r in results),inherited_pilot_frames=len(prior),
        numerical_screen_all_pass=None if probe else sum(all(c['stable'] for c in r['comparisons']['dense_contour']) for r in results)))
    print('STABILITY_'+command.upper()+'_COMPLETE',len(results),round(time.monotonic()-start,2),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=('prepare','resume-prepare','tests','pilot','infer','probe'))
    parser.add_argument('--workers',type=int,default=7);args=parser.parse_args()
    if args.command in ('prepare','resume-prepare'):prepare(resume=args.command=='resume-prepare')
    elif args.command=='tests':tests()
    else:calculate(args.command,args.workers)


if __name__=='__main__':main()
