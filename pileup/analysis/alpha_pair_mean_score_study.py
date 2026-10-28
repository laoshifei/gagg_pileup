"""Frozen mean-score development diagnostics, with an explicit advancement gate."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
import importlib.metadata
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import unittest

import numpy as np

from .alpha_pair_information_study import load,save,digest,rows,serializable
from .alpha_pair_joint_likelihood import ConservingResponse,JointNoise
from .alpha_pair_mean_score import DEFAULTS,fit_estimators,oracle_diagnostics,METHODS

PROJECT=Path(os.environ.get('GAGG_MEANSCORE_PROJECT',Path(__file__).resolve().parents[2])).resolve()
NAME='alpha_pair_mean_score_20261028'
EXPERIMENT=PROJECT/'pileup/experiments'/NAME
OUTPUT=PROJECT/'jupyter/artifacts'/NAME
RECORD=PROJECT/'jupyter/run_records'/NAME
PREVIOUS=PROJECT/'pileup/experiments/alpha_pair_selection_20261027'
INPUT=PROJECT/'jupyter/artifacts/alpha_pair_selection_20261027'
OLD_RECORD=PROJECT/'jupyter/run_records/alpha_pair_selection_20261027'
NOTEBOOK=PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
TASKS=[dict(index=i,forced_K=2) for i in list(range(320,384,4))+list(range(384,448,4))]
TASKS += [dict(index=i,forced_K=1) for i in list(range(64,68))+list(range(256,260))]
TASKS += [dict(index=i,forced_K=0) for i in range(4)]


def now():return datetime.now(timezone.utc).isoformat()


def copy_chunks(source,target):
    """Explicit offsets avoid trusting broken shared-folder write counts.

    Full source/destination SHA256 equality is required before accepting a copy.
    The destination is new; failed copies remain for diagnosis.
    """
    source=Path(source);target=Path(target);offset=0;anomalies=[]
    fd=os.open(target,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o644)
    try:
        with source.open('rb') as stream:
            while block:=stream.read(256*1024):
                n=os.pwrite(fd,block,offset)
                if n!=len(block):anomalies.append(dict(offset=offset,expected=len(block),returned=n))
                offset+=len(block)
        os.ftruncate(fd,offset);os.fsync(fd)
    finally:os.close(fd)
    assert digest(source)==digest(target),(source,target,'copy hash mismatch')
    if anomalies:print('SHARED_COPY_VERIFIED_AFTER_COUNT_ANOMALY',str(target),anomalies,flush=True)


def prepare():
    previous=load(OLD_RECORD/'final-validation.json');assert previous['passed']
    assert digest(NOTEBOOK)==previous['notebook_sha256']
    for n,h in previous['source_sha256'].items():assert digest(PROJECT/n)==h,n
    for n,h in previous['documentation_sha256'].items():assert digest(PROJECT/n)==h,n
    for folder in (EXPERIMENT,OUTPUT,RECORD):folder.mkdir(parents=True,exist_ok=False)
    for src,dst in ((NOTEBOOK,'notebook.before.ipynb'),(PROJECT/'README.md','README.before.md'),(PROJECT/'docs/literature-use.md','literature.before.md')):
        copy_chunks(src,RECORD/dst)
    shutil.copytree(PREVIOUS/'source',EXPERIMENT/'source')
    for name in ('alpha_pair_mean_score.py','alpha_pair_mean_score_study.py','test_alpha_pair_mean_score.py'):
        copy_chunks(PROJECT/'pileup/analysis'/name,EXPERIMENT/'source/pileup/analysis'/name)
    for name in ('joint-response.npz','noise.npz','gates.json'):copy_chunks(PREVIOUS/name,EXPERIMENT/name)
    protocol=dict(created=now(),question='Separate finite-search failures,charge/time information and mean-only estimating equations before a new post-selection confirmation.',
        scope='Purealpha,fixed95/351ns,complete1024x2ns analog voltage,ample dynamic range; unchanged calibrated mean and JointNoise covariance.',
        development='Previously exposed selection confirmation data;systematic every4th pair:16 each40/200ns. Four each800/18800pe singles and four empties. This is development,not independent confirmation.',
        tasks=TASKS,pilot_indices=[320,384,64,0],settings=DEFAULTS,
        known_count='forced_K explicitly provided for development benchmarks; no claim of blind count recovery. Workers block labels and use only observed starts. K0 unchanged.',
        algorithms='Same coordinates,bounds,observed seed bank and L-BFGS-B budgets:wide full NLL;one fixed-C GLS;damped iterative GLS with C rebuilt only between solves. Baseline included. No logdet-dropping with varying C.',
        common_seeds='Observed K1 anchor +/-4/8/16/32/64/128/256/512ns;chronological shares .01/.03/.1/.3/.5/.7/.9/.97/.99. Union of3 diverse seeds by each initial objective plus baseline. K1 uses +/-8ns and .9/1/1.1 charge.',
        time_domain='Observed K1 trigger[-96,+512]ns,expanded only to contain cached observed candidate,clipped0..2046ns. Ordered times. Each Q50..50000pe.',
        irgls='Start frozenGLS;up to12 outer steps,each current/baseline/wideNLL starts,150iteration solves. Try damping .5/.25/.125/.0625 with fixed-C descent. Stop when undamped proposal changes all lnQ<=.002 and times<=.02ns,optimizer converged. Failures retained.',
        oracle='Only after primary predictions hash locks:fixed true times with observed Q starts;fixed true charges with observed time starts;truth-initialized free NLL;unfitted true parameters. Separate output,no feedback to primary.',
        advancement=dict(candidate_methods=['wide_nll','frozen_gls','iterated_gls'],minimum_40ns_strict_gain=2,maximum_200ns_strict_loss=1,minimum_pair_convergence_fraction=.8,
            maximum_isolated_strict_loss=0,rank='Eligible method with most strict successes in32pairs,then smaller all-pair median absolute ordered logQ2 error,then method order. Development choice only.',
            failure='If no method qualifies,report failure and do not generate a fresh confirmation sample this round.'),
        future_confirmation='Only after advancement passes:freeze chosen postfit before new112 Gaussian draws:16each empty,800/3000/8000/18800 singles,and18000/800 pairs40/200ns;t1=200.3ns. SeedSequence[202610280901,2,index,996].',
        confirmation_selection='Original and matched .025 gates from previous study remain fixed and act on unchanged baseline candidates. Only accepted K2 estimates are postfit. Counts retained even when revisedQ<500,flag separately;no truth in selection.',
        metrics='Strict bothQ errors<=20%,both times<=10ns,no extras;all known-K conditions marked. Also matched errors,sigma68,missed/extra/split/merge,boundary and convergence. Q is equivalentpe,not deposited energy;gap is not injection rate.',
        limitations='Finite optimizers,small exposed development,fixed conditional Gaussian parameter points. No nominal interval coverage,native transfer,free tau,alpha-gamma orK>=3. No default change.',
        workers_limit=8,workers=2,blas_threads=1,new_ROOTs=0,new_training=0,seed='No new randomness in PhaseA;conditional fresh design above only if advancement passes.')
    save(EXPERIMENT/'protocol.json',protocol)
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for key in ('MEANSCORE','SELECTION','CONTOUR','STABILITY','EVIDENCE','GATE','RISK','INFORMATION'):os.environ['GAGG_'+key+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_mean_score_study import main
    main()
''')
    files={str(p.relative_to(EXPERIMENT)):digest(p) for p in (EXPERIMENT/'source').rglob('*.py')}
    for n in ('run.py','protocol.json','joint-response.npz','noise.npz','gates.json'):files[n]=digest(EXPERIMENT/n)
    save(EXPERIMENT/'freeze.json',dict(created=now(),files=files))
    names=['test/observations.npz','test/predictions.jsonl','test/labels.truth.npz','test/labels.cases.json','report.json']
    save(RECORD/'input-provenance.json',dict(previous_validation=previous,previous_validation_sha256=digest(OLD_RECORD/'final-validation.json'),
        previous_freeze_sha256=digest(PREVIOUS/'freeze.json'),inputs={n:digest(INPUT/n) for n in names}))
    status=subprocess.check_output(['git','-c','core.fileMode=false','status','--short','--untracked-files=no','--','README.md','docs/literature-use.md','pileup/analysis'],cwd=PROJECT,text=True,timeout=60)
    save(RECORD/'environment.json',dict(created=now(),python=sys.version,executable=sys.executable,
        dependencies={n:importlib.metadata.version(n) for n in ('numpy','scipy','matplotlib','pandas','nbformat')},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),scoped_git_status=status,
        uncommitted='Existing changes retained;actual source/model/config snapshots and Notebook/document backups. No commit or push.'))
    print('MEAN_SCORE_PREPARED',len(files),flush=True)


def check_freeze():
    for n,h in load(EXPERIMENT/'freeze.json')['files'].items():assert digest(EXPERIMENT/n)==h,n


def models():
    return ConservingResponse.load(EXPERIMENT/'joint-response.npz'),JointNoise.load(EXPERIMENT/'noise.npz')


def tests():
    check_freeze();suite=unittest.TestSuite()
    for name in ('test_alpha_pair_mean_score','test_alpha_pair_model_sampler'):
        suite.addTests(unittest.defaultTestLoader.loadTestsFromName('pileup.analysis.'+name))
    result=unittest.TextTestRunner(verbosity=2).run(suite);assert result.wasSuccessful()
    save(RECORD/'tests.json',dict(passed=True,tests=result.testsRun,created=now()))


def initialize(probe=False):
    global RESPONSE,NOISE,OBS,BASE,CONFIG,PRIMARY
    from .alpha_pair_charge_allocation import forbid_labels
    RESPONSE,NOISE=models();sys.addaudithook(forbid_labels)
    with np.load(INPUT/'test/observations.npz') as f:OBS=f['voltage_mV']
    BASE={r['index']:r['candidates'] for r in rows(INPUT/'test/predictions.jsonl')}
    CONFIG=load(EXPERIMENT/'protocol.json')['settings']
    PRIMARY={r['index']:r for r in rows(OUTPUT/'predictions.jsonl')} if probe else {}


def infer_one(task):
    began=time.monotonic();i=task['index'];k=task['forced_K'];base=BASE[i]
    if k==0:
        return dict(**task,fits={m:base[0] for m in METHODS},known_count=True,seconds=time.monotonic()-began,cost_replay_error=0.,truth_inputs=[])
    result=fit_estimators(RESPONSE,NOISE,OBS[i],base[1],base[k],CONFIG)
    return dict(**task,**result,known_count=True,seconds=time.monotonic()-began)


def oracle_one(task):
    began=time.monotonic();i=task['index'];k=len(task['true_charges'])
    result=oracle_diagnostics(RESPONSE,NOISE,OBS[i],BASE[i][1],BASE[i][k],PRIMARY[i],task['true_times'],task['true_charges'],CONFIG)
    return dict(index=i,**result,seconds=time.monotonic()-began)


def calculate(command,workers):
    check_freeze();assert load(RECORD/'tests.json')['passed'] and 1<=workers<=8
    for n,h in load(RECORD/'input-provenance.json')['inputs'].items():assert digest(INPUT/n)==h,n
    tasks=load(EXPERIMENT/'protocol.json')['tasks'];prior=[];oracle=command=='oracle'
    if command=='pilot':tasks=[t for t in tasks if t['index'] in load(EXPERIMENT/'protocol.json')['pilot_indices']]
    elif command=='infer':
        v=load(RECORD/'pilot-validation.json');assert v['passed'] and digest(OUTPUT/'pilot.jsonl')==v['sha256']
        prior=rows(OUTPUT/'pilot.jsonl');done={r['index'] for r in prior};tasks=[t for t in tasks if t['index'] not in done]
    elif oracle:
        v=load(RECORD/'inference-validation.json');assert v['passed'] and digest(OUTPUT/'predictions.jsonl')==v['sha256']
        with np.load(INPUT/'test/labels.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
        tasks=[dict(index=c['index'],true_charges=q[c['index'],:c['forced_K']].tolist(),true_times=t[c['index'],:c['forced_K']].tolist()) for c in tasks if c['forced_K']]
    target=OUTPUT/('oracle-diagnostics.jsonl' if oracle else 'pilot.jsonl' if command=='pilot' else 'predictions.jsonl')
    began=time.monotonic();result=list(prior)
    with target.open('x') as stream,ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize,initargs=(oracle,)) as pool:
        for row in prior:stream.write(json.dumps(row,default=serializable,allow_nan=False)+'\n')
        stream.flush()
        futures=[pool.submit(oracle_one if oracle else infer_one,t) for t in tasks]
        for f in as_completed(futures):
            row=f.result();result.append(row);stream.write(json.dumps(row,default=serializable,allow_nan=False)+'\n');stream.flush()
            print('MEAN_SCORE_'+command.upper(),len(result),len(prior)+len(tasks),row['index'],round(row['seconds'],2),flush=True)
    expected=4 if command=='pilot' else 40 if oracle else 44
    assert len(result)==expected and len({r['index'] for r in result})==expected
    validation=dict(passed=True,created=now(),frames=expected,sha256=digest(target),seconds=time.monotonic()-began,workers=workers,
        inherited_pilot=len(prior),primary_truth_inputs=[],known_count=True,label_files_blocked=True,
        max_cost_replay_error=max(r.get('cost_replay_error',0.) for r in result),oracle_truth_assistance=oracle,
        primary_sha256=digest(OUTPUT/'predictions.jsonl') if oracle else None)
    save(RECORD/('pilot-validation.json' if command=='pilot' else 'oracle-validation.json' if oracle else 'inference-validation.json'),validation)
    print('MEAN_SCORE_COMPLETE',command,expected,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('prepare','tests','pilot','infer','oracle'));p.add_argument('--workers',type=int,default=2);a=p.parse_args()
    if a.command=='prepare':prepare()
    elif a.command=='tests':tests()
    else:calculate(a.command,a.workers)


if __name__=='__main__':main()
