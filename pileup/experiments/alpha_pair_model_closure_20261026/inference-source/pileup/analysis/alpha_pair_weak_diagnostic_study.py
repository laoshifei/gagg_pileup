"""Frozen, resumable development-only diagnosis of exposed weak-pair data."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

from .alpha_pair_resolution_data import digest, load, save

PROJECT = Path(os.environ.get('GAGG_WEAK_PROJECT', Path(__file__).resolve().parents[2])).resolve()
NAME = 'alpha_pair_weak_diagnostics_20261026'
EXPERIMENT = PROJECT/'pileup/experiments'/NAME
OUTPUT = PROJECT/'jupyter/artifacts'/NAME
RECORD = PROJECT/'jupyter/run_records'/NAME
PARENT = PROJECT/'pileup/experiments/alpha_pair_null_gate_20261023'
INPUT = PROJECT/'jupyter/artifacts/alpha_pair_null_gate_20261023/test'


def prepare():
    for path in (EXPERIMENT, OUTPUT, RECORD):
        path.mkdir(parents=True, exist_ok=False)
    notebook = PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    shutil.copyfile(notebook, RECORD/'notebook.before.ipynb')
    for name in ('response.npz', 'settings.json', 'calibration.json', 'model-null-calibration.json'):
        shutil.copyfile(PARENT/name, EXPERIMENT/name)
    shutil.copytree(PARENT/'source', EXPERIMENT/'source')
    cases = load(INPUT/'labels.cases.json')
    indices = [i for i, case in enumerate(cases) if case['kind'] != 'pair' or case['gap_ns'] in (8.,40.,200.)]
    assert len(indices) == 1600
    plan = dict(created=datetime.now().astimezone().isoformat(), indices=indices,
                gaps_ns=[8.,40.,200.], frames=len(indices), parent_frames=3200,
                optical_parents=320, new_native_simulation=False, neural_training=False,
                input_role='Previously published independent test, now exposed development data',
                question='Separate search failures from time-charge coupling and waveform fluctuations',
                primary_focus='high_weak at40ns; retain all ten charge combinations and native controls',
                truth_diagnostics=['Known K and times, fit Q with five observed-total charge starts',
                                   'True K/Q/t initialization, joint local refit',
                                   'Same two diagnostics with OLS',
                                   'Exact other-component subtraction with true neighbor Q/t and true weak start'],
                bounds=dict(charge_pe=[50.,50000.], local_radius_ns=24., max_nfev=120),
                selection='All pairs at8/40/200ns and every empty/native-single frame; no success-based selection',
                scope='Pure alpha, fixed95/351ns, 2ns samples, 2.048us observed record, sufficient dynamic range',
                uncertainty='Parents reused across gaps and conditions; not pooled independent exposure',
                inference_claim='Truth-assisted diagnoses are not blind recovery or a theoretical bound',
                examples=[dict(condition='high_weak',gap_ns=gap,replicate=0) for gap in (40.,200.)])
    save(EXPERIMENT/'protocol.json', plan)
    inputs = [INPUT/name for name in ('observations.npz','labels.truth.npz','labels.cases.json',
                                       'predictions.model.jsonl','predictions.jsonl','report.json')]
    inputs += [PARENT/'freeze.json', PARENT/'response.npz']
    save(RECORD/'input-provenance.json', dict(files={str(p.relative_to(PROJECT)):digest(p) for p in inputs},
                                            notebook_before_sha256=digest(notebook)))
    wrapper = '''from pathlib import Path
import os, sys
folder=Path(__file__).resolve().parent
for name in ('WEAK','NULL','CONDITIONAL','RESIDUAL','RESOLUTION','PAIR','COMPONENT'):
    os.environ['GAGG_'+name+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_weak_diagnostic_study import main
    main()
'''
    (EXPERIMENT/'run.py').write_text(wrapper)
    snapshot = EXPERIMENT/'inference-source'
    for namespace in ('analysis','pileup/analysis'):
        destination=snapshot/namespace; destination.mkdir(parents=True)
        for source in (PROJECT/namespace).glob('*.py'):
            shutil.copyfile(source,destination/source.name)
    if (PROJECT/'pileup/__init__.py').exists():
        shutil.copyfile(PROJECT/'pileup/__init__.py',snapshot/'pileup/__init__.py')
    paths=[p for folder in ('source','inference-source') for p in (EXPERIMENT/folder).rglob('*') if p.is_file()]
    paths += [EXPERIMENT/name for name in ('response.npz','settings.json','calibration.json',
                                         'model-null-calibration.json','protocol.json','run.py')]
    save(EXPERIMENT/'freeze.json',dict(frozen_at=datetime.now().astimezone().isoformat(),
        files={str(p.relative_to(EXPERIMENT)):digest(p) for p in paths},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        data_role='Exposed development; no new independent test yet',frozen_before_diagnostics=True))
    save(RECORD/'environment.json',dict(python=sys.version,executable=sys.executable,numpy=np.__version__,
        cpu_affinity=sorted(os.sched_getaffinity(0)),workers=4,
        git_changes=subprocess.check_output(['git','--no-optional-locks','-c','core.filemode=false','status',
            '--short','--','pileup/analysis/alpha_pair_weak_diagnostics.py',
            'pileup/analysis/alpha_pair_weak_diagnostic_study.py',
            'pileup/analysis/test_alpha_pair_weak_diagnostics.py'],cwd=PROJECT,text=True)))
    print('DIAGNOSTIC_PREPARED',len(paths),'frozen files',len(indices),'frames',flush=True)


def check_freeze():
    for name,sha in load(EXPERIMENT/'freeze.json')['files'].items():
        assert digest(EXPERIMENT/name)==sha,name


def predictions(filename):
    rows=[json.loads(line) for line in (INPUT/filename).read_text().splitlines()]
    result={r['index']:r['candidates'] for r in rows}
    assert len(result)==len(rows)==3200
    return result


def initialize_worker():
    global RESPONSE
    from .alpha_pair_conditional import ConditionalResponse
    RESPONSE=ConditionalResponse.load(EXPERIMENT/'response.npz')


def one(task):
    from .alpha_pair_weak_diagnostics import diagnose
    i,y,t,q,components,candidates=task
    started=time.monotonic()
    return dict(index=i,diagnostics=diagnose(RESPONSE,y,t,q,components,candidates),
                seconds=time.monotonic()-started)


def run(pilot=False,workers=4):
    check_freeze()
    indices=load(EXPERIMENT/'protocol.json')['indices'];cases=load(INPUT/'labels.cases.json')
    if pilot:
        indices=[i for i in indices if cases[i]['replicate']==0]
    with np.load(INPUT/'observations.npz') as f:y=f['voltage_mV']
    with np.load(INPUT/'labels.truth.npz') as f:
        q=f['charge_pe'];t=f['times_ns'];components=f['components_mV']
    old=predictions('predictions.model.jsonl')
    path=OUTPUT/'diagnostics.jsonl'
    previous=[json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []
    done={r['index'] for r in previous};assert len(done)==len(previous)
    tasks=[(i,y[i],t[i],q[i],components[i],old[i]) for i in indices if i not in done]
    started=time.monotonic();count=0
    with path.open('a') as stream:
        with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),
                                 initializer=initialize_worker) as pool:
            futures=[pool.submit(one,task) for task in tasks]
            for future in as_completed(futures):
                row=future.result()
                stream.write(json.dumps(row,allow_nan=False)+'\n');stream.flush()
                count+=1
                if count%20==0:
                    os.fsync(stream.fileno())
                    print('DIAGNOSED',count,'/',len(tasks),'seconds',round(time.monotonic()-started,1),flush=True)
        os.fsync(stream.fileno())
    save(RECORD/('pilot-run.json' if pilot else 'run.json'),dict(completed=True,requested=len(indices),
        newly_evaluated=count,seconds=time.monotonic()-started,workers=workers,
        explicit_truth_inputs=['K','times','charges for initialization/conditioning','component ancestry'],
        previous_blind_predictions_unchanged=True))
    print('DIAGNOSTIC_COMPLETE','pilot' if pilot else 'full',len(indices),flush=True)


def report():
    from .alpha_pair_conditional_study import evaluate_fit
    from .alpha_pair_null_report import summarize
    from .alpha_pair_null_gate import select
    from .alpha_pair_residual import choose
    cases=load(INPUT/'labels.cases.json')
    with np.load(INPUT/'labels.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
    original=predictions('predictions.jsonl');model=predictions('predictions.model.jsonl')
    calibration=load(EXPERIMENT/'model-null-calibration.json')['relative_residual']
    diagnostics=sorted([json.loads(line) for line in (OUTPUT/'diagnostics.jsonl').read_text().splitlines()],key=lambda r:r['index'])
    groups={};scored=[]
    for row in diagnostics:
        i=row['index'];case=cases[i];d=row['diagnostics']
        fits=dict(observed_previous=choose(original[i],.001),blind_previous=select(model[i],calibration))
        fits.update({k:d[k] for k in ('known_count_previous','fixed_true_times','fixed_true_times_ols',
            'truth_start_joint','truth_start_joint_ols','oracle_subtract_other') if k in d})
        for method,fit in fits.items():
            score=dict(evaluate_fit(fit,t[i],q[i],case),index=i,method=method,kind=case['kind'],
                condition=case['condition'],gap_ns=case['gap_ns'],true_q=q[i].tolist(),true_t=t[i].tolist(),
                known_count=method not in ('observed_previous','blind_previous'),
                diagnostic_truth=fit.get('diagnostic_truth', ['K'] if method=='known_count_previous' else []))
            scored.append(score)
            groups.setdefault((method,case['kind'],case['condition'],case['gap_ns']),[]).append(score)
    summary=[summarize(rows) for rows in groups.values()]
    save(OUTPUT/'report.json',dict(data_role='Exposed development, truth-assisted diagnosis',
        frames=len(diagnostics),score_rows=len(scored),summary=summary,
        complete=len(diagnostics)==len(load(EXPERIMENT/'protocol.json')['indices']),
        maximum_metric_replay_error=max(r['diagnostics']['metric_replay_relative_error'] for r in diagnostics),
        cost_comparison='Only common-metric full-frame GLS costs can be compared; OLS and oracle-subtracted costs differ',
        no_independent_performance_claim=True))
    with (OUTPUT/'scores.jsonl').open('w') as stream:
        for row in scored:stream.write(json.dumps(row,allow_nan=False)+'\n')
    for row in summary:
        if row['condition']=='high_weak':
            print('HIGH_WEAK',row['gap_ns'],row['method'],row['success'],'/',row['eligible_frames'],flush=True)
    print('DIAGNOSTIC_REPORT',len(diagnostics),len(scored),flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('command',choices=['prepare','run','report'])
    parser.add_argument('--pilot',action='store_true')
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    if args.command=='prepare':prepare()
    elif args.command=='run':run(args.pilot,args.workers)
    else:report()


if __name__=='__main__':
    main()
