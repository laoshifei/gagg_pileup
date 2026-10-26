"""Develop, freeze, and independently test precision regularization and births."""
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
import tempfile
import time

import numpy as np

from . import alpha_pair_conditional_study as previous
from .alpha_pair_conditional import ConditionalResponse
from .alpha_pair_resolution_data import load, save, digest, MODEL_NAMES
from .alpha_pair_residual import fit_candidates, choose

PROJECT=Path(os.environ.get('GAGG_RESIDUAL_PROJECT',Path(__file__).resolve().parents[2])).resolve()
NAME='alpha_pair_residual_20261023'
EXPERIMENT=PROJECT/'pileup/experiments'/NAME
OUTPUT=PROJECT/'jupyter/artifacts'/NAME
RECORD=PROJECT/'jupyter/run_records'/NAME
PARENT=PROJECT/'pileup/experiments/alpha_pair_conditional_20261022'
PARENT_OUTPUT=PROJECT/'jupyter/artifacts/alpha_pair_conditional_20261022'
POLICIES=dict(original=dict(fraction=0.,rho=None),rho90=dict(fraction=0.,rho=.9),rho80=dict(fraction=0.,rho=.8))


def cache():
    path=Path(load(RECORD/'local-cache.json')['directory'])
    if not path.exists():
        path=Path(tempfile.mkdtemp(prefix='gagg-residual-restored-'))
        for pattern in ('validation.*','test.*'):
            for source in OUTPUT.glob(pattern):shutil.copyfile(source,path/source.name)
        save(RECORD/'local-cache.json',dict(directory=str(path),restored_from=str(OUTPUT)))
    return path


def prepare():
    RECORD.mkdir(parents=True,exist_ok=False);OUTPUT.mkdir(parents=True,exist_ok=True)
    c=Path(tempfile.mkdtemp(prefix='gagg-residual-'))
    save(RECORD/'local-cache.json',dict(directory=str(c)))
    for name in ('codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb','Run_Catalog.ipynb'):
        shutil.copyfile(PROJECT/'jupyter'/name,RECORD/(name+'.before'))
    shutil.copyfile(PROJECT/'docs/literature-use.md',RECORD/'literature-use.before.md')
    for name in ('settings.json','response.npz','calibration.json','selection.json','fitting-selection.json'):
        shutil.copyfile(PARENT/name,EXPERIMENT/name)
    old_cache=Path(load(PROJECT/'jupyter/run_records/alpha_pair_conditional_20261022/local-cache.json')['directory'])
    source=old_cache if old_cache.exists() else PARENT_OUTPUT
    cases=load(source/'test.cases.json')
    # All native single/empty controls; first20 pairs at every condition/gap.
    # Previous test parents have been exposed and are explicitly development.
    indices=[i for i,r in enumerate(cases) if r['kind']!='pair' or r['replicate']<20]
    with np.load(source/'test.waveforms.npz') as f:
        np.savez_compressed(c/'validation.waveforms.npz',voltage_mV=f['voltage_mV'][indices])
    with np.load(source/'test.truth.npz') as f:
        np.savez_compressed(c/'validation.truth.npz',charge_pe=f['charge_pe'][indices],times_ns=f['times_ns'][indices])
    save(c/'validation.cases.json',[cases[i] for i in indices])
    old_seeds=load(source/'test.seeds.json')
    save(c/'validation.seeds.json',{k:[v[i] for i in indices] for k,v in old_seeds.items()})
    old_pred={r['index']:r for r in map(json.loads,(source/'test.predictions.jsonl').read_text().splitlines())}
    with (c/'validation.previous.jsonl').open('w') as stream:
        for j,i in enumerate(indices):stream.write(json.dumps(dict(old_pred[i],index=j))+'\n')
    for p in c.glob('validation.*'):shutil.copyfile(p,OUTPUT/p.name)
    save(EXPERIMENT/'development.json',dict(previous_test_now_development=True,
        calibration_parent=str((PROJECT/'pileup/experiments/alpha_pair_resolution_20261022').relative_to(PROJECT)),
        validation_parent=str(PARENT.relative_to(PROJECT)),validation_indices=indices,
        validation_frames=len(indices),validation_optical_parents=240,pairs_per_point=20,
        controls_per_group=80,policies=POLICIES,
        acceptance=dict(extra_output_each_control_max=.05,single_success_each_min=.95,
                        single_median_absolute_relative_bias_max=.05),
        selection='Maximize equal-weight pair success, subject to all control and bias constraints; tie prefer larger count penalty',
        initial_probe='104 training-parent frames; preliminary LS blend and rho shrinkage; precision numeric scale subsequently increased to1e6',
        files={p.name:digest(p) for p in c.glob('validation.*')}))
    (RECORD/'git-status-before.txt').write_text(subprocess.check_output(
        ['git','--no-optional-locks','-c','core.fileMode=false','status','--short'],cwd=PROJECT,text=True))
    save(RECORD/'environment.json',dict(python=sys.version,executable=sys.executable,numpy=np.__version__,
        dependency_path=os.environ.get('PYTHONPATH'),cpu_affinity=sorted(os.sched_getaffinity(0)),
        head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip()))
    print('PREPARED',len(indices),'validation frames',c,flush=True)


_RESPONSE=None
_POLICIES=None


def initialize(response_path,policies):
    global _RESPONSE,_POLICIES
    sys.addaudithook(previous.forbid_labels)
    _RESPONSE=ConditionalResponse.load(response_path);_POLICIES=policies
    previous._RESPONSE=_RESPONSE;previous._BASE=_RESPONSE.base


def infer_one(task):
    i,y,seeds,reference=task;started=time.monotonic()
    result=dict(index=i,policies={name:fit_candidates(_RESPONSE,y,seeds,**params) for name,params in _POLICIES.items()})
    if reference:result['previous']=previous.infer_one((i,y,seeds))
    result['seconds']=time.monotonic()-started
    return result


def infer(dataset,workers=4,limit=None):
    if dataset=='test':check_freeze()
    c=cache();policies=POLICIES
    if dataset=='test':
        selected=load(EXPERIMENT/'residual-selection.json')['chosen_policy']
        policies={name:POLICIES[name] for name in ('original',selected)}
    with np.load(c/(dataset+'.waveforms.npz')) as f:
        assert f.files==['voltage_mV'];y=f['voltage_mV']
    seeds=load(c/(dataset+'.seeds.json'))
    dest=c/(dataset+'.predictions.jsonl');done={}
    if dest.exists():
        for line in dest.read_text().splitlines():
            r=json.loads(line);assert r['index'] not in done;done[r['index']]=r
    indices=list(range(len(y))) if limit is None else list(range(min(limit,len(y))))
    enabled=[True];opened=set()
    def guard(event,args):
        if not enabled[0]:return
        previous.forbid_labels(event,args)
        if event=='open' and isinstance(args[0],str) and args[0].endswith(('.json','.npz','.pt')):opened.add(args[0])
    sys.addaudithook(guard)
    try:
        with ProcessPoolExecutor(max_workers=workers,initializer=initialize,
                initargs=(EXPERIMENT/'response.npz',policies),mp_context=mp.get_context('spawn')) as pool:
            futures=[pool.submit(infer_one,(i,y[i],{k:v[i] for k,v in seeds.items()},dataset=='test')) for i in indices if i not in done]
            with dest.open('a') as stream:
                for future in as_completed(futures):
                    row=future.result();done[row['index']]=row
                    stream.write(json.dumps(row,allow_nan=False)+'\n');stream.flush()
                    if len(done)%50==0 or len(done)==len(indices):print('INFER',dataset,len(done),'/',len(indices),flush=True)
    finally:enabled[0]=False
    shutil.copyfile(dest,OUTPUT/dest.name)
    save(RECORD/(dataset+'-inference.json'),dict(passed=True,events=len(done),requested=len(indices),
        policies=policies,workers=workers,truth_inputs=[],label_access_blocked=True,
        observed_files=[str(c/(dataset+'.waveforms.npz')),str(c/(dataset+'.seeds.json'))],
        opened_while_inference=sorted(opened),prediction_sha256=digest(dest),
        waveform_sha256=digest(c/(dataset+'.waveforms.npz')),mean_seconds=float(np.mean([r['seconds'] for r in done.values()]))))


def select():
    c=cache();cases=load(c/'validation.cases.json')
    predictions=sorted(map(json.loads,(c/'validation.predictions.jsonl').read_text().splitlines()),key=lambda r:r['index'])
    assert len(predictions)==len(cases)
    with np.load(c/'validation.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
    table=[]
    for policy in POLICIES:
        for variant in ('expanded','existing','waveform_only'):
            for penalty in (.00001,.00003,.0001,.0003,.001,.003,.01,.03,.1):
                scores=[previous.evaluate_fit(choose(p['policies'][policy][variant],penalty),t[i],q[i],cases[i]) for i,p in enumerate(predictions)]
                pair=[s['success'] for case,s in zip(cases,scores) if case['kind']=='pair']
                control={}
                for label in ('all','low','mid','high'):
                    ss=[s for case,s in zip(cases,scores) if case['kind']!='pair' and case['condition']==label]
                    errors=[e for s in ss for e in s['q_errors']]
                    control[label]=dict(frames=len(ss),success=float(np.mean([s['success'] for s in ss])),
                        extra_fraction=float(np.mean([s['extra']>0 for s in ss])),
                        q_median_relative_error=float(np.median(errors)) if errors else None)
                eligible=all(r['extra_fraction']<=.05 for r in control.values()) and all(
                    r['success']>=.95 and abs(r['q_median_relative_error'])<=.05 for k,r in control.items() if k!='all')
                table.append(dict(policy=policy,variant=variant,penalty=penalty,pair_success=float(np.mean(pair)),
                                  eligible=eligible,controls=control))
    eligible=[r for r in table if r['eligible'] and r['variant']=='expanded']
    assert eligible,'No policy passed the predeclared control and bias requirements'
    chosen=max(eligible,key=lambda r:(r['pair_success'],r['penalty'],r['policy']))
    penalties={}
    for policy in ('original',chosen['policy']):
        penalties[policy]={}
        for variant in ('expanded','existing','waveform_only'):
            rows=[r for r in table if r['policy']==policy and r['variant']==variant]
            # The original precision is a declared ablation even when its low-Q
            # bias fails. Its penalty still respects the false-output budget.
            allowed=[r for r in rows if all(s['extra_fraction']<=.05 for s in r['controls'].values())]
            best=max(allowed,key=lambda r:(r['pair_success'],r['penalty']))
            if policy==chosen['policy']:
                passed=[r for r in allowed if r['eligible']]
                if passed:best=max(passed,key=lambda r:(r['pair_success'],r['penalty']))
            penalties[policy][variant]=best
    save(EXPERIMENT/'residual-selection.json',dict(chosen_policy=chosen['policy'],chosen=chosen,
        penalties=penalties,grid=table,selection_uses='Exposed development parents only',
        inference_truth_inputs=[],models_retrained=False))
    print('SELECTED',chosen,flush=True)


def check_freeze():
    frozen=load(EXPERIMENT/'freeze.json')
    for name,sha in frozen['files'].items():assert digest(EXPERIMENT/name)==sha,name
    return frozen


def adapt():
    previous.PROJECT=PROJECT;previous.EXPERIMENT=EXPERIMENT;previous.OUTPUT=OUTPUT;previous.RECORD=RECORD
    previous.cache=cache;previous.check_freeze=check_freeze
    return previous


def freeze():
    assert not (EXPERIMENT/'freeze.json').exists()
    choice=load(EXPERIMENT/'residual-selection.json')
    for name in ('source','models','bin'):shutil.copytree(PARENT/name,EXPERIMENT/name)
    for namespace in ('analysis','pileup/analysis'):
        dest=EXPERIMENT/'inference-source'/namespace;dest.mkdir(parents=True)
        for path in (PROJECT/namespace).glob('*.py'):shutil.copyfile(path,dest/path.name)
    (EXPERIMENT/'inference-source/pileup/__init__.py').write_text('')
    for name in ('macros','logs'):(EXPERIMENT/name).mkdir()
    protocol=load(PARENT/'protocol.json')
    protocol.update(created=datetime.now().astimezone().isoformat(),source_seed_start=263500000,
        digitizer_seed_start=263510000,purpose='Regularized precision and residual count expansion',
        validation_selection=str((EXPERIMENT/'residual-selection.json').relative_to(PROJECT)),
        selected_policy=choice['chosen_policy'],policies=POLICIES,
        new_timing_search_radius_ns=24.,original_baseline_radius_ns=16.,test_not_started=True,
        examples_predeclared=[dict(condition='high_low',replicate=0,gap_ns=gap) for gap in (40.,200.)],
        development_input_sha256=digest(EXPERIMENT/'development.json'),
        development_predictions_sha256=digest(OUTPUT/'validation.predictions.jsonl'),
        neural_runtime_preflight_sha256=digest(RECORD/'neural-runtime-preflight.json'),
        time_definition='First photon inside the SiPM active plane; pair gap is between these arrivals, not emission times or maxima',
        single_crop_jitter='Previous deterministic index pattern, not optical or SiPM RNG reuse')
    save(EXPERIMENT/'protocol.json',protocol)
    plans=[]
    for energy in protocol['energies_MeV']:
        for rep in range(80):
            plans.append(dict(label=f'alpha_{energy:g}_{rep:03d}',particle='alpha',energy_MeV=energy,
                x_mm=0.,y_mm=0.,replicate=rep,role='independent_residual_test',gamma_model='brems',
                seed=protocol['source_seed_start']+len(plans)))
    seeds={r['seed'] for r in plans};audit=[]
    for path in sorted((PROJECT/'pileup/experiments').glob('*/manifest.json')):
        old=load(path);used={r.get('seed') for key in ('runs','planned_runs') for r in old.get(key,[])}
        assert not seeds&used,str(path)
        audit.append(dict(path=str(path.relative_to(PROJECT)),sha256=digest(path),source_seeds=len(used-{None})))
    save(RECORD/'parent-independence-audit.json',dict(passed=True,new_seeds=240,manifests=audit,collisions=[]))
    # Check actual old digitization inputs as well as G4 source manifests.
    # Gap scans intentionally share one digitizer base seed within a pair.
    bases={protocol['digitizer_seed_start']+ci*10000+rep for ci in range(5) for rep in range(40)}
    bases.update(protocol['digitizer_seed_start']+100000+rep for rep in range(80))
    new_streams={seed+offset for seed in bases for offset in (0,100000,200000)}
    old_case_audit=[]
    for parent in ('alpha_pair_resolution_20261022','alpha_pair_conditional_20261022'):
        for path in sorted((PROJECT/'jupyter/artifacts'/parent).glob('*.cases.json')):
            cases=load(path)
            old_bases={r['seed'] for r in cases if 'seed' in r}
            old_streams={seed+offset for seed in old_bases for offset in (0,100000,200000)}
            assert not new_streams&old_streams,str(path)
            old_case_audit.append(dict(path=str(path.relative_to(PROJECT)),sha256=digest(path),frames=len(cases)))
    save(RECORD/'digitizer-independence-audit.json',dict(passed=True,new_base_seeds=len(bases),
        comparisons=old_case_audit,collisions=[],
        note='Common random numbers across gaps retained; new base/avalanche/readout seed values disjoint from old cases'))
    old=load(PARENT/'manifest.json')
    save(EXPERIMENT/'manifest.json',dict(created=protocol['created'],planned_runs=plans,runs=[],
        preflight_validated=False,max_workers=4,window_us=20.,environment=old['environment'],
        binary_sha256=digest(EXPERIMENT/'bin/gagg-brems')))
    files={}
    for folder in ('source','inference-source','models','bin'):
        for path in (EXPERIMENT/folder).rglob('*'):
            if path.is_file():files[str(path.relative_to(EXPERIMENT))]=digest(path)
    for name in ('settings.json','response.npz','calibration.json','selection.json','fitting-selection.json',
                 'residual-selection.json','protocol.json','development.json','run.py','report.py',
                 'validate.py','deliver.py','prefill_photons.py','production.py'):
        files[name]=digest(EXPERIMENT/name)
    save(EXPERIMENT/'freeze.json',dict(frozen_at=datetime.now().astimezone().isoformat(),files=files,
        new_test_generation_started=False,git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip()))
    (RECORD/'git-status-at-freeze.txt').write_text(subprocess.check_output(
        ['git','--no-optional-locks','-c','core.fileMode=false','status','--short'],cwd=PROJECT,text=True))
    print('FROZEN',len(files),'files',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=('prepare','infer','select','freeze','native','physical','seeds'))
    p.add_argument('--dataset',default='validation');p.add_argument('--workers',type=int,default=4);p.add_argument('--limit',type=int)
    a=p.parse_args()
    if a.stage=='prepare':prepare()
    elif a.stage=='infer':infer(a.dataset,a.workers,a.limit)
    elif a.stage=='select':select()
    elif a.stage=='freeze':freeze()
    elif a.stage=='native':adapt().native(a.workers)
    elif a.stage=='physical':adapt().physical(a.workers)
    elif a.stage=='seeds':adapt().test_seeds()


if __name__=='__main__':main()
