"""Frozen independent optical-parent confirmation of the joint likelihood."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime
import importlib
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

from .alpha_pair_resolution_data import load,save,digest

PROJECT=Path(os.environ.get('GAGG_CONFIRM_PROJECT',Path(__file__).resolve().parents[2])).resolve()
NAME='alpha_pair_joint_confirmation_20261026'
EXPERIMENT=PROJECT/'pileup/experiments'/NAME
OUTPUT=PROJECT/'jupyter/artifacts'/NAME
RECORD=PROJECT/'jupyter/run_records'/NAME
PARENT=PROJECT/'pileup/experiments/alpha_pair_joint_likelihood_20261026'
PHYSICS=PROJECT/'pileup/experiments/alpha_pair_null_gate_20261023'
CONDITIONS=['low_low','mid_mid','high_high','low_high','high_low','low_weak','weak_low','high_weak','weak_high','weak_weak']


def check_freeze():
    for p,sha in load(EXPERIMENT/'freeze.json')['files'].items():assert digest(EXPERIMENT/p)==sha,p


def prepare():
    selection=load(PARENT/'development-selection.json')
    assert selection['chosen']=='joint_conserved_nll'
    assert load(PROJECT/'jupyter/run_records'/PARENT.name/'science-validation.json')['passed']
    for folder in (EXPERIMENT,OUTPUT,RECORD):folder.mkdir(parents=True,exist_ok=False)
    for part in ('source','bin'):shutil.copytree(PHYSICS/part,EXPERIMENT/part)
    for part in ('macros','logs'):(EXPERIMENT/part).mkdir()
    for name in ('response.npz','joint-response.npz','settings.json','joint_conserved_nll-noise.npz','gate-calibrations.json','development-selection.json'):
        shutil.copyfile(PARENT/name,EXPERIMENT/name)
    shutil.copyfile(PHYSICS/'model-null-calibration.json',EXPERIMENT/'legacy-gate-calibrations.json')
    source_seed=276000000;digitizer_seed=276500000
    plans=[dict(label=f'alpha_{energy:g}_{rep:03d}',particle='alpha',energy_MeV=energy,x_mm=0.,y_mm=0.,replicate=rep,
                role='independent_joint_confirmation',gamma_model='brems') for energy in (1.625,2.5,5.486,8.) for rep in range(80)]
    for i,r in enumerate(plans):r['seed']=source_seed+i
    # Audit primary, avalanche, readout and joint-digitization seeds before production.
    bases={r['seed'] for r in plans}|{digitizer_seed+ci*1000+rep for ci in range(10) for rep in range(40)}|{digitizer_seed+50000+rep for rep in range(80)}
    streams={s+o for s in bases for o in (0,100000,200000)};audit=[]
    paths=list((PROJECT/'pileup/experiments').glob('*/manifest.json'))+list((PROJECT/'pileup/experiments').glob('*/*/manifest.json'))
    for path in paths:
        data=load(path)
        if not isinstance(data,dict):continue
        old={r['seed'] for key in ('planned_runs','runs') for r in data.get(key,[]) if 'seed' in r}
        assert not streams&{s+o for s in old for o in (0,100000,200000)},str(path)
        audit.append(dict(path=str(path.relative_to(PROJECT)),sha256=digest(path)))
    case_paths=list((PROJECT/'jupyter/artifacts').glob('*/*.cases.json'))+list((PROJECT/'jupyter/artifacts').glob('*/*/labels.cases.json'))
    for path in case_paths:
        data=load(path)
        if isinstance(data,dict):data=data.get('cases',[])
        old={r['seed'] for r in data if isinstance(r,dict) and isinstance(r.get('seed'),int)}
        assert not streams&{s+o for s in old for o in (0,100000,200000)},str(path)
    save(RECORD/'seed-independence.json',dict(passed=True,source_seeds=[r['seed'] for r in plans],new_streams=len(streams),
        manifests=audit,case_files=[str(p.relative_to(PROJECT)) for p in case_paths]))
    protocol=dict(created=datetime.now().astimezone().isoformat(),source_energies_MeV=[1.625,2.5,5.486,8.],
        native_roots=320,per_energy=80,gaps_ns=[8.,40.,200.],conditions=CONDITIONS,pairs_per_point=40,empty_frames=80,
        expected_frames=1600,source_seed_start=source_seed,digitizer_seed_start=digitizer_seed,
        selected_method='joint_conserved_nll',baselines=['observed_penalty','model_relative_residual','model_matched_gate'],
        fixed_taus_ns=[95.,351.],new_test_generation_started=False,gate_refit_on_test=False,
        metrics='Both Q>=500pe; each Q within20%, each t within min(10ns,gap/4); no extra output',
        scope='Pure central alpha, fixed double tau, sufficient analog range, 2ns x1024 local samples, unknown K<=2',
        selection=load(PARENT/'method-plan.json')['selection'],
        examples=[dict(condition='high_weak',gap_ns=gap,replicate=0) for gap in (40.,200.)],
        independence='Fresh optical parents and seeds after method freeze; parents reused across conditions/gaps; point estimates descriptive')
    save(EXPERIMENT/'protocol.json',protocol)
    previous=load(PHYSICS/'test/manifest.json')
    save(EXPERIMENT/'manifest.json',dict(created=protocol['created'],planned_runs=plans,runs=[],preflight_validated=False,
        max_workers=4,window_us=20.,environment=previous['environment'],binary_sha256=digest(EXPERIMENT/'bin/gagg-brems')))
    for namespace in ('analysis','pileup/analysis'):
        dest=EXPERIMENT/'inference-source'/namespace;dest.mkdir(parents=True)
        for p in (PROJECT/namespace).glob('*.py'):shutil.copyfile(p,dest/p.name)
    (EXPERIMENT/'inference-source/pileup/__init__.py').write_text('')
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for name in ('CONFIRM','JOINT','NULL','CONDITIONAL','RESOLUTION','PAIR','COMPONENT'):
    os.environ['GAGG_'+name+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_joint_confirmation import main
    main()
''')
    files=[p for part in ('source','bin','inference-source') for p in (EXPERIMENT/part).rglob('*') if p.is_file()]
    files += [EXPERIMENT/n for n in ('response.npz','joint-response.npz','settings.json','joint_conserved_nll-noise.npz',
        'gate-calibrations.json','legacy-gate-calibrations.json','development-selection.json','protocol.json','run.py')]
    save(EXPERIMENT/'freeze.json',dict(files={str(p.relative_to(EXPERIMENT)):digest(p) for p in files},
        frozen_at=datetime.now().astimezone().isoformat(),git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        original_development_freeze_sha256=digest(PARENT/'development-freeze-v2.json'),new_test_generation_started=False,
        uncommitted_source='Actual current analysis and unchanged physics copied; Git HEAD alone is not the source version'))
    shutil.copyfile(PROJECT/'jupyter/Run_Catalog.ipynb',RECORD/'Run_Catalog.before.ipynb')
    print('CONFIRMATION_FROZEN',len(files),'files',len(plans),'planned ROOTs',flush=True)


def data_module():
    from . import alpha_pair_resolution_data as data
    data.PROJECT=PROJECT;data.EXPERIMENT=EXPERIMENT;data.OUTPUT=OUTPUT;data.RECORD=RECORD;data.check_freeze=check_freeze
    return data


def native(workers=4):data_module().native(workers)


def physical(workers=4):
    from .alpha_pair_null_study import initialize_physical,physical_one
    from .alpha_pair_resolution_data import joint_case
    from . import sipm_response
    from .alpha_pair_digitizer_optimized import simulate_avalanches
    check_freeze();(OUTPUT/'photons').mkdir(exist_ok=True)
    records=data_module().donor_cache();assert len(records)==320
    local=Path(tempfile.mkdtemp(prefix='gagg-joint-confirm-'));(local/'photons').mkdir()
    for i in range(len(records)):shutil.copyfile(OUTPUT/'photons'/f'{i:04}.npz',local/'photons'/f'{i:04}.npz')
    save(local/'sensor-config.json',load(PROJECT/records[0]['metadata'])['config'])
    save(RECORD/'digitizer-cache.json',dict(directory=str(local)))
    p=load(EXPERIMENT/'protocol.json');base=p['digitizer_seed_start'];energies=dict(weak=1.625,low=2.5,mid=5.486,high=8.)
    groups={name:[i for i,r in enumerate(records) if r['energy_MeV']==e] for name,e in energies.items()}
    assert all(len(v)==80 for v in groups.values());cases=[]
    for ci,condition in enumerate(CONDITIONS):
        first,second=condition.split('_')
        pairs=(list(zip(groups[first][::2],groups[second][1::2])) if first==second else list(zip(groups[first][:40],groups[second][:40])))
        for gap in p['gaps_ns']:
            for rep,parents in enumerate(pairs):
                seed=base+ci*1000+rep
                cases.append(dict(kind='pair',condition=condition,gap_ns=gap,replicate=rep,parents=list(parents),seed=seed,
                    first_ns=float(np.random.default_rng(seed).uniform(96.,256.)),cluster='parents-'+'-'.join(map(str,sorted(parents)))))
    for rep in range(80):cases.append(dict(kind='empty',condition='all',gap_ns=None,replicate=rep,parents=[],seed=base+50000+rep,first_ns=128.,cluster=f'empty-{rep}'))
    (OUTPUT/'physical-cases').mkdir(exist_ok=True)
    save(OUTPUT/'generation-plan.json',dict(cases=cases,digitizer_base=base,common_random_numbers_across_gaps=True,photon_scaling=False,direct_voltage_addition=False))
    save(OUTPUT/'physical-parents.json',dict(records=records))
    initialize_physical(local,len(records));_,first=physical_one((0,cases[0],OUTPUT))
    importlib.reload(sipm_response)
    try:reference=joint_case(cases[0])
    finally:sipm_response.simulate_avalanches=simulate_avalanches
    for actual,expected in zip(first[:4],reference[:4]):np.testing.assert_array_equal(actual,expected)
    save(RECORD/'physical-preflight.json',dict(passed=True,scalar_digitizer_exact=True,closure_error_mV=first[4]['component_sum_max_error_mV']))
    result=[None]*len(cases);result[0]=first
    with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize_physical,initargs=(local,len(records))) as pool:
        futures=[pool.submit(physical_one,(i,c,OUTPUT)) for i,c in enumerate(cases) if i]
        for n,future in enumerate(as_completed(futures),2):
            i,r=future.result();result[i]=r
            if n%100==0:print('CONFIRM_PHYSICAL',n,'/',len(cases),flush=True)
    for i,r in enumerate(records):
        with np.load(PROJECT/r['waveform']) as f:grid=f['time_ns'];y=f['analog_mV']
        with np.load(PROJECT/r['truth']) as f:q=float(f['charge_pe'][0]);first=float(f['first_active_ns'][0])
        pre=float(np.random.default_rng(base+80000+i).uniform(96.,256.));start=int(np.floor((5000.-pre-grid[0])/2.));origin=float(grid[start]-.25)
        case=dict(kind='single',condition=next(k for k,v in energies.items() if r['energy_MeV']==v),gap_ns=None,replicate=r['replicate'],
            parents=[i],seed=r['seed'],cluster=f'parent-{i}',crop_first_index=start,origin_ns=origin,crop='Nominal source trigger plus independent jitter; not true onset')
        result.append((y[start:start+1024],np.array([q,0.]),np.array([first-origin,0.]),np.zeros((2,1024)),case))
    assert len(result)==p['expected_frames']
    np.savez_compressed(OUTPUT/'observations.npz',voltage_mV=np.asarray([r[0] for r in result],np.float32))
    np.savez_compressed(OUTPUT/'labels.truth.npz',charge_pe=[r[1] for r in result],times_ns=[r[2] for r in result],components_mV=np.asarray([r[3] for r in result],np.float32))
    save(OUTPUT/'labels.cases.json',[r[4] for r in result])
    save(RECORD/'physical-validation.json',dict(passed=True,frames=len(result),native_parents=len(records),workers=workers,
        closure_max_mV=max(r[4].get('component_sum_max_error_mV',0.) for r in result),observation_sha256=digest(OUTPUT/'observations.npz')))
    print('CONFIRM_PHYSICAL_COMPLETE',len(result),flush=True)


def initialize(method):
    global RESPONSE,NOISE,METHOD
    from .alpha_pair_conditional import ConditionalResponse
    from .alpha_pair_joint_likelihood import ConservingResponse,JointNoise
    from .alpha_pair_conditional_study import forbid_labels
    METHOD=method;NOISE=None
    if method=='joint':
        RESPONSE=ConservingResponse.load(EXPERIMENT/'joint-response.npz');NOISE=JointNoise.load(EXPERIMENT/'joint_conserved_nll-noise.npz')
    else:RESPONSE=ConditionalResponse.load(EXPERIMENT/'response.npz')
    sys.addaudithook(forbid_labels)


def one(task):
    i,y=task;start=time.monotonic()
    if METHOD=='joint':
        from .alpha_pair_joint_likelihood import waveform_candidates
        candidates=waveform_candidates(RESPONSE,NOISE,y)
    else:
        from .alpha_pair_null_gate import waveform_candidates
        candidates=waveform_candidates(RESPONSE,y,pilot=METHOD)
    return dict(index=i,candidates=candidates,seconds=time.monotonic()-start)


def infer(method,workers=4):
    check_freeze()
    with np.load(OUTPUT/'observations.npz') as f:y=f['voltage_mV']
    path=OUTPUT/f'predictions.{method}.jsonl'
    old=[json.loads(s) for s in path.read_text().splitlines()] if path.exists() else [];done={r['index'] for r in old}
    assert len(old)==len(done);started=time.monotonic()
    with path.open('a') as stream:
        with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize,initargs=(method,)) as pool:
            futures=[pool.submit(one,(i,row)) for i,row in enumerate(y) if i not in done]
            for n,future in enumerate(as_completed(futures),1):
                r=future.result();stream.write(json.dumps(r,allow_nan=False)+'\n');stream.flush()
                if n%50==0:os.fsync(stream.fileno());print('CONFIRM_INFER',method,n,'/',len(futures),flush=True)
        os.fsync(stream.fileno())
    save(RECORD/f'{method}-inference.json',dict(passed=True,frames=len(y),new_frames=len(futures),seconds=time.monotonic()-started,
        workers=workers,worker_labels_blocked=True,truth_inputs=[],observation_sha256=digest(OUTPUT/'observations.npz'),predictions_sha256=digest(path)))
    print('CONFIRM_INFERENCE_COMPLETE',method,len(y),flush=True)


def report():
    from .alpha_pair_joint_likelihood import select as joint_select
    from .alpha_pair_null_gate import select as old_select
    from .alpha_pair_residual import choose
    from .alpha_pair_conditional_study import evaluate_fit
    from .alpha_pair_null_report import summarize
    check_freeze();cases=load(OUTPUT/'labels.cases.json')
    with np.load(OUTPUT/'labels.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
    predictions={name:{r['index']:r['candidates'] for r in map(json.loads,(OUTPUT/f'predictions.{name}.jsonl').read_text().splitlines())} for name in ('joint','observed','model')}
    assert all(len(p)==len(cases)==1600 for p in predictions.values())
    gates=load(EXPERIMENT/'gate-calibrations.json');legacy=load(EXPERIMENT/'legacy-gate-calibrations.json')['relative_residual']
    summary=[];all_rows=[];results={}
    for method in ('observed_penalty','model_relative_residual','model_matched_gate','joint_conserved_nll'):
        groups={}
        for i,c in enumerate(cases):
            if method=='observed_penalty':fit=choose(predictions['observed'][i],.001)
            elif method=='model_relative_residual':fit=old_select(predictions['model'][i],legacy)
            elif method=='model_matched_gate':fit=old_select(predictions['model'][i],gates['old_ar'])
            else:fit=joint_select(predictions['joint'][i],gates['joint_conserved_nll'])
            row=dict(evaluate_fit(fit,t[i],q[i],c),index=i,method=method,kind=c['kind'],condition=c['condition'],gap_ns=c['gap_ns'],
                true_q=q[i].tolist(),true_t=t[i].tolist(),fit=fit)
            all_rows.append(row);groups.setdefault((c['kind'],c['condition'],c['gap_ns']),[]).append(row)
        s=[summarize(rr) for rr in groups.values()];summary+=s;controls=[r for r in s if r['kind']!='pair']
        passes=all(r['extra_fraction']<=.05 for r in controls) and all(r['efficiency']>=.95 and abs(r['q_error_p16_p50_p84'][1])<=.05 for r in controls if r['kind']=='single')
        results[method]=dict(weak_mean=float(np.mean([r['efficiency'] for r in s if r['kind']=='pair' and 'weak' in r['condition']])),
            strong_mean=float(np.mean([r['efficiency'] for r in s if r['kind']=='pair' and 'weak' not in r['condition']])),controls_pass=passes,
            controls={r['condition']:dict(extra=r['extra_frames'],success=r['success'],frames=r['frames'],median_q_bias=None if r['q_error_p16_p50_p84'] is None else r['q_error_p16_p50_p84'][1]) for r in controls})
        print('CONFIRM_RESULT',method,json.dumps(results[method]),flush=True)
    primary=results['joint_conserved_nll'];confirmed=(primary['controls_pass'] and primary['weak_mean']>max(results[k]['weak_mean'] for k in ('observed_penalty','model_relative_residual')) and primary['strong_mean']>=results['observed_penalty']['strong_mean']-.02)
    save(OUTPUT/'report.json',dict(results=results,summary=summary,chosen_before_test='joint_conserved_nll',test_requirements_passed=confirmed,
        independent_optical_parents=True,gate_recalibrated=False,interpretation='Descriptive independent confirmation; reused parents across points, finite controls, no population guarantee'))
    with (OUTPUT/'scores.jsonl').open('w') as stream:
        for row in all_rows:stream.write(json.dumps(row)+'\n')
    print('CONFIRMATION_PASSED',confirmed,flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['prepare','native','physical','infer','report'])
    parser.add_argument('--workers',type=int,default=4);parser.add_argument('--method',choices=['joint','observed','model'],default='joint');args=parser.parse_args()
    if args.command=='infer':infer(args.method,args.workers)
    elif args.command in ('native','physical'):globals()[args.command](args.workers)
    else:globals()[args.command]()


if __name__=='__main__':main()
