"""Development of empirical single-pulse covariance, with disjoint gate parents."""
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

PROJECT=Path(os.environ.get('GAGG_EMPIRICAL_PROJECT',Path(__file__).resolve().parents[2])).resolve()
NAME='alpha_pair_empirical_covariance_20261026'
EXPERIMENT=PROJECT/'pileup/experiments'/NAME
OUTPUT=PROJECT/'jupyter/artifacts'/NAME
RECORD=PROJECT/'jupyter/run_records'/NAME
PARENT=PROJECT/'pileup/experiments/alpha_pair_null_gate_20261023'
OLD=PROJECT/'jupyter/artifacts/alpha_pair_null_gate_20261023'
VARIANTS=dict(diagonal=dict(rank=0,fraction=0.),rank16_half=dict(rank=16,fraction=.5),
              rank16_full=dict(rank=16,fraction=.9),rank32_full=dict(rank=32,fraction=.9))


def prepare():
    from .alpha_pair_conditional import ConditionalResponse
    from .alpha_pair_empirical_noise import calibrate
    for p in (EXPERIMENT,OUTPUT,RECORD):p.mkdir(parents=True,exist_ok=False)
    for name in ('response.npz','settings.json','calibration.json','model-null-calibration.json'):
        shutil.copyfile(PARENT/name,EXPERIMENT/name)
    shutil.copytree(PARENT/'source',EXPERIMENT/'source')
    shutil.copytree(PARENT/'bin',EXPERIMENT/'bin')
    records=[]
    for name,part in [('alpha_pair_conditional_20261022',''),('alpha_pair_residual_20261023',''),
                      ('alpha_pair_null_gate_20261023','development')]:
        records+=load(PROJECT/'jupyter/artifacts'/name/part/'native-records.json')['test_single']
    assert len(records)==560 and len({r['seed'] for r in records})==560
    training=[r for r in records if r['replicate']%2==0]
    gate_seeds={r['seed'] for r in records if r['replicate']%2==1}
    assert len(training)==len(gate_seeds)==280
    validation_records=load(OLD/'test/native-records.json')['test_single']
    assert not {r['seed'] for r in records}&{r['seed'] for r in validation_records}
    plan=dict(created=datetime.now().astimezone().isoformat(),variants=VARIANTS,
        covariance_training_parents=280,gate_calibration_parents=280,
        development_validation_parents=320,development_frames=1600,
        covariance_split='Even native replicate; disjoint odd replicates for gate ranks',
        response_unchanged=True,neural_training=False,fixed_taus_ns=[95.,351.],
        relative_covariance_grid_ns=[-64.,2112.,2.],maximum_rank=32,
        covariance='Centered native single residuals; regularized diagonal plus empirical modes',
        interpolation='Linear 2ns alignment; electronic variance floor0.01mV^2; no empirical mean added to template',
        covariance_at_inference='Fixed from observation-only OLS single pilot Q/t for every candidate in a frame',
        gate=dict(score='relative_residual',alpha=.025,edges_pe=[1500.,5000.,12500.],require_birth_penalty=False),
        selection=dict(extra_fraction_max=.05,single_success_min=.95,abs_median_charge_bias_max=.05,
            weak_mean_must_exceed_both_prior_blind_methods=True,strong_max_drop_percentage_points=2.),
        development_data='Previous null-gate independent test is exposed development; not independent evidence',
        control_scope='Control acceptance uses point estimates; not a simultaneous population guarantee',
        scope='Pure alpha, fixed lifetimes, local unknown K<=2; no high-rate or energy-spectrum claim')
    save(EXPERIMENT/'development-plan.json',plan)
    response=ConditionalResponse.load(EXPERIMENT/'response.npz')
    grid=np.arange(-64.,2112.01,2.);residuals=[];charges=[];provenance=[]
    for j,r in enumerate(training):
        wave=PROJECT/r['waveform'];truth=PROJECT/r['truth']
        with np.load(wave) as f:sample_time=f['time_ns'];voltage=f['analog_mV']
        with np.load(truth) as f:t=float(f['first_active_ns'][0]);q=float(f['charge_pe'][0])
        assert grid[0]+t>=sample_time[0] and grid[-1]+t<=sample_time[-1]
        aligned=np.interp(grid+t,sample_time,voltage)-response.baseline
        residuals.append(aligned-q*response.unit(grid,q));charges.append(q)
        provenance.append(dict(seed=r['seed'],root=r['root'],root_sha256=r['root_original_sha256'],
            waveform=r['waveform'],waveform_sha256=digest(wave),truth=r['truth'],truth_sha256=digest(truth)))
        if (j+1)%40==0:print('CALIBRATION_PARENTS',j+1,'/',len(training),flush=True)
    noise=calibrate(grid,np.asarray(residuals),np.asarray(charges),rank=32)
    noise.save(EXPERIMENT/'empirical-noise.npz')
    np.savez_compressed(OUTPUT/'aligned-single-residuals.truth.npz',residuals=residuals,charge_pe=charges,grid=grid)
    save(RECORD/'covariance-calibration.json',dict(parents=provenance,counts=noise.counts.tolist(),
        gate_seeds=sorted(gate_seeds),training_seeds=sorted(r['seed'] for r in training),
        parents_disjoint=True,grid_fully_observed=True,mean_response_changed=False))
    # The held-out single records are already cropped using nominal triggers.
    dev=OLD/'model-development';cases=load(dev/'labels.cases.json')
    gate_indices=[i for i,c in enumerate(cases) if c['kind']=='single' and c['seed'] in gate_seeds]
    assert len(gate_indices)==280 and {cases[i]['seed'] for i in gate_indices}==gate_seeds
    diagnostic=load(PROJECT/'pileup/experiments/alpha_pair_weak_diagnostics_20261026/protocol.json')
    validation_indices=diagnostic['indices'];assert len(validation_indices)==1600
    all_y=[];all_q=[];all_t=[];all_cases=[];old_predictions=[]
    for role,folder,indices in [('calibration',dev,gate_indices),('development',OLD/'test',validation_indices)]:
        source_cases=load(folder/'labels.cases.json')
        with np.load(folder/'observations.npz') as f:y=f['voltage_mV']
        with np.load(folder/'labels.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
        preds={r['index']:r['candidates'] for r in map(json.loads,(folder/'predictions.model.jsonl').read_text().splitlines())}
        for i in indices:
            all_y.append(y[i]);all_q.append(q[i]);all_t.append(t[i])
            all_cases.append(dict(source_cases[i],role=role,source_index=i))
            old_predictions.append(dict(index=len(all_cases)-1,candidates=preds[i]))
    np.savez_compressed(OUTPUT/'observations.npz',voltage_mV=np.asarray(all_y))
    np.savez_compressed(OUTPUT/'labels.truth.npz',charge_pe=all_q,times_ns=all_t)
    save(OUTPUT/'labels.cases.json',all_cases)
    with (OUTPUT/'predictions.old_ar.jsonl').open('w') as stream:
        for row in old_predictions:stream.write(json.dumps(row)+'\n')
    save(RECORD/'data-provenance.json',dict(exposed_development=True,frames=len(all_cases),
        development_indices=validation_indices,gate_indices=gate_indices,
        input_files={str((folder/file).relative_to(PROJECT)):digest(folder/file)
            for folder in (dev,OLD/'test') for file in ('observations.npz','labels.truth.npz','labels.cases.json','predictions.model.jsonl')}))
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os, sys
folder=Path(__file__).resolve().parent
for name in ('EMPIRICAL','WEAK','NULL','CONDITIONAL','RESIDUAL','RESOLUTION','PAIR','COMPONENT'):
    os.environ['GAGG_'+name+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'development-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_empirical_study import main
    main()
''')
    for namespace in ('analysis','pileup/analysis'):
        destination=EXPERIMENT/'development-source'/namespace;destination.mkdir(parents=True)
        for p in (PROJECT/namespace).glob('*.py'):shutil.copyfile(p,destination/p.name)
    if (PROJECT/'pileup/__init__.py').exists():
        shutil.copyfile(PROJECT/'pileup/__init__.py',EXPERIMENT/'development-source/pileup/__init__.py')
    frozen=[p for name in ('source','bin','development-source') for p in (EXPERIMENT/name).rglob('*') if p.is_file()]
    frozen += [EXPERIMENT/name for name in ('response.npz','settings.json','calibration.json',
        'model-null-calibration.json','empirical-noise.npz','development-plan.json','run.py')]
    save(EXPERIMENT/'development-freeze.json',dict(created=datetime.now().astimezone().isoformat(),
        files={str(p.relative_to(EXPERIMENT)):digest(p) for p in frozen},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        independent_test_not_started=True))
    save(RECORD/'environment.json',dict(python=sys.version,executable=sys.executable,numpy=np.__version__,workers=4,
        git_changes=subprocess.check_output(['git','--no-optional-locks','-c','core.filemode=false','status','--short','--',
            'pileup/analysis/alpha_pair_empirical_noise.py','pileup/analysis/alpha_pair_empirical_study.py',
            'pileup/analysis/test_alpha_pair_empirical_noise.py'],cwd=PROJECT,text=True)))
    print('EMPIRICAL_PREPARED',len(all_cases),'inference frames',len(frozen),'frozen files',flush=True)


def initialize(variant):
    global RESPONSE,NOISE,OPTIONS
    from .alpha_pair_conditional import ConditionalResponse
    from .alpha_pair_empirical_noise import EmpiricalNoise
    from .alpha_pair_conditional_study import forbid_labels
    RESPONSE=ConditionalResponse.load(EXPERIMENT/'response.npz')
    NOISE=EmpiricalNoise.load(EXPERIMENT/'empirical-noise.npz');OPTIONS=VARIANTS[variant]
    sys.addaudithook(forbid_labels)


def one(task):
    from .alpha_pair_empirical_noise import waveform_candidates
    i,y=task;started=time.monotonic()
    return dict(index=i,candidates=waveform_candidates(RESPONSE,NOISE,y,**OPTIONS),seconds=time.monotonic()-started)


def infer(variant,workers=4,pilot=False):
    for name,sha in load(EXPERIMENT/'development-freeze.json')['files'].items():
        assert digest(EXPERIMENT/name)==sha,name
    with np.load(OUTPUT/'observations.npz') as f:y=f['voltage_mV']
    indices=list(range(len(y)))
    if pilot:
        cases=load(OUTPUT/'labels.cases.json')
        indices=[i for i,c in enumerate(cases) if c['replicate']<=3 and
                 (c['role']=='calibration' or c['replicate']==0)]
    path=OUTPUT/('predictions.'+variant+'.jsonl')
    previous=[json.loads(s) for s in path.read_text().splitlines()] if path.exists() else []
    done={r['index'] for r in previous};assert len(done)==len(previous)
    started=time.monotonic();count=0
    with path.open('a') as stream:
        with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize,initargs=(variant,)) as pool:
            futures=[pool.submit(one,(i,y[i])) for i in indices if i not in done]
            for future in as_completed(futures):
                result=future.result();stream.write(json.dumps(result,allow_nan=False)+'\n');stream.flush();count+=1
                if count%50==0:
                    os.fsync(stream.fileno());print('EMPIRICAL_INFER',variant,count,'/',len(futures),flush=True)
        os.fsync(stream.fileno())
    save(RECORD/(variant+('.pilot' if pilot else '')+'-inference.json'),dict(completed=True,requested=len(indices),new_frames=count,
        seconds=time.monotonic()-started,workers=workers,truth_inputs=[],worker_labels_blocked=True,
        waveform_sha256=digest(OUTPUT/'observations.npz'),prediction_sha256=digest(path)))
    print('EMPIRICAL_INFERENCE_COMPLETE',variant,len(indices),flush=True)


def select():
    from .alpha_pair_null_gate import calibrate as calibrate_gate,select as gate_select
    from .alpha_pair_conditional_study import evaluate_fit
    from .alpha_pair_null_report import summarize
    cases=load(OUTPUT/'labels.cases.json')
    with np.load(OUTPUT/'labels.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
    previous=load(OLD/'test/report.json')['summary']
    groups_weak=lambda rows: [r['efficiency'] for r in rows if r['kind']=='pair' and 'weak' in r['condition']]
    groups_strong=lambda rows: [r['efficiency'] for r in rows if r['kind']=='pair' and 'weak' not in r['condition']]
    baselines={}
    for method in ('observed_penalty','model_relative_residual'):
        rows=[r for r in previous if r['method']==method and (r['kind']!='pair' or r['gap_ns'] in (8.,40.,200.))]
        baselines[method]=dict(weak_mean=float(np.mean(groups_weak(rows))),strong_mean=float(np.mean(groups_strong(rows))))
    results={};calibrations={};all_summary=[]
    for name in [*VARIANTS,'old_ar']:
        rows=[json.loads(s) for s in (OUTPUT/('predictions.'+name+'.jsonl')).read_text().splitlines()]
        predictions={r['index']:r['candidates'] for r in rows}
        assert len(predictions)==len(cases)==1880,name
        singles=[predictions[i] for i,c in enumerate(cases) if c['role']=='calibration']
        calibration=calibrate_gate(singles,edges=[1500.,5000.,12500.],mode='relative_residual',alpha=.025,
                                   penalty=.001,require_birth_penalty=False)
        calibrations[name]=calibration;groups={}
        for i,c in enumerate(cases):
            if c['role']!='development':continue
            fit=gate_select(predictions[i],calibration)
            r=dict(evaluate_fit(fit,t[i],q[i],c),method=name,kind=c['kind'],condition=c['condition'],gap_ns=c['gap_ns'],
                true_q=q[i].tolist(),true_t=t[i].tolist())
            groups.setdefault((c['kind'],c['condition'],c['gap_ns']),[]).append(r)
        summary=[summarize(rows) for rows in groups.values()];all_summary+=summary
        controls=[r for r in summary if r['kind']!='pair']
        controls_pass=all(r['extra_fraction']<=.05 for r in controls) and all(
            r['efficiency']>=.95 and abs(r['q_error_p16_p50_p84'][1])<=.05 for r in controls if r['kind']=='single')
        weak_mean=float(np.mean(groups_weak(summary)));strong_mean=float(np.mean(groups_strong(summary)))
        eligible=(controls_pass and weak_mean>max(b['weak_mean'] for b in baselines.values())
                  and strong_mean>=baselines['observed_penalty']['strong_mean']-.02 and name in VARIANTS)
        results[name]=dict(controls_pass=controls_pass,eligible=eligible,weak_mean=weak_mean,strong_mean=strong_mean,
            controls={r['condition']:dict(extra=r['extra_frames'],success=r['success'],frames=r['frames'],
                median_q_bias=None if r['q_error_p16_p50_p84'] is None else r['q_error_p16_p50_p84'][1]) for r in controls})
        print('EMPIRICAL_SELECTION',name,json.dumps(results[name]),flush=True)
    eligible=[name for name in VARIANTS if results[name]['eligible']]
    chosen=max(eligible,key=lambda name:results[name]['weak_mean']+results[name]['strong_mean']) if eligible else None
    save(EXPERIMENT/'development-selection.json',dict(chosen=chosen,results=results,baselines=baselines,
        independent_test_not_started=True,calibration_role='280 disjoint odd-replicate native parents',
        validation_role='320 previously exposed parents; development only'))
    save(EXPERIMENT/'gate-calibrations.json',calibrations)
    save(OUTPUT/'development-report.json',dict(chosen=chosen,summary=all_summary,results=results,baselines=baselines,
        independent_performance=False,paired_parent_reuse_across_conditions=True))
    print('EMPIRICAL_CHOSEN',chosen,flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['prepare','infer','select'])
    parser.add_argument('--variant',choices=list(VARIANTS),default='rank16_full')
    parser.add_argument('--workers',type=int,default=4);parser.add_argument('--pilot',action='store_true')
    args=parser.parse_args()
    if args.command=='prepare':prepare()
    elif args.command=='infer':infer(args.variant,args.workers,args.pilot)
    else:select()


if __name__=='__main__':main()
