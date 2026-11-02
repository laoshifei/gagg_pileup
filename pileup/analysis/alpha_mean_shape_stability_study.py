"""Frozen, training-pool-only cross-fit and refit-bootstrap diagnostic."""
from __future__ import annotations

import argparse
import ast
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

from .alpha_mean_shape_correction import MeanShapeCorrection,positivity_scale
from .alpha_mean_shape_stability import (mother_folds,bootstrap_mothers,fit_update,
    many_voltage,cached_positivity_scale,crossfit,projected_variance)
from .alpha_pair_energy_calibrated import CalibratedEnergyResponse
from .alpha_pair_recovery_energy import RecoveryEnergyResponse
from .alpha_pair_spatial_recovery import SpatialRecovery

NAME='alpha_mean_shape_stability_20261102'
ROOT=Path(os.environ.get('GAGG_SHAPE_STABILITY_ROOT',Path(__file__).resolve().parents[2])).resolve()
EXP=ROOT/'pileup/experiments'/NAME
ART=ROOT/'jupyter/artifacts'/NAME
REC=ROOT/'jupyter/run_records'/NAME
PARENT=ROOT/'pileup/experiments/alpha_mean_shape_20261102'
PARENT_ART=ROOT/'jupyter/artifacts/alpha_mean_shape_20261102'
WORKERS=8
BOOTSTRAPS=1000
PARTITIONS=10


def read(p):return json.loads(Path(p).read_text())


def save(p,value):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    temporary=p.with_suffix(p.suffix+'.partial')
    temporary.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n');temporary.replace(p)


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        while block:=f.read(1024*1024):h.update(block)
    return h.hexdigest()


def load_arrays(p):
    with np.load(p) as f:return {k:f[k].copy() for k in f.files}


def baseline():
    if not hasattr(np,'trapezoid'):np.trapezoid=np.trapz
    anchor=CalibratedEnergyResponse.load(RecoveryEnergyResponse.load(EXP),EXP/'mean-correction.npz')
    return SpatialRecovery(anchor,read(EXP/'sensor-config.json'),EXP/'bin/age.so',
        load_arrays(EXP/'charge-map-spatial12.npz'),load_arrays(EXP/'spatial-profile.npz')['profile_12'])


def snapshot():
    target=EXP/'source/pileup/analysis';target.mkdir(parents=True,exist_ok=True)
    (target.parent/'__init__.py').write_text('');(target/'__init__.py').write_text('')
    pending=['alpha_mean_shape_stability_study','test_alpha_mean_shape_stability'];copied={}
    while pending:
        name=pending.pop()
        if name in copied:continue
        old=PARENT/'source/pileup/analysis'/f'{name}.py'
        source=old if old.exists() else ROOT/'pileup/analysis'/f'{name}.py'
        shutil.copyfile(source,target/source.name);copied[name]=dict(source=str(source.relative_to(ROOT)),sha256=sha(source))
        for node in ast.walk(ast.parse(source.read_text())):
            if isinstance(node,ast.ImportFrom) and node.level==1:
                pending.extend([node.module] if node.module else [a.name for a in node.names])
    for name in ('sipm_age_mean.cpp','sipm_age_mean_jacobian.cpp'):
        shutil.copyfile(PARENT/'source/pileup/analysis'/name,target/name)
    (EXP/'bin').mkdir(exist_ok=True)
    command=['g++','-O3','-std=c++17','-fPIC','-shared',str((target/'sipm_age_mean.cpp').relative_to(ROOT)),
             str((target/'sipm_age_mean_jacobian.cpp').relative_to(ROOT)),'-o',str((EXP/'bin/age.so').relative_to(ROOT))]
    subprocess.run(command,cwd=ROOT,check=True)
    save(EXP/'build.json',dict(command=command,binary_sha256=sha(EXP/'bin/age.so'),
         compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0]))
    save(EXP/'source-origins.json',copied)
    (EXP/'run.py').write_text("from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\n"
        "os.environ['GAGG_SHAPE_STABILITY_ROOT']=str(f.parents[2])\nos.environ['OPENBLAS_NUM_THREADS']='1'\nos.environ['OMP_NUM_THREADS']='1'\n"
        "sys.path.insert(0,str(f/'source'))\nif __name__=='__main__':\n    from pileup.analysis.alpha_mean_shape_stability_study import main\n    main()\n")


def prepare(resume=False):
    if resume:
        assert EXP.exists() and not (EXP/'design-freeze.json').exists()
        assert not (EXP/'results').exists() and not (EXP/'pilot-gate.json').exists()
    for p in (EXP,ART,REC):p.mkdir(parents=True,exist_ok=resume)
    for p in ('original-metadata','original-runinfo'):(EXP/p).mkdir(exist_ok=resume)
    (ART/'figures').mkdir(exist_ok=resume)
    for name in ('template.npz','energy-model.json','response-parameters.json','mean-correction.npz',
                 'sensor-config.json','charge-map-spatial12.npz','spatial-profile.npz','nominal-influences.npz','literature-evidence.json'):
        shutil.copyfile(PARENT/name,EXP/name)
    shutil.copyfile(PARENT/'shape-correction.npz',EXP/'original-shape.npz')
    mothers=[r for r in read(PARENT/'mothers.json') if r['role']=='train']
    assert len(mothers)==128 and [r['id'] for r in mothers]==list(range(128))
    save(EXP/'mothers.json',mothers)
    train_plan=read(PARENT/'train-plan.json');save(EXP/'original-digitization-plan.json',train_plan)
    for r in mothers:
        for folder in ('original-metadata','original-runinfo'):
            name=f"{r['id']:04d}.json";shutil.copyfile(PARENT/folder/name,EXP/folder/name)
    groups=np.array([r['group'] for r in mothers]);identities=[r['root'] for r in mothers]
    folds=np.array([mother_folds(identities,groups,737000+i) for i in range(PARTITIONS)])
    save(EXP/'folds.json',dict(seeds=list(range(737000,737000+PARTITIONS)),assignments=folds.tolist()))
    save(EXP/'protocol.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),
        base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        git_status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),
        scope='Conditional stability of incremental isolated mean correction;truth-labelled diagnostics,not energy/efficiency validation',
        data='Only128 previous mean/profile-training mothers,32per source group;reuse1024 digitizations (8each) at phase0. No64 confirmation mothers or32noise-fold strong mothers used in fitting/scoring',
        exposure='Old baseline/prior were already trained on these mothers;only the new correction is cross-fitted. Method chosen afterL93 failure;not independent confirmation',
        recipe='ExactlyL93:4ns smoothing,2ns causal guard,8ns rise,160ns end taper to1600ns,support shrink,zero complete area,logq0 interpolation,global positivity scale. No hyperparameter search',
        cross_validation='4folds stratified by4source groups,8mothers/group/fold;96train/32heldout. Ten fixed path-hash partitions,first primary;all8repeats of a mother stay together',
        bootstrap='1000 whole-mother resamples within each primary group/fold stratum;duplicates retain fold. Refit all4corrections each resample. Percentile intervals are conditional sensitivity estimates,not guaranteed population coverage',
        bootstrap_seed_base=737100000,partitions=PARTITIONS,bootstraps=BOOTSTRAPS,
        metrics='Strong-mother8repeat residual along old fixed nominal8/40/200ns weak-logE operators;mean and per-mother squared projection. Also4group mean-waveform and individual-mother waveform MSE',
        uncertainty='Refit bootstrap includes variation in sampled mother means and learned correction,conditional on fixed old baseline,operator,ADC phase and finite8repeat records. No independent within-mother noise resampling or microphysical systematics',
        decision='Development stability only:8/40ns primary projection MSEratio<1 and refit percentile upper<1;atleast8of10partitions have ratio<1;all4primary group-mean waveform MSEratios<=1.1;allobserved scales>=.5. Passing only motivates new isolated confirmation,never direct pair expansion',
        variance='ReportW=mean within-mother sample variance,B=var(mother means)-W/8 without clipping. Planning normal95%halfwidth10 usesceil(1.96^2*var(mother means)/100);infinite-repeat floor clipsBonlyforplanning',
        repeat_count=8,formal_new_digitizations=0,new_ROOT=0,energy_fits=0,
        nominal_true_energy_MeV=[7.5,.25],projection_units='100delta ln E,not relative energy percent',
        phases_ns=[0.],positivity_phases_ns=[0.,.5,1.,1.5],positivity_energy_points=81,
        workers=WORKERS,BLAS_threads=1,pilot_bootstraps=8,cost_limit_seconds=900.,
        literature='ReuseL91/Fixsen methodology boundaries andL92/L93project evidence;resampling design is project analysis,not a cited detector-performance claim',
        default_changed=False,remote_upload_authorized=False,resumed_before_freeze=resume))
    for name in ('README.md','docs/literature-use.md','jupyter/codex_0007_SiPM_Energy_Calibration.ipynb',
                 'jupyter/Research_Plan_20260921.ipynb','jupyter/Run_Catalog.ipynb'):
        p=ROOT/name;shutil.copyfile(p,REC/(p.name+'.before'))
    snapshot()
    # Verify the exact old training records, then build a self-contained input.
    old_inputs=read(PARENT/'train-input-freeze.json');old_cal=read(PARENT/'model-freeze.json')
    data=load_arrays(PARENT_ART/'calibration.npz');repeats=[];input_paths=[PARENT_ART/'calibration.npz']
    assert sha(PARENT_ART/'calibration.npz')==old_cal[str((PARENT_ART/'calibration.npz').relative_to(ROOT))]
    for i,r in enumerate(mothers):
        p=PARENT_ART/'truth/train'/f'{i:04d}.npz'
        assert sha(p)==old_inputs[str(p.relative_to(ROOT))]
        a=load_arrays(p)['components_mV'][:,0]
        np.testing.assert_array_equal(a.mean(0),data['means_mV'][i]);repeats.append(a);input_paths.append(p)
        assert sha(ROOT/r['root'])==r['root_original_sha256']
    model=baseline();time_ns=model.readout.adc_centers
    data.update(signal_repeats_mV=np.array(repeats),time_ns=time_ns,operator=load_arrays(EXP/'nominal-influences.npz')['high'][1:],
                energy_MeV=np.array([r['edep_MeV'] for r in mothers]))
    np.testing.assert_array_equal(data['groups'],groups)
    np.savez_compressed(ART/'inputs.npz',**data)
    energy=np.geomspace(*model.energy_bounds,81);phases=np.array([0.,.5,1.,1.5])
    cached=np.array([[model.isolated([128.+phase],[e])[0][0] for e in energy] for phase in phases])
    np.savez_compressed(ART/'positivity-cache.npz',time_ns=time_ns,energy_MeV=energy,
        phases_ns=phases,q0=model.anchor.q_isolated(energy),old_mV=cached,
        nominal_high_q0=model.anchor.q_isolated([7.5])[0])
    save(EXP/'parent-input-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in input_paths})
    import scipy
    save(EXP/'environment.json',dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,
        workers=WORKERS,BLAS_threads=1))
    paths=[p for p in EXP.rglob('*') if p.is_file() and 'bin' not in p.relative_to(EXP).parts]
    paths+=[ART/'inputs.npz',ART/'positivity-cache.npz']
    save(EXP/'design-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})
    print('Frozen128 development mothers;10 grouped partitions and1000 refit bootstraps;no new digitization',flush=True)


def verify():
    paths=read(EXP/'design-freeze.json')
    for p,d in paths.items():assert sha(ROOT/p)==d,p
    return dict(design_files=len(paths))


def initialize(guard=True):
    global DATA,CACHE,FOLDS
    DATA=load_arrays(ART/'inputs.npz');CACHE=load_arrays(ART/'positivity-cache.npz')
    FOLDS=np.array(read(EXP/'folds.json')['assignments'])
    if guard:
        def audit(event,args):
            if event!='open' or not isinstance(args[0],(str,bytes,os.PathLike)):return
            p=Path(os.fsdecode(args[0])).resolve()
            if p.suffix=='.root' or p.is_relative_to(PARENT) or p.is_relative_to(PARENT_ART):
                raise PermissionError('Crossfit worker cannot access previous confirmation data: '+str(p))
        sys.addaudithook(audit)
        for p in (PARENT/'single-result.json',PARENT_ART/'truth/single/0000.npz'):
            try:p.open('rb')
            except PermissionError:pass
            else:raise AssertionError('Old confirmation read guard failed')


def controls():
    verify();initialize(False)
    c=fit_update(DATA,np.arange(128),CACHE);old=MeanShapeCorrection.load(EXP/'original-shape.npz')
    np.testing.assert_array_equal(c.profiles,old.profiles);np.testing.assert_array_equal(c.log_charge_knots,old.log_charge_knots)
    assert c.scale==old.scale
    model=baseline();scalar=positivity_scale(model,c,CACHE['energy_MeV'],CACHE['phases_ns'])
    assert abs(scalar-cached_positivity_scale(c,CACHE))<1e-14
    errors=[]
    for phase in CACHE['phases_ns']:
        expected=np.array([c.voltage(DATA['time_ns'],128.+phase,q) for q in DATA['q0']])
        errors.append(float(np.max(abs(many_voltage(c,DATA['q0'],DATA['time_ns'],128.+phase)-expected))))
    assert max(errors)<1e-10
    save(EXP/'numerical.json',dict(passed=True,full_recipe_exactly_reproduces_L93=True,
        cached_positivity_matches_scalar=True,max_vectorized_difference_mV=max(errors),
        original_model_unchanged=True))
    print('Full128-mother recipe exactly reproduced;cached positivity and interpolation agree',flush=True)


def bootstrap_job(i):
    start=time.monotonic();sample=bootstrap_mothers(DATA['groups'],FOLDS[0],737100000+i)
    summary=crossfit(DATA,CACHE,FOLDS[0],sample=sample)
    return dict(id=i,seconds=time.monotonic()-start,result=summary)


def partition_job(i):
    start=time.monotonic();summary,arrays=crossfit(DATA,CACHE,FOLDS[i],keep_arrays=True)
    return dict(id=i,seconds=time.monotonic()-start,result=summary),arrays


def pilot():
    verify();assert read(EXP/'numerical.json')['passed'] and not (EXP/'pilot-gate.json').exists()
    start=time.monotonic()
    with ProcessPoolExecutor(max_workers=WORKERS,initializer=initialize) as pool:rows=list(pool.map(bootstrap_job,range(8)))
    estimate=max(r['seconds'] for r in rows)*(BOOTSTRAPS+PARTITIONS)/WORKERS*1.5
    save(EXP/'pilot-gate.json',dict(projected_seconds=estimate,limit_seconds=900.,passed=bool(estimate<900.)))
    save(REC/'pilot.json',dict(seconds=time.monotonic()-start,rows=rows,calibrations=32,new_digitizations=0))
    assert estimate<900.,estimate
    print('Pilot passed,projected seconds',estimate,flush=True)


def analyze(label):
    verify();assert read(EXP/'pilot-gate.json')['passed']
    out=EXP/'results'/label;out.mkdir(parents=True,exist_ok=False)
    array_out=ART/'results'/label;array_out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic();partitions=[];bootstraps=[]
    with ProcessPoolExecutor(max_workers=WORKERS,initializer=initialize) as pool:
        for row,arrays in pool.map(partition_job,range(PARTITIONS)):
            partitions.append(row);np.savez_compressed(array_out/f"partition{row['id']:02d}.npz",**arrays)
        for row in pool.map(bootstrap_job,range(BOOTSTRAPS)):
            bootstraps.append(row)
            if len(bootstraps)%100==0:print(label,'refit bootstrap',len(bootstraps),'/',BOOTSTRAPS,flush=True)
    for pilot_row,row in zip(read(REC/'pilot.json')['rows'],bootstraps[:8]):assert pilot_row['result']==row['result']
    summary=summarize(partitions,bootstraps)
    save(out/'summary.json',summary)
    save(out/'partitions.json',[r['result'] for r in partitions])
    save(out/'bootstrap.json',[r['result'] for r in bootstraps])
    save(REC/(label+'.json'),dict(seconds=time.monotonic()-start,calibrations=4*(PARTITIONS+BOOTSTRAPS),
        partitions=partitions,bootstrap_seconds=[r['seconds'] for r in bootstraps],new_digitizations=0,
        energy_fits=0,pilot_equal=True,worker_confirmation_read_guard=True))
    save(out/'output-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in [out/'summary.json',out/'partitions.json',out/'bootstrap.json',*array_out.glob('*.npz')]})
    print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)


def summarize(partitions,bootstraps):
    data=load_arrays(ART/'inputs.npz');primary=partitions[0]['result'];projection=[];failures=[]
    for j,row in enumerate(primary['projection']):
        draws=[r['result']['projection'][j] for r in bootstraps]
        intervals={k:np.percentile([r[k] for r in draws],[2.5,97.5]).tolist()
                   for k in ('old_mean','corrected_mean','mse_ratio')}
        ratios=[r['result']['projection'][j]['mse_ratio'] for r in partitions]
        improvements=sum(x<1 for x in ratios)
        projection.append(dict(row,refit_bootstrap_ci95=intervals,improving_partitions=improvements,
            partition_mse_ratio_range=[min(ratios),max(ratios)]))
        if row['gap_ns'] in (8.,40.) and not (row['mse_ratio']<1 and intervals['mse_ratio'][1]<1 and improvements>=8):
            failures.append(f"gap{row['gap_ns']:g}:insufficient development stability")
    for row in primary['waveform']:
        if row['mean_mse_ratio']>1.1:failures.append(f"group{row['group']}:mean waveform MSE")
    scales=[m['scale'] for r in [*partitions,*bootstraps] for m in r['result']['models']]
    if min(scales)<.5:failures.append('positivity scale below .5')
    variance=projected_variance(data['signal_repeats_mV'][data['groups']==3],
                                data['old_predictions_mV'][data['groups']==3],data['operator'])
    return dict(primary=primary,projection=projection,variance=variance,
        minimum_scale=min(scales),subhalf_scale_count=sum(s<.5 for s in scales),
        partitions=PARTITIONS,refit_bootstraps=BOOTSTRAPS,development_stability_passed=not failures,
        failures=failures,independent_confirmation=False,energy_fits=0,new_digitizations=0,
        next='Only consider new isolated confirmation;no pair promotion' if not failures else
             'Retain failed correction;do not tune on confirmation;separate training uncertainty and mother variation')


def replay():
    verify();assert not (EXP/'replay.json').exists()
    original=EXP/'results/primary'
    for p,d in read(original/'output-freeze.json').items():assert sha(ROOT/p)==d,p
    old_partitions=read(original/'partitions.json');old_bootstrap=read(original/'bootstrap.json')
    directory=ART/'replay';directory.mkdir()
    arrays_checked=0;rows=[]
    with ProcessPoolExecutor(max_workers=WORKERS,initializer=initialize) as pool:
        for row,arrays in pool.map(partition_job,[0,9]):
            i=row['id'];assert row['result']==old_partitions[i]
            with np.load(ART/'results/primary'/f'partition{i:02d}.npz') as f:
                assert set(f.files)==set(arrays)
                for k in f.files:np.testing.assert_array_equal(f[k],arrays[k]);arrays_checked+=1
            np.savez_compressed(directory/f'partition{i:02d}.npz',**arrays);rows.append(row)
        for row in pool.map(bootstrap_job,[8,511,999]):
            assert row['result']==old_bootstrap[row['id']];rows.append(row)
    save(REC/'replay.json',dict(rows=rows,calibrations=20,new_digitizations=0))
    save(EXP/'replay.json',dict(passed=True,partitions=[0,9],bootstrap_ids=[8,511,999],
        exact_arrays=arrays_checked,all_selected_summaries_equal=True,calibrations=20))
    print('Replay passed:',arrays_checked,'exact arrays and5selected summaries',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['prepare','controls','pilot','analyze','replay','verify'])
    p.add_argument('--label',default='primary');p.add_argument('--resume-prepare',action='store_true');args=p.parse_args()
    if not args.label.replace('_','').isalnum():p.error('Simple alphanumeric label required')
    if args.stage=='prepare':prepare(args.resume_prepare)
    elif args.stage=='controls':controls()
    elif args.stage=='pilot':pilot()
    elif args.stage=='analyze':analyze(args.label)
    elif args.stage=='replay':replay()
    else:print(json.dumps(verify()))


if __name__=='__main__':main()
