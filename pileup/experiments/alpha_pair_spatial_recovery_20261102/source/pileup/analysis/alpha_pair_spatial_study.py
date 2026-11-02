"""Frozen development-prior spatial recovery controls; no pair energy fit."""
from __future__ import annotations

import argparse
import ast
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

from . import alpha_pair_causal_controls as single_control
from . import alpha_pair_causal_shared_study as pair_control
from .alpha_pair_causal_controls import compatibility, read, save, sha
from .alpha_pair_causal_recovery import CausalRecovery, DifferentiableAge, light_bin_probabilities
from .alpha_pair_energy_calibrated import CalibratedEnergyResponse
from .alpha_pair_recovery_energy import RecoveryEnergyResponse
from .alpha_pair_spatial_recovery import SpatialAge, SpatialRecovery, symmetric_profile
from .sipm_readout import ReadoutConfig
from .sipm_response import SensorConfig, HC_EV_NM, map_cells, read_photon_hits

NAME='alpha_pair_spatial_recovery_20261102'
ROOT=Path(os.environ.get('GAGG_SPATIAL_ROOT',Path(__file__).resolve().parents[2])).resolve()
EXP=ROOT/'pileup/experiments'/NAME
ART=ROOT/'jupyter/artifacts'/NAME
REC=ROOT/'jupyter/run_records'/NAME
VARIANTS={'spatial6':(6,1.),'spatial12':(12,1.),'spatial12_fine':(12,.5)}
PRIMARY='spatial12'
PHASES=(0.,.5,1.,1.5)


def extract_development(job):
    index,row,configuration=job
    sensor=SensorConfig(**configuration['sensor'])
    path=ROOT/row['root']
    assert sha(path)==row['root_original_sha256'],row['root']
    metadata=read(ROOT/row['metadata_file'])
    assert metadata['config']['sensor']==configuration['sensor']
    assert metadata['config']['readout']==configuration['readout']
    hits=read_photon_hits(path)
    cells=map_cells(hits,sensor)
    active=cells>=0
    assert active.any()
    start=float(hits.time_ns[active].min())
    photons=dict(time_ns=hits.time_ns-start,x_mm=hits.x_mm,y_mm=hits.y_mm,energy_eV=hits.energy_eV)
    dest=ART/'development-photons'/f'{index:04d}.npz'
    np.savez_compressed(dest,**photons)
    save(ART/'development-runinfo'/f'{index:04d}.json',hits.metadata)
    weight=np.interp(HC_EV_NM/hits.energy_eV[active],sensor.pde_wavelength_nm,sensor.pde_probabilities[-1])
    counts=np.bincount(cells[active],weights=weight,minlength=sensor.n_cells)
    probability=counts/counts.sum()
    return dict(row,profile_id=index,photon_sha256=sha(dest),active_photons=int(active.sum()),
                full_pde_equivalent_light=float(weight.sum()),first_active_ns=start),probability


def snapshot():
    source=EXP/'source/pileup/analysis'
    source.mkdir(parents=True)
    (source.parent/'__init__.py').write_text('')
    (source/'__init__.py').write_text('')
    pending=['alpha_pair_spatial_study'];copied=set()
    while pending:
        name=pending.pop()
        if name in copied:continue
        p=ROOT/'pileup/analysis'/(name+'.py')
        tree=ast.parse(p.read_text());shutil.copyfile(p,source/p.name);copied.add(name)
        for node in ast.walk(tree):
            if isinstance(node,ast.ImportFrom) and node.level==1:
                pending.extend([node.module] if node.module else [item.name for item in node.names])
    for name in ('sipm_age_mean.cpp','sipm_age_mean_jacobian.cpp'):
        shutil.copyfile(ROOT/'pileup/analysis'/name,source/name)
    (EXP/'run.py').write_text("from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\n"
        "os.environ['GAGG_SPATIAL_ROOT']=str(f.parents[2])\n"
        "os.environ['OPENBLAS_NUM_THREADS']='1'\nos.environ['OMP_NUM_THREADS']='1'\n"
        "sys.path.insert(0,str(f/'source'))\n"
        "if __name__=='__main__':\n    from pileup.analysis.alpha_pair_spatial_study import main\n    main()\n")
    (EXP/'bin').mkdir()
    command=['g++','-O3','-std=c++17','-fPIC','-shared',
        str((source/'sipm_age_mean.cpp').relative_to(ROOT)),str((source/'sipm_age_mean_jacobian.cpp').relative_to(ROOT)),
        '-o',str((EXP/'bin/age.so').relative_to(ROOT))]
    subprocess.run(command,cwd=ROOT,check=True)
    save(EXP/'build.json',dict(command=command,compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0],binary_sha256=sha(EXP/'bin/age.so')))


def prepare():
    compatibility()
    for p in (EXP,ART,REC):p.mkdir(parents=True,exist_ok=False)
    for folder in ('development-photons','development-runinfo','photons'):(ART/folder).mkdir()
    previous=ROOT/'pileup/experiments/alpha_pair_causal_shared_20261101'
    for name in ('template.npz','energy-model.json','response-parameters.json','mean-correction.npz',
                 'noise-calibration.json','sensor-config.json','charge-map-1.0.npz'):
        shutil.copyfile(previous/name,EXP/name)
    config=read(EXP/'sensor-config.json')
    pool=read(ROOT/'pileup/experiments/alpha_pair_energy_calibrated_20261030/mother-provenance.json')['development']
    development=[]
    for source_energy in (1.625,2.5,5.486,8.):
        rows=sorted((r for r in pool if r['fold']==0 and r['source_energy_MeV']==source_energy),key=lambda r:r['root'])[:32]
        assert len(rows)==32
        development.extend(rows)
    prior=read(ROOT/'pileup/experiments/alpha_pair_causal_controls_20261101/mother-provenance.json')
    previous_roots={r['source']['root'] for r in prior}
    history=ROOT/'jupyter/artifacts/sipm_energy_history_20261029/inputs'
    parents=read(history/'parents.json')['records']
    cases=read(history/'cases.json')
    with np.load(history/'labels.npz') as f:energies=f['energy_MeV']
    truth={c['parents'][0]:float(energies[i,0]) for i,c in enumerate(cases) if c['kind']=='single'}
    evaluation=[];single_tasks=[];pair_tasks=[]
    for group,condition in enumerate(('weak','low','mid','high')):
        for block in (5,15,25,35):
            parent=group*80+block
            p=ROOT/'jupyter/artifacts/alpha_pair_joint_confirmation_20261026/photons'/f'{parent:04d}.npz'
            shutil.copyfile(p,ART/'photons'/p.name)
            evaluation.append(dict(parent=parent,block=block,condition=condition,energy_dep_MeV=truth[parent],
                                   source=parents[parent],photon_sha256=sha(p)))
            for phase in PHASES:
                single_tasks.append(dict(id=len(single_tasks),kind='single',parent=parent,block=block,condition=condition,
                    phase_ns=phase,time_ns=128.+phase,repeats=8,seed_base=715000000+parent*100,window_ns=4096.))
    for block in (5,15,25,35):
        for condition,gap in (('high_weak',8.),('high_weak',40.),('high_weak',200.),('high_weak',1200.),('weak_high',40.)):
            ps=[block+240,block] if condition=='high_weak' else [block,block+240]
            pair_tasks.append(dict(id=len(pair_tasks),kind='pair',block=block,condition=condition,gap_ns=gap,parents=ps,
                                   times_ns=[128.,128.+gap],repeats=8,seed_base=715100000+block*1000,window_ns=4096.))
    dev_roots={r['root'] for r in development}
    eval_roots={r['source']['root'] for r in evaluation}
    assert len(dev_roots)==128 and len(eval_roots)==16
    assert not dev_roots & eval_roots and not (dev_roots|eval_roots)&previous_roots
    save(EXP/'mother-split.json',dict(development=sorted(dev_roots),evaluation=sorted(eval_roots),
         prior_controls=sorted(previous_roots),disjoint=True,evaluation_historically_exposed=True))
    save(EXP/'evaluation-mothers.json',evaluation)
    save(EXP/'single-plan.json',single_tasks);save(EXP/'pair-plan.json',pair_tasks)
    save(EXP/'protocol.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        git_status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),
        scope='Pure alpha, known-input response diagnostics; single E fits know t, pairs use true E/t for forward scoring only. No pair energy fit or blind recovery claim.',
        development='128 historical fold0 mothers,32 per source-energy group, sorted by ROOT path before outputs. Equal mother weighting of PDE-weighted cell probabilities; no evaluation positions used to estimate profile.',
        evaluation='16 different previously exposed20261026 mothers, four per energy group, block5/15/25/35; disjoint from prior16 causal controls and from128 profile mothers. Retrospective parameter holdout, not fresh optical confirmation.',
        spatial='Energy-independent D4-symmetric illumination prior;6x6 and12x12 exact cell tiles. Equal-illumination tiles combined without altering their independent age equations. No effective-cell or strength fit.',
        approximation='Same95/351ns fraction0.3 light shape,29ns recovery,14400cells and regionwise prompt-XT mean field; no explicit neighbour transport between regions or dark occupancy in expectation.',
        calibration='81 charge-to-light knots250..40000full-PDE equivalents integrated6144ns; invert the SAME spatial operator used for pairs, target inherited isolated Q(E_dep).',
        primary=PRIMARY,variants=VARIANTS,primary_dt_ns=1.,secondary_dt_ns=.5,
        single='64 parent/phase controls x8 repeats=512 digitizations. Old isolated-voltage anchor retained exactly by construction; mean-E fit from voltage only and frozen metric. Add raw truth-time spatial forward single response to exercise the spatial operator.',
        pairs='20cases x8repeats x(shared+2isolated) x(recovery on/off)=960 digitizations, only after single gates pass.',
        raw_diagnostic='Actual accepted photon times and full-PDE weights, but frozen development spatial prior; diagnostic truth inputs, not reconstruction regressors.',
        gates=dict(single_group_mean_energy_percent=5.,single_phase_span_pp=1.,raw_single_group_area_error_percent=5.,
                   shared_weak_area_error_pp=5.,recovery_off_abs_pp=3.,spatial_grid_change_pp=1.,time_grid_change_pp=1.,
                   pre_event_mV=1e-10,numerical_uniform_abs=1e-9,gradient_relative=1e-5,
                   stage_projected_seconds=900.),
        decision='No pair expansion if single gate fails. If pair response gates pass, plan separate E-fit study; do not claim new inference from response agreement.',
        inference_boundary='Primary component predictor receives E/t and fixed profiles/models only. Evaluation photon truth is confined to digitization and separately labeled raw diagnostics.',
        workers=4,BLAS_threads=1,no_new_ROOT=True,default_changed=False,
        literature='Reuse project L87 limitations; configured shared library unavailable, no new external paper assertion.'))
    snapshot()
    start=time.monotonic();records=[];probabilities=[]
    with ProcessPoolExecutor(max_workers=4) as executor:
        jobs=[(i,row,config) for i,row in enumerate(development)]
        for row,p in executor.map(extract_development,jobs):
            records.append(row);probabilities.append(p)
            if len(records)%32==0:print('profile inputs',len(records),'/128',flush=True)
    probability=np.mean(probabilities,axis=0).reshape(120,120)
    profiles={f'profile_{n}':symmetric_profile(probability,n) for n in (6,12)}
    np.savez_compressed(EXP/'spatial-profile.npz',cell_probability=probability,**profiles)
    np.savez_compressed(ART/'development-cell-probabilities.npz',probabilities=np.asarray(probabilities))
    save(EXP/'development-mothers.json',records)
    save(EXP/'profile-summary.json',dict(mothers=128,weighting='Equal per mother;32 per source-energy group; full-recovery PDE weighted within each mother',
         elapsed_seconds=time.monotonic()-start,cell_overlap=float(14400*(probability**2).sum()),
         coarse_probability_sums={k:float(v.sum()) for k,v in profiles.items()},
         grouped_states={k:int(n*(n+2)/8) for k,(n,dt) in VARIANTS.items()}))
    import scipy
    save(EXP/'environment.json',dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,numpy_trapezoid_adapter=compatibility()))
    save(EXP/'source-config-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in EXP.rglob('*') if p.is_file() and 'bin' not in p.relative_to(EXP).parts})
    save(EXP/'input-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in ART.rglob('*') if p.is_file()})
    print('Prepared128 development/16 evaluation mothers and frozen spatial prior',flush=True)


def verify(models=False):
    count=0
    for name in ('source-config-freeze.json','input-freeze.json')+ (('model-freeze.json',) if models else ()):
        for path,expected in read(EXP/name).items():
            assert sha(ROOT/path)==expected,path
            count+=1
    return count


def initialize(models=True):
    compatibility()
    global CONFIG,ANCHOR,PROFILES,MODELS,DONORS
    CONFIG=read(EXP/'sensor-config.json')
    ANCHOR=CalibratedEnergyResponse.load(RecoveryEnergyResponse.load(EXP),EXP/'mean-correction.npz')
    with np.load(EXP/'spatial-profile.npz') as f:PROFILES={k:f[k].copy() for k in f.files}
    with np.load(EXP/'charge-map-1.0.npz') as f:mapping={k:f[k].copy() for k in f.files}
    MODELS={'uniform':CausalRecovery(ANCHOR,CONFIG,EXP/'bin/age.so',mapping)}
    if models:
        for name,(n,dt) in VARIANTS.items():
            with np.load(EXP/f'charge-map-{name}.npz') as f:mapping={k:f[k].copy() for k in f.files}
            MODELS[name]=SpatialRecovery(ANCHOR,CONFIG,EXP/'bin/age.so',mapping,PROFILES[f'profile_{n}'])
    DONORS={}
    for p in (ART/'photons').glob('*.npz'):
        with np.load(p) as f:DONORS[int(p.stem)]={k:f[k].copy() for k in f.files}
    sensor=SensorConfig(**CONFIG['sensor']);readout=ReadoutConfig(**CONFIG['readout'])
    single_control.SENSOR=sensor;single_control.READOUT=readout;single_control.ANCHOR=ANCHOR
    single_control.MODEL=MODELS.get(PRIMARY,MODELS['uniform']);single_control.NOISE=read(EXP/'noise-calibration.json')
    single_control.DONORS=DONORS
    pair_control.SENSOR=sensor;pair_control.READOUT=readout;pair_control.DONORS=DONORS


def calibration_initialize():
    compatibility()
    global CALIBRATION_AGES,CALIBRATION_SHAPES
    sensor=read(EXP/'sensor-config.json')['sensor']
    with np.load(EXP/'spatial-profile.npz') as f:profiles={k:f[k].copy() for k in f.files}
    CALIBRATION_AGES={name:SpatialAge(sensor,EXP/'bin/age.so',profiles[f'profile_{n}'],dt) for name,(n,dt) in VARIANTS.items()}
    CALIBRATION_SHAPES={name:light_bin_probabilities([0.],dt,6144.,.3)[0] for name,(n,dt) in VARIANTS.items()}


def map_point(job):
    name,amplitude=job;start=time.monotonic()
    charge=float(CALIBRATION_AGES[name].forward(amplitude*CALIBRATION_SHAPES[name]).sum())
    return dict(variant=name,light_pe=amplitude,charge_pe=charge,seconds=time.monotonic()-start)


def calibrate():
    verify();assert not (EXP/'model-freeze.json').exists()
    pilot_jobs=[(name,a) for name in VARIANTS for a in (1000.,20000.)]
    with ProcessPoolExecutor(max_workers=4,initializer=calibration_initialize) as executor:
        pilot=list(executor.map(map_point,pilot_jobs))
    estimate=sum(np.mean([r['seconds'] for r in pilot if r['variant']==name])*81 for name in VARIANTS)/4*1.5
    save(REC/'mapping-pilot.json',dict(rows=pilot,projected_seconds=float(estimate),limit_seconds=900.,passed=estimate<900.))
    assert estimate<900.,estimate
    light=np.geomspace(250.,40000.,81);start=time.monotonic()
    rows=[]
    with ProcessPoolExecutor(max_workers=4,initializer=calibration_initialize) as executor:
        jobs=[(name,float(a)) for name in VARIANTS for a in light]
        for row in executor.map(map_point,jobs):
            rows.append(row)
            if len(rows)%27==0:print('spatial charge map',len(rows),'/243',round(time.monotonic()-start,1),flush=True)
    anchor=CalibratedEnergyResponse.load(RecoveryEnergyResponse.load(EXP),EXP/'mean-correction.npz')
    target=np.array([anchor.components([0.],[e])[1][0] for e in np.linspace(*anchor.energy_bounds,81)])
    files=[]
    for name,(n,dt) in VARIANTS.items():
        q=np.array([r['charge_pe'] for r in rows if r['variant']==name])
        assert np.all(np.diff(q)>0) and min(target)>min(q) and max(target)<max(q)
        p=EXP/f'charge-map-{name}.npz'
        np.savez_compressed(p,light_pe=light,charge_pe=q,dt_ns=dt,fast_fraction=.3,duration_ns=6144.,divisions=n)
        files.append(p)
    save(EXP/'calibration-report.json',dict(rows=rows,seconds=time.monotonic()-start,
         target_charge_range=[float(target.min()),float(target.max())],tuned_parameters=0,
         input_boundary='Only frozen development spatial prior, inherited isolated calibration and fixed light shape; no evaluation photon arrays read by calibration workers.'))
    files.append(EXP/'calibration-report.json')
    save(EXP/'model-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in files})


def numerical():
    verify(True);initialize()
    sensor=CONFIG['sensor'];rng=np.random.default_rng(715010001)
    x=rng.uniform(.1,400.,160);directions=np.array([x,x*np.linspace(-.3,.3,len(x))])
    uniform_errors=[];gradient_errors=[]
    for dt in (1.,.5):
        uniform=SpatialAge(sensor,EXP/'bin/age.so',np.full((12,12),1/144),dt)
        old=DifferentiableAge(sensor,EXP/'bin/age.so',dt)
        uniform_errors.append(float(abs(uniform.forward(x)-old.forward(x)).max()))
        actual=SpatialAge(sensor,EXP/'bin/age.so',PROFILES['profile_12'],dt)
        y,dy=actual.with_gradient(x,directions)
        assert np.allclose(y,actual.forward(x),rtol=1e-12,atol=1e-12)
        for d,g in zip(directions,dy):
            finite=(actual.forward(x+1e-5*d)-actual.forward(x-1e-5*d))/2e-5
            gradient_errors.append(float(np.linalg.norm(g-finite)/np.linalg.norm(finite)))
    profile=PROFILES['profile_6'];explicit=np.zeros(32)
    engine=DifferentiableAge(sensor,EXP/'bin/age.so',1.);engine.cells=400.
    for p in profile.ravel():explicit+=engine.forward(x[:32]*p)
    grouped=SpatialAge(sensor,EXP/'bin/age.so',profile,1.).forward(x[:32])
    grouping_error=float(abs(explicit-grouped).max())
    single_error=0.;pre_event=0.
    for name,model in MODELS.items():
        for phase in PHASES:
            t=128.+phase
            for energy in (.22,1.5,4.2,7.5):
                curve=model.components([t],[energy])[0]
                old=MODELS['uniform'].components([t],[energy])[0]
                single_error=max(single_error,float(abs(curve-old).max()))
                pre_event=max(pre_event,float(abs(curve[model.readout.adc_centers+.25<=t]).max()))
        joint=model.components([128.,168.],[7.5,.22])
        swapped=model.components([168.,128.],[.22,7.5])
        assert np.allclose(joint,swapped[::-1],atol=1e-10,rtol=1e-12)
        assert float(joint.min())>=-1e-10
    # Validate the interpolation used in the primary charge-to-light inverse.
    model=MODELS[PRIMARY];mapping=model.mapping;relative=[]
    shape=light_bin_probabilities([0.],1.,6144.,.3)[0]
    for energy in (.22,1.5,4.2,7.5):
        target=model.isolated([0.],[energy])[1][0]
        amplitude=np.interp(target,mapping['charge_pe'],mapping['light_pe'])
        relative.append(float(model.age.forward(amplitude*shape).sum()/target-1))
    result=dict(uniform_limit_max_abs=max(uniform_errors),gradient_max_relative=max(gradient_errors),
         grouped_vs_explicit_tiles_max_abs=grouping_error,isolated_anchor_max_difference_mV=single_error,
         pre_event_max_mV=pre_event,charge_map_relative_errors=relative,
         passing_limitations='Anchor invariance is by construction; raw single response still requires physical checking.')
    result['passed']=bool(max(uniform_errors)<1e-9 and max(gradient_errors)<1e-5 and grouping_error<1e-9 and
                          single_error<1e-10 and pre_event<1e-10 and max(abs(v) for v in relative)<.001)
    save(REC/'numerical.json',result);print(json.dumps(result),flush=True)
    assert result['passed']


def prediction_one(job):
    kind,task,energies=job;curves={};row=dict(kind=kind,id=task['id'],predictions={});start=time.monotonic()
    if kind=='pair':
        weak=task['parents'].index(task['block'])
        for name,model in MODELS.items():
            solo=model.components(task['times_ns'],energies,interaction=False)
            joint=model.components(task['times_ns'],energies)
            row['predictions'][name]=float(100*(joint[weak].sum()/solo[weak].sum()-1))
            curves[name+'_solo']=solo;curves[name+'_joint']=joint
    else:
        donor=DONORS[task['parent']];t=donor['time_ns']+128.
        sensor=SensorConfig(**CONFIG['sensor'])
        # Only footprint membership is used; realized spatial cells do not
        # determine this raw diagnostic's frozen illumination prior.
        keep=((donor['x_mm']>=sensor.x_min_mm)&(donor['x_mm']<=sensor.x_max_mm)&
              (donor['y_mm']>=sensor.y_min_mm)&(donor['y_mm']<=sensor.y_max_mm)&(t>=0)&(t<2048.))
        weights=np.interp(HC_EV_NM/donor['energy_eV'][keep],sensor.pde_wavelength_nm,sensor.pde_probabilities[-1])
        for name,model in MODELS.items():
            light=np.bincount(np.floor(t[keep]/model.dt).astype(int),weights=weights,minlength=int(2048/model.dt))
            assert np.isclose(light.sum(),weights.sum(),rtol=1e-13)
            curves[name+'_raw']=model.readout(model.age.forward(light))
            row['predictions'][name]=float(curves[name+'_raw'].sum())
        row['parent']=task['parent']
    row['seconds']=time.monotonic()-start
    return row,curves


def predict():
    verify(True);assert read(REC/'numerical.json')['passed']
    truth={r['parent']:r['energy_dep_MeV'] for r in read(EXP/'evaluation-mothers.json')}
    jobs=[('single',t,[truth[t['parent']]]) for t in read(EXP/'single-plan.json') if t['phase_ns']==0]
    jobs += [('pair',t,[truth[p] for p in t['parents']]) for t in read(EXP/'pair-plan.json')]
    folder=ART/'predictions';folder.mkdir(exist_ok=False)
    rows=[];start=time.monotonic()
    with ProcessPoolExecutor(max_workers=4,initializer=initialize) as executor:
        for row,curves in executor.map(prediction_one,jobs):
            np.savez_compressed(folder/f"{row['kind']}-{row['id']:04d}.npz",**curves)
            rows.append(row)
            if len(rows)%4==0:print('frozen predictions',len(rows),'/36',round(time.monotonic()-start,1),flush=True)
    save(EXP/'predictions.json',dict(rows=rows,seconds=time.monotonic()-start,
         inputs='Primary anchored models take E/t only. Raw single model additionally takes accepted truth photons time/energy, never individual spatial cells. No fit to control outputs.'))
    save(EXP/'prediction-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in [EXP/'predictions.json',*folder.glob('*.npz')]})


def physical_one(task):
    if task['kind']=='single':
        assert 'energy_dep_MeV' not in task
        row,arrays=single_control.single(task)
    else:row,arrays=pair_control.one(task)
    row['kind']=task['kind'];return row,arrays


def run(mode,output=None,kind=None,ids=None):
    verify(True);assert read(REC/'numerical.json')['passed']
    assert (EXP/'prediction-freeze.json').exists()
    for p,h in read(EXP/'prediction-freeze.json').items():assert sha(ROOT/p)==h,p
    if mode=='pilot':
        tasks=[dict(read(EXP/'single-plan.json')[0],repeats=2),dict(read(EXP/'pair-plan.json')[1],repeats=2)]
        folder=ART/'pilot'
    elif mode=='replay':
        if not output or not kind or not ids:raise ValueError('Replay needs --output, --kind and --ids')
        tasks=[t for t in read(EXP/f'{kind}-plan.json') if t['id'] in ids]
        assert {t['id'] for t in tasks}==set(ids)
        folder=Path(output)
    else:
        assert read(REC/'pilot.json')['passed']
        if mode=='pair':assert read(EXP/'single-gate.json')['passed'],'Single gate stops pair expansion'
        tasks=read(EXP/f'{mode}-plan.json');folder=ART/mode
    folder.mkdir(parents=True,exist_ok=False)
    for part in {t['kind'] for t in tasks}:(folder/part).mkdir()
    start=time.monotonic();rows=[]
    with ProcessPoolExecutor(max_workers=4,initializer=initialize) as executor:
        for row,arrays in executor.map(physical_one,tasks):
            sub=folder/row['kind'];stem=f"{row['id']:04d}"
            if mode=='replay':
                original=ART/row['kind']/row['kind']/(stem+'.npz')
                with np.load(original) as f:row['physical_arrays_exact']=all(np.array_equal(arrays[k],f[k]) for k in arrays)
                assert row['physical_arrays_exact']
            np.savez_compressed(sub/(stem+'.npz'),**arrays);save(sub/(stem+'.json'),row);rows.append(row)
            if len(rows)%4==0 or mode!='single' and mode!='pair':
                print(mode,len(rows),'/',len(tasks),round(time.monotonic()-start,1),flush=True)
    record=dict(tasks=len(rows),seconds=time.monotonic()-start,workers=4,
                digitizations=sum(t['repeats']*(1 if t['kind']=='single' else 6) for t in tasks))
    if mode=='pilot':
        times={r['kind']:r['seconds'] for r in rows}
        projected=times['single']*4*64/4*1.5+times['pair']*4*20/4*1.5
        record.update(projected_seconds=projected,passed=projected<900.,rows=rows)
        save(REC/'pilot.json',record);assert record['passed']
    elif mode=='replay':
        record.update(checks=rows,passed=all(r['physical_arrays_exact'] for r in rows));save(folder/'replay.json',record)
    else:save(REC/(mode+'-run.json'),record)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=('prepare','calibrate','numerical','predict','pilot','single','pair','replay'))
    parser.add_argument('--output');parser.add_argument('--kind',choices=('single','pair'));parser.add_argument('--ids',type=int,nargs='+')
    args=parser.parse_args()
    if args.mode in ('prepare','calibrate','numerical','predict'):globals()[args.mode]()
    else:run(args.mode,args.output,args.kind,args.ids)


if __name__=='__main__':main()
