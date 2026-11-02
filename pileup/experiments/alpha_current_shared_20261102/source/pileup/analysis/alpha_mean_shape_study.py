"""Frozen area-preserving mean calibration and staged response validation."""
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

from .alpha_mean_shape_correction import (MeanShapeCorrection,MeanShapeRecovery,
    calibrate_profiles,positivity_scale)
from .alpha_pair_digitizer_optimized import simulate_avalanches
from .alpha_pair_energy_calibrated import CalibratedEnergyResponse
from .alpha_pair_recovery_energy import RecoveryEnergyResponse
from .alpha_pair_spatial_recovery import SpatialRecovery
from .alpha_pair_spatial_energy import common_inputs,fit_energy
from .sipm_response import SensorConfig,PhotonHits,AvalancheResult,read_photon_hits
from .sipm_readout import ReadoutConfig,synthesize_readout

NAME='alpha_mean_shape_20261102'
ROOT=Path(os.environ.get('GAGG_MEAN_SHAPE_ROOT',Path(__file__).resolve().parents[2])).resolve()
EXP=ROOT/'pileup/experiments'/NAME
ART=ROOT/'jupyter/artifacts'/NAME
REC=ROOT/'jupyter/run_records'/NAME
PARENT=ROOT/'pileup/experiments/alpha_single_residual_sources_20261102'
SPATIAL=ROOT/'pileup/experiments/alpha_pair_spatial_recovery_20261102'
SPATIAL_ART=ROOT/'jupyter/artifacts/alpha_pair_spatial_recovery_20261102'
CONFIRM=ROOT/'jupyter/artifacts/alpha_pair_joint_confirmation_20261026'
PHASES=(0.,.5,1.,1.5)
GROUPS=('weak','low','mid','high')
WORKERS=8


def read(p):return json.loads(Path(p).read_text())


def save(p,value):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    temporary=p.with_suffix(p.suffix+'.partial')
    temporary.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    temporary.replace(p)


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        while block:=f.read(1024*1024):h.update(block)
    return h.hexdigest()


def snapshot():
    target=EXP/'source/pileup/analysis';target.mkdir(parents=True)
    (target.parent/'__init__.py').write_text('');(target/'__init__.py').write_text('')
    pending=['alpha_mean_shape_study','test_alpha_mean_shape_correction'];copied=set()
    while pending:
        name=pending.pop()
        if name in copied:continue
        source=ROOT/'pileup/analysis'/f'{name}.py'
        tree=ast.parse(source.read_text());shutil.copyfile(source,target/source.name);copied.add(name)
        for node in ast.walk(tree):
            if isinstance(node,ast.ImportFrom) and node.level==1:
                pending.extend([node.module] if node.module else [a.name for a in node.names])
    for name in ('sipm_age_mean.cpp','sipm_age_mean_jacobian.cpp'):
        shutil.copyfile(ROOT/'pileup/analysis'/name,target/name)
    (EXP/'bin').mkdir()
    command=['g++','-O3','-std=c++17','-fPIC','-shared',
        str((target/'sipm_age_mean.cpp').relative_to(ROOT)),str((target/'sipm_age_mean_jacobian.cpp').relative_to(ROOT)),
        '-o',str((EXP/'bin/age.so').relative_to(ROOT))]
    subprocess.run(command,cwd=ROOT,check=True)
    save(EXP/'build.json',dict(command=command,binary_sha256=sha(EXP/'bin/age.so'),
        compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0]))
    (EXP/'run.py').write_text("from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\n"
        "os.environ['GAGG_MEAN_SHAPE_ROOT']=str(f.parents[2])\nos.environ['OPENBLAS_NUM_THREADS']='1'\nos.environ['OMP_NUM_THREADS']='1'\n"
        "sys.path.insert(0,str(f/'source'))\nif __name__=='__main__':\n    from pileup.analysis.alpha_mean_shape_study import main\n    main()\n")


def prepare_one(row):
    assert sha(ROOT/row['root'])==row['root_original_sha256']
    meta=read(ROOT/row['metadata_file']);config=read(EXP/'sensor-config.json')
    assert meta['config']['sensor']==config['sensor'] and meta['config']['readout']==config['readout']
    with np.load(ROOT/row['truth_file']) as f:
        assert len(f['edep_MeV'])==1 and abs(float(f['edep_MeV'][0])-row['edep_MeV'])<1e-10
    destination=ART/'photons'/f"{row['id']:04d}.npz"
    shutil.copyfile(ROOT/row['photon_file'],destination)
    shutil.copyfile(ROOT/row['metadata_file'],EXP/'original-metadata'/f"{row['id']:04d}.json")
    hits=read_photon_hits(ROOT/row['root'])
    assert hits.metadata['run_info']['random_seed']==row['seed']
    save(EXP/'original-runinfo'/f"{row['id']:04d}.json",hits.metadata)
    return dict(id=row['id'],photon_sha256=sha(destination),original_runinfo_saved=True,truth_energy_matches=True)


def prepare():
    for p in (EXP,ART,REC):p.mkdir(parents=True,exist_ok=False)
    for p in ('photons','voltage/train','voltage/single','voltage/pair','truth/train','truth/single','truth/pair',
              'fits','figures','pilot','replay'):(ART/p).mkdir(parents=True)
    for p in ('original-metadata','original-runinfo'):(EXP/p).mkdir()
    for name in ('template.npz','energy-model.json','response-parameters.json','mean-correction.npz',
                 'noise-calibration.json','sensor-config.json','charge-map-spatial12.npz','spatial-profile.npz'):
        shutil.copyfile(PARENT/name,EXP/name)
    mothers=[];train_plan=[];single_plan=[]
    source_levels=(1.625,2.5,5.486,8.)
    for row in read(SPATIAL/'development-mothers.json'):
        group=source_levels.index(row['source_energy_MeV']);i=len(mothers)
        mothers.append(dict(row,id=i,group=group,condition=GROUPS[group],role='train',
            photon_file=str((SPATIAL_ART/'development-photons'/f"{row['profile_id']:04d}.npz").relative_to(ROOT))))
        train_plan.append(dict(id=len(train_plan),kind='train',parents=[i],times_ns=[128.],phase_ns=0.,
            seed_base=736000000+i*100,repeats=8))
    assert len(mothers)==128
    parents=read(CONFIRM/'physical-parents.json')['records']
    excluded={r['source']['root'] for r in read(SPATIAL/'evaluation-mothers.json')}
    excluded|={r['source']['root'] for r in read(ROOT/'pileup/experiments/alpha_pair_causal_controls_20261101/mother-provenance.json')}
    for group,condition in enumerate(GROUPS):
        candidates=[(p,r) for p,r in enumerate(parents) if p//80==group and r['root'] not in excluded]
        candidates.sort(key=lambda x:hashlib.sha256(('736000:'+x[1]['root']).encode()).hexdigest())
        for block,(parent,row) in enumerate(candidates[:16]):
            with np.load(ROOT/row['truth']) as f:energy=float(f['edep_MeV'][0])
            i=len(mothers)
            mothers.append(dict(row,id=i,group=group,condition=condition,role='validation',block=block,
                source_energy_MeV=row['energy_MeV'],edep_MeV=energy,
                metadata_file=row['metadata'],truth_file=row['truth'],original_parent=parent,
                photon_file=str((CONFIRM/'photons'/f'{parent:04d}.npz').relative_to(ROOT))))
            for phase in PHASES:
                single_plan.append(dict(id=len(single_plan),kind='single',parents=[i],times_ns=[128.+phase],
                    phase_ns=phase,seed_base=736100000+i*100,repeats=8))
    train_roots={r['root'] for r in mothers if r['role']=='train'}
    eval_roots={r['root'] for r in mothers if r['role']=='validation'}
    assert len(train_roots)==128 and len(eval_roots)==64 and not train_roots&eval_roots
    assert not eval_roots&excluded
    pair_plan=[]
    for block in range(16):
        weak=128+block;high=176+block
        for condition,gap in [('high_weak',8.),('high_weak',40.),('high_weak',200.),('high_weak',1200.),('weak_high',40.)]:
            pair_plan.append(dict(id=len(pair_plan),kind='pair',block=block,condition=condition,gap_ns=gap,
                parents=[high,weak] if condition=='high_weak' else [weak,high],times_ns=[128.,128.+gap],
                repeats=8,seed_base=736300000+block*1000))
    save(EXP/'mothers.json',mothers)
    for name,plan in [('train',train_plan),('single',single_plan),('pair',pair_plan)]:save(EXP/(name+'-plan.json'),plan)
    save(EXP/'split.json',dict(train=sorted(train_roots),validation=sorted(eval_roots),excluded_recent=sorted(excluded),
        disjoint=True,exposure='Training is128 former mean/profile-development mothers;64 validation mothers belong to the historical confirmation pool and have older exposure,not fresh transport'))
    save(EXP/'protocol.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        git_status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),
        preceding_push='origin/main6d8146a before this new work',
        scope='Conditional mean SHAPE validation only; knownK/t. Single-energy fits and truth-assisted forward pairs; no pair-energy fit or blind reconstruction claim',
        calibration='128 existing profile/mean-training mothers,32 each of4source energies;8new independent digitizations each;one phase0',
        validation='64 different historical confirmation mothers,16per group,SHA256 path selection excluding32recent causal/spatial control mothers;4phasesx8repeats. Never train/select on validation',
        correction='Normalize each mother mean residual by frozen q0(true E_dep);equal mother mean per energy group;Gaussian4ns smoothing;2ns causal guard then8ns onset ramp,160ns end taper to1600ns;old-template support weight reference/(reference+.01peak)',
        area='Subtract old positive-template-weighted area mode;zero full2ns-grid sum for every knot. Linear interpolation in logq0,endpoint clamping;time interpolation preserves full area at fractional phases. No posthoc clipping',
        energy='Old q0(E),isolated returned Q(E),charge-to-light mapping and spatial shared-minus-isolated correction unchanged. Shape-only empirical anchor is a project hypothesis,not altered microscopic light-time simulation',
        positivity='One global scale<=1 set only by old model and training correction on81logE points x4phases;99%of admissible bound ifneeded. Minimum allowedscale=.5;no validation tuning',
        fitting='Same voltage-only initial E and fixed baseline GLS weights for old/corrected;same scalar bounded solver. Each8individuals plus8repeat mean;knownonset/K1;truth files forbidden in fit workers',
        projection='Use frozen nominal old-model influences fromL92;diagnostic at true E,100delta ln E units;phase0strong mother means;not nonlinear pair-energy error',
        single_gates=dict(mean_energy_bias_percent=5.,single_rmse_ratio=1.1,phase_bias_span_pp=1.,
            mean_waveform_mse_ratio=1.1,weak_projection_abs_mean_ratio=.7,weak_projection_floor_100logE=10.,
            max_nonconverged=0),
        shared_gates=dict(weak_area_error_pp=5.,mean_component_mse_ratio=1.1),
        gate_scope='Prespecified engineering point estimates;not literature thresholds or independent confirmation',
        stages=['prepare','pilot','train','calibrate','numerical','single','fit','single_report','pair_pilot onlyif singles pass','pair','pair_report','replay','report'],
        observation='4096ns digitization;first2048ns/1024samples;2nsperiod,.5nsaperture;first active photon at128+phase ns;actual photon times/positions/energies retained',
        phases_ns=list(PHASES),repeats=8,training_digitizations=1024,validation_single_digitizations=2048,
        conditional_pair_digitizations=1920,workers=WORKERS,BLAS_threads=1,cost_limit_seconds=900.,
        randomness='Independent repeat seeds;common seeds across phases/gaps change consumption,not identical avalanches. Shared/isolated truth components use separate digitizations',
        source_rate='One/two controlled optical mothers;no Poisson scan. Historical ROOT RunInfo and macros preserved',
        literature='ReuseL91/Fixsen and L92 project evidence;no new external paper claims',
        no_new_ROOT=True,default_changed=False))
    for name in ('README.md','docs/literature-use.md','jupyter/codex_0007_SiPM_Energy_Calibration.ipynb',
                 'jupyter/Research_Plan_20260921.ipynb','jupyter/Run_Catalog.ipynb'):
        p=ROOT/name;shutil.copyfile(p,REC/(p.name+'.before'))
    shutil.copyfile(PARENT/'literature-evidence.json',EXP/'literature-evidence.json')
    snapshot()
    save(EXP/'design-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in EXP.rglob('*') if p.is_file() and 'bin' not in p.relative_to(EXP).parts})
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:provenance=list(pool.map(prepare_one,mothers))
    save(EXP/'provenance.json',provenance)
    paths=[ROOT/r[k] for r in mothers for k in ('root','truth_file','metadata_file')]
    paths+=list((ART/'photons').glob('*.npz'))+list((EXP/'original-metadata').glob('*.json'))+list((EXP/'original-runinfo').glob('*.json'))
    previous_influence=ROOT/'jupyter/artifacts/alpha_single_residual_sources_20261102/influences.npz'
    shutil.copyfile(previous_influence,EXP/'nominal-influences.npz');paths.append(EXP/'nominal-influences.npz')
    save(EXP/'input-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})
    import scipy
    save(EXP/'environment.json',dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,workers=WORKERS,BLAS_threads=1))
    print('Frozen128 calibration/64 validation mothers;1024+2048 singles;pair expansion conditional',flush=True)


def verify(model=False):
    counts={}
    for name in ('design-freeze.json','input-freeze.json')+(('model-freeze.json',) if model else ()):
        items=read(EXP/name)
        for p,d in items.items():assert sha(ROOT/p)==d,p
        counts[name]=len(items)
    return counts


def initialize(corrected=False):
    if not hasattr(np,'trapezoid'):np.trapezoid=np.trapz
    global SENSOR,READOUT,OLD,NEW,NOISE,CONFIG
    CONFIG=read(EXP/'sensor-config.json');SENSOR=SensorConfig(**CONFIG['sensor']);READOUT=ReadoutConfig(**CONFIG['readout'])
    anchor=CalibratedEnergyResponse.load(RecoveryEnergyResponse.load(EXP),EXP/'mean-correction.npz')
    with np.load(EXP/'charge-map-spatial12.npz') as f:mapping={k:f[k].copy() for k in f.files}
    with np.load(EXP/'spatial-profile.npz') as f:profile=f['profile_12'].copy()
    OLD=SpatialRecovery(anchor,CONFIG,EXP/'bin/age.so',mapping,profile)
    NEW=MeanShapeRecovery(anchor,CONFIG,EXP/'bin/age.so',mapping,profile,MeanShapeCorrection.load(EXP/'shape-correction.npz')) if corrected else None
    NOISE=read(EXP/'noise-calibration.json')


def signal(parents,times,seed):
    fields={k:[] for k in ('time_ns','x_mm','y_mm','energy_eV')};owners=[]
    for j,(parent,t0) in enumerate(zip(parents,times)):
        with np.load(ART/'photons'/f'{parent:04d}.npz') as f:
            t=f['time_ns']+t0;keep=t<4096.
            for key in fields:fields[key].append((t if key=='time_ns' else f[key])[keep])
            owners.append(np.full(int(keep.sum()),j))
    arrays={k:np.concatenate(v) for k,v in fields.items()};owner=np.concatenate(owners)
    order=np.argsort(arrays['time_ns'],kind='stable');owner=owner[order]
    hits=PhotonHits(**{k:v[order] for k,v in arrays.items()},window_ns=4096.)
    av=simulate_avalanches(hits,SENSOR,seed)
    voltage=synthesize_readout(av,READOUT,seed+10000000).arrays['analog_mV'][:1024]
    rec=av.records;assignment=np.full(len(rec),-1,int);direct=rec['kind']==0
    assignment[direct]=owner[rec['photon_index'][direct]]
    for j in np.flatnonzero(rec['parent_id']>=0):assignment[j]=assignment[rec['parent_id'][j]]
    curves=[];charges=[]
    for j in range(len(parents)):
        sub=rec[assignment==j];event=AvalancheResult(sub,4096.,av.full_charge_C,np.array([]),{})
        curves.append(synthesize_readout(event,READOUT,0).arrays['signal_mV'][:1024]);charges.append(float(sub['charge_C'].sum()/av.full_charge_C))
    return voltage,np.array(curves),np.array(charges)


def simulate(job):
    task,output=job;start=time.monotonic();voltage=[];components=[];charge=[];isolated=[]
    for repeat in range(task['repeats']):
        seed=task['seed_base']+repeat
        v,s,q=signal(task['parents'],task['times_ns'],seed)
        voltage.append(v);components.append(s);charge.append(q)
        if task['kind']=='pair':
            isolated.append(np.concatenate([signal([p],[t],seed)[1] for p,t in zip(task['parents'],task['times_ns'])]))
    prefix='' if output=='formal' else output+'/'
    vp=ART/(prefix+'voltage')/task['kind'];tp=ART/(prefix+'truth')/task['kind']
    vp.mkdir(parents=True,exist_ok=True);tp.mkdir(parents=True,exist_ok=True)
    voltage=np.asarray(voltage);components=np.asarray(components)
    assert np.isfinite(voltage).all() and np.isfinite(components).all()
    for j,t in enumerate(task['times_ns']):
        assert np.max(abs(components[:,j,OLD.readout.adc_centers+.25<=t]),initial=0)<1e-10
    np.savez_compressed(vp/f"{task['id']:04d}.npz",voltage_mV=voltage)
    arrays=dict(components_mV=components,charge_pe=charge)
    if isolated:arrays['isolated_mV']=isolated
    np.savez_compressed(tp/f"{task['id']:04d}.npz",**arrays)
    result=dict(id=task['id'],kind=task['kind'],seconds=time.monotonic()-start,
        digitizations=task['repeats']*(3 if task['kind']=='pair' else 1),finite=True,pre_event_zero=True)
    save(REC/output/task['kind']/f"{task['id']:04d}.json",result)
    return result


def generate(stage):
    verify();output='formal'
    if stage=='pilot':
        tasks=[dict(read(EXP/'train-plan.json')[i],repeats=2) for i in (0,96)];output='pilot'
    elif stage=='pair_pilot':
        assert read(EXP/'single-gate.json')['passed']
        tasks=[dict(read(EXP/'pair-plan.json')[1],repeats=2)];output='pair_pilot'
    elif stage=='replay':
        tasks=[dict(read(EXP/'train-plan.json')[1],repeats=2),dict(read(EXP/'single-plan.json')[197],repeats=2)]
        if (EXP/'pair-input-freeze.json').exists():tasks.append(dict(read(EXP/'pair-plan.json')[6],repeats=2))
        output='replay'
    else:
        assert read(EXP/'pilot-gate.json')['passed']
        if stage=='single':assert read(EXP/'numerical.json')['passed']
        if stage=='pair':
            assert read(EXP/'single-gate.json')['passed'] and read(EXP/'pair-pilot-gate.json')['passed']
        tasks=read(EXP/(stage+'-plan.json'))
    for task in tasks:
        prefix='' if output=='formal' else output+'/'
        assert not (ART/(prefix+'voltage')/task['kind']/f"{task['id']:04d}.npz").exists()
    start=time.monotonic()
    with ProcessPoolExecutor(max_workers=WORKERS,initializer=initialize) as pool:
        rows=[]
        for r in pool.map(simulate,[(t,output) for t in tasks]):
            rows.append(r)
            if len(rows)%16==0:print(stage,len(rows),'/',len(tasks),flush=True)
    save(REC/(stage+'.json'),dict(rows=rows,seconds=time.monotonic()-start,digitizations=sum(r['digitizations'] for r in rows)))
    if stage in ('pilot','pair_pilot'):
        estimate=max(r['seconds']/r['digitizations'] for r in rows)*(3072 if stage=='pilot' else 1920)/WORKERS*1.5
        save(EXP/('pilot-gate.json' if stage=='pilot' else 'pair-pilot-gate.json'),dict(projected_seconds=estimate,limit_seconds=900.,passed=bool(estimate<900.)))
        assert estimate<900.,estimate
    elif stage=='replay':
        for task in tasks:
            for folder in ('voltage','truth'):
                with np.load(ART/folder/task['kind']/f"{task['id']:04d}.npz") as a,np.load(ART/'replay'/folder/task['kind']/f"{task['id']:04d}.npz") as b:
                    assert a.files==b.files
                    for k in a.files:np.testing.assert_array_equal(a[k][:2],b[k])
        save(EXP/'replay-validation.json',dict(passed=True,tasks=len(tasks),exact_arrays=True))
    else:
        paths=[p for folder in ('voltage','truth') for p in (ART/folder/stage).glob('*.npz')]
        save(EXP/(stage+'-input-freeze.json'),{str(p.relative_to(ROOT)):sha(p) for p in paths})
    print(stage,'complete',sum(r['digitizations'] for r in rows),'digitizations',flush=True)


def calibrate():
    verify();initialize();assert not (EXP/'shape-correction.npz').exists()
    for p,d in read(EXP/'train-input-freeze.json').items():assert sha(ROOT/p)==d,p
    mothers=read(EXP/'mothers.json');means=[];predicted=[];q=[];groups=[]
    for task in read(EXP/'train-plan.json'):
        row=mothers[task['parents'][0]];e=row['edep_MeV']
        with np.load(ART/'truth/train'/f"{task['id']:04d}.npz") as f:means.append(f['components_mV'].mean(0)[0])
        predicted.append(OLD.isolated([128.],[e])[0][0]);q.append(OLD.anchor.q_isolated([e])[0]);groups.append(row['group'])
    correction,rows=calibrate_profiles(OLD.readout.adc_centers,means,predicted,q,groups)
    correction.scale=positivity_scale(OLD,correction,np.geomspace(*OLD.energy_bounds,81),PHASES)
    correction.save(EXP/'shape-correction.npz')
    save(EXP/'calibration.json',dict(groups=rows,scale=correction.scale,smoothing_ns=4.,
        minimum_allowed_scale=.5,passed=bool(correction.scale>=.5),true_E_used_only_for_training=True,
        area_preserving=True,shared_charge_map_unchanged=True))
    np.savez_compressed(ART/'calibration.npz',means_mV=means,old_predictions_mV=predicted,q0=q,groups=groups)
    save(EXP/'model-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in [EXP/'shape-correction.npz',EXP/'calibration.json',ART/'calibration.npz']})
    assert correction.scale>=.5,'Positivity constraint prevents the proposed update; do not inspect validation'
    print('Calibration frozen, global positivity scale',correction.scale,flush=True)


def numerical():
    verify(True);initialize(True)
    full_grid=np.arange(3072)*2.+.25;areas=[];pre=[];minimum=0.
    for e in np.geomspace(*OLD.energy_bounds,81):
        q=OLD.anchor.q_isolated([e])[0]
        for phase in PHASES:
            t=128.+phase
            correction=NEW.correction.voltage(full_grid,t,q)
            areas.append(float(abs(2*correction.sum())))
            pre.append(float(np.max(abs(correction[full_grid+.25<=t]),initial=0)))
            minimum=min(minimum,float(NEW.isolated([t],[e])[0].min()))
    rows=[]
    for gap in (8.,40.,200.,1200.):
        times=[128.,128.+gap];energy=[7.5,.25]
        old_a,old_q=OLD.isolated(times,energy);new_a,new_q=NEW.isolated(times,energy)
        old_mu=OLD.components(times,energy);new_mu=NEW.components(times,energy)
        difference=float(np.max(abs((new_mu-new_a)-(old_mu-old_a))))
        np.testing.assert_array_equal(old_q,new_q)
        rows.append(dict(gap_ns=gap,shared_minus_isolated_difference_mV=difference,
            returned_charge_equal=True,model_min_component_mV=float(new_mu.min())))
    t=[128.5,168.5];loge=np.log([7.2,.28]);prob=NEW.prepare(t)
    _,gradient=NEW.mean_and_gradient(t,loge,prob,strength=1.)
    errors=[]
    for j in range(2):
        step=np.zeros(2);step[j]=1e-5
        up=NEW.components(t,np.exp(loge+step)).sum(0);down=NEW.components(t,np.exp(loge-step)).sum(0)
        fd=(up-down)/2e-5
        errors.append(float(np.linalg.norm(fd-gradient[j])/max(np.linalg.norm(fd),1e-12)))
    scale=NEW.correction.scale
    NEW.correction.scale=0.
    identity=float(np.max(abs(NEW.components(t,np.exp(loge))-OLD.components(t,np.exp(loge)))))
    NEW.correction.scale=scale
    result=dict(max_complete_area_change_mV_ns=max(areas),max_pre_event_mV=max(pre),
        minimum_isolated_mV=minimum,shared_controls=rows,gradient_relative_errors=errors,
        zero_correction_identity_mV=identity,
        passed=bool(max(areas)<1e-8 and max(pre)==0 and minimum>=-1e-8 and identity<1e-10 and
            max(errors)<1e-4 and max(r['shared_minus_isolated_difference_mV'] for r in rows)<1e-9))
    save(EXP/'numerical.json',result);assert result['passed'],result
    print('Area,phase,positivity,charge/shared-invariance and gradient controls passed',flush=True)


def fit_initialize():
    initialize(True)
    def audit(event,args):
        if event!='open' or not isinstance(args[0],(str,bytes,os.PathLike)):return
        p=Path(os.fsdecode(args[0])).resolve()
        forbidden=p.suffix=='.root' or p.name=='mothers.json'
        forbidden|=p.is_relative_to(ART/'truth') or p.is_relative_to(ART/'photons')
        forbidden|=p.is_relative_to(EXP/'original-metadata') or p.is_relative_to(EXP/'original-runinfo')
        if forbidden:raise PermissionError('Single-energy fit cannot read event truth: '+str(p))
    sys.addaudithook(audit)
    for p in (EXP/'mothers.json',ART/'truth/single/0000.npz',ART/'photons/0000.npz'):
        try:p.open('rb')
        except PermissionError:pass
        else:raise AssertionError('Truth-read guard failed')


def fit_one(job):
    i,t,output=job;start=time.monotonic()
    with np.load(ART/'voltage/single'/f'{i:04d}.npz') as f:y=f['voltage_mV'].copy()
    rows=[]
    for repeat,voltage in [*enumerate(y),(-1,y.mean(0))]:
        common=common_inputs(OLD,voltage,[t],NOISE)
        signature=hashlib.sha256(common['metric'].factor.tobytes()+common['metric'].modes.tobytes()).hexdigest()
        for method,model in [('old',OLD),('corrected',NEW)]:
            fit=fit_energy(model,voltage,[t],common)
            rows.append(dict(method=method,repeat=repeat,known_time_ns=t,metric_signature=signature,**fit))
    p=ART/output/f'{i:04d}.json';save(p,dict(id=i,rows=rows,seconds=time.monotonic()-start,
        true_E_input=False,truth_guard=True,known_K1_t=True))
    return dict(id=i,seconds=time.monotonic()-start,fits=len(rows))


def fit_singles():
    verify(True)
    for p,d in read(EXP/'single-input-freeze.json').items():assert sha(ROOT/p)==d,p
    tasks=read(EXP/'single-plan.json')
    assert not list((ART/'fits').glob('*.json'))
    with ProcessPoolExecutor(max_workers=WORKERS,initializer=fit_initialize) as pool:
        pilot=list(pool.map(fit_one,[(i,tasks[i]['times_ns'][0],'fit-pilot') for i in (0,192)]))
    estimate=max(r['seconds'] for r in pilot)*len(tasks)/WORKERS*1.5
    save(EXP/'fit-pilot-gate.json',dict(rows=pilot,projected_seconds=estimate,passed=bool(estimate<900.)))
    assert estimate<900.,estimate
    start=time.monotonic()
    with ProcessPoolExecutor(max_workers=WORKERS,initializer=fit_initialize) as pool:
        rows=[]
        for result in pool.map(fit_one,[(t['id'],t['times_ns'][0],'fits') for t in tasks]):
            rows.append(result)
            if len(rows)%32==0:print('single fits',len(rows),'/',len(tasks),flush=True)
    for i in (0,192):
        assert read(ART/'fits'/f'{i:04d}.json')['rows']==read(ART/'fit-pilot'/f'{i:04d}.json')['rows'] or all(
            {k:v for k,v in a.items() if k!='seconds'}=={k:v for k,v in b.items() if k!='seconds'}
            for a,b in zip(read(ART/'fits'/f'{i:04d}.json')['rows'],read(ART/'fit-pilot'/f'{i:04d}.json')['rows']))
    save(REC/'fits.json',dict(rows=rows,fits=sum(r['fits'] for r in rows),seconds=time.monotonic()-start,
        pilot_fits=sum(r['fits'] for r in pilot),pilot_results_equal_except_timing=True))
    save(EXP/'fit-output-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in (ART/'fits').glob('*.json')})


def single_report():
    verify(True);initialize(True)
    for name in ('single-input-freeze.json','fit-output-freeze.json'):
        for p,d in read(EXP/name).items():assert sha(ROOT/p)==d,p
    mothers=read(EXP/'mothers.json');tasks=read(EXP/'single-plan.json')
    details=[];waveforms={g:[] for g in range(4)};old_res={g:[] for g in range(4)};new_res={g:[] for g in range(4)}
    for task in tasks:
        row=mothers[task['parents'][0]];energy=row['edep_MeV']
        fit=read(ART/'fits'/f"{task['id']:04d}.json")
        assert fit['truth_guard'] and not fit['true_E_input']
        for f in fit['rows']:
            details.append(dict(case=task['id'],parent=row['id'],group=row['group'],condition=row['condition'],
                phase_ns=task['phase_ns'],method=f['method'],repeat=f['repeat'],
                relative_error_percent=100*(f['energy_MeV'][0]/energy-1),converged=f['converged'],boundary=f['boundary']))
        if task['phase_ns']==0:
            with np.load(ART/'truth/single'/f"{task['id']:04d}.npz") as f:actual=f['components_mV'].mean(0)[0]
            old=OLD.isolated(task['times_ns'],[energy])[0][0];new=NEW.isolated(task['times_ns'],[energy])[0][0]
            waveforms[row['group']].append(actual);old_res[row['group']].append(actual-old);new_res[row['group']].append(actual-new)
    energy_rows=[];failures=[]
    for group in range(4):
        for method in ('old','corrected'):
            for phase in PHASES:
                for kind in ('individual','mean'):
                    rows=[r for r in details if r['group']==group and r['method']==method and r['phase_ns']==phase and (r['repeat']>=0 if kind=='individual' else r['repeat']==-1)]
                    err=np.array([r['relative_error_percent'] for r in rows])
                    energy_rows.append(dict(group=group,condition=GROUPS[group],method=method,phase_ns=phase,kind=kind,
                        mothers=16,rows=len(rows),bias_percent=float(err.mean()),rmse_percent=float(np.sqrt(np.mean(err**2))),
                        error_sd_percent=float(err.std(ddof=1)),strict10_count=int(np.sum(abs(err)<=10)),
                        boundary=sum(r['boundary'] for r in rows),nonconverged=sum(not r['converged'] for r in rows)))
        for phase in PHASES:
            old=next(r for r in energy_rows if r['group']==group and r['method']=='old' and r['phase_ns']==phase and r['kind']=='individual')
            new=next(r for r in energy_rows if r['group']==group and r['method']=='corrected' and r['phase_ns']==phase and r['kind']=='individual')
            mean=next(r for r in energy_rows if r['group']==group and r['method']=='corrected' and r['phase_ns']==phase and r['kind']=='mean')
            if abs(mean['bias_percent'])>5:failures.append(f'{GROUPS[group]} phase{phase}:mean energy bias')
            if new['rmse_percent']>1.1*old['rmse_percent']:failures.append(f'{GROUPS[group]} phase{phase}:individual energy RMSE')
        span=np.ptp([r['bias_percent'] for r in energy_rows if r['group']==group and r['method']=='corrected' and r['kind']=='mean'])
        if span>1:failures.append(GROUPS[group]+':phase bias span')
    if any(r['nonconverged'] for r in energy_rows if r['method']=='corrected'):failures.append('nonconverged energy fit')
    wave_rows=[]
    arrays={}
    for group in range(4):
        old=np.asarray(old_res[group]);new=np.asarray(new_res[group])
        before=float(np.mean(old.mean(0)**2));after=float(np.mean(new.mean(0)**2))
        ratio=after/before
        wave_rows.append(dict(group=group,condition=GROUPS[group],mothers=16,old_mean_residual_mse_mV2=before,
            corrected_mean_residual_mse_mV2=after,ratio=ratio))
        if ratio>1.1:failures.append(GROUPS[group]+':mean waveform MSE')
        arrays[GROUPS[group]+'_old_residual_mV']=old;arrays[GROUPS[group]+'_corrected_residual_mV']=new
        arrays[GROUPS[group]+'_actual_mean_mV']=waveforms[group]
    with np.load(EXP/'nominal-influences.npz') as f:operator=f['high'][1:]
    before=100*np.asarray(old_res[3])@operator.T;after=100*np.asarray(new_res[3])@operator.T
    rng=np.random.default_rng(736900);idx=rng.integers(16,size=(4000,16))
    projections=[]
    for j,gap in enumerate((8.,40.,200.)):
        b=float(before[:,j].mean());a=float(after[:,j].mean());limit=max(10.,.7*abs(b))
        ci_b=np.percentile(before[idx,j].mean(1),[2.5,97.5]);ci_a=np.percentile(after[idx,j].mean(1),[2.5,97.5])
        ci_change=np.percentile((after-before)[idx,j].mean(1),[2.5,97.5])
        passed=abs(a)<=limit
        projections.append(dict(gap_ns=gap,mothers=16,old_mean=b,corrected_mean=a,old_ci95=ci_b.tolist(),
            corrected_ci95=ci_a.tolist(),paired_change_ci95=ci_change.tolist(),absolute_mean_limit=limit,
            passed=bool(passed),units='100delta ln E;fixed old nominal operator,not pair energy error'))
        if gap in (8.,40.) and not passed:failures.append(f'weak projection gap{int(gap)}')
    arrays['old_projections_100logE']=before;arrays['corrected_projections_100logE']=after
    arrays['time_ns']=OLD.readout.adc_centers
    np.savez_compressed(ART/'single-summary.npz',**arrays)
    save(EXP/'single-result.json',dict(energy=energy_rows,waveform=wave_rows,projection=projections))
    save(EXP/'single-fit-details.json',details)
    save(EXP/'single-gate.json',dict(passed=not failures,failures=failures,point_estimate_engineering_gate=True,
        next='Perform frozen shared-response controls' if not failures else 'Stop before shared expansion; retain isolated failures'))
    print(json.dumps(dict(passed=not failures,failures=failures,projection=projections),indent=2),flush=True)


def pair_report():
    verify(True);initialize(True)
    for p,d in read(EXP/'pair-input-freeze.json').items():assert sha(ROOT/p)==d,p
    mothers=read(EXP/'mothers.json');details=[];examples={}
    for task in read(EXP/'pair-plan.json'):
        energies=[mothers[i]['edep_MeV'] for i in task['parents']];weak=int(np.argmin(energies))
        with np.load(ART/'truth/pair'/f"{task['id']:04d}.npz") as f:
            observed=f['components_mV'].mean(0);isolated=f['isolated_mV'].mean(0)
        actual_ratio=100*(observed[weak].sum()/isolated[weak].sum()-1)
        for method,model in [('old',OLD),('corrected',NEW)]:
            prediction=model.components(task['times_ns'],energies);solo=model.isolated(task['times_ns'],energies)[0]
            ratio=100*(prediction[weak].sum()/solo[weak].sum()-1)
            residual=observed-prediction
            details.append(dict(case=task['id'],block=task['block'],condition=task['condition'],gap_ns=task['gap_ns'],
                method=method,actual_weak_area_change_percent=float(actual_ratio),predicted_weak_area_change_percent=float(ratio),
                error_pp=float(ratio-actual_ratio),weak_mse_mV2=float(np.mean(residual[weak]**2)),
                strong_mse_mV2=float(np.mean(residual[1-weak]**2)),model_min_mV=float(prediction.min())))
            if task['block']==0:
                examples[f"case{task['id']}_actual"]=observed
                examples[f"case{task['id']}_{method}"]=prediction
    groups=[];failures=[]
    for condition,gap in [('high_weak',8.),('high_weak',40.),('high_weak',200.),('high_weak',1200.),('weak_high',40.)]:
        matched={m:[r for r in details if r['condition']==condition and r['gap_ns']==gap and r['method']==m] for m in ('old','corrected')}
        rows={m:{k:float(np.mean([r[k] for r in v])) for k in ('actual_weak_area_change_percent','predicted_weak_area_change_percent','error_pp','weak_mse_mV2','strong_mse_mV2')} for m,v in matched.items()}
        ratios={k:rows['corrected'][k]/rows['old'][k] for k in ('weak_mse_mV2','strong_mse_mV2')}
        passed=abs(rows['corrected']['error_pp'])<=5 and max(ratios.values())<=1.1
        if not passed:failures.append(f'{condition}_{gap}')
        groups.append(dict(condition=condition,gap_ns=gap,mother_pairs=16,old=rows['old'],corrected=rows['corrected'],
            mse_ratios=ratios,passed=bool(passed)))
    save(EXP/'pair-result.json',dict(groups=groups,details=details))
    save(EXP/'pair-gate.json',dict(passed=not failures,failures=failures,
        energy_reconstruction_not_tested=True,default_changed=False))
    np.savez_compressed(ART/'pair-examples.npz',**examples)
    print(json.dumps(dict(passed=not failures,groups=groups),indent=2),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['prepare','pilot','train','calibrate','numerical','single','fit',
        'single_report','pair_pilot','pair','pair_report','replay','verify'])
    stage=p.parse_args().stage
    if stage=='prepare':prepare()
    elif stage in ('pilot','train','single','pair_pilot','pair','replay'):generate(stage)
    elif stage=='calibrate':calibrate()
    elif stage=='numerical':numerical()
    elif stage=='fit':fit_singles()
    elif stage=='single_report':single_report()
    elif stage=='pair_report':pair_report()
    else:print(json.dumps(verify(model=(EXP/'model-freeze.json').exists())))


if __name__=='__main__':main()
