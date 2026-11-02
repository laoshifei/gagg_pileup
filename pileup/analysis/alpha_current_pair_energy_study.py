"""Frozen known-K/time energy pilot on L97 current-environment voltages."""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
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

from . import alpha_mean_shape_study as base
from .alpha_pair_spatial_energy import common_inputs,fit_energy
from .alpha_current_pair_energy import CONDITIONS,selected_blocks,summarize,paired_statistics

NAME='alpha_current_pair_energy_20261102'
ROOT=Path(os.environ.get('GAGG_CURRENT_PAIR_ENERGY_ROOT',Path(__file__).resolve().parents[2])).resolve()
EXP=ROOT/'pileup/experiments'/NAME
ART=ROOT/'jupyter/artifacts'/NAME
REC=ROOT/'jupyter/run_records'/NAME
PARENT=ROOT/'pileup/experiments/alpha_current_shared_20261102'
PARENT_ART=ROOT/'jupyter/artifacts/alpha_current_shared_20261102'
read,save,sha=base.read,base.save,base.sha


def prepare():
    assert read(PARENT/'single-gate.json')['passed'] and read(PARENT/'pair-gate.json')['passed']
    cost=read(REC/'cost-preflight.json');assert len(cost['rows'])==6
    estimate=max(r['seconds'] for r in cost['rows'])*112/8*1.5
    assert estimate<3600,estimate
    EXP.mkdir(parents=True,exist_ok=False);ART.mkdir(parents=True,exist_ok=False)
    for part in ('voltage','fits','numerical','replay','figures'):(ART/part).mkdir()
    for part in ('bin','macros','original-runinfo'):(EXP/part).mkdir()
    shutil.copytree(PARENT/'source',EXP/'source',ignore=shutil.ignore_patterns('__pycache__'))
    for name in ('alpha_current_pair_energy.py','alpha_current_pair_energy_study.py','test_alpha_current_pair_energy.py'):
        shutil.copyfile(ROOT/'pileup/analysis'/name,EXP/'source/pileup/analysis'/name)
    models=('template.npz','energy-model.json','response-parameters.json','mean-correction.npz',
            'noise-calibration.json','sensor-config.json','charge-map-spatial12.npz','spatial-profile.npz',
            'shape-correction.npz','nominal-influences.npz')
    for name in models:shutil.copyfile(PARENT/name,EXP/name)
    save(EXP/'model-freeze.json',{str((EXP/n).relative_to(ROOT)):sha(EXP/n) for n in models})
    shutil.copyfile(PARENT/'bin/age.so',EXP/'bin/age.so')
    for name in ('build.json','analysis-environment.json','transport-environment.json'):
        shutil.copyfile(PARENT/name,EXP/name)
    for name in ('literature-evidence.json','cost-design.json','cost-preflight.json','cost_preflight.py'):
        shutil.copyfile(REC/name,EXP/name)
    (EXP/'pilot-fits').mkdir()
    for p in REC.glob('cost-fit-*.json'):shutil.copyfile(p,EXP/'pilot-fits'/p.name)
    mothers=read(PARENT/'mothers.json');plans=read(PARENT/'pair-plan.json');blocks=selected_blocks()
    cases=[]
    for b in blocks:
        for condition,gap in CONDITIONS:
            task=next(t for t in plans if t['block']==b and t['condition']==condition and t['gap_ns']==gap and t['phase_ns']==0)
            cases.append(dict(id=len(cases),kind='pair',block=b,condition=condition,gap_ns=gap,
                source_id=task['id'],parents=task['parents'],times_ns=task['times_ns'],
                repeats=[0,-1] if condition=='high_weak' and gap in (8.,40.) else [0],
                energy_dep_MeV=[mothers[i]['edep_MeV'] for i in task['parents']],
                source_voltage=str((PARENT_ART/'voltage/pair'/f"{task['id']:04d}.npz").relative_to(ROOT))))
    for g in range(4):
        for b in blocks[:2]:
            mother=8+32*g+b;sid=mother*4
            cases.append(dict(id=len(cases),kind='single',block=b,condition=f'single_group{g}',group=g,gap_ns=None,
                source_id=sid,parents=[mother],times_ns=[128.],repeats=[0,-1],energy_dep_MeV=[mothers[mother]['edep_MeV']],
                source_voltage=str((PARENT_ART/'voltage/single'/f'{sid:04d}.npz').relative_to(ROOT))))
    for c in cases:
        target=ART/'voltage'/f"{c['id']:04d}.npz";shutil.copyfile(ROOT/c['source_voltage'],target)
        c['voltage']=str(target.relative_to(ROOT));c['voltage_sha256']=sha(target)
    save(EXP/'cases.json',cases)
    used=sorted({p for c in cases for p in c['parents']});catalog=[]
    for i in used:
        row=mothers[i];macro=EXP/'macros'/Path(row['macro']).name
        shutil.copyfile(ROOT/row['macro'],macro)
        shutil.copyfile(PARENT/'original-runinfo'/f'{i:04d}.json',EXP/'original-runinfo'/f'{i:04d}.json')
        catalog.append(dict(row,snapshot_macro=str(macro.relative_to(ROOT))))
    save(EXP/'input-mothers.json',catalog)
    tasks=[]
    for c in cases:
        for repeat in c['repeats']:
            for method in ('old','corrected'):
                tasks.append(dict(id=len(tasks),case=c['id'],kind=c['kind'],times_ns=c['times_ns'],
                    repeat=repeat,method=method,voltage=c['voltage'],voltage_sha256=c['voltage_sha256']))
    assert len(tasks)==144 and sum(t['kind']=='pair' for t in tasks)==112
    save(EXP/'fit-plan.json',tasks)
    save(EXP/'protocol.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),
        base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        git_status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),
        scope='Known true K/time energy pilot; same L97 waveforms already exposed to response diagnostics, not new independent confirmation',
        cohort=dict(blocks=blocks,selection='Lowest8 SHA256(current-pair-energy-20261102:block) ranks; no outcome selection',available=32,used=8),
        pair_observations='40individual waveforms:8pairs x5conditions xrepeat0.16eight-repeat mean waveforms:8pairs xHW8/40 only.112paired-method fits.',
        single_controls='8phase0 cases:4source groups xfirst2 selected blocks;repeat0 and8repeatmean,2methods=32fits;exactly compare saved L97 fits',
        phase_ns=0.,known_time_definition='First active photon aligned128ns; declared second time=128+gap, not a blind timestamp estimate',
        previous_repeats='Every copied input contains8saved repeats; individual inference uses only index0; means use all8, overlap with individual, not independent',
        new_ROOTs=0,new_digitizations=0,new_training=0,fixed_models='Unchanged L93 old spatial12 and mean correction; original shared-minus-isolated term, Q map, noise model and physical parameters',
        search='Inherited fit_energy unchanged:common voltage-only NNLS+charge-shares0.1/0.9 starts,logE bounds,analytic gradients,max_nfev80/start,ftol/xtol1e-8,gtol1e-5;minimum cost including unconverged candidates',
        metric='Inherited common learned_metric computed from old voltage initial E and known times; identical weights for both models, fixed during fitting, not a newly calibrated full shared likelihood',
        isolation='Worker jobs contain only waveform path/hash,known times,kind,repeat,method and IDs. Audit hooks reject roots,labels,mothers,photon/truth/prediction arrays; truth E only in numerical synthetic controls and scoring',
        numerical=dict(gradient_relative_limit=3e-4,component_mean_limit_mV=1e-8,inverse_relative_limit=.001,synthetic_fits=4),
        metrics='Both E_dep bias/RMSE/errorSD;strict pair both relativeE within10% AND converged;boundaries/nonconvergence retained;no blind efficiency/miss/split/merge measurements',
        uncertainty='10000paired whole-mother resamples,seed741900000;all conditions/methods and the individual/mean observations travel with the same block;conditional descriptive intervals with only8mothers,not simultaneous coverage',
        gates=dict(short_strict_gain_pp=10.,short_pooled_weak_RMSE_ratio=.95,each_short_weak_RMSE_ratio=1.1,
            each_short_mean_weak_abs_bias_percent=10.,mean_strong_abs_bias_percent=5.,nonconverged=0),
        gate_origin='Retain L89 short-gap advancement thresholds;add explicit strong-mean and convergence controls;project decisions,not literature guarantees',
        decision='Only if all gates pass plan a larger frozen confirmation;otherwise retain failures and diagnose on development data. No default replacement,unknownK/t or high-rate expansion in this stage',
        cost_preflight=cost,projected_pair_seconds=estimate,cost_limit_seconds=3600.,workers=8,BLAS_threads=1,
        literature='Read-only paper-search:existing Fixsen2004/Zhang2024 notes and original review context;no new PDF review,web retrieval or human approval'))
    for name in ('README.md','docs/literature-use.md','jupyter/codex_0007_SiPM_Energy_Calibration.ipynb',
                 'jupyter/Research_Plan_20260921.ipynb','jupyter/Run_Catalog.ipynb'):
        p=ROOT/name;shutil.copyfile(p,REC/(p.name+'.before'))
    (EXP/'run.py').write_text("from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\n"
        "os.environ['GAGG_CURRENT_PAIR_ENERGY_ROOT']=str(f.parents[2])\nos.environ['OPENBLAS_NUM_THREADS']='1'\nos.environ['OMP_NUM_THREADS']='1'\n"
        "sys.path.insert(0,str(f/'source'))\nif __name__=='__main__':\n    from pileup.analysis.alpha_current_pair_energy_study import main\n    main()\n")
    paths=[p for p in EXP.rglob('*') if p.is_file() and 'bin' not in p.relative_to(EXP).parts]
    save(EXP/'design-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})
    save(EXP/'input-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in (ART/'voltage').glob('*.npz')})
    print('Frozen blocks',blocks,'112pair+32single fits;projected',round(estimate),'seconds',flush=True)


def verify():
    counts={}
    for name in ('design-freeze.json','model-freeze.json','input-freeze.json'):
        data=read(EXP/name)
        for path,digest in data.items():assert sha(ROOT/path)==digest,path
        counts[name]=len(data)
    assert sha(EXP/'bin/age.so')==read(EXP/'build.json')['age_binary_sha256']
    return counts


def initialize(guard=False):
    base.ROOT,base.EXP,base.ART,base.REC=ROOT,EXP,ART,REC;base.initialize(True)
    if not guard:return
    def audit(event,args):
        if event!='open' or not isinstance(args[0],(str,bytes,os.PathLike)):return
        p=Path(os.fsdecode(args[0])).resolve()
        forbidden=p.suffix=='.root' or p.name in ('cases.json','mothers.json','input-mothers.json')
        forbidden|=p.suffix=='.npz' and not p.is_relative_to(ART/'voltage')
        forbidden|=p.is_relative_to(EXP/'original-runinfo')
        if forbidden:raise PermissionError('Energy inference cannot read event truth: '+str(p))
    sys.addaudithook(audit)
    for p in (EXP/'cases.json',EXP/'input-mothers.json',ART/'truth/0000.npz',ROOT/'g4/data/forbidden.root'):
        try:p.open('rb')
        except PermissionError:pass
        else:raise AssertionError('Truth guard failed')


def fit_one(job):
    task,output=job
    assert set(task)=={'id','case','kind','times_ns','repeat','method','voltage','voltage_sha256'}
    path=ROOT/task['voltage'];assert sha(path)==task['voltage_sha256']
    with np.load(path) as data:all_y=data['voltage_mV']
    y=all_y[task['repeat']] if task['repeat']>=0 else all_y.mean(0)
    common=common_inputs(base.OLD,y,task['times_ns'],base.NOISE)
    signature=hashlib.sha256(common['metric'].factor.tobytes()+common['metric'].modes.tobytes()).hexdigest()
    model=base.OLD if task['method']=='old' else base.NEW
    result=fit_energy(model,y,task['times_ns'],common)
    result.update(task,metric_signature=signature,voltage_array_sha256=hashlib.sha256(y.tobytes()).hexdigest(),
        known_K=len(task['times_ns']),true_E_input=False,truth_guard=True)
    save(ART/output/f"{task['id']:04d}.json",result)
    return dict(id=task['id'],kind=task['kind'],seconds=result['seconds'],fits=1)


def inverse_one(job):
    i,method,times,energy=job;model=base.OLD if method=='old' else base.NEW
    y=model.components(times,energy).sum(0)+model.background
    common=common_inputs(base.OLD,y,times,base.NOISE)
    result=fit_energy(model,y,times,common)
    result.update(id=i,method=method,times_ns=times,synthetic_energy_MeV=energy,
                  max_relative_error=float(np.max(abs(np.asarray(result['energy_MeV'])/energy-1))))
    save(ART/'numerical'/f'{i:04d}.json',result);return result


def numerical():
    verify();assert not (EXP/'numerical.json').exists();initialize()
    rows=[];jobs=[]
    for method,model in [('old',base.OLD),('corrected',base.NEW)]:
        for times,energy in [([128.,136.],[7.2,.25]),([128.,168.],[.25,7.2])]:
            loge=np.log(energy);p=model.prepare(times)
            mean,grad=model.mean_and_gradient(times,loge,p,strength=1.)
            component_error=float(np.max(abs(mean-model.components(times,energy).sum(0))))
            h=2e-5;direction=np.array([.7,-.3])
            finite=(model.components(times,np.exp(loge+h*direction)).sum(0)-model.components(times,np.exp(loge-h*direction)).sum(0))/(2*h)
            error=float(np.linalg.norm(grad.T@direction-finite)/np.linalg.norm(finite))
            rows.append(dict(method=method,times_ns=times,gradient_relative_error=error,component_mean_max_error_mV=component_error))
            assert component_error<1e-8 and error<3e-4,rows[-1]
            jobs.append((len(jobs),method,times,energy))
    with ProcessPoolExecutor(max_workers=4,initializer=initialize) as pool:inverse=list(pool.map(inverse_one,jobs))
    passed=all(r['converged'] and r['max_relative_error']<.001 for r in inverse)
    save(EXP/'numerical.json',dict(passed=passed,gradient_rows=rows,inverse_rows=inverse,fits=4));assert passed
    print('Gradient,component identity and4synthetic inversions passed',flush=True)


def run_fits(kind):
    verify();assert read(EXP/'numerical.json')['passed']
    if kind=='pair':assert read(EXP/'single-control.json')['passed']
    tasks=[t for t in read(EXP/'fit-plan.json') if t['kind']==kind]
    assert not (EXP/f'{kind}-fits.json').exists()
    for t in tasks:assert not (ART/'fits'/f"{t['id']:04d}.json").exists()
    started=time.monotonic();rows=[]
    with ProcessPoolExecutor(max_workers=8,initializer=initialize,initargs=(True,)) as pool:
        pending={pool.submit(fit_one,(t,'fits')):t['id'] for t in tasks}
        for future in as_completed(pending):
            rows.append(future.result())
            if len(rows)%8==0:print(kind,'fits',len(rows),'/',len(tasks),flush=True)
    save(EXP/f'{kind}-fits.json',dict(rows=rows,fits=len(rows),seconds=time.monotonic()-started))
    save(EXP/f'{kind}-fit-freeze.json',{str((ART/'fits'/f"{t['id']:04d}.json").relative_to(ROOT)):sha(ART/'fits'/f"{t['id']:04d}.json") for t in tasks})
    if kind=='single':
        cases=read(EXP/'cases.json');checked=0
        for task in tasks:
            current=read(ART/'fits'/f"{task['id']:04d}.json")
            previous=read(PARENT_ART/'fits'/f"{cases[task['case']]['source_id']:04d}.json")
            row=next(r for r in previous['rows'] if r['repeat']==task['repeat'] and r['method']==task['method'])
            for key,value in row.items():
                if key not in ('seconds','known_time_ns'):assert current[key]==value,(key,task)
            checked+=1
        save(EXP/'single-control.json',dict(passed=True,exact_fits_except_timing=checked,reused_previous_isolated_gates=True))


def score():
    verify();initialize()
    for name in ('single-fit-freeze.json','pair-fit-freeze.json'):
        for path,digest in read(EXP/name).items():assert sha(ROOT/path)==digest,path
    cases=read(EXP/'cases.json');tasks=read(EXP/'fit-plan.json');details=[]
    for task in tasks:
        row=read(ART/'fits'/f"{task['id']:04d}.json");c=cases[task['case']]
        assert row['truth_guard'] and not row['true_E_input']
        if c['kind']=='single':continue
        weak=int(np.argmin(c['energy_dep_MeV']));order=[weak,1-weak]
        error=100*(np.asarray(row['energy_MeV'])/c['energy_dep_MeV']-1)
        details.append(dict(id=task['id'],case=c['id'],block=c['block'],condition=c['condition'],gap_ns=c['gap_ns'],
            method=task['method'],kind='individual' if task['repeat']>=0 else 'mean',
            weak_error_percent=float(error[weak]),strong_error_percent=float(error[1-weak]),
            weak_true_MeV=c['energy_dep_MeV'][weak],strong_true_MeV=c['energy_dep_MeV'][1-weak],
            weak_fit_MeV=row['energy_MeV'][weak],strong_fit_MeV=row['energy_MeV'][1-weak],
            strict=bool(np.all(abs(error)<=10) and row['converged']),converged=row['converged'],boundary=row['boundary'],
            selected_start=row['selected_start'],evaluations=row['evaluations'],cost=row['cost'],seconds=row['seconds']))
    blocks=selected_blocks();index=np.random.default_rng(741900000).integers(8,size=(10000,8))
    groups=[]
    def selected(method,kind,condition,gap):
        selected_rows=[r for r in details if r['method']==method and r['kind']==kind and r['condition']==condition and r['gap_ns']==gap]
        selected_rows.sort(key=lambda r:blocks.index(r['block']));assert len(selected_rows)==8
        return selected_rows
    def metrics(rows):
        return summarize([[r['weak_error_percent'],r['strong_error_percent']] for r in rows],
            [r['strict'] for r in rows],[r['boundary'] for r in rows],[r['converged'] for r in rows])
    for condition,gap in CONDITIONS:
        for kind in (['individual','mean'] if condition=='high_weak' and gap in (8.,40.) else ['individual']):
            before=selected('old',kind,condition,gap);after=selected('corrected',kind,condition,gap)
            comparison=paired_statistics(np.array([r['weak_error_percent'] for r in before])[:,None],
                np.array([r['weak_error_percent'] for r in after])[:,None],
                np.array([r['strict'] for r in before])[:,None],np.array([r['strict'] for r in after])[:,None],index)
            groups.append(dict(condition=condition,gap_ns=gap,kind=kind,mother_pairs=8,
                               old=metrics(before),corrected=metrics(after),comparison=comparison))
    pooled=[]
    for kind in ('individual','mean'):
        rows={m:[selected(m,kind,'high_weak',gap) for gap in (8.,40.)] for m in ('old','corrected')}
        errors={m:np.array([[r['weak_error_percent'] for r in condition] for condition in conditions]).T for m,conditions in rows.items()}
        success={m:np.array([[r['strict'] for r in condition] for condition in conditions]).T for m,conditions in rows.items()}
        comparison=paired_statistics(errors['old'],errors['corrected'],success['old'],success['corrected'],index)
        pooled.append(dict(kind=kind,old=metrics(sum(rows['old'],[])),corrected=metrics(sum(rows['corrected'],[])),comparison=comparison,mother_pairs=8))
    primary=pooled[0]['comparison'];failures=[]
    if primary['strict_gain_pp']<10:failures.append('short strict-pair gain <10pp')
    if primary['weak_rmse_ratio']>.95:failures.append('short pooled weak RMSEratio >0.95')
    for row in groups:
        if row['condition']!='high_weak' or row['gap_ns'] not in (8.,40.):continue
        if row['kind']=='individual' and row['comparison']['weak_rmse_ratio']>1.1:
            failures.append(f"{row['gap_ns']:g}ns individual weak RMSEratio >1.1")
        if row['kind']=='mean':
            if abs(row['corrected']['weak_bias_percent'])>10:failures.append(f"{row['gap_ns']:g}ns mean weak bias outside10percent")
            if abs(row['corrected']['strong_bias_percent'])>5:failures.append(f"{row['gap_ns']:g}ns mean strong bias outside5percent")
    if any(not r['converged'] for r in details):failures.append('nonconverged pair fit')
    save(EXP/'fit-details.json',details);save(EXP/'result-summary.json',dict(groups=groups,short_pooled=pooled,blocks=blocks,
        pair_fits=len(details),single_fits=32,noisy_pair_waveforms=40,mean_pair_waveforms=16,new_ROOTs=0,new_digitizations=0))
    save(EXP/'gate.json',dict(passed=not failures,failures=failures,default_changed=False,
        next='Freeze a larger confirmation on reserved mother pairs before changing inference scope' if not failures else
             'Stop expansion;retain these energy failures and diagnose on development data without retuning on the pilot'))
    examples={'time_ns':base.OLD.readout.adc_centers,'background_mV':base.OLD.background}
    for cid in (0,1):
        c=cases[cid]
        with np.load(ART/'voltage'/f'{cid:04d}.npz') as data:examples[f'case{cid}_voltage_mV']=data['voltage_mV'][0]
        with np.load(PARENT_ART/'truth/pair'/f"{c['source_id']:04d}.npz") as data:
            examples[f'case{cid}_truth_components_mV']=data['components_mV'][0]
        for method,model in [('old',base.OLD),('corrected',base.NEW)]:
            task=next(t for t in tasks if t['case']==cid and t['repeat']==0 and t['method']==method)
            fit=read(ART/'fits'/f"{task['id']:04d}.json")
            examples[f'case{cid}_{method}_components_mV']=model.components(c['times_ns'],fit['energy_MeV'])
    np.savez_compressed(ART/'examples.npz',**examples)
    paths=[EXP/'fit-details.json',EXP/'result-summary.json',EXP/'gate.json',ART/'examples.npz']
    save(EXP/'result-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})
    print(json.dumps(read(EXP/'gate.json'),ensure_ascii=False,indent=2),flush=True)


def replay():
    verify();assert not (EXP/'replay.json').exists()
    tasks=read(EXP/'fit-plan.json');chosen=[tasks[i] for i in (0,1,112,113,114,115)]
    with ProcessPoolExecutor(max_workers=6,initializer=initialize,initargs=(True,)) as pool:
        rows=list(pool.map(fit_one,[(t,'replay') for t in chosen]))
    for t in chosen:
        a=read(ART/'fits'/f"{t['id']:04d}.json");b=read(ART/'replay'/f"{t['id']:04d}.json")
        assert {k:v for k,v in a.items() if k!='seconds'}=={k:v for k,v in b.items() if k!='seconds'},t
    save(EXP/'replay.json',dict(passed=True,fits=6,pair_fits=2,single_fits=4,rows=rows,exact_except_timing=True))
    print('Replay passed:2pair+4single fits identical except timing',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['prepare','verify','numerical','single','pair','score','replay'])
    stage=parser.parse_args().stage
    if stage=='prepare':prepare()
    elif stage=='verify':print(json.dumps(verify()))
    elif stage=='numerical':numerical()
    elif stage in ('single','pair'):run_fits(stage)
    elif stage=='score':score()
    elif stage=='replay':replay()


if __name__=='__main__':main()
