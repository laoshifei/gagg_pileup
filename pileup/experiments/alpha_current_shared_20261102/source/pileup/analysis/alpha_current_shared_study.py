"""Fresh four-energy/phase and staged shared-cell response controls."""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import replace
from datetime import datetime,timezone
import os
from pathlib import Path
import shutil
import subprocess
import threading
import time
import json
import numpy as np

from . import alpha_mean_shape_study as frozen
from .alpha_current_response import energy_metrics
from .alpha_current_shared import PHASES,ENERGIES,CONDITIONS,make_plans,anchored_pair_predictions,paired_mean_ci
from .alpha_transport_environment import relocate_macro
from .alpha_transport_environment_study import load_root
from .sipm_response import SensorConfig,read_photon_hits,map_cells

NAME='alpha_current_shared_20261102'
ROOT=Path(os.environ.get('GAGG_CURRENT_SHARED_ROOT',Path(__file__).resolve().parents[2])).resolve()
EXP=ROOT/'pileup/experiments'/NAME
ART=ROOT/'jupyter/artifacts'/NAME
REC=ROOT/'jupyter/run_records'/NAME
PARENT=ROOT/'pileup/experiments/alpha_current_confirmation_20261102'
read,save,sha=frozen.read,frozen.save,frozen.sha


def initialize(mode='normal'):
    frozen.ROOT,frozen.EXP,frozen.ART,frozen.REC=ROOT,EXP,ART,REC
    if mode=='fit':frozen.fit_initialize()
    else:
        frozen.initialize(True)
        if mode=='off':frozen.SENSOR=replace(frozen.SENSOR,recovery_enabled=False)


def prepare():
    assert read(PARENT/'gate.json')['passed']
    for p in (EXP,ART,REC):p.mkdir(parents=True,exist_ok=False)
    for part in ('bin','build','macros','original-runinfo'):(EXP/part).mkdir()
    for part in ('photons','voltage/single','truth/single','voltage/pair','truth/pair','fits','predictions','figures'):
        (ART/part).mkdir(parents=True)
    for part in ('logs','attempts'):(REC/part).mkdir()
    shutil.copytree(PARENT/'source',EXP/'source',ignore=shutil.ignore_patterns('__pycache__'))
    for name in ('alpha_current_shared.py','alpha_current_shared_study.py','test_alpha_current_shared.py'):
        shutil.copyfile(ROOT/'pileup/analysis'/name,EXP/'source/pileup/analysis'/name)
    models=('template.npz','energy-model.json','response-parameters.json','mean-correction.npz',
            'noise-calibration.json','sensor-config.json','charge-map-spatial12.npz','spatial-profile.npz',
            'shape-correction.npz','nominal-influences.npz')
    for name in models:shutil.copyfile(PARENT/name,EXP/name)
    for name in ('gagg','age.so'):shutil.copyfile(PARENT/'bin'/name,EXP/'bin'/name)
    (EXP/'bin/gagg').chmod(0o755)
    build=read(PARENT/'build.json')
    for path,digest in build['native_build']['source_sha256'].items():
        assert sha(ROOT/path)==digest and sha(EXP/'source'/path)==digest
    assert sha(EXP/'bin/gagg')==build['native_binary_sha256']
    assert sha(EXP/'bin/age.so')==build['age_binary_sha256']
    save(EXP/'build.json',build)
    save(EXP/'model-freeze.json',{str((EXP/name).relative_to(ROOT)):sha(EXP/name) for name in models})
    for name in ('transport-environment.json','analysis-environment.json'):
        shutil.copyfile(PARENT/name,EXP/name)
    save(EXP/'preceding-gate.json',read(PARENT/'gate.json'))
    save(EXP/'preceding-result.json',read(PARENT/'result-summary.json'))
    shutil.copyfile(PARENT/'macro-template.mac',EXP/'macro-template.mac')
    native,single,pairs=make_plans()
    for name,plan in [('native',native),('single',single),('pair',pairs)]:save(EXP/(name+'-plan.json'),plan)
    pilot=[]
    for gap in (8.,40.,1200.):
        for phase in (0.,1.5):
            pilot.append(dict(id=len(pilot),kind='pair',block=-1,condition='high_weak',gap_ns=gap,phase_ns=phase,
                parents=[6,0],times_ns=[128.+phase,128.+phase+gap],seed_base=740700000,repeats=2))
    save(EXP/'pair-pilot-plan.json',pilot)
    streams={r[k] for r in native for k in ('seed','sipm_seed','readout_seed')};audit=[]
    for p in sorted((ROOT/'pileup/experiments').glob('*/manifest.json')):
        data=read(p)
        if not isinstance(data,dict):continue
        old={r['seed']+offset for k in ('runs','planned_runs') for r in data.get(k,[])
             if isinstance(r,dict) and isinstance(r.get('seed'),int) for offset in (0,100000,200000)}
        assert not old&streams,p;audit.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
    save(EXP/'seed-audit.json',dict(passed=True,native_streams=len(streams),manifests=audit,
         scope='Native manifest seeds and standard avalanche/readout offsets; new offline7403/7405/7407 streams separated'))
    save(EXP/'protocol.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),
        base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        git_status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),
        purpose='Fresh current-environment four-energy/four-phase singles followed by shared-cell forward controls; no pair-energy fitting',
        native='136new roots:2pilot perenergy=8;32formal perenergy=128;source1.625/2.5/5.486/8MeV; no old mothers reused',
        physical='Unchanged39G4 files,95/351ns fastfraction.3,29ns recoveryproxy,14400cells;Geant4 11.3.2/ROOT6.34.10',
        native_window='One20us event perROOT,onecentral normal alpha at5us;no Poisson rate or background;G4 serial,max8jobs,>=1.05s launch spacing',
        output='Only g4/data/YYYYMMDD_HHhMMmSSs.root; no overwrite; actual macro/config/RunInfo saved',
        fixed_models='All L93 models copied unchanged;no new training,shape smoothing,noise learning or per-event nuisance parameter',
        single='128formal mothers x4phases x8repeats=4096 digitizations;each8individual+mean x2models=9216fits;knownK1/time,voltage-only inference',
        single_pilot='8mothers x4phases x8=256digitizations;fourphase0groups x18=72fits;only computational/structure cost gate',
        single_gates=dict(absolute_group_mean_energy_bias_percent=5.,corrected_to_old_individual_RMSE=1.1,
            corrected_phase_bias_span_pp=1.,group_mean_waveform_MSE_ratio=1.1,max_nonconverged=0),
        single_scope='All groups/phases; mean-waveform MSE gate uses square of group-mean forward residual. Strong weak-direction projection was independently checked in L96; not a new small-N gate here',
        shared='32disjoint weak/strong motherpairs x5conditions x4phases=640cases;8repeats of(shared+2isolated)=15360digitizations',
        conditions=[dict(order=order,gap_ns=gap) for order,gap in CONDITIONS],phases_ns=list(PHASES),
        recovery_off='Same32pairs and5conditions atphase0 only;160cases x8 x3=3840digitizations;diagnostic sensor.recovery_enabled True->False, not model tuning',
        shared_gates=dict(abs_mean_weak_window_area_error_pp=5.,corrected_to_old_per_mother_component_MSE=1.1,
                          recovery_off_abs_mean_weak_window_area_change_pp=3.),
        area_definition='Ratio of mean component voltage integrals inside observed0..2048ns;not complete charge. Late1200ns weak pulse is partly outside observation',
        mse_definition='Eachmotherpair8repeatmean residual squared over1024samples,then average32pairs;gate both weak/strong components at each condition/phase',
        confidence='10000whole-mother-pair paired draws percondition/phase,seed740900000+conditionID*4+phaseID; pointwise descriptive intervals, gates use inherited point estimates; no simultaneous coverage claim',
        pair_prediction='True E/time forwarded as explicit diagnostic inputs;all other predictor inputs frozen. Reuse old shared-minus-isolated term algebraically only after direct corrected-model equality on pilot inputs',
        observation='4096ns offline;first2048ns1024samples;2nsperiod,.5nsaperture;firstactivephoton at128+phase. All photons retained before sensor footprint rejection',
        randomness='Single740300000+motherID*100+repeat;pair740500000+block*1000+repeat;pilotpair740700000+repeat;readout+10000000. Common seeds across phase/gap do not imply identical avalanches',
        pilot='Nativepreflight ID6;nativepilot remaining7;all8readout audited;shared pilot6cases x2repeats x3=36, direct model comparisons; no efficacy selection on pilot',
        cost_limit_seconds=dict(native=900,single=900,fit=900,pair=1800,prediction=900),workers=8,BLAS_threads=1,
        stages=['native preflight','native pilot','isolated costpilot','native formal','singles','singlefit','single gate',
                'shared pilot only if singlepass','fixed predictions','shared on/off','shared gate','replay'],
        decision='Stop shared expansion if any single gate fails;no pair energy fit or default replacement even if response passes',
        literature='Reuse L91/L93/L96 recorded method conditions; no new literature claim or PDF review',default_changed=False))
    for name in ('README.md','docs/literature-use.md','jupyter/codex_0007_SiPM_Energy_Calibration.ipynb',
                 'jupyter/Research_Plan_20260921.ipynb','jupyter/Run_Catalog.ipynb'):
        p=ROOT/name;shutil.copyfile(p,REC/(p.name+'.before'))
    (EXP/'run.py').write_text("from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\n"
        "os.environ['GAGG_CURRENT_SHARED_ROOT']=str(f.parents[2])\nos.environ['OPENBLAS_NUM_THREADS']='1'\nos.environ['OMP_NUM_THREADS']='1'\n"
        "sys.path.insert(0,str(f/'source'))\nif __name__=='__main__':\n    from pileup.analysis.alpha_current_shared_study import main\n    main()\n")
    paths=[p for p in EXP.rglob('*') if p.is_file() and not {'build','bin'}&set(p.relative_to(EXP).parts)]
    save(EXP/'design-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})
    save(EXP/'manifest.json',dict(planned_runs=native,runs=[],failed_attempts=[]))
    print('Frozen136new mothers,4096formal singles and conditional19200shared/off digitizations',flush=True)


def verify():
    counts={}
    for name in ('design-freeze.json','model-freeze.json'):
        data=read(EXP/name)
        for path,digest in data.items():assert sha(ROOT/path)==digest,path
        counts[name]=len(data)
    build=read(EXP/'build.json')
    assert sha(EXP/'bin/gagg')==build['native_binary_sha256']
    assert sha(EXP/'bin/age.so')==build['age_binary_sha256']
    return counts


def run_native(stage):
    verify();manifest=read(EXP/'manifest.json')
    if stage=='native_pilot':assert read(EXP/'native-cost-gate.json')['passed']
    if stage=='native_formal':assert read(EXP/'pilot-gate.json')['passed']
    done={r['id'] for r in manifest['runs']}
    tasks=[r for r in manifest['planned_runs'] if r['id'] not in done and
           (r['id']==6 if stage=='preflight' else r['role']==('pilot' if stage=='native_pilot' else 'formal'))]
    assert tasks,'Stage already completed'
    env=os.environ.copy();env.update(read(EXP/'transport-environment.json')['current_environment']['g4_dataset_environment'])
    sensor=SensorConfig(**read(EXP/'sensor-config.json')['sensor'])
    last=[max([r['started_epoch'] for r in manifest['runs']]+[0.])];lock=threading.Lock()
    def one(task):
        with lock:
            while True:
                delay=1.05-(time.time()-last[0])
                if delay>0:time.sleep(delay)
                root=ROOT/'g4/data'/(datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss')+'.root')
                if not root.exists():break
                time.sleep(1.05)
            macro=EXP/'macros'/f"mother{task['id']:04d}_{root.stem}.mac"
            body=relocate_macro((EXP/'macro-template.mac').read_text(),root,
                                EXP/'source/g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json')
            replacements={'/gagg/run/seed':str(task['seed']),'/gagg/sipm/seed':str(task['sipm_seed']),
                          '/gagg/sipm/readoutSeed':str(task['readout_seed']),
                          '/gun/energy':f"{task['source_energy_MeV']:.8g} MeV"}
            found=dict.fromkeys(replacements,0);lines=[]
            for line in body.splitlines():
                command=line.strip().split(maxsplit=1)[0] if line.strip() else ''
                if command in replacements:found[command]+=1;line=command+' '+replacements[command]
                lines.append(line)
            assert set(found.values())=={1}
            body='\n'.join(lines)+'\n';macro.write_text(body)
            log=REC/'logs'/f'{root.stem}.log';stream=log.open('w');started=time.monotonic()
            proc=subprocess.Popen([str(EXP/'bin/gagg'),str(macro)],cwd=EXP/'build',env=env,stdout=stream,stderr=subprocess.STDOUT)
            last[0]=time.time()
            row=dict(task,root=str(root.relative_to(ROOT)),macro=str(macro.relative_to(ROOT)),macro_sha256=sha(macro),
                     log=str(log.relative_to(ROOT)),started_epoch=last[0],started_utc=datetime.now(timezone.utc).isoformat())
        code=proc.wait();stream.close();row.update(returncode=code,seconds=time.monotonic()-started)
        save(REC/'attempts'/f'{root.stem}.json',row);assert code==0,row
        info,primary,photons,voltage,scalars=load_root(root)
        assert info['random_seed']==task['seed'] and info['sipm_avalanche_base_seed']==task['sipm_seed']
        assert info['sipm_readout_base_seed']==task['readout_seed'] and info['macro_contents']==body
        assert info['t_length_ps']==20000000 and primary['primary_time_ps']==5000000
        assert np.isclose(primary['primary_energy_MeV'],task['source_energy_MeV'],rtol=2e-7,atol=1e-9)
        assert 0<primary['edep_total_MeV']<=task['source_energy_MeV']
        assert json.loads(info['sipm_config_json'])==read(EXP/'source/g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json')
        assert abs(scalars['charge_balance_error_C'])<1e-18
        assert all(np.isfinite(a).all() for a in voltage.values())
        hits=read_photon_hits(root);active=map_cells(hits,sensor)>=0;assert active.any()
        first=float(hits.time_ns[active].min())
        np.savez_compressed(ART/'photons'/f"{task['id']:04d}.npz",time_ns=hits.time_ns-first,
                            x_mm=hits.x_mm,y_mm=hits.y_mm,energy_eV=hits.energy_eV)
        save(EXP/'original-runinfo'/f"{task['id']:04d}.json",info)
        row.update(root_sha256=sha(root),root_bytes=root.stat().st_size,status='completed',
                   edep_MeV=primary['edep_total_MeV'],first_active_ns=first,photons=len(hits.time_ns),active_photons=int(active.sum()))
        save(REC/'attempts'/f'{root.stem}.json',row);return row
    failures=[]
    with ThreadPoolExecutor(max_workers=1 if stage=='preflight' else 8) as pool:
        futures={pool.submit(one,t):t for t in tasks}
        for future in as_completed(futures):
            try:row=future.result()
            except Exception as exc:
                error=dict(task=futures[future],error=repr(exc));manifest['failed_attempts'].append(error);failures.append(error)
            else:
                manifest['runs'].append(row)
                if len(manifest['runs'])<=8 or len(manifest['runs'])%16==0:
                    print(stage,len(manifest['runs']),'/136',round(row['seconds'],2),'s',flush=True)
            save(EXP/'manifest.json',manifest)
    save(EXP/'mothers.json',sorted(manifest['runs'],key=lambda r:r['id']))
    assert not failures,failures
    if stage=='preflight':
        cost=manifest['runs'][0]['seconds'];estimate=136*max(1.05,cost/8)*1.5
        save(EXP/'native-cost-gate.json',dict(passed=estimate<900,one_seconds=cost,projected_seconds=estimate));assert estimate<900
    if stage=='native_formal':
        assert len(manifest['runs'])==136
        paths=[ROOT/r['root'] for r in manifest['runs']]+list((ART/'photons').glob('*.npz'))
        paths += [EXP/'manifest.json',EXP/'mothers.json',*sorted((EXP/'macros').glob('*.mac')),*sorted((EXP/'original-runinfo').glob('*.json'))]
        save(EXP/'native-output-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})


def pilot():
    verify();mothers=read(EXP/'mothers.json');assert len(mothers)==8
    assert not (EXP/'pilot-gate.json').exists()
    from analysis.single_particle_sipm import audit_response,load_root as native_load
    audits=[dict(id=r['id'],**audit_response(native_load(ROOT/r['root']))) for r in mothers]
    tasks=[t for t in read(EXP/'single-plan.json') if t['role']=='pilot']
    with ProcessPoolExecutor(max_workers=8,initializer=initialize) as pool:
        digits=list(pool.map(frozen.simulate,[(task,'formal') for task in tasks]))
    with ProcessPoolExecutor(max_workers=4,initializer=initialize,initargs=('fit',)) as pool:
        fits=list(pool.map(frozen.fit_one,[(i,128.,'pilot-fits') for i in (0,8,16,24)]))
    estimates=dict(digitization_seconds=max(r['seconds'] for r in digits)*512/8*1.5,
                   fit_seconds=max(r['seconds'] for r in fits)*512/8*1.5)
    passed=all(v<900 for v in estimates.values())
    save(EXP/'pilot-gate.json',dict(passed=passed,readout_audits=audits,digitization_rows=digits,fit_rows=fits,
        projected_seconds=estimates,pilot_mothers=8,digitizations=256,fits=72,efficacy_outcomes_not_used=True))
    print(json.dumps(dict(passed=passed,**estimates)),flush=True);assert passed


def singles():
    verify();assert len(read(EXP/'mothers.json'))==136
    assert not (EXP/'single-generation.json').exists()
    tasks=[t for t in read(EXP/'single-plan.json') if t['role']=='formal'];started=time.monotonic()
    for task in tasks:assert not (ART/'voltage/single'/f"{task['id']:04d}.npz").exists()
    with ProcessPoolExecutor(max_workers=8,initializer=initialize) as pool:
        rows=[]
        for row in pool.map(frozen.simulate,[(task,'formal') for task in tasks]):
            rows.append(row)
            if len(rows)%64==0:print('single',len(rows),'/512',flush=True)
    save(EXP/'single-generation.json',dict(rows=rows,digitizations=4096,seconds=time.monotonic()-started))
    paths=[*sorted((ART/'voltage/single').glob('*.npz')),*sorted((ART/'truth/single').glob('*.npz'))]
    save(EXP/'single-output-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})


def fit_singles():
    verify()
    for path,digest in read(EXP/'single-output-freeze.json').items():assert sha(ROOT/path)==digest
    assert not list((ART/'fits').glob('*.json'))
    tasks=[t for t in read(EXP/'single-plan.json') if t['role']=='formal'];started=time.monotonic()
    with ProcessPoolExecutor(max_workers=8,initializer=initialize,initargs=('fit',)) as pool:
        rows=[]
        for row in pool.map(frozen.fit_one,[(t['id'],t['times_ns'][0],'fits') for t in tasks]):
            rows.append(row)
            if len(rows)%64==0:print('singlefit',len(rows),'/512',flush=True)
    save(EXP/'single-fits.json',dict(rows=rows,fits=9216,seconds=time.monotonic()-started))
    save(EXP/'single-fit-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in (ART/'fits').glob('*.json')})


def score_singles():
    verify();initialize()
    for name in ('native-output-freeze.json','single-output-freeze.json','single-fit-freeze.json'):
        for path,digest in read(EXP/name).items():assert sha(ROOT/path)==digest
    mothers=read(EXP/'mothers.json');tasks=[t for t in read(EXP/'single-plan.json') if t['role']=='formal']
    details=[];residuals={};actual={};failures=[]
    for task in tasks:
        mother=mothers[task['parents'][0]];fit=read(ART/'fits'/f"{task['id']:04d}.json")
        assert fit['truth_guard'] and not fit['true_E_input'] and fit['known_K1_t']
        for row in fit['rows']:
            details.append(dict(id=task['id'],mother=mother['id'],group=mother['group'],phase_ns=task['phase_ns'],
                method=row['method'],repeat=row['repeat'],relative_error_percent=100*(row['energy_MeV'][0]/mother['edep_MeV']-1),
                converged=row['converged'],boundary=row['boundary']))
        with np.load(ART/'truth/single'/f"{task['id']:04d}.npz") as data:signal=data['components_mV'].mean(0)[0]
        group=mother['group'];phase=task['phase_ns'];actual.setdefault((group,phase),[]).append(signal)
        for method,model in [('old',frozen.OLD),('corrected',frozen.NEW)]:
            pred=model.isolated(task['times_ns'],[mother['edep_MeV']])[0][0]
            residuals.setdefault((group,phase,method),[]).append(signal-pred)
    energy=[];wave=[];arrays={'time_ns':frozen.OLD.readout.adc_centers}
    for group in range(4):
        for phase in PHASES:
            for method in ('old','corrected'):
                for kind in ('individual','mean'):
                    rows=[r for r in details if r['group']==group and r['phase_ns']==phase and r['method']==method and
                          (r['repeat']>=0 if kind=='individual' else r['repeat']==-1)]
                    energy.append(dict(group=group,phase_ns=phase,method=method,kind=kind,mothers=32,
                        **energy_metrics([r['relative_error_percent'] for r in rows],[r['boundary'] for r in rows],[r['converged'] for r in rows])))
            def e(method,kind):return next(r for r in energy if r['group']==group and r['phase_ns']==phase and r['method']==method and r['kind']==kind)
            for method in ('old','corrected'):
                if abs(e(method,'mean')['bias_percent'])>5:failures.append(f'{group}/{phase}/{method}:meanE bias')
            if e('corrected','individual')['rmse_percent']>1.1*e('old','individual')['rmse_percent']:
                failures.append(f'{group}/{phase}:individual E RMSEratio')
            b,a=(np.asarray(residuals[(group,phase,m)]) for m in ('old','corrected'))
            old=float(np.mean(b.mean(0)**2));new=float(np.mean(a.mean(0)**2));ratio=new/old
            wave.append(dict(group=group,phase_ns=phase,old_group_mean_mse_mV2=old,corrected_group_mean_mse_mV2=new,
                group_mean_mse_ratio=ratio,old_per_mother_mse_mV2=float(np.mean(b*b)),corrected_per_mother_mse_mV2=float(np.mean(a*a))))
            if ratio>1.1:failures.append(f'{group}/{phase}:group-mean waveform MSEratio')
            arrays[f'g{group}_p{phase}_actual_mV']=actual[(group,phase)]
            arrays[f'g{group}_p{phase}_old_residual_mV']=b;arrays[f'g{group}_p{phase}_corrected_residual_mV']=a
        span=np.ptp([r['bias_percent'] for r in energy if r['group']==group and r['method']=='corrected' and r['kind']=='mean'])
        if span>1:failures.append(f'{group}:phase bias span')
    if any(r['nonconverged'] for r in energy):failures.append('unconverged fits')
    save(EXP/'single-details.json',details);save(EXP/'single-result.json',dict(energy=energy,waveform=wave))
    save(EXP/'single-gate.json',dict(passed=not failures,failures=failures,
        next='Frozen shared-response controls' if not failures else 'Stop before shared response; retain all isolated outcomes'))
    np.savez_compressed(ART/'single-summary.npz',**arrays)
    paths=[EXP/'single-result.json',EXP/'single-gate.json',EXP/'single-details.json',ART/'single-summary.npz']
    save(EXP/'single-result-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})
    print(json.dumps(read(EXP/'single-gate.json'),ensure_ascii=False,indent=2),flush=True)


def predict_one(job):
    task,energies,label,direct=job;started=time.monotonic()
    predictions=anchored_pair_predictions(frozen.OLD,frozen.NEW,task['times_ns'],energies)
    error=0.
    if direct:
        full=frozen.NEW.components(task['times_ns'],energies)
        error=float(np.max(abs(predictions['corrected_joint']-full)))
        assert error<1e-10,error
    destination=ART/label;destination.mkdir(exist_ok=True)
    np.savez_compressed(destination/f"{task['id']:04d}.npz",**predictions)
    return dict(id=task['id'],seconds=time.monotonic()-started,direct_checked=direct,max_direct_difference_mV=error)


def pair_pilot():
    verify();assert read(EXP/'single-gate.json')['passed']
    assert not (EXP/'pair-pilot-gate.json').exists()
    mothers=read(EXP/'mothers.json');tasks=read(EXP/'pair-pilot-plan.json')
    with ProcessPoolExecutor(max_workers=6,initializer=initialize) as pool:
        digits=list(pool.map(frozen.simulate,[(task,'pair-pilot') for task in tasks]))
    with ProcessPoolExecutor(max_workers=6,initializer=initialize) as pool:
        prediction=list(pool.map(predict_one,[(task,[mothers[i]['edep_MeV'] for i in task['parents']],
                                            'pilot-predictions',True) for task in tasks]))
    # Direct verification costs more than the frozen reuse route; conservative cost estimate.
    estimates=dict(digitization_seconds=max(r['seconds']/r['digitizations'] for r in digits)*19200/8*1.5,
                   prediction_seconds=max(r['seconds'] for r in prediction)*640/8*1.5)
    passed=estimates['digitization_seconds']<1800 and estimates['prediction_seconds']<900
    save(EXP/'pair-pilot-gate.json',dict(passed=passed,digitization_rows=digits,prediction_rows=prediction,
        estimates=estimates,digitizations=36,pilot_mothers=[0,6],no_formal_sample_used=True))
    print(json.dumps(dict(passed=passed,**estimates)),flush=True);assert passed


def predictions():
    verify();assert read(EXP/'pair-pilot-gate.json')['passed']
    assert not list((ART/'predictions').glob('*.npz'))
    mothers=read(EXP/'mothers.json');tasks=read(EXP/'pair-plan.json');started=time.monotonic()
    jobs=[(t,[mothers[i]['edep_MeV'] for i in t['parents']],'predictions',False) for t in tasks]
    with ProcessPoolExecutor(max_workers=8,initializer=initialize) as pool:
        rows=[]
        for row in pool.map(predict_one,jobs):
            rows.append(row)
            if len(rows)%80==0:print('predictions',len(rows),'/640',flush=True)
    save(EXP/'predictions.json',dict(rows=rows,seconds=time.monotonic()-started,
        inputs='Declared true deposited E/t plus frozen models only; no measured shared waveform inputs'))
    paths=[EXP/'predictions.json',*sorted((ART/'predictions').glob('*.npz'))]
    save(EXP/'prediction-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})


def pairs(mode):
    verify();assert read(EXP/'single-gate.json')['passed'] and read(EXP/'pair-pilot-gate.json')['passed']
    for path,digest in read(EXP/'prediction-freeze.json').items():assert sha(ROOT/path)==digest
    tasks=read(EXP/'pair-plan.json')
    if mode=='off':tasks=[t for t in tasks if t['phase_ns']==0]
    output='formal' if mode=='on' else 'off'
    assert not (EXP/f'pair-{mode}-generation.json').exists()
    prefix='' if mode=='on' else 'off/'
    for task in tasks:assert not (ART/(prefix+'voltage/pair')/f"{task['id']:04d}.npz").exists()
    started=time.monotonic()
    with ProcessPoolExecutor(max_workers=8,initializer=initialize,initargs=('normal' if mode=='on' else 'off',)) as pool:
        rows=[]
        for row in pool.map(frozen.simulate,[(task,output) for task in tasks]):
            rows.append(row)
            if len(rows)%40==0:print('shared recovery',mode,len(rows),'/',len(tasks),flush=True)
    save(EXP/f'pair-{mode}-generation.json',dict(rows=rows,digitizations=sum(r['digitizations'] for r in rows),seconds=time.monotonic()-started))
    paths=[p for part in ('voltage/pair','truth/pair') for p in (ART/(prefix+part)).glob('*.npz')]
    save(EXP/f'pair-{mode}-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})


def score_pairs():
    verify()
    for name in ('pair-on-freeze.json','pair-off-freeze.json','prediction-freeze.json'):
        for path,digest in read(EXP/name).items():assert sha(ROOT/path)==digest
    tasks=read(EXP/'pair-plan.json');mothers=read(EXP/'mothers.json');details=[];off_details=[];examples={}
    for task in tasks:
        i=task['id'];weak=int(np.argmin([mothers[p]['edep_MeV'] for p in task['parents']]))
        with np.load(ART/'truth/pair'/f'{i:04d}.npz') as data:
            observed=data['components_mV'].mean(0);solo=data['isolated_mV'].mean(0)
        actual_ratio=100*(observed[weak].sum()/solo[weak].sum()-1)
        with np.load(ART/'predictions'/f'{i:04d}.npz') as pred:
            for method in ('old','corrected'):
                predicted=pred[method+'_joint'];single=pred[method+'_solo']
                ratio=100*(predicted[weak].sum()/single[weak].sum()-1)
                residual=observed-predicted
                details.append(dict(id=i,block=task['block'],condition=task['condition'],gap_ns=task['gap_ns'],
                    phase_ns=task['phase_ns'],method=method,actual_weak_area_change_percent=float(actual_ratio),
                    predicted_weak_area_change_percent=float(ratio),weak_area_error_pp=float(ratio-actual_ratio),
                    weak_mse_mV2=float(np.mean(residual[weak]**2)),strong_mse_mV2=float(np.mean(residual[1-weak]**2)),
                    total_mse_mV2=float(np.mean(residual.sum(0)**2)),model_min_mV=float(predicted.min())))
                if task['block']==0 and task['phase_ns']==0:examples[f'case{i}_{method}_mV']=predicted
        if task['block']==0 and task['phase_ns']==0:
            examples[f'case{i}_actual_mV']=observed;examples[f'case{i}_isolated_mV']=solo
        if task['phase_ns']==0:
            with np.load(ART/'off/truth/pair'/f'{i:04d}.npz') as data:
                observed=data['components_mV'].mean(0);solo=data['isolated_mV'].mean(0)
            off_details.append(dict(id=i,block=task['block'],condition=task['condition'],gap_ns=task['gap_ns'],
                weak_area_change_percent=float(100*(observed[weak].sum()/solo[weak].sum()-1))))
    groups=[];failures=[]
    for condition_id,(condition,gap) in enumerate(CONDITIONS):
        for phase_id,phase in enumerate(PHASES):
            chosen={m:[r for r in details if r['method']==m and r['condition']==condition and r['gap_ns']==gap and r['phase_ns']==phase]
                    for m in ('old','corrected')}
            assert all(len(v)==32 for v in chosen.values())
            assert [r['block'] for r in chosen['old']]==list(range(32))
            index=np.random.default_rng(740900000+condition_id*4+phase_id).integers(32,size=(10000,32))
            summaries={m:{key:float(np.mean([r[key] for r in rows])) for key in ('actual_weak_area_change_percent',
                'predicted_weak_area_change_percent','weak_area_error_pp','weak_mse_mV2','strong_mse_mV2','total_mse_mV2')} for m,rows in chosen.items()}
            ratios={};intervals={}
            for key in ('weak_mse_mV2','strong_mse_mV2','total_mse_mV2'):
                before=np.array([r[key] for r in chosen['old']]);after=np.array([r[key] for r in chosen['corrected']])
                ratios[key]=float(after.mean()/before.mean())
                intervals[key]=np.percentile(after[index].mean(1)/before[index].mean(1),[2.5,97.5]).tolist()
            error=np.array([r['weak_area_error_pp'] for r in chosen['corrected']])
            passed=abs(error.mean())<=5 and max(ratios['weak_mse_mV2'],ratios['strong_mse_mV2'])<=1.1
            if not passed:failures.append(f'{condition}/{gap}/phase{phase}')
            groups.append(dict(condition=condition,gap_ns=gap,phase_ns=phase,mother_pairs=32,
                old=summaries['old'],corrected=summaries['corrected'],mse_ratios=ratios,mse_ratio_ci95=intervals,
                corrected_area_error_ci95=paired_mean_ci(error,index),passed=bool(passed)))
    off_groups=[]
    for j,(condition,gap) in enumerate(CONDITIONS):
        values=np.array([r['weak_area_change_percent'] for r in off_details if r['condition']==condition and r['gap_ns']==gap])
        assert len(values)==32
        index=np.random.default_rng(740910000+j).integers(32,size=(10000,32));passed=abs(values.mean())<=3
        off_groups.append(dict(condition=condition,gap_ns=gap,phase_ns=0.,mother_pairs=32,
            mean_weak_area_change_percent=float(values.mean()),ci95=paired_mean_ci(values,index),passed=bool(passed)))
        if not passed:failures.append(f'recovery off:{condition}/{gap}')
    save(EXP/'pair-details.json',dict(on=details,off=off_details))
    save(EXP/'pair-result.json',dict(groups=groups,recovery_off=off_groups,minimum_model_mV=min(r['model_min_mV'] for r in details)))
    save(EXP/'pair-gate.json',dict(passed=not failures,failures=failures,no_pair_energy_fits=True,default_changed=False,
        next='Freeze a paired known-K/time energy-fit comparison on these unchanged waveforms' if not failures else
             'Stop before pair-energy fitting; retain response mismatch and diagnose on a separate development design'))
    np.savez_compressed(ART/'pair-examples.npz',**examples)
    paths=[EXP/'pair-details.json',EXP/'pair-result.json',EXP/'pair-gate.json',ART/'pair-examples.npz']
    save(EXP/'pair-result-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})
    print(json.dumps(read(EXP/'pair-gate.json'),ensure_ascii=False,indent=2),flush=True)


def replay():
    verify();assert not (EXP/'replay.json').exists()
    singles=read(EXP/'single-plan.json');single_tasks=[dict(singles[i],repeats=2) for i in (32,543)]
    pair_tasks=[]
    if (EXP/'pair-result.json').exists():
        pairs=read(EXP/'pair-plan.json');pair_tasks=[dict(pairs[i],repeats=2) for i in (0,639)]
    tasks=single_tasks+pair_tasks
    with ProcessPoolExecutor(max_workers=4,initializer=initialize) as pool:
        rows=list(pool.map(frozen.simulate,[(t,'replay') for t in tasks]))
    exact=0
    for task in tasks:
        for part in ('voltage','truth'):
            path=f"{part}/{task['kind']}/{task['id']:04d}.npz"
            with np.load(ART/path) as a,np.load(ART/'replay'/path) as b:
                for key in b.files:assert np.array_equal(a[key][:2],b[key]),path;exact+=1
    with ProcessPoolExecutor(max_workers=2,initializer=initialize,initargs=('fit',)) as pool:
        fits=list(pool.map(frozen.fit_one,[(t['id'],t['times_ns'][0],'fit-replay') for t in single_tasks]))
    for task in single_tasks:
        a=read(ART/'fits'/f"{task['id']:04d}.json");b=read(ART/'fit-replay'/f"{task['id']:04d}.json")
        assert [{k:v for k,v in r.items() if k!='seconds'} for r in a['rows']]==[{k:v for k,v in r.items() if k!='seconds'} for r in b['rows']]
    off=[]
    if pair_tasks:
        off_task=dict(read(EXP/'pair-plan.json')[636],repeats=2)
        with ProcessPoolExecutor(max_workers=1,initializer=initialize,initargs=('off',)) as pool:
            off=list(pool.map(frozen.simulate,[(off_task,'off-replay')]))
        for part in ('voltage','truth'):
            path=f'{part}/pair/0636.npz'
            with np.load(ART/'off'/path) as a,np.load(ART/'off-replay'/path) as b:
                for key in b.files:assert np.array_equal(a[key][:2],b[key]),path;exact+=1
    save(EXP/'replay.json',dict(passed=True,exact_arrays=exact,digitizations=sum(r['digitizations'] for r in rows+off),
        single_energy_fits=sum(r['fits'] for r in fits),digital_rows=rows,off_rows=off,fit_rows=fits,
        outputs_equal_except_fit_timing=True))
    print('Replay passed',exact,'arrays',sum(r['fits'] for r in fits),'single fits',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['prepare','verify','preflight','native_pilot','pilot','native_formal',
        'single','fit','single_score','pair_pilot','predict','pair','off','pair_score','replay'])
    stage=parser.parse_args().stage
    if stage=='prepare':prepare()
    elif stage=='verify':print(json.dumps(verify()))
    elif stage in ('preflight','native_pilot','native_formal'):run_native(stage)
    elif stage=='pilot':pilot()
    elif stage=='single':singles()
    elif stage=='fit':fit_singles()
    elif stage=='single_score':score_singles()
    elif stage=='pair_pilot':pair_pilot()
    elif stage=='predict':predictions()
    elif stage in ('pair','off'):pairs('on' if stage=='pair' else 'off')
    elif stage=='pair_score':score_pairs()
    elif stage=='replay':replay()


if __name__=='__main__':main()
