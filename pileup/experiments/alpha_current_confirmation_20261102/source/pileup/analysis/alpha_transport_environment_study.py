"""Frozen same-seed migration check before releasing new optical mothers."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import threading
import time

import awkward as ak
import numpy as np
import uproot
from uproot.source.file import MemmapSource

from .alpha_transport_environment import relocate_macro,macro_semantics,compare_photons,event_statistics

NAME='alpha_transport_environment_20261102'
ROOT=Path(os.environ.get('GAGG_TRANSPORT_ENV_ROOT',Path(__file__).resolve().parents[2])).resolve()
EXP=ROOT/'pileup/experiments'/NAME
ART=ROOT/'jupyter/artifacts'/NAME
REC=ROOT/'jupyter/run_records'/NAME
PARENT=ROOT/'pileup/experiments/alpha_mean_shape_20261102'


def read(p):return json.loads(Path(p).read_text())


def save(p,value):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.partial')
    tmp.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n');tmp.replace(p)


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        while block:=f.read(1024*1024):h.update(block)
    return h.hexdigest()


def load_root(path):
    with uproot.open(path,handler=MemmapSource,num_workers=1,decompression_executor=uproot.TrivialExecutor(),
                     interpretation_executor=uproot.TrivialExecutor()) as f:
        assert all(f[k].num_entries==1 for k in ('RunInfo','WaveformEvent','PrimaryPhoton','SiPMEvent'))
        info=ak.to_list(f['RunInfo'].arrays(library='ak'))[0]
        primary=ak.to_list(f['PrimaryPhoton'].arrays(['primary_particle_name','primary_energy_MeV','primary_time_ps',
            'edep_total_MeV','n_scint_photons_generated'],library='ak'))[0]
        waveform=f['WaveformEvent'].arrays(library='ak');assert int(waveform['n_primary'][0])==1
        photons={new:ak.to_numpy(waveform[old][0]) for old,new in
                 [('all_photon_arrival_t_ps','time_ps'),('all_photon_arrival_x_mm','x_mm'),
                  ('all_photon_arrival_y_mm','y_mm'),('photon_energy_eV','energy_eV')]}
        assert len(photons['time_ps'])==int(waveform['n_optical_photons_arrived_total'][0])
        sipm=f['SiPMEvent'].arrays(['sample_time_ns','signal_mV','analog_mV','charge_total_C','charge_balance_error_C',
            'n_photon_avalanches','n_dark_avalanches','n_crosstalk_avalanches'],library='ak')
        voltage={k:ak.to_numpy(sipm[k][0]) for k in ('sample_time_ns','signal_mV','analog_mV')}
        scalars={k:ak.to_list(sipm[k])[0] for k in ('charge_total_C','charge_balance_error_C',
            'n_photon_avalanches','n_dark_avalanches','n_crosstalk_avalanches')}
    assert primary['primary_particle_name']=='alpha' and info['n_events']==1
    return info,primary,photons,voltage,scalars


def prepare():
    assert EXP.exists() and not (EXP/'design-freeze.json').exists()
    for part in ('macros','runinfo','reference-runinfo','reference-metadata','bin','analysis-source/pileup/analysis'):
        (EXP/part).mkdir(parents=True,exist_ok=True)
    for part in ('logs','attempts'):(REC/part).mkdir(exist_ok=True)
    for part in ('comparisons','figures'):(ART/part).mkdir(exist_ok=True)
    for p,d in read(EXP/'source-freeze.json').items():assert sha(ROOT/p)==d,p
    shutil.copyfile(EXP/'build/gagg',EXP/'bin/gagg');(EXP/'bin/gagg').chmod(0o755)
    tests=subprocess.run(['ctest','--test-dir',str(EXP/'build'),'--output-on-failure'],capture_output=True,text=True)
    (REC/'ctest.log').write_text(tests.stdout+tests.stderr);assert tests.returncode==0
    save(EXP/'build.json',dict(binary_sha256=sha(EXP/'bin/gagg'),ctest_passed=True,
        configure=['cmake','-S',str((EXP/'source/g4').relative_to(ROOT)),'-B',str((EXP/'build').relative_to(ROOT)),
            '-DCMAKE_BUILD_TYPE=Release','-DGeant4_DIR=/opt/geant4/geant4-v1132/lib/cmake/Geant4','-DROOT_DIR=/opt/root/root63410/cmake'],
        build=['cmake','--build',str((EXP/'build').relative_to(ROOT)),'-j','4'],
        compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0]))
    candidates=[r for r in read(PARENT/'mothers.json') if r['role']=='train'];plan=[]
    for group in range(4):
        rows=[r for r in candidates if r['group']==group]
        rows.sort(key=lambda r:hashlib.sha256(('738000:'+r['root']).encode()).hexdigest())
        for repeat,r in enumerate(rows[:4]):
            assert sha(ROOT/r['root'])==r['root_original_sha256']
            reference=r['id'];i=len(plan)
            shutil.copyfile(PARENT/'original-runinfo'/f'{reference:04d}.json',EXP/'reference-runinfo'/f'{i:02d}.json')
            shutil.copyfile(PARENT/'original-metadata'/f'{reference:04d}.json',EXP/'reference-metadata'/f'{i:02d}.json')
            info=read(EXP/'reference-runinfo'/f'{i:02d}.json')['run_info']
            assert info['random_seed']==r['seed'] and info['sipm_enabled']
            plan.append(dict(id=i,group=group,condition=r['condition'],source_energy_MeV=r['source_energy_MeV'],
                replicate=repeat,seed=r['seed'],reference_root=r['root'],reference_sha256=r['root_original_sha256'],
                reference_mother_id=reference,role='same_seed_environment_control',independent_new_mother=False))
    save(EXP/'plan.json',plan)
    save(EXP/'protocol.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),
        base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        git_status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),
        purpose='Resolve software-environment transfer before independent strong-mother confirmation',
        data='16 exact old seeds,4per source energy1.625/2.5/5.486/8MeV;deterministic path-hash selection from128L93training mothers',
        old_environment='Geant4 11.4.2,ROOT6.40.04,Fedora;old executable cannot load on current libc/Qt',
        current_environment='Geant4 11.3.2,ROOT6.34.10,Ubuntu/WSL;all39physics-source files exactly preserved',
        changes='Only dependency/build environment and macro output/config paths;no physical parameter or algorithm adjustment',
        windows='One20us WaveformEvent and SiPMEvent per ROOT;one central normal alpha at5us;actual original macro retained',
        primary_gate='Every seed: initial RNGstate equal,E_dep and generated scintillation count exactly equal,and sorted full photon records (integerps/float32bitpatterns) exactly equal. A conservative same-seed reproducibility check,not a statistical equivalence test',
        secondary='Ordered photon equality;native SiPM voltage/charge comparison. Native electronics mismatch alone does not identify optical transport cause',
        stop='Do not release336fresh strong mothers if primary transport gate fails;dependency changes cannot be mistaken for an independent same-environment confirmation',
        inferential_scope='16paired engineering controls,not16independent new optical mothers;different environments may decorrelate RNGpaths without population bias. No detector-efficiency or energy-reconstruction claim',
        preflight_id=12,cost_limit_seconds=900.,max_workers=4,launch_spacing_seconds=1.05,
        output='Only g4/data/YYYYMMDD_HHhMMmSSs.root;no overwrite;logs and analysis remain under jupyter',
        no_new_offline_digitization=True,mean_model_frozen=True,default_changed=False))
    future=dict(status='frozen_design_not_released',purpose='Independent strong-only confirmation of the unchangedL93mean correction',
        pilot_mothers=8,formal_mothers=336,source_energy_MeV=8.,native_seed_base=738000000,
        native_SiPM_seed_offset=100000,native_readout_seed_offset=200000,offline_seed_base=738300000,
        offline_seed_formula='738300000+motherID*100+repeatID',offline_readout_seed_offset=10000000,
        repeats=8,phase_ns=0.,known_K=1,known_time='First active-plane photon aligned to128ns',
        primary_endpoints='Fixed old8/40ns weak-logE projected per-mother MSE ratios<1 with paired mother-bootstrap95%upper<1;absolute corrected mean<=max(10,.7absold)',
        secondary_controls='Known-onset isolated E mean bias<=5%,individual RMSEratio<=1.1,zero unconverged fits;200ns projection secondary',
        source_scope='New central source8MeV only;truthE_dep remains separate;no Q/E cuts or replacement seeds based on results',
        model_paths_and_hashes=read(PARENT/'model-freeze.json'),
        hold='Transport environment must be reconciled first. No fresh ROOT from this design generated in the migration experiment',
        limitations='336rounds theL94approximation334;not guaranteed precision. Fresh strong-only pass would still require other energies/phases/shared controls before pair-energy work')
    streams={738000000+i+offset for i in range(344) for offset in (0,100000,200000)}
    audit=[]
    for p in sorted((ROOT/'pileup/experiments').glob('*/manifest.json')):
        m=read(p)
        if not isinstance(m,dict):continue
        old={int(r['seed'])+off for key in ('runs','planned_runs') for r in m.get(key,[]) if isinstance(r,dict)
             and isinstance(r.get('seed'),int) for off in (0,100000,200000)}
        assert not streams&old,p;audit.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
    save(EXP/'confirmation-design.json',future);save(EXP/'future-seed-audit.json',dict(passed=True,manifests=audit,streams=len(streams),new_mothers_run=0))
    source=EXP/'analysis-source/pileup/analysis'
    (source.parent/'__init__.py').write_text('');(source/'__init__.py').write_text('')
    for name in ('alpha_transport_environment.py','alpha_transport_environment_study.py','test_alpha_transport_environment.py'):
        shutil.copyfile(ROOT/'pileup/analysis'/name,source/name)
    (EXP/'run.py').write_text("from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\n"
        "os.environ['GAGG_TRANSPORT_ENV_ROOT']=str(f.parents[2])\nos.environ['OPENBLAS_NUM_THREADS']='1'\nos.environ['OMP_NUM_THREADS']='1'\n"
        "sys.path.insert(0,str(f/'analysis-source'))\nif __name__=='__main__':\n    from pileup.analysis.alpha_transport_environment_study import main\n    main()\n")
    paths=[p for p in EXP.rglob('*') if p.is_file() and not {'build','bin'}&set(p.relative_to(EXP).parts)]
    save(EXP/'design-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})
    save(EXP/'manifest.json',dict(planned_runs=plan,runs=[],failed_attempts=[],new_independent_mothers=0))
    print('Frozen16same-seed controls;336fresh-mother design held pending environment checks',flush=True)


def verify():
    paths=read(EXP/'design-freeze.json')
    for p,d in paths.items():assert sha(ROOT/p)==d,p
    assert sha(EXP/'bin/gagg')==read(EXP/'build.json')['binary_sha256']
    return dict(files=len(paths),binary_matches=True)


def run(stage):
    verify();manifest=read(EXP/'manifest.json');done={r['id'] for r in manifest['runs']}
    if stage=='complete':assert read(EXP/'preflight.json')['passed']
    plan=[r for r in manifest['planned_runs'] if r['id'] not in done and (stage=='complete' or r['id']==12)]
    assert plan,'Requested stage already complete;use saved outputs'
    env=os.environ.copy();env.update(read(EXP/'environment.json')['current_environment']['g4_dataset_environment'])
    lock=threading.Lock();last=[max([r['started_epoch'] for r in manifest['runs']]+[0.])]
    def one(task):
        with lock:
            while True:
                delay=1.05-(time.time()-last[0])
                if delay>0:time.sleep(delay)
                output=ROOT/'g4/data'/(datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss')+'.root')
                if not output.exists():break
                time.sleep(1.05)
            suffix=output.stem;macro=EXP/'macros'/f"control{task['id']:02d}_{suffix}.mac"
            old=read(EXP/'reference-runinfo'/f"{task['id']:02d}.json")['run_info']
            body=relocate_macro(old['macro_contents'],output,EXP/'source/g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json')
            assert macro_semantics(body)==macro_semantics(old['macro_contents']);macro.write_text(body)
            log=REC/'logs'/f'{suffix}.log';stream=log.open('w');start=time.monotonic()
            row=dict(task,root=str(output.relative_to(ROOT)),macro=str(macro.relative_to(ROOT)),macro_sha256=sha(macro),
                log=str(log.relative_to(ROOT)),started_epoch=time.time(),started_utc=datetime.now(timezone.utc).isoformat())
            proc=subprocess.Popen([str(EXP/'bin/gagg'),str(macro)],cwd=ROOT/'g4/build',env=env,
                                  stdout=stream,stderr=subprocess.STDOUT)
            last[0]=time.time();row['started_epoch']=last[0]
        code=proc.wait();stream.close();row.update(returncode=code,seconds=time.monotonic()-start)
        if code:
            row['status']='failed';save(REC/'attempts'/f'{suffix}.json',row)
            raise RuntimeError(f"Native control exited{code};preserved {row['log']}")
        info,primary,photons,voltage,scalars=load_root(output)
        assert info['random_seed']==task['seed'] and info['macro_contents']==body
        assert info['t_length_ps']==old['t_length_ps']==20000000
        assert np.isclose(primary['primary_energy_MeV'],task['source_energy_MeV'],rtol=2e-7,atol=1e-9)
        assert json.loads(info['sipm_config_json'])==json.loads(old['sipm_config_json'])
        for key in ('gagg_max_step_um','mylar_front_max_step_um','mylar_side_max_step_um',
                    'sipm_avalanche_base_seed','sipm_readout_base_seed'):
            assert info[key]==old[key],key
        assert abs(scalars['charge_balance_error_C'])<1e-18
        save(EXP/'runinfo'/f"{task['id']:02d}.json",info)
        row.update(status='completed',root_sha256=sha(output),root_bytes=output.stat().st_size,
            structure_macro_config_verified=True,edep_MeV=primary['edep_total_MeV'])
        save(REC/'attempts'/f'{suffix}.json',row);return row
    failures=[]
    with ThreadPoolExecutor(max_workers=1 if stage=='preflight' else 4) as pool:
        futures={pool.submit(one,t):t for t in plan}
        for f in as_completed(futures):
            try:r=f.result()
            except Exception as exc:
                error=dict(task=futures[f],error=repr(exc));failures.append(error);manifest['failed_attempts'].append(error)
            else:
                manifest['runs'].append(r);print(stage,len(manifest['runs']),'/16',r['id'],round(r['seconds'],2),'s',flush=True)
            save(EXP/'manifest.json',manifest)
    if failures:raise RuntimeError(failures)
    if stage=='preflight':
        r=next(x for x in manifest['runs'] if x['id']==12);estimate=16*max(1.05,r['seconds']/4)*1.5
        save(EXP/'preflight.json',dict(passed=estimate<900,projected_seconds=estimate,structure_verified=True,
             physics_matching_not_used_for_cost_gate=True,root=r['root']))
        assert estimate<900


def compare():
    verify();manifest=read(EXP/'manifest.json');assert len(manifest['runs'])==16
    rows=[];hist_edges=np.arange(0.,1604.,4.);examples={}
    for task in sorted(manifest['runs'],key=lambda r:r['id']):
        assert sha(ROOT/task['root'])==task['root_sha256']
        assert sha(ROOT/task['reference_root'])==task['reference_sha256']
        a=load_root(ROOT/task['reference_root']);b=load_root(ROOT/task['root'])
        ai,ap,aph,av,asc=a;bi,bp,bph,bv,bsc=b
        light=compare_photons(aph,bph)
        old=event_statistics(aph,ap['edep_total_MeV']);new=event_statistics(bph,bp['edep_total_MeV'])
        energy_equal=ap['edep_total_MeV']==bp['edep_total_MeV']
        generated_equal=ap['n_scint_photons_generated']==bp['n_scint_photons_generated']
        state_equal=ai['random_engine_state']==bi['random_engine_state']
        time_equal=np.array_equal(av['sample_time_ns'],bv['sample_time_ns'])
        assert time_equal
        row=dict(id=task['id'],group=task['group'],source_energy_MeV=task['source_energy_MeV'],seed=task['seed'],
            old=old,current=new,photons=light,initial_random_state_equal=state_equal,edep_exact_equal=energy_equal,
            generated_exact_equal=generated_equal,old_generated=ap['n_scint_photons_generated'],current_generated=bp['n_scint_photons_generated'],
            transport_gate_passed=bool(state_equal and energy_equal and generated_equal and light['sorted_exact_equal']),
            native_sample_grid_equal=time_equal,native_signal_equal=bool(np.array_equal(av['signal_mV'],bv['signal_mV'])),
            native_analog_equal=bool(np.array_equal(av['analog_mV'],bv['analog_mV'])),
            native_signal_rmse_mV=float(np.sqrt(np.mean((av['signal_mV']-bv['signal_mV'])**2))),
            old_charge_C=asc['charge_total_C'],current_charge_C=bsc['charge_total_C'])
        rows.append(row)
        for label,ph,v in [('old',aph,av),('current',bph,bv)]:
            t=(ph['time_ps'].astype(float)-5000000)/1000
            examples[f"case{task['id']}_{label}_photon_density"]=np.histogram(t,hist_edges)[0]/len(t)/4
            examples[f"case{task['id']}_{label}_signal_mV"]=v['signal_mV']
        examples['native_time_ns']=av['sample_time_ns'];examples['photon_edges_ns']=hist_edges
    np.savez_compressed(ART/'comparisons/events.npz',**examples)
    passed=all(r['transport_gate_passed'] for r in rows)
    summary=dict(controls=16,independent_new_mothers=0,transport_matching_passed=passed,
        exact_edep=sum(r['edep_exact_equal'] for r in rows),exact_generated=sum(r['generated_exact_equal'] for r in rows),
        exact_sorted_photons=sum(r['photons']['sorted_exact_equal'] for r in rows),
        exact_initial_rng=sum(r['initial_random_state_equal'] for r in rows),
        exact_native_signal=sum(r['native_signal_equal'] for r in rows),
        fresh_confirmation_released=False,
        next='Same-seed controls passed;review dependency provenance before a separate fresh campaign' if passed else
             'Do not run336fresh-mother confirmation under mismatched dependencies;restore matching environment or separately rebaseline current version',rows=rows)
    save(EXP/'result-summary.json',summary)
    save(EXP/'output-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in
        [EXP/'manifest.json',EXP/'result-summary.json',ART/'comparisons/events.npz',*sorted((EXP/'macros').glob('*.mac')),
         *sorted((EXP/'runinfo').glob('*.json'))]})
    print(json.dumps({k:v for k,v in summary.items() if k!='rows'},indent=2),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','preflight','complete','compare','verify']);args=p.parse_args()
    if args.stage=='prepare':prepare()
    elif args.stage in ('preflight','complete'):run(args.stage)
    elif args.stage=='compare':compare()
    else:print(json.dumps(verify()))


if __name__=='__main__':main()
