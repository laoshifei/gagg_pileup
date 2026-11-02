"""New optical mothers in the current environment, with unchanged L93 models."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import threading
import time

import numpy as np

from . import alpha_mean_shape_study as frozen
from .alpha_current_response import energy_metrics
from .alpha_transport_environment import relocate_macro
from .alpha_transport_environment_study import load_root
from .sipm_response import SensorConfig, read_photon_hits, map_cells

NAME = 'alpha_current_confirmation_20261102'
ROOT = Path(os.environ.get('GAGG_CURRENT_CONFIRMATION_ROOT', Path(__file__).resolve().parents[2])).resolve()
EXP = ROOT / 'pileup/experiments' / NAME
ART = ROOT / 'jupyter/artifacts' / NAME
REC = ROOT / 'jupyter/run_records' / NAME
CONTROL = ROOT / 'pileup/experiments/alpha_current_response_20261102'
TRANSPORT = ROOT / 'pileup/experiments/alpha_transport_environment_20261102'
read, save, sha = frozen.read, frozen.save, frozen.sha


def initialize(fitting=False):
    frozen.ROOT, frozen.EXP, frozen.ART, frozen.REC = ROOT, EXP, ART, REC
    if fitting:
        frozen.fit_initialize()
    else:
        frozen.initialize(True)


def prepare():
    assert read(CONTROL / 'gate.json')['passed']
    for p in (EXP, ART, REC): p.mkdir(parents=True, exist_ok=False)
    for part in ('bin', 'build', 'macros', 'original-runinfo'): (EXP / part).mkdir()
    for part in ('photons', 'voltage/single', 'truth/single', 'fits', 'figures'): (ART / part).mkdir(parents=True)
    for part in ('logs', 'attempts'): (REC / part).mkdir()
    shutil.copytree(CONTROL / 'source', EXP / 'source', ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copyfile(ROOT / 'pileup/analysis/alpha_current_confirmation_study.py',
                    EXP / 'source/pileup/analysis/alpha_current_confirmation_study.py')
    models = ('template.npz', 'energy-model.json', 'response-parameters.json', 'mean-correction.npz',
              'noise-calibration.json', 'sensor-config.json', 'charge-map-spatial12.npz', 'spatial-profile.npz',
              'shape-correction.npz', 'nominal-influences.npz')
    for name in models: shutil.copyfile(CONTROL / name, EXP / name)
    shutil.copyfile(CONTROL / 'bin/age.so', EXP / 'bin/age.so')
    shutil.copyfile(ROOT / 'g4/build/gagg', EXP / 'bin/gagg'); (EXP / 'bin/gagg').chmod(0o755)
    native = read(TRANSPORT / 'main-build.json')
    for path, digest in native['source_sha256'].items():
        assert sha(ROOT / path) == digest
        assert sha(EXP / 'source' / path) == digest
    assert sha(EXP / 'bin/gagg') == native['binary_sha256']
    save(EXP / 'build.json', dict(native_binary_sha256=sha(EXP / 'bin/gagg'),
        age_binary_sha256=sha(EXP / 'bin/age.so'), native_build=native,
        age_rebuild=read(CONTROL / 'build.json'), binaries_reused_from_current_validated_build=True))
    save(EXP / 'model-freeze.json', {str((EXP / name).relative_to(ROOT)): sha(EXP / name) for name in models})
    save(EXP / 'control-gate.json', read(CONTROL / 'gate.json'))
    save(EXP / 'control-result.json', read(CONTROL / 'result-summary.json'))
    shutil.copyfile(TRANSPORT / 'confirmation-design.json', EXP / 'prior-design.json')
    shutil.copyfile(TRANSPORT / 'environment.json', EXP / 'transport-environment.json')
    shutil.copyfile(CONTROL / 'analysis-environment.json', EXP / 'analysis-environment.json')
    template = read(TRANSPORT / 'runinfo/12.json')['macro_contents']
    (EXP / 'macro-template.mac').write_text(template)
    plan = [dict(id=i, role='pilot' if i < 8 else 'formal', seed=738000000+i,
                 sipm_seed=738100000+i, readout_seed=738200000+i, source_energy_MeV=8.) for i in range(344)]
    streams = {r[key] for r in plan for key in ('seed', 'sipm_seed', 'readout_seed')}
    audited = []
    for path in sorted((ROOT / 'pileup/experiments').glob('*/manifest.json')):
        data = read(path)
        if not isinstance(data, dict): continue
        historical = {r['seed']+offset for key in ('runs', 'planned_runs') for r in data.get(key, [])
                      if isinstance(r, dict) and isinstance(r.get('seed'), int) for offset in (0,100000,200000)}
        assert not streams & historical, path
        audited.append(dict(path=str(path.relative_to(ROOT)), sha256=sha(path)))
    save(EXP / 'seed-audit.json', dict(passed=True, native_streams=len(streams), manifests=audited,
        scope='Historical manifest native seed fields and established +100000/+200000 streams',
        offline='738300000+motherID*100+repeat; readout+10000000; distinct from L93 736... and controls739...'))
    save(EXP / 'native-plan.json', plan)
    tasks = [dict(id=r['id'], kind='single', parents=[r['id']], times_ns=[128.], phase_ns=0.,
                  seed_base=738300000+r['id']*100, repeats=8, role=r['role']) for r in plan]
    save(EXP / 'single-plan.json', tasks)
    save(EXP / 'protocol.json', dict(created_utc=datetime.now(timezone.utc).isoformat(),
        base_commit=subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip(),
        git_status=subprocess.check_output(['git','status','--short'], cwd=ROOT, text=True),
        purpose='Independent current-environment transfer test of frozen L93 old and corrected means; not same-environment historical confirmation',
        prerequisite='Current-response four-energy/four-phase engineering and native readout controls passed before fresh generation',
        sample='8 pilot mothers IDs0..7 excluded from inference;336 formal mothers IDs8..343;source8MeV central normal alpha; no Q/E cuts',
        independence='Native G4 seeds are unused in audited historical manifests; no fresh event can train or change fixed model/thresholds',
        physics='Unchanged39G4 source/config files;Geant4 11.3.2 ROOT6.34.10;95/351ns fastfraction.3,29ns recovery proxy,14400cells',
        window='One20us native WaveformEvent/SiPMEvent perROOT;onealpha at5us; no Poisson source rate or background',
        output='g4/data/YYYYMMDD_HHhMMmSSs.root; no overwrite; launches>=1.05s; serial G4, at most8 concurrent tasks',
        formal_digitizations=2688, pilot_digitizations=64, formal_fits=6048, pilot_fits=36,
        observation='4096ns offline, first2048ns/1024samples;2nsperiod,.5nsaperture;first active photon aligned128ns, phase0;all photons retained',
        primary='8and40ns fixed nominal weak-logE projected per-mother MSEratio<1 AND pairedmother bootstrap95%upper<1 AND abs(correctedmean)<=max(10,.7absoldmean)',
        secondary='200ns projection;single-E corrected mean bias<=5%,individual RMSEratio<=1.1,zero unconverged fits; report boundaries and conditional±10%counts',
        interpretation='Projection units100delta ln E are not actual pair energy errors. Known K1/t isolated fits are not blind efficiency; no double-pulse fitting or high-rate expansion',
        bootstrap='10000 paired whole-mother draws,seed738900000; conditional on frozen models and this current environment; no model retraining uncertainty',
        pilot_gate='Only structure,config/seeds,finitearrays,causality,readoutaudit and projected cost<900s; no efficacy threshold using pilot outcomes',
        stages=['preflight nativeID0','native pilot IDs1..7','pilot digital/readout/fit checks','native formal','formal offline','formal fit','score','local replay'],
        pilot_fit_ids=[0,1], workers=8, BLAS_threads=1, no_new_model=True, default_changed=False,
        limitations='Strong-source only;336 rounds old approx334 and does not guarantee precision; environment shift and finite training uncertainty remain'))
    (EXP / 'run.py').write_text("from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\n"
        "os.environ['GAGG_CURRENT_CONFIRMATION_ROOT']=str(f.parents[2])\nos.environ['OPENBLAS_NUM_THREADS']='1'\nos.environ['OMP_NUM_THREADS']='1'\n"
        "sys.path.insert(0,str(f/'source'))\nif __name__=='__main__':\n    from pileup.analysis.alpha_current_confirmation_study import main\n    main()\n")
    paths = [p for p in EXP.rglob('*') if p.is_file() and not {'bin','build'} & set(p.relative_to(EXP).parts)]
    save(EXP / 'design-freeze.json', {str(p.relative_to(ROOT)): sha(p) for p in paths})
    save(EXP / 'manifest.json', dict(planned_runs=plan, runs=[], failed_attempts=[]))
    print('Frozen fresh current-environment campaign:8pilot+336formal mothers', flush=True)


def verify():
    counts = {}
    for name in ('design-freeze.json', 'model-freeze.json'):
        data = read(EXP / name)
        for path, digest in data.items(): assert sha(ROOT / path) == digest, path
        counts[name] = len(data)
    build = read(EXP / 'build.json')
    assert sha(EXP / 'bin/gagg') == build['native_binary_sha256']
    assert sha(EXP / 'bin/age.so') == build['age_binary_sha256']
    return counts


def run_native(stage):
    verify(); manifest = read(EXP / 'manifest.json')
    if stage == 'native_pilot': assert read(EXP / 'native-cost-gate.json')['passed']
    if stage == 'native_formal': assert read(EXP / 'pilot-gate.json')['passed']
    done = {r['id'] for r in manifest['runs']}
    tasks = [r for r in manifest['planned_runs'] if r['id'] not in done and
             (r['id'] == 0 if stage == 'preflight' else r['role'] == ('pilot' if stage == 'native_pilot' else 'formal'))]
    assert tasks, 'Stage already complete; preserve saved results'
    env = os.environ.copy(); env.update(read(EXP / 'transport-environment.json')['current_environment']['g4_dataset_environment'])
    sensor = SensorConfig(**read(EXP / 'sensor-config.json')['sensor'])
    last = [max([r['started_epoch'] for r in manifest['runs']] + [0.])]; lock = threading.Lock()
    def one(task):
        with lock:
            while True:
                delay = 1.05-(time.time()-last[0])
                if delay > 0: time.sleep(delay)
                root = ROOT / 'g4/data' / (datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss')+'.root')
                if not root.exists(): break
                time.sleep(1.05)
            macro = EXP / 'macros' / f"mother{task['id']:04d}_{root.stem}.mac"
            body = relocate_macro((EXP / 'macro-template.mac').read_text(), root,
                                  EXP / 'source/g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json')
            replacements = {'/gagg/run/seed':task['seed'], '/gagg/sipm/seed':task['sipm_seed'],
                            '/gagg/sipm/readoutSeed':task['readout_seed']}
            found = dict.fromkeys(replacements, 0); lines=[]
            for line in body.splitlines():
                command = line.strip().split(maxsplit=1)[0] if line.strip() else ''
                if command in replacements:
                    found[command] += 1; line=command+' '+str(replacements[command])
                lines.append(line)
            assert set(found.values()) == {1}
            body='\n'.join(lines)+'\n'; macro.write_text(body)
            log=REC/'logs'/f'{root.stem}.log'; stream=log.open('w')
            started=time.monotonic()
            proc=subprocess.Popen([str(EXP/'bin/gagg'),str(macro)],cwd=EXP/'build',env=env,
                                  stdout=stream,stderr=subprocess.STDOUT)
            last[0]=time.time()
            row=dict(task,root=str(root.relative_to(ROOT)),macro=str(macro.relative_to(ROOT)),macro_sha256=sha(macro),
                     log=str(log.relative_to(ROOT)),started_epoch=last[0],started_utc=datetime.now(timezone.utc).isoformat())
        code=proc.wait();stream.close()
        row.update(returncode=code,seconds=time.monotonic()-started)
        save(REC/'attempts'/f'{root.stem}.json',row)
        assert code == 0, row
        info,primary,photons,voltage,scalars=load_root(root)
        assert info['random_seed']==task['seed'] and info['sipm_avalanche_base_seed']==task['sipm_seed']
        assert info['sipm_readout_base_seed']==task['readout_seed'] and info['macro_contents']==body
        assert info['t_length_ps']==20000000 and primary['primary_energy_MeV']==8.
        assert primary['primary_time_ps']==5000000 and 0<primary['edep_total_MeV']<=8.
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
        save(REC/'attempts'/f'{root.stem}.json',row)
        return row
    failures=[]
    with ThreadPoolExecutor(max_workers=1 if stage=='preflight' else 8) as pool:
        futures={pool.submit(one,task):task for task in tasks}
        for future in as_completed(futures):
            try: row=future.result()
            except Exception as exc:
                error=dict(task=futures[future],error=repr(exc));manifest['failed_attempts'].append(error);failures.append(error)
            else:
                manifest['runs'].append(row)
                if len(manifest['runs'])<=8 or len(manifest['runs'])%16==0:
                    print(stage,len(manifest['runs']),'/344',round(row['seconds'],2),'s',flush=True)
            save(EXP/'manifest.json',manifest)
    save(EXP/'mothers.json',sorted(manifest['runs'],key=lambda r:r['id']))
    assert not failures, failures
    if stage=='preflight':
        cost=manifest['runs'][0]['seconds'];estimate=344*max(1.05,cost/8)*1.5
        save(EXP/'native-cost-gate.json',dict(passed=estimate<900,one_seconds=cost,projected_seconds=estimate))
        assert estimate<900
    if stage=='native_formal':
        assert len(manifest['runs'])==344
        paths=[ROOT/r['root'] for r in manifest['runs']]+list((ART/'photons').glob('*.npz'))
        paths += [EXP/'manifest.json',EXP/'mothers.json',*sorted((EXP/'macros').glob('*.mac')),
                  *sorted((EXP/'original-runinfo').glob('*.json'))]
        save(EXP/'native-output-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})


def pilot():
    verify();mothers=read(EXP/'mothers.json');assert len(mothers)==8
    assert not (EXP/'pilot-gate.json').exists()
    from analysis.single_particle_sipm import audit_response,load_root as native_load
    audits=[dict(id=r['id'],**audit_response(native_load(ROOT/r['root']))) for r in mothers]
    tasks=read(EXP/'single-plan.json')[:8]
    with ProcessPoolExecutor(max_workers=8,initializer=initialize) as pool:
        digits=list(pool.map(frozen.simulate,[(task,'formal') for task in tasks]))
    with ProcessPoolExecutor(max_workers=2,initializer=initialize,initargs=(True,)) as pool:
        fits=list(pool.map(frozen.fit_one,[(i,128.,'pilot-fits') for i in (0,1)]))
    estimates=dict(digitization_seconds=max(r['seconds'] for r in digits)*336/8*1.5,
                   fit_seconds=max(r['seconds'] for r in fits)*336/8*1.5)
    assert all(np.isfinite(list(estimates.values())))
    passed=all(v<900 for v in estimates.values())
    save(EXP/'pilot-gate.json',dict(passed=passed,readout_audits=audits,digitization_rows=digits,fit_rows=fits,
         formal_cost_estimates=estimates,pilot_mothers=8,digitizations=64,fits=36,
         efficacy_outcomes_not_used_for_gate=True))
    print(json.dumps(dict(passed=passed,**estimates)),flush=True);assert passed


def generate():
    verify();assert len(read(EXP/'manifest.json')['runs'])==344
    assert not (EXP/'formal-generation.json').exists()
    tasks=read(EXP/'single-plan.json')[8:];started=time.monotonic()
    with ProcessPoolExecutor(max_workers=8,initializer=initialize) as pool:
        rows=[]
        for row in pool.map(frozen.simulate,[(task,'formal') for task in tasks]):
            rows.append(row)
            if len(rows)%32==0:print('formal digits',len(rows),'/336',flush=True)
    save(EXP/'formal-generation.json',dict(rows=rows,digitizations=2688,seconds=time.monotonic()-started))
    paths=[*sorted((ART/'voltage/single').glob('*.npz')),*sorted((ART/'truth/single').glob('*.npz'))]
    save(EXP/'digitization-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})


def fit():
    verify()
    for path,digest in read(EXP/'digitization-freeze.json').items():assert sha(ROOT/path)==digest
    assert not list((ART/'fits').glob('*.json'))
    started=time.monotonic()
    with ProcessPoolExecutor(max_workers=8,initializer=initialize,initargs=(True,)) as pool:
        rows=[]
        for row in pool.map(frozen.fit_one,[(i,128.,'fits') for i in range(8,344)]):
            rows.append(row)
            if len(rows)%32==0:print('formal fits',len(rows),'/336',flush=True)
    save(EXP/'fits.json',dict(rows=rows,fits=6048,seconds=time.monotonic()-started))
    save(EXP/'fit-output-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in (ART/'fits').glob('*.json')})


def score():
    verify();initialize()
    for name in ('native-output-freeze.json','digitization-freeze.json','fit-output-freeze.json'):
        for path,digest in read(EXP/name).items():assert sha(ROOT/path)==digest
    mothers=read(EXP/'mothers.json')[8:];assert len(mothers)==336 and all(r['role']=='formal' for r in mothers)
    details=[];residuals={'old':[],'corrected':[]};observed=[]
    for mother in mothers:
        i=mother['id'];fit=read(ART/'fits'/f'{i:04d}.json')
        assert fit['truth_guard'] and not fit['true_E_input'] and fit['known_K1_t']
        for row in fit['rows']:
            details.append(dict(id=i,method=row['method'],repeat=row['repeat'],
                relative_error_percent=100*(row['energy_MeV'][0]/mother['edep_MeV']-1),
                boundary=row['boundary'],converged=row['converged']))
        with np.load(ART/'truth/single'/f'{i:04d}.npz') as data:actual=data['components_mV'].mean(0)[0]
        observed.append(actual)
        for method,model in [('old',frozen.OLD),('corrected',frozen.NEW)]:
            residuals[method].append(actual-model.isolated([128.],[mother['edep_MeV']])[0][0])
    with np.load(EXP/'nominal-influences.npz') as data:operator=data['high'][1:]
    before=100*np.asarray(residuals['old'])@operator.T;after=100*np.asarray(residuals['corrected'])@operator.T
    index=np.random.default_rng(738900000).integers(336,size=(10000,336))
    projections=[];failures=[]
    for j,gap in enumerate((8.,40.,200.)):
        b,a=before[:,j],after[:,j];mse_b=float(np.mean(b*b));mse_a=float(np.mean(a*a))
        bootstrap_ratio=np.mean(a[index]**2,axis=1)/np.mean(b[index]**2,axis=1)
        ci=np.percentile(bootstrap_ratio,[2.5,97.5]);limit=max(10.,.7*abs(b.mean()))
        passed=mse_a/mse_b<1 and ci[1]<1 and abs(a.mean())<=limit
        row=dict(gap_ns=gap,mothers=336,old_mean=float(b.mean()),corrected_mean=float(a.mean()),
            old_mse=mse_b,corrected_mse=mse_a,mse_ratio=mse_a/mse_b,mse_ratio_ci95=ci.tolist(),
            old_mean_ci95=np.percentile(b[index].mean(1),[2.5,97.5]).tolist(),
            corrected_mean_ci95=np.percentile(a[index].mean(1),[2.5,97.5]).tolist(),
            change_mean_ci95=np.percentile((a-b)[index].mean(1),[2.5,97.5]).tolist(),
            absolute_mean_limit=float(limit),passed=bool(passed),primary=gap in (8.,40.))
        projections.append(row)
        if row['primary']:
            if not row['mse_ratio']<1:failures.append(f'{int(gap)}ns: projected MSE ratio >=1')
            if not ci[1]<1:failures.append(f'{int(gap)}ns: projected MSE ratio 95%upper >=1')
            if not abs(a.mean())<=limit:failures.append(f'{int(gap)}ns: corrected absolute mean exceeds frozen limit')
    energies=[]
    for method in ('old','corrected'):
        for kind in ('individual','mean'):
            chosen=[r for r in details if r['method']==method and (r['repeat']>=0 if kind=='individual' else r['repeat']==-1)]
            energies.append(dict(method=method,kind=kind,mothers=336,
                **energy_metrics([r['relative_error_percent'] for r in chosen],[r['boundary'] for r in chosen],
                                  [r['converged'] for r in chosen])))
    def e(method,kind):return next(r for r in energies if r['method']==method and r['kind']==kind)
    if abs(e('corrected','mean')['bias_percent'])>5:failures.append('Isolated corrected mean E bias >5%')
    if e('corrected','individual')['rmse_percent']>1.1*e('old','individual')['rmse_percent']:
        failures.append('Isolated corrected individual E RMSE >1.1 baseline')
    if any(r['nonconverged'] for r in energies):failures.append('Unconverged isolated E fit')
    wave=[dict(method=method,mean_residual_mse_mV2=float(np.mean(np.mean(value,axis=0)**2)),
               per_mother_residual_mse_mV2=float(np.mean(np.asarray(value)**2))) for method,value in residuals.items()]
    np.savez_compressed(ART/'summary.npz',old_projection_100logE=before,corrected_projection_100logE=after,
        old_residual_mV=residuals['old'],corrected_residual_mV=residuals['corrected'],actual_mV=observed,
        time_ns=frozen.OLD.readout.adc_centers,edep_MeV=[r['edep_MeV'] for r in mothers],mother_ids=[r['id'] for r in mothers])
    save(EXP/'fit-details.json',details)
    save(EXP/'result-summary.json',dict(formal_mothers=336,pilot_mothers=8,new_ROOTs=344,projection=projections,
        energy=energies,waveform=wave,formal_digitizations=2688,pilot_digitizations=64,formal_fits=6048,pilot_fits=36,
        source_energy_MeV=8.,edep_range_MeV=[min(r['edep_MeV'] for r in mothers),max(r['edep_MeV'] for r in mothers)]))
    save(EXP/'gate.json',dict(passed=not failures,failures=failures,
        next='Independent other-energy/phase/shared controls before any pair-energy study' if not failures else
             'Retain negative transfer result; do not tune correction or expand pair-energy/high-rate work on these data'))
    save(EXP/'result-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in
         [EXP/'result-summary.json',EXP/'fit-details.json',EXP/'gate.json',ART/'summary.npz']})
    print(json.dumps(read(EXP/'gate.json'),ensure_ascii=False,indent=2),flush=True)


def replay():
    verify();assert not (EXP/'replay.json').exists()
    tasks=read(EXP/'single-plan.json');tasks=[dict(tasks[i],repeats=2) for i in (8,343)]
    with ProcessPoolExecutor(max_workers=2,initializer=initialize) as pool:
        digits=list(pool.map(frozen.simulate,[(task,'replay') for task in tasks]))
    exact=0
    for i in (8,343):
        for part in ('voltage','truth'):
            with np.load(ART/f'{part}/single/{i:04d}.npz') as a,np.load(ART/f'replay/{part}/single/{i:04d}.npz') as b:
                for key in b.files:assert np.array_equal(a[key][:2],b[key]);exact+=1
    with ProcessPoolExecutor(max_workers=1,initializer=initialize,initargs=(True,)) as pool:
        fits=list(pool.map(frozen.fit_one,[(8,128.,'fit-replay')]))
    a,b=read(ART/'fits/0008.json'),read(ART/'fit-replay/0008.json')
    assert [{k:v for k,v in r.items() if k!='seconds'} for r in a['rows']]==[
        {k:v for k,v in r.items() if k!='seconds'} for r in b['rows']]
    save(EXP/'replay.json',dict(passed=True,exact_arrays=exact,digitizations=4,fits=18,
        digitization_rows=digits,fit_rows=fits,all_fit_outputs_equal_except_timing=True))
    print('Replay passed:6 arrays and18 fits exactly equal, timing excluded',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['prepare','verify','preflight','native_pilot','pilot','native_formal','single','fit','score','replay'])
    stage=parser.parse_args().stage
    if stage=='prepare':prepare()
    elif stage=='verify':print(json.dumps(verify()))
    elif stage in ('preflight','native_pilot','native_formal'):run_native(stage)
    elif stage=='pilot':pilot()
    elif stage=='single':generate()
    elif stage=='fit':fit()
    elif stage=='score':score()
    elif stage=='replay':replay()


if __name__=='__main__':main()
