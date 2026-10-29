"""Native isolated-alpha guard-band campaign; one serial G4 event per ROOT.

Freeze before pilot, run pilot before development, freeze energy mapping before
new confirmation. All stages preserve prior attempts and historical datasets.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import threading
import time

import awkward as ak
import numpy as np

from .alpha_rate_campaign import digest
from .alpha_rate_data import export_window, open_root
from .alpha_gamma_regularized_campaign import macro

PROJECT=Path(os.environ.get('GAGG_ENERGY_GUARD_PROJECT',Path(__file__).resolve().parents[2])).resolve()
NAME='sipm_energy_guardband_20261029'
EXPERIMENT=PROJECT/'pileup/experiments'/NAME
OUTPUT=PROJECT/'jupyter/artifacts'/NAME
RECORD=PROJECT/'jupyter/run_records'/NAME
PHYSICS=PROJECT/'pileup/experiments/alpha_pair_joint_confirmation_20261026'
OLD=PROJECT/'jupyter/artifacts/sipm_energy_calibration_v2_20261028'


def load(path):return json.loads(Path(path).read_text())
def save(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def now():return datetime.now().astimezone().isoformat()


def check_freeze():
    for path,sha in load(EXPERIMENT/'freeze.json')['files'].items():assert digest(EXPERIMENT/path)==sha,path


def prepare():
    for p in (EXPERIMENT,OUTPUT,RECORD):p.mkdir(parents=True,exist_ok=False)
    for part in ('source','bin'):shutil.copytree(PHYSICS/part,EXPERIMENT/part)
    for part in ('macros','logs'):(EXPERIMENT/part).mkdir()
    dest=EXPERIMENT/'analysis-source/pileup/analysis';dest.mkdir(parents=True)
    (dest.parent/'__init__.py').write_text('');(dest/'__init__.py').write_text('')
    for name in ('sipm_energy_guard_data.py','sipm_energy_calibration.py','alpha_rate_data.py',
                 'alpha_rate_campaign.py','alpha_gamma_regularized_campaign.py'):
        shutil.copyfile(PROJECT/'pileup/analysis'/name,dest/name)
    # Preserve the exact parent frozen waveform models for any later fixed-fit replay.
    for name in ('response.npz','joint-response.npz','joint_conserved_nll-noise.npz','settings.json',
                 'gate-calibrations.json','legacy-gate-calibrations.json'):
        shutil.copyfile(PHYSICS/name,EXPERIMENT/name)
    save(EXPERIMENT/'sensor-config.json',load(PROJECT/'pileup/experiments/sipm_energy_calibration_v2_20261028/sensor-config.json'))
    shutil.copyfile(OLD/'calibration-records.json',OUTPUT/'old-calibration-records.json')
    plan=[]
    stages=[('pilot',[1.5,2.,3.5,9.],4,690000000),
            ('fit',[1.5,2.,3.5,4.5,6.5,9.],16,690010000),
            ('select',[2.,3.5,4.5,6.5],16,690020000),
            ('confirmation',[1.625,1.8,2.,2.5,3.5,4.5,5.486,6.5,7.,8.],32,690030000)]
    for stage,energies,reps,base in stages:
        for ei,energy in enumerate(energies):
            for rep in range(reps):
                plan.append(dict(label=f'{stage}_alpha_{energy:g}_{rep:03d}',stage=stage,role=stage,particle='alpha',
                    gamma_model='brems',energy_MeV=energy,x_mm=0.,y_mm=0.,replicate=rep,seed=base+ei*reps+rep))
    streams={r['seed']+offset for r in plan for offset in (0,100000,200000)}
    assert len(streams)==3*len(plan)
    audit=[]
    paths=list((PROJECT/'pileup/experiments').glob('*/manifest.json'))+list((PROJECT/'pileup/experiments').glob('*/*/manifest.json'))
    for path in paths:
        value=load(path)
        if not isinstance(value,dict):continue
        oldseeds={int(r['seed'])+offset for k in ('runs','planned_runs') for r in value.get(k,[])
                  if isinstance(r,dict) and isinstance(r.get('seed'),int) for offset in (0,100000,200000)}
        assert not streams&oldseeds,path
        audit.append(dict(path=str(path.relative_to(PROJECT)),sha256=digest(path)))
    save(RECORD/'seed-independence.json',dict(passed=True,streams=len(streams),manifests=audit))
    protocol=dict(created=now(),purpose='Repair calibration coverage and separate internal target from wider support',
        physics_unchanged=True,fixed_taus_ns=[95.,351.],sampling_ns=2.,native_window_us=20.,adequate_dynamic_range=True,
        particles='Pure central alpha; source energy is not deposited energy',
        target_source_energies_MeV=stages[-1][1],target_nominal_edep_MeV=[.2,7.7],
        source_groups_are_design_strata=True,target_selection='All events in prespecified source groups; no cuts on measured Q or true E within a group',
        guard_source_energies_MeV=[1.5,9.],intermediate_fit_source_energies_MeV=[2.,3.5,4.5,6.5],
        interpolation_test_only_source_energies_MeV=[1.8,7.],
        phases={stage:dict(energies_MeV=energies,repetitions=reps,seed_start=base) for stage,energies,reps,base in stages},
        total_new_roots=len(plan),calibration='Old280 fold0 + new96; select old280 fold1 + new64; no pilot used for fit/selection',
        new_confirmation=320,confirmation_after_energy_freeze=True,
        pilot_gate='All16 exports/config/RunInfo correct; lower guard positive deposited energy and Q with max E<0.20MeV; upper guard min E>8.0MeV; intermediate means increasing; halt on failure',
        mapping_candidates=['linear_energy','yield_scale','yield_saturation','empirical_pchip'],
        selection='Equal mean of relative-MSE across eight INTERNAL source groups, domain rejection penalty1; no guard groups in selection loss',
        fitting='Equal source-group weights; all positive native calibration Q/E; no truth-Q threshold cut',
        inference_domain='Observed fit E range expanded by physical guard data; no inverse extrapolation or endpoint clipping',
        confirmation_control='Each of ten source groups >=95% finite inverse, abs median relative energy error<=5%; report all 320 including failures',
        target_energy='Total per-primary E_dep in MeV',charge='Primary-descendant realized avalanche charge/q0, excluding dark ancestry',
        no_extra_Birks=True,
        evaluation='First new isolated truth-Q energy calibration; then fixed blind fits on historical shared-cell pairs; keep truth-assisted and blind results separate',
        limitations=['No new mixed alpha/gamma','No free tau','No new high-rate stream','No hardware calibration'],
        max_workers=4)
    save(EXPERIMENT/'protocol.json',protocol)
    old=load(PHYSICS/'manifest.json')
    environment=dict(old['environment'],python=sys.version,platform=platform.platform(),
        **{k:version(k) for k in ('numpy','scipy','uproot','awkward','matplotlib')})
    save(EXPERIMENT/'manifest.json',dict(created=now(),planned_runs=plan,runs=[],failed_attempts=[],
        preflight_validated=False,max_workers=4,window_us=20.,environment=environment,
        binary_sha256=digest(EXPERIMENT/'bin/gagg-brems')))
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
os.environ['GAGG_ENERGY_GUARD_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'analysis-source'))
if __name__=='__main__':
    from pileup.analysis.sipm_energy_guard_data import main
    main()
''')
    files=[p for part in ('source','bin','analysis-source') for p in (EXPERIMENT/part).rglob('*') if p.is_file()]
    files += [p for p in EXPERIMENT.iterdir() if p.is_file() and p.name!='manifest.json']
    save(EXPERIMENT/'freeze.json',dict(created=now(),files={str(p.relative_to(EXPERIMENT)):digest(p) for p in files},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        uncommitted_source='Actual unchanged physics and current analysis copied; prior project changes preserved',
        parent_physics_freeze_sha256=digest(PHYSICS/'freeze.json'),new_root_generation_started=False))
    for name in ('Run_Catalog.ipynb','codex_0007_SiPM_Energy_Calibration.ipynb','Research_Plan_20260921.ipynb'):
        shutil.copyfile(PROJECT/'jupyter'/name,RECORD/(name+'.before'))
    save(RECORD/'environment.json',dict(environment,executable=sys.executable,threads=1))
    print('GUARD_PREPARED',len(plan),'planned ROOTs',flush=True)


def native(stage,workers=4):
    check_freeze();assert 1<=workers<=4
    if stage in ('fit','select'):assert load(RECORD/'pilot-validation.json')['passed']
    if stage=='confirmation':
        for p,sha in load(EXPERIMENT/'energy-freeze.json')['files'].items():assert digest(EXPERIMENT/p)==sha,p
    manifest=load(EXPERIMENT/'manifest.json');done={r['label'] for r in manifest['runs']}
    plans=[p for p in manifest['planned_runs'] if p['stage']==stage and p['label'] not in done]
    if not plans:print('GUARD_NATIVE_ALREADY_COMPLETE',stage,flush=True);return
    cache=Path(tempfile.mkdtemp(prefix='gagg-energy-guard-'));binary=cache/'gagg-brems'
    shutil.copyfile(EXPERIMENT/'bin/gagg-brems',binary);binary.chmod(0o755)
    save(RECORD/(stage+'-cache.json'),dict(path=str(cache),binary_sha256=digest(binary)))
    env=os.environ.copy();env.update(manifest['environment']['g4_dataset_environment'])
    lock=threading.Lock();last=[0.];config=load(EXPERIMENT/'sensor-config.json')
    def one(plan):
        with lock:
            while True:
                wait=1.05-(time.monotonic()-last[0])
                if wait>0:time.sleep(wait)
                root=PROJECT/'g4/data'/(datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss')+'.root')
                if not root.exists():break
                time.sleep(1.05)
            attempts=list((EXPERIMENT/'macros').glob(plan['label']+'*.mac'))
            suffix='' if not attempts else f'_attempt{len(attempts)+1}'
            path=EXPERIMENT/'macros'/(plan['label']+suffix+'.mac');log=EXPERIMENT/'logs'/(plan['label']+suffix+'.log')
            row=dict(plan,root=str(root.relative_to(PROJECT)),started=now(),sipm_seed=plan['seed']+100000,
                     readout_seed=plan['seed']+200000,macro=str(path.relative_to(PROJECT)),log=str(log.relative_to(PROJECT)))
            body=macro(EXPERIMENT,row,root);path.write_text(body);row['macro_sha256']=digest(path)
            stream=log.open('w');started=time.monotonic()
            proc=subprocess.Popen([str(binary),str(path)],cwd=PROJECT/'g4/build',env=env,stdout=stream,stderr=subprocess.STDOUT)
            last[0]=time.monotonic()
        code=proc.wait();stream.close();row.update(returncode=code,elapsed_s=time.monotonic()-started,finished=now())
        if code:
            row['status']='failed';save(RECORD/'attempts'/(root.stem+'.json'),row)
            raise RuntimeError(f'G4 exited {code}: {row["log"]}')
        local=cache/root.name;shutil.copyfile(root,local);sha=digest(root);assert digest(local)==sha
        meta=export_window(local,OUTPUT/'native-inputs',allowed_particles=('alpha',))
        with open_root(local) as f:
            info=ak.to_list(f['RunInfo'].arrays(library='ak'))[0]
            assert f['WaveformEvent'].num_entries==f['SiPMEvent'].num_entries==1
            assert info['random_seed']==row['seed'] and info['macro_contents']==body
            assert info['sipm_avalanche_base_seed']==row['sipm_seed'] and info['sipm_readout_base_seed']==row['readout_seed']
            assert info['t_length_ps']==20000000
        assert meta['primaries']==1 and meta['config']==config and meta['root_sha256']==sha
        meta.update(root=row['root'],io_note='Verified regular temporary copy; no links')
        save(OUTPUT/'native-inputs/metadata'/(root.stem+'.json'),meta)
        with np.load(OUTPUT/'native-inputs/truth'/(root.stem+'.npz')) as f:
            assert len(f['edep_MeV'])==1 and str(f['particle_name'][0])=='alpha'
            assert np.isclose(float(f['source_energy_MeV'][0]),row['energy_MeV'],rtol=2e-7,atol=1e-9)
        row.update(status='completed',root_sha256=sha,root_bytes=root.stat().st_size,
            validation=dict(runinfo_macro_verified=True,seed_verified=True,config_verified=True,
                            charge_balance_error_C=meta['charge_balance_error_C']))
        save(RECORD/'attempts'/(root.stem+'.json'),row)
        return row
    if not manifest['preflight_validated']:
        row=one(plans.pop(0));manifest['runs'].append(row);manifest['preflight_validated']=True
        save(EXPERIMENT/'manifest.json',manifest);print('GUARD_PREFLIGHT',row['elapsed_s'],flush=True)
    failures=[]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(one,p):p for p in plans}
        for future in as_completed(futures):
            try:row=future.result()
            except Exception as exc:
                failure=dict(plan=futures[future],error=repr(exc));failures.append(failure)
                manifest['failed_attempts'].append(failure)
                print('GUARD_FAILED',failure,flush=True)
            else:
                manifest['runs'].append(row)
                count=sum(r['stage']==stage for r in manifest['runs'])
                print('GUARD_NATIVE',stage,count,row['label'],round(row['elapsed_s'],2),flush=True)
            save(EXPERIMENT/'manifest.json',manifest)
    if failures:raise RuntimeError(f'{len(failures)} failed attempts preserved; do not discard successful outputs')
    records=[]
    for r in sorted(manifest['runs'],key=lambda x:x['seed']):
        stem=Path(r['root']).stem
        records.append(dict(r,root_original_sha256=r['root_sha256'],source_manifest=str((EXPERIMENT/'manifest.json').relative_to(PROJECT)),
            **{key:str((OUTPUT/'native-inputs'/folder/(stem+suffix)).relative_to(PROJECT)) for key,folder,suffix in
                [('truth','truth','.npz'),('metadata','metadata','.json'),('waveform','waveforms','.npz')]}))
    save(OUTPUT/'native-records.json',dict(test_single=records))
    print('GUARD_NATIVE_COMPLETE',stage,len([r for r in records if r['stage']==stage]),flush=True)


def scalar_rows(stage):
    rows=[]
    for r in load(OUTPUT/'native-records.json')['test_single']:
        if r['stage']!=stage:continue
        with np.load(PROJECT/r['truth']) as f:
            rows.append(dict(root=r['root'],root_sha256=r['root_sha256'],seed=r['seed'],replicate=r['replicate'],
                stage=stage,source_energy_MeV=r['energy_MeV'],root_source_energy_MeV=float(f['source_energy_MeV'][0]),
                edep_MeV=float(f['edep_MeV'][0]),charge_pe=float(f['charge_pe'][0]),
                truth_file=r['truth'],truth_sha256=digest(PROJECT/r['truth']),metadata_file=r['metadata'],
                metadata_sha256=digest(PROJECT/r['metadata'])))
    return rows


def pilot_check():
    check_freeze();rows=scalar_rows('pilot');assert len(rows)==16
    groups={}
    for energy in [1.5,2.,3.5,9.]:
        rr=[r for r in rows if r['source_energy_MeV']==energy]
        groups[str(energy)]=dict(n=len(rr),edep_min=min(r['edep_MeV'] for r in rr),
            edep_max=max(r['edep_MeV'] for r in rr),charge_min=min(r['charge_pe'] for r in rr),charge_max=max(r['charge_pe'] for r in rr))
    passed=all(r['edep_MeV']>0 and r['charge_pe']>0 for r in rows)
    passed=passed and groups['1.5']['edep_max']<.2 and groups['9.0']['edep_min']>8.
    passed=passed and all(groups[str(a)]['edep_max']<groups[str(b)]['edep_min'] for a,b in [(1.5,2.),(2.,3.5),(3.5,9.)])
    save(OUTPUT/'pilot.records.json',rows)
    save(RECORD/'pilot-validation.json',dict(passed=bool(passed),created=now(),groups=groups,
        root_count=16,config_and_runinfo_verified=True,used_for_calibration=False))
    print('GUARD_PILOT_CHECK',passed,groups,flush=True)
    assert passed,'Pilot failed: preserve outputs and redesign in a new experiment'


def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['prepare','pilot','fit','select','confirmation','pilot-check'])
    parser.add_argument('--workers',type=int,default=4);a=parser.parse_args()
    if a.command=='prepare':prepare()
    elif a.command=='pilot-check':pilot_check()
    else:native(a.command,a.workers)


if __name__=='__main__':main()
