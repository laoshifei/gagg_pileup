"""Truth-assisted response diagnostics; no new blind fits or G4 transports.

Replay an immutable optical sample with shared versus separate cell histories.
The recovery-disabled configuration is a counterfactual model control only.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
from datetime import datetime
import hashlib
import importlib.metadata
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

from .sipm_response import PhotonHits, SensorConfig
from .alpha_pair_digitizer_optimized import simulate_avalanches
from .sipm_energy_calibration import invert_charge

PROJECT=Path(os.environ.get('GAGG_ENERGY_HISTORY_PROJECT',Path(__file__).resolve().parents[2])).resolve()
NAME='sipm_energy_history_20261029'
EXPERIMENT=PROJECT/'pileup/experiments'/NAME
OUTPUT=PROJECT/'jupyter/artifacts'/NAME
RECORD=PROJECT/'jupyter/run_records'/NAME
HIST=PROJECT/'jupyter/artifacts/sipm_energy_calibration_v2_20261028'
GUARD=PROJECT/'jupyter/artifacts/sipm_energy_guardband_20261029'
PARENT=PROJECT/'jupyter/artifacts/alpha_pair_joint_confirmation_20261026'
METHODS=['model_relative_residual','joint_conserved_nll']
GAPS=[8.,40.,200.]
BLOCKS=[0,5,10,15,20,25,30,35]


def load(path):return json.loads(Path(path).read_text())
def lines(path):return [json.loads(s) for s in Path(path).read_text().splitlines()]
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def now():return datetime.now().astimezone().isoformat()
def save(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.partial')
    tmp.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n');tmp.replace(path)
def write_lines(path,rows):
    with Path(path).open('w') as f:
        for r in rows:f.write(json.dumps(r,ensure_ascii=False,allow_nan=False)+'\n')


def planned_tasks():
    tasks=[]
    for slot,block in enumerate(BLOCKS):
        parents={'weak':block,'high':240+block}
        for mc in range(4):
            seed=691000000+slot*100+mc
            for recovery in [True,False]:
                base=dict(block=block,slot=slot,mc=mc,avalanche_seed=seed,recovery=recovery,window_ns=4096.)
                prefix=f'p{slot}_m{mc}_r{int(recovery)}'
                for label,parent in parents.items():
                    for offset in [0.,*GAPS]:
                        tasks.append(dict(base,id=f'{prefix}_single_{label}_{int(offset)}',kind='single',
                            parents=[parent],times_ns=[128.+offset],condition=label,offset_ns=offset,gap_ns=None))
                for condition in ['high_weak','weak_high']:
                    for gap in GAPS:
                        tasks.append(dict(base,id=f'{prefix}_{condition}_{int(gap)}',kind='pair',
                            parents=[parents[k] for k in condition.split('_')],times_ns=[128.,128.+gap],
                            condition=condition,gap_ns=gap))
            if mc==0:
                for label,parent in parents.items():
                    tasks.append(dict(id=f'p{slot}_long_{label}',block=block,slot=slot,mc=mc,
                        avalanche_seed=seed,recovery=True,window_ns=20000.,kind='long_single',
                        parents=[parent],times_ns=[128.],condition=label,offset_ns=0.,gap_ns=None))
    assert len(tasks)==912 and len({t['id'] for t in tasks})==912
    return tasks


def prepare():
    for folder in (EXPERIMENT,OUTPUT,RECORD):folder.mkdir(parents=True,exist_ok=False)
    inputs=OUTPUT/'inputs';inputs.mkdir();(inputs/'photons').mkdir()
    sources={'cases.json':HIST/'evaluation.cases.json','labels.npz':HIST/'evaluation.labels.npz',
             'fits.jsonl':GUARD/'historical-energy.predictions.jsonl','scores.jsonl':GUARD/'historical-energy.scores.jsonl',
             'parents.json':PARENT/'physical-parents.json'}
    for target,source in sources.items():shutil.copyfile(source,inputs/target)
    guards=PROJECT/'pileup/experiments/sipm_energy_guardband_20261029'
    for name in ['energy-model.json','sensor-config.json']:shutil.copyfile(guards/name,EXPERIMENT/name)
    off=load(EXPERIMENT/'sensor-config.json');assert off['sensor']['recovery_enabled']
    off['sensor']['recovery_enabled']=False
    off['counterfactual_note']='Diagnostic only: no recovery or occupancy loss; not the default or a measured device.'
    save(EXPERIMENT/'sensor-config-no-recovery.json',off)
    parents=load(inputs/'parents.json')['records'];assert len(parents)==320
    photon_index=[]
    for parent in BLOCKS+[240+b for b in BLOCKS]:
        source=PARENT/'photons'/f'{parent:04}.npz';dest=inputs/'photons'/source.name
        shutil.copyfile(source,dest);assert digest(source)==digest(dest)
        photon_index.append(dict(parent=parent,root=parents[parent]['root'],root_sha256=parents[parent]['root_original_sha256'],
            optical_seed=parents[parent]['seed'],source_energy_MeV=parents[parent]['energy_MeV'],
            photon_file=str(source.relative_to(PROJECT)),sha256=digest(source)))
    save(EXPERIMENT/'photon-inputs.json',photon_index)
    tasks=planned_tasks();oldcases=load(inputs/'cases.json')
    assert not {t['avalanche_seed'] for t in tasks}&{c['seed']+100000 for c in oldcases}
    save(EXPERIMENT/'simulation-plan.json',tasks)
    protocol=dict(created=now(),kind='Retrospective mechanistic diagnostic, not new blind confirmation',
        committed_baseline=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        historical_frames=1600,historical_pairs=1200,historical_optical_parents=320,
        decomposition='(g(Qfit)-E)/E = (g(Qnative_single)-E)/E + (g(Qshared)-g(Qnative_single))/E + (g(Qfit)-g(Qshared))/E',
        historical_caveat='Native singles use C++/20us; shared pairs use offline/4.096us; different RNG. Historical response term alone is not a causal recovery effect.',
        decomposition_summary='Means of the SAME finite, time-matched subset; medians of terms must not be added; missed/extra and matching denominators retained',
        paired_counterfactual='Same optical parents and placements, joint shared-cell versus separately digitized single injections',
        optical_blocks=BLOCKS,optical_pairs=8,optical_parents=16,digitizer_replicates=4,
        orders=['high_weak','weak_high'],source_energies_MeV=[1.625,8.],gaps_ns=GAPS,
        first_ns=128.,window_ns=4096.,observed_historical_ns=2048.,long_window_ns=20000.,
        simulations=912,shared_simulations=384,isolated_simulations=512,long_window_controls=16,max_workers=2,
        avalanche_seed_base=691000000,
        rng='Same seed across order/gap/single/joint/recovery controls; changing photons changes RNG consumption, so not eventwise identical avalanche draws; four seeds per optical pair',
        recovery_ablation={'original':True,'counterfactual':False,'all_other_sensor_parameters_unchanged':True},
        disabled_recovery='No energy mapping or blind-performance claim for this altered detector; charge interaction control only',
        charges='Primary-ancestor avalanche charge/q0; dark ancestry excluded, its occupancy and crosstalk remain simulated',
        uncertainty='Paired nested bootstrap across eight optical pairs and four digitizer replicates; descriptive for this selected finite sample',
        bootstrap_seed=202610295701,pilot_gate='Reproduce one historical joint charge vector within 1e-10 relative/1e-8pe absolute; charge closure; projected two-worker run <=1800s',
        physics_default_changed=False,new_ROOTs=False,new_training=False,new_fits=False,
        restrictions=['Pure central alpha, fixed95/351ns transport','No gamma or free tau','No rate limit or universal separation claim','Existing optical mothers reused, no fresh transport test'])
    save(EXPERIMENT/'protocol.json',protocol)
    src=EXPERIMENT/'source/pileup/analysis';src.mkdir(parents=True)
    (src.parent/'__init__.py').write_text('');(src/'__init__.py').write_text('')
    for name in ['sipm_energy_history_study.py','sipm_energy_calibration.py','sipm_response.py','sipm_readout.py','alpha_pair_digitizer_optimized.py']:
        shutil.copyfile(PROJECT/'pileup/analysis'/name,src/name)
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
os.environ['GAGG_ENERGY_HISTORY_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'source'))
if __name__=='__main__':
    from pileup.analysis.sipm_energy_history_study import main
    main()
''')
    save(EXPERIMENT/'input-hashes.json',{str(p.relative_to(OUTPUT)):digest(p) for p in inputs.rglob('*') if p.is_file()})
    files=[p for p in EXPERIMENT.rglob('*') if p.is_file()]
    save(EXPERIMENT/'freeze.json',dict(created=now(),files={str(p.relative_to(EXPERIMENT)):digest(p) for p in files},
        source='Actual new diagnostic and unchanged digitizer snapshots; historical data already exposed'))
    for name in ['codex_0007_SiPM_Energy_Calibration.ipynb','Research_Plan_20260921.ipynb','Run_Catalog.ipynb']:
        shutil.copyfile(PROJECT/'jupyter'/name,RECORD/(name+'.before'))
    save(RECORD/'environment.json',dict(python=sys.version,executable=sys.executable,threads=1,
        packages={n:importlib.metadata.version(n) for n in ['numpy','scipy','matplotlib','uproot','awkward','nbformat']}))
    print('HISTORY_PREPARED',len(tasks),'offline simulations; no new G4',flush=True)


def check_freeze():
    for name,sha in load(EXPERIMENT/'freeze.json')['files'].items():assert digest(EXPERIMENT/name)==sha,name


def initialize(folder):
    global DONORS,CONFIG
    folder=Path(folder);DONORS={}
    for p in (folder/'photons').glob('*.npz'):
        with np.load(p) as f:DONORS[int(p.stem)]={k:f[k].copy() for k in f.files}
    CONFIG=SensorConfig(**load(EXPERIMENT/'sensor-config.json')['sensor'])


def digitize(task):
    started=time.monotonic();columns={k:[] for k in ['time_ns','x_mm','y_mm','energy_eV']};owners=[];dropped=[]
    for j,(parent,first) in enumerate(zip(task['parents'],task['times_ns'])):
        donor=DONORS[parent];t=donor['time_ns']+first;keep=(t>=0)&(t<task['window_ns'])
        dropped.append(int((~keep).sum()))
        for k in columns:columns[k].append((t if k=='time_ns' else donor[k])[keep])
        owners.extend([j]*int(keep.sum()))
    arrays={k:np.concatenate(v) for k,v in columns.items()};owner=np.asarray(owners,int)
    order=np.argsort(arrays['time_ns'],kind='stable');arrays={k:v[order] for k,v in arrays.items()};owner=owner[order]
    hits=PhotonHits(**arrays,window_ns=task['window_ns'])
    result=simulate_avalanches(hits,replace(CONFIG,recovery_enabled=task['recovery']),task['avalanche_seed'])
    rec=result.records;assigned=np.full(len(rec),-1,int);direct=rec['kind']==0
    assigned[direct]=owner[rec['photon_index'][direct]]
    for i in np.flatnonzero(rec['parent_id']>=0):
        assert rec['parent_id'][i]<i
        assigned[i]=assigned[rec['parent_id'][i]]
    components=[]
    for j,parent in enumerate(task['parents']):
        block=rec[assigned==j];q=block['charge_C']/result.full_charge_C
        components.append(dict(parent=parent,charge_pe=float(q.sum()),q_before_2048_pe=float(q[block['time_ns']<2048.].sum()),
            q_before_4096_pe=float(q[block['time_ns']<4096.].sum()),direct_avalanches=int(np.sum(block['kind']==0)),
            secondary_avalanches=int(np.sum(block['kind']>=2)),mean_recovery=float(block['recovery_fraction'].mean()) if len(block) else None,
            photons_dropped=dropped[j]))
    closure=abs(sum(c['charge_pe'] for c in components)+float(rec['charge_C'][assigned<0].sum()/result.full_charge_C)-float(rec['charge_C'].sum()/result.full_charge_C))
    assert closure<1e-8 and all(c['charge_pe']>0 for c in components)
    return dict(task=task,components=components,charge_closure_error_pe=closure,seconds=time.monotonic()-started,
                counters=result.metadata['counters'])


def pilot():
    check_freeze();initialize(OUTPUT/'inputs')
    cases=load(OUTPUT/'inputs/cases.json');index=next(i for i,c in enumerate(cases) if c['kind']=='pair' and c['condition']=='high_weak' and c['gap_ns']==8. and c['replicate']==0)
    c=cases[index];task=dict(id='historical_reproduction',parents=c['parents'],times_ns=[c['first_ns'],c['first_ns']+8.],
        avalanche_seed=c['seed']+100000,recovery=True,window_ns=4096.)
    replay=digitize(task)
    with np.load(OUTPUT/'inputs/labels.npz') as f:q=f['charge_pe'][index]
    np.testing.assert_allclose([r['charge_pe'] for r in replay['components']],q,rtol=1e-10,atol=1e-8)
    tasks=load(EXPERIMENT/'simulation-plan.json');selected=[next(t for t in tasks if t['kind']=='pair' and t['condition']=='high_weak' and t['gap_ns']==8. and t['slot']==0 and t['mc']==0 and t['recovery']==r) for r in [True,False]]
    samples=[digitize(t) for t in selected];seconds=max(r['seconds'] for r in samples)
    projected=seconds*len(tasks)/2
    passed=projected<=1800.
    save(RECORD/'pilot-validation.json',dict(passed=passed,historical_case=index,historical_charge_reproduced=True,
        max_simulation_seconds=seconds,projected_two_worker_seconds=projected,inputs_sha256=load(EXPERIMENT/'input-hashes.json')))
    save(OUTPUT/'pilot-results.json',dict(reproduction=replay,samples=samples))
    print('HISTORY_PILOT',passed,'max seconds',round(seconds,3),'projected',round(projected,1),flush=True)
    assert passed,'Cost gate failed; preserve pilot and revise a separate plan'


def historical():
    check_freeze();dest=OUTPUT/'historical-components.jsonl';assert not dest.exists()
    cases=load(OUTPUT/'inputs/cases.json');model=load(EXPERIMENT/'energy-model.json')
    with np.load(OUTPUT/'inputs/labels.npz') as f:q=f['charge_pe'];e=f['energy_MeV'];t=f['times_ns']
    singles={c['parents'][0]:i for i,c in enumerate(cases) if c['kind']=='single'};assert len(singles)==320
    predictions={(r['method'],r['index']):r for r in lines(OUTPUT/'inputs/fits.jsonl')}
    scores={(r['method'],r['index']):r for r in lines(OUTPUT/'inputs/scores.jsonl') if r['method'] in METHODS}
    components=[];fitted=[];closure=[]
    for i,c in enumerate(cases):
        if c['kind']!='pair':continue
        assert np.all(q[i]>=500.) and np.all(e[i]>0)
        for j,parent in enumerate(c['parents']):
            si=singles[parent];assert np.isclose(e[i,j],e[si,0],rtol=0,atol=1e-12)
            native=float(q[si,0]);shared=float(q[i,j]);energy=float(e[i,j]);en=float(invert_charge(model,native));es=float(invert_charge(model,shared))
            assert np.isfinite([en,es]).all()
            base=dict(index=i,component=j,parent=parent,condition=c['condition'],gap_ns=c['gap_ns'],replicate=c['replicate'],
                position='first' if j==0 else 'second',species=c['condition'].split('_')[j],edep_MeV=energy,
                native_charge_pe=native,shared_charge_pe=shared,charge_change_relative=shared/native-1.,
                native_energy_relative=en/energy-1.,shared_energy_relative=es/energy-1.,response_energy_change=(es-en)/energy)
            components.append(base)
            for method in METHODS:
                p=predictions[method,i];s=scores[method,i];matches={m['target_index']:m for m in s['matches']}
                m=matches.get(j);row=dict(base,method=method,matched=m is not None,frame_success=s['frame_success_10'],missed_in_frame=s['missed'],extra_in_frame=s['extra'])
                if m is None:row.update(fit_charge_pe=None,fit_energy_relative=None,fit_energy_change=None,time_error_ns=None,cancellation=None)
                else:
                    k=m['prediction_index'];ef=p['energy_MeV'][k];qf=p['charge_pe'][k]
                    if ef is None:row.update(fit_charge_pe=qf,fit_energy_relative=None,fit_energy_change=None,time_error_ns=m['time_error_ns'],cancellation=None)
                    else:
                        total=(ef-energy)/energy;fit_change=(ef-es)/energy
                        residual=total-(base['native_energy_relative']+base['response_energy_change']+fit_change)
                        closure.append(abs(residual))
                        row.update(fit_charge_pe=qf,fit_energy_relative=total,fit_energy_change=fit_change,time_error_ns=m['time_error_ns'],
                            cancellation=bool(base['shared_energy_relative']*fit_change<0 and abs(total)<abs(base['shared_energy_relative'])))
                fitted.append(row)
    assert len(components)==2400 and len(fitted)==4800 and max(closure)<1e-12
    write_lines(dest,components);write_lines(OUTPUT/'historical-fitted-components.jsonl',fitted)
    save(RECORD/'historical-validation.json',dict(passed=True,pairs=1200,components=2400,fitted_component_rows=4800,
        single_parent_references=320,energy_identity_checked=True,additive_error_closure_max=max(closure),new_fits=False,
        baseline_window_ns=20000.,shared_charge_window_ns=4096.,observed_waveform_ns=2048.,pure_causal_interpretation=False))
    print('HISTORY_DECOMPOSITION',len(components),len(fitted),'closure',max(closure),flush=True)


def simulate(workers=2):
    assert 1<=workers<=2;check_freeze();assert load(RECORD/'pilot-validation.json')['passed']
    for name,sha in load(EXPERIMENT/'input-hashes.json').items():assert digest(OUTPUT/name)==sha,name
    tasks=load(EXPERIMENT/'simulation-plan.json');folder=OUTPUT/'simulations';folder.mkdir(exist_ok=True)
    cache=Path(tempfile.mkdtemp(prefix='gagg-energy-history-'));shutil.copytree(OUTPUT/'inputs/photons',cache/'photons')
    for p in (cache/'photons').glob('*.npz'):assert digest(p)==digest(OUTPUT/'inputs/photons'/p.name)
    save(RECORD/'simulation-cache.json',dict(path=str(cache),links_used=False,workers=workers))
    pending=[]
    for task in tasks:
        p=folder/(task['id']+'.json')
        if p.exists():assert load(p)['task']==task
        else:pending.append(task)
    start=time.monotonic();done=len(tasks)-len(pending)
    with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize,initargs=(cache,)) as pool:
        futures={pool.submit(digitize,t):t for t in pending}
        for future in as_completed(futures):
            result=future.result();assert result['task']==futures[future]
            save(folder/(result['task']['id']+'.json'),result);done+=1
            if done%32==0 or done==len(tasks):print('HISTORY_SIMULATION',done,'/',len(tasks),flush=True)
    results=[load(folder/(t['id']+'.json')) for t in tasks]
    assert all(r['charge_closure_error_pe']<1e-8 for r in results)
    save(RECORD/'simulation-validation.json',dict(passed=True,simulations=len(results),new_runs=len(pending),seconds=time.monotonic()-start,
        workers=workers,new_ROOTs=0,source_photons_preserved=True,config_difference_only_recovery_switch=True,
        maximum_charge_closure_pe=max(r['charge_closure_error_pe'] for r in results)))
    print('HISTORY_SIMULATION_COMPLETE',len(results),flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','pilot','historical','simulate']);p.add_argument('--workers',type=int,default=2)
    a=p.parse_args()
    if a.command=='simulate':simulate(a.workers)
    else:globals()[a.command]()


if __name__=='__main__':main()
