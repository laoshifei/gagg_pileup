"""Development, freezing and independent replay of conditional alpha-pair fits."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from functools import wraps
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

from .alpha_pair_conditional import ConditionalResponse, calibrate, select_count, selective_refit
from .alpha_pair_fitting import PairFitter
from .alpha_pair_resolution_data import load, save, digest, MODEL_NAMES
from .alpha_pileup_fit import score

PROJECT=Path(os.environ.get('GAGG_CONDITIONAL_PROJECT',Path(__file__).resolve().parents[2])).resolve()
NAME='alpha_pair_conditional_20261022'
EXPERIMENT=PROJECT/'pileup/experiments'/NAME
OUTPUT=PROJECT/'jupyter/artifacts'/NAME
RECORD=PROJECT/'jupyter/run_records'/NAME
PARENT=PROJECT/'pileup/experiments/alpha_pair_resolution_20261022'
OLD=PROJECT/'jupyter/artifacts/alpha_pair_resolution_20261022'
GAPS=(4.,8.,12.,20.,30.,40.,60.,80.,120.,200.)


def prepare():
    for p in (EXPERIMENT,OUTPUT,RECORD):p.mkdir(exist_ok=False,parents=True)
    cache=Path(tempfile.mkdtemp(prefix='gagg-conditional-'))
    save(RECORD/'local-cache.json',dict(directory=str(cache)))
    for name in ('codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb','Run_Catalog.ipynb'):
        shutil.copyfile(PROJECT/'jupyter'/name,RECORD/(name+'.before'))
    shutil.copyfile(PROJECT/'docs/literature-use.md',RECORD/'literature-use.before.md')
    shutil.copyfile(PARENT/'settings.json',EXPERIMENT/'settings.json')
    oldcases=load(OLD/'assets/shared.json')['cases']
    with np.load(OLD/'assets/shared.waveforms.npz') as f:y=f['voltage_mV']
    with np.load(OLD/'assets/shared.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
    with np.load(OLD/'assets/shared.pair-components.npz') as f:pair=f['components_mV']
    with np.load(OLD/'native-pool.npz') as f:pool={k:f[k] for k in f.files}
    rows=[];components=[];cases=[];train=[];pi=0
    for i,c in enumerate(oldcases):
        curves=pair[pi] if c['kind']=='pair' else np.zeros((2,1024))
        if c['kind']=='pair':pi+=1
        if c['kind']=='pair' and c['gap_ns'] not in GAPS:continue
        parents=c['parents']
        folds={(p//2)%2 for p in parents}
        assert len(folds)<=1,'Parent leakage between folds'
        fold=next(iter(folds)) if folds else c['replicate']%2
        if c['kind']=='single':
            p=parents[0];grid=np.arange(1024)*2.+.25
            curves[0]=np.interp(grid-t[i,0]+pool['times_ns'][p],
                                np.arange(1280)*2.+.25,pool['component'][p],left=0.,right=0.)
        rows.append(i);components.append(curves);cases.append(c);train.append(fold==0)
    rows=np.asarray(rows)
    np.savez_compressed(cache/'development.waveforms.npz',voltage_mV=y[rows])
    np.savez_compressed(cache/'development.truth.npz',charge_pe=q[rows],times_ns=t[rows],
                        components_mV=np.asarray(components),train=np.asarray(train))
    save(cache/'development.cases.json',cases)
    candidates={r['index']:r['candidates'] for r in map(json.loads,(OLD/'predictions/shared/candidates.jsonl').read_text().splitlines())}
    seed=dict(candidates=[candidates[i] for i in rows])
    for j,name in enumerate(MODEL_NAMES,1):
        old=load(OLD/'predictions/shared'/(name+'.json'))
        seed[f'ml_{j}']=[old['fits'][i] for i in rows]
        seed[f'old_hybrid_{j}']=[old['refined'][i]['fit'] for i in rows]
    save(cache/'development.seeds.json',seed)
    a={p for c,b in zip(cases,train) if b for p in c['parents']}
    b={p for c,b in zip(cases,train) if not b for p in c['parents']}
    assert not a&b and len(a)==len(b)==120
    save(EXPERIMENT/'development.json',dict(parent=str(OLD.relative_to(PROJECT)),
        previous_test_now_development=True,train_parents=sorted(a),validation_parents=sorted(b),
        parent_folds_disjoint=True,frames=len(rows),training_frames=int(sum(train)),
        gap_ns=list(GAPS),old_indices=rows.tolist(),
        input_sha256={p.name:digest(p) for p in cache.glob('development.*')}))
    (RECORD/'git-status-before.txt').write_text(subprocess.check_output(
        ['git','--no-optional-locks','-c','core.fileMode=false','status','--short'],cwd=PROJECT,text=True))
    save(RECORD/'environment.json',dict(python=sys.version,executable=sys.executable,
        numpy=np.__version__,dependency_path=os.environ.get('PYTHONPATH'),
        head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip()))
    print('PREPARED',len(rows),'frames',len(a),len(b),'parents',cache,flush=True)


def cache():
    path=Path(load(RECORD/'local-cache.json')['directory'])
    if not path.exists():
        path=Path(tempfile.mkdtemp(prefix='gagg-conditional-restored-'))
        for pattern in ('development.*','test.*'):
            for source in OUTPUT.glob(pattern):shutil.copyfile(source,path/source.name)
        shutil.copyfile(EXPERIMENT/'response.npz',path/'response.npz')
        save(RECORD/'local-cache.json',dict(directory=str(path),restored_from=str(OUTPUT)))
    return path


def calibration():
    c=cache()
    with np.load(c/'development.waveforms.npz') as f:y=f['voltage_mV']
    with np.load(c/'development.truth.npz') as f:
        train=f['train']&(f['charge_pe'].sum(1)>0.)
        response,info=calibrate(load(EXPERIMENT/'settings.json'),f['components_mV'][train],
            f['times_ns'][train],f['charge_pe'][train],y[train])
    response.save(EXPERIMENT/'response.npz');shutil.copyfile(EXPERIMENT/'response.npz',c/'response.npz')
    save(EXPERIMENT/'calibration.json',info);print('CALIBRATED',info,flush=True)


_RESPONSE=None
_BASE=None


def forbid_labels(event,args):
    if event=='open':
        path=str(args[0])
        if path.endswith(('.truth.npz','.cases.json','.root','native-records.json')) or '/truth/' in path or '/native-inputs/' in path:
            raise PermissionError('Inference attempted to open labels: '+path)


def observation_only(function):
    @wraps(function)
    def wrapped(*args,**kwargs):
        enabled=[True];paths=set()
        def hook(event,values):
            if not enabled[0]:return
            forbid_labels(event,values)
            if event=='open' and isinstance(values[0],(str,bytes)):
                value=str(values[0])
                if value.endswith(('.npz','.json','.jsonl','.pt')):paths.add(value)
        sys.addaudithook(hook)
        try:
            result=function(*args,**kwargs)
        finally:
            enabled[0]=False
        dataset=str(args[0]) if function.__name__=='infer' and args else 'test'
        save(RECORD/(dataset+'-'+function.__name__+'-access-audit.json'),dict(
            passed=True,opened_paths=sorted(paths),truth_inputs=[],label_access_blocked=True,
            workers_use_independent_label_guard=function.__name__=='infer'))
        return result
    return wrapped


def initialize(path):
    global _RESPONSE,_BASE
    sys.addaudithook(forbid_labels)
    _RESPONSE=ConditionalResponse.load(path);_BASE=_RESPONSE.base


def infer_one(task):
    i,y,seeds=task;start=time.monotonic()
    raw=seeds['candidates'] if seeds else _BASE.candidates(y)
    raw=[{k:v for k,v in c.items() if k!='model_mV'} for c in raw]
    metric=_RESPONSE.metric(y)
    candidates=[]
    for c in raw:
        candidates.append(_RESPONSE.refine(y,c,metric=metric))
    fitting_only=[dict(fit) for fit in candidates]
    hybrids={};ols={}
    for j in (1,2):
        initial=seeds[f'ml_{j}']
        hybrids[str(j)]=_RESPONSE.refine(y,initial,metric=metric)
        ols[str(j)]=_RESPONSE.refine(y,initial)
    # Two observed network starts can help the blind K=2 search. Neither K nor
    # source condition is supplied by a label; compare under the same metric.
    for fit in hybrids.values():
        k=len(fit['times_ns'])
        if k<len(candidates) and fit['cost']<candidates[k]['cost']:candidates[k]=fit
    return dict(index=i,candidates=candidates,original_candidates=raw,
                fitting_only_candidates=fitting_only,
                correlated_hybrids=hybrids,mean_only_hybrids=ols,seconds=time.monotonic()-start)


@observation_only
def infer(dataset,workers=4,limit=None):
    if dataset=='test':check_freeze()
    c=cache()
    with np.load(c/(dataset+'.waveforms.npz'),allow_pickle=False) as f:
        assert f.files==['voltage_mV'];y=f['voltage_mV']
    seeds=load(c/(dataset+'.seeds.json'))
    dest=c/(dataset+'.predictions.jsonl');done={}
    if dest.exists():
        for line in dest.read_text().splitlines():
            row=json.loads(line);assert row['index'] not in done;done[row['index']]=row
    indices=list(range(len(y))) if limit is None else list(range(min(len(y),limit)))
    with ProcessPoolExecutor(max_workers=workers,initializer=initialize,initargs=(c/'response.npz',),
                             mp_context=mp.get_context('spawn')) as pool:
        futures={pool.submit(infer_one,(i,y[i],{k:v[i] for k,v in seeds.items()})):i for i in indices if i not in done}
        with dest.open('a') as stream:
            for future in as_completed(futures):
                row=future.result();done[row['index']]=row;stream.write(json.dumps(row,allow_nan=False)+'\n');stream.flush()
                if len(done)%100==0 or len(done)==len(indices):print('INFER',dataset,len(done),'/',len(indices),flush=True)
    save(RECORD/(dataset+'-inference.json'),dict(events=len(done),requested=len(indices),truth_inputs=[],
        observation_sha256=digest(c/(dataset+'.waveforms.npz')),prediction_sha256=digest(dest),
        mean_seconds=float(np.mean([r['seconds'] for r in done.values()])),workers=workers))
    shutil.copyfile(dest,OUTPUT/dest.name)


def evaluate_fit(fit,t,q,case):
    ids=np.flatnonzero(q>=500.);tt=t[ids];qq=q[ids]
    truth=dict(first_active_ns=tt,emission_ns=tt,charge_pe=qq)
    normal=score(fit,truth,[0.,2048.])
    tol=min(10.,case['gap_ns']/4.) if case['kind']=='pair' else 10.
    strict=score(fit,truth,[0.,2048.],tolerance_ns=tol)
    return dict(success=strict['recovered_count']==len(ids) and normal['unmatched_predictions']==0,
        target=len(ids),recovered=normal['recovered_count'],strict_recovered=strict['recovered_count'],
        matched=normal['matched_count'],extra=normal['unmatched_predictions'],
        missed=normal['missed_count'],split=normal['split_flags'],merge=normal['merge_flags'],
        q_errors=normal['relative_charge_errors'],t_errors=normal['time_errors_ns'],
        count=len(fit['times_ns']),accepted=fit.get('accepted_refit'))


def choose_rules():
    c=cache();cases=load(c/'development.cases.json');seed=load(c/'development.seeds.json')
    pred={r['index']:r for r in map(json.loads,(c/'development.predictions.jsonl').read_text().splitlines())}
    with np.load(c/'development.truth.npz') as f:q=f['charge_pe'];t=f['times_ns'];valid=~f['train']
    assert len(pred)==len(cases)
    indices=np.flatnonzero(valid)
    def utility(fits):
        scores=[evaluate_fit(fits[i],t[i],q[i],cases[i]) for i in indices]
        groups={}
        for i,r in zip(indices,scores):
            if cases[i]['kind']!='pair':groups.setdefault(cases[i]['condition'],[]).append(r['extra']>0)
        extra={k:float(np.mean(v)) for k,v in groups.items()}
        pair=[r['success'] for i,r in zip(indices,scores) if cases[i]['kind']=='pair']
        return dict(success_fraction=float(np.mean(pair)),false_fractions=extra,
                    eligible=all(v<=.05 for v in extra.values()))
    count_selection={}
    for variant in ('candidates','fitting_only_candidates'):
        choices=[]
        for penalty in (0.,.00001,.00003,.0001,.0003,.001,.003,.01,.03,.1):
            fits=[select_count(pred[i][variant],penalty) for i in range(len(cases))]
            choices.append(dict(penalty=penalty,**utility(fits)))
        eligible=[r for r in choices if r['eligible']];assert eligible
        best=max(eligible,key=lambda r:(r['success_fraction'],r['penalty']))
        count_selection[variant]=dict(chosen=best,grid=choices)
    gates={}
    for j in (1,2):
        grid=[]
        for mode in ('mean_only_hybrids','correlated_hybrids'):
            for gap in (0.,8.,16.,24.,40.,80.,1e9):
                for change in (.1,.25,.5,1.,10.):
                    fits=[selective_refit(seed[f'ml_{j}'][i],pred[i][mode][str(j)],
                        gap_threshold_ns=gap,max_log_charge_change=change) for i in range(len(cases))]
                    grid.append(dict(mode=mode,gap_threshold_ns=gap,max_log_charge_change=change,**utility(fits)))
        eligible=[r for r in grid if r['eligible']];assert eligible
        chosen=max(eligible,key=lambda r:(r['success_fraction'],r['gap_threshold_ns'],-r['max_log_charge_change']))
        gates[str(j)]=dict(chosen=chosen,grid=grid)
    save(EXPERIMENT/'selection.json',dict(penalty=count_selection['candidates']['chosen']['penalty'],
        fitting_only_penalty=count_selection['fitting_only_candidates']['chosen']['penalty'],
        count_selection=count_selection,gates=gates,
        selection='Validation parents only; mean pair success equally weights five conditions and ten gaps',
        false_output_budget=.05,truth_inputs_at_inference=[],
        frozen_at=datetime.now().astimezone().isoformat()))
    print('SELECTED',{k:v['chosen'] for k,v in count_selection.items()},{j:g['chosen'] for j,g in gates.items()},flush=True)


def freeze():
    assert not (EXPERIMENT/'freeze.json').exists()
    assert (EXPERIMENT/'selection.json').exists()
    shutil.copytree(PARENT/'source',EXPERIMENT/'source')
    (EXPERIMENT/'bin').mkdir();shutil.copyfile(PARENT/'bin/gagg-brems',EXPERIMENT/'bin/gagg-brems')
    (EXPERIMENT/'bin/gagg-brems').chmod(0o755)
    shutil.copytree(PARENT/'models',EXPERIMENT/'models')
    shutil.copyfile(PARENT/'fitting-selection.json',EXPERIMENT/'fitting-selection.json')
    for namespace in ('analysis','pileup/analysis'):
        dest=EXPERIMENT/'inference-source'/namespace;dest.mkdir(parents=True)
        for p in (PROJECT/namespace).glob('*.py'):shutil.copyfile(p,dest/p.name)
    (EXPERIMENT/'inference-source/pileup/__init__.py').write_text('')
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for key in ('GAGG_CONDITIONAL_PROJECT','GAGG_RESOLUTION_PROJECT','GAGG_PAIR_PROJECT','GAGG_COMPONENT_PROJECT'):
    os.environ[key]=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_conditional_study import main
    main()
''')
    (EXPERIMENT/'report.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for key in ('GAGG_CONDITIONAL_PROJECT','GAGG_RESOLUTION_PROJECT','GAGG_PAIR_PROJECT','GAGG_COMPONENT_PROJECT'):
    os.environ[key]=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_conditional_report import report,figures
    if sys.argv[1]=='figures':figures()
    else:report(sys.argv[1])
''')
    p=dict(created=datetime.now().astimezone().isoformat(),gaps_ns=list(GAPS),
        energies_MeV=[2.5,5.486,8.],roots_per_energy=80,pairs_per_condition_gap=40,
        empty_frames=80,source_seed_start=262300000,digitizer_seed_start=262310000,
        conditions=['low_low','mid_mid','high_high','low_high','high_low'],
        observation='Analog mV, 2 ns x1024 samples, known baseline, adequate dynamic range',
        inference_truth_inputs=[],fixed_taus_ns=[95.,351.],
        parent_independence='Fresh optical ROOTs; parents reused across gaps/conditions, not independent exposures',
        metrics='Both Q within20%, both t within min(10ns,gap/4), no extra output; Q>=500pe',
        q_units='Realized avalanche-equivalent pe, not energy',false_output_budget=.05,
        count='Unknown local0/1/2; no source-rate or true count input',
        validation_selection=str((EXPERIMENT/'selection.json').relative_to(PROJECT)),
        test_not_started=True,limits='Pure alpha, fixed lifetimes; no mixed alpha/gamma or K>=3 claim')
    save(EXPERIMENT/'protocol.json',p)
    plans=[]
    for energy in p['energies_MeV']:
        for rep in range(80):
            plans.append(dict(label=f'alpha_{energy:g}_{rep:03d}',particle='alpha',energy_MeV=energy,
                x_mm=0.,y_mm=0.,replicate=rep,role='independent_conditional_test',gamma_model='brems',
                seed=p['source_seed_start']+len(plans)))
    oldseeds={r['seed'] for r in load(PARENT/'manifest.json')['planned_runs']}
    assert not oldseeds&{r['seed'] for r in plans}
    for name in ('macros','logs'):(EXPERIMENT/name).mkdir()
    save(EXPERIMENT/'manifest.json',dict(created=p['created'],planned_runs=plans,runs=[],
        preflight_validated=False,max_workers=4,window_us=20.,
        environment=load(PARENT/'manifest.json')['environment'],binary_sha256=digest(EXPERIMENT/'bin/gagg-brems')))
    (RECORD/'git-status-at-freeze.txt').write_text(subprocess.check_output(
        ['git','--no-optional-locks','-c','core.fileMode=false','status','--short'],cwd=PROJECT,text=True))
    files={}
    for folder in ('source','inference-source','models','bin'):
        for path in (EXPERIMENT/folder).rglob('*'):
            if path.is_file():files[str(path.relative_to(EXPERIMENT))]=digest(path)
    for name in ('settings.json','response.npz','calibration.json','selection.json','protocol.json','fitting-selection.json','run.py','report.py'):
        files[name]=digest(EXPERIMENT/name)
    save(EXPERIMENT/'freeze.json',dict(frozen_at=p['created'],files=files,
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        working_tree_status_sha256=digest(RECORD/'git-status-at-freeze.txt'),new_test_generation_started=False))
    print('FROZEN',len(files),'files',flush=True)


def check_freeze():
    frozen=load(EXPERIMENT/'freeze.json')
    for name,sha in frozen['files'].items():assert digest(EXPERIMENT/name)==sha,name


def adapt_data():
    from . import alpha_pair_resolution_data as data
    data.PROJECT=PROJECT;data.EXPERIMENT=EXPERIMENT;data.OUTPUT=OUTPUT;data.RECORD=RECORD
    data.check_freeze=check_freeze
    return data


def native(workers):
    data=adapt_data();data.native(workers)


def init_physical(local):
    from . import alpha_pair_resolution_data as data
    from . import sipm_response
    from .alpha_pair_digitizer_optimized import simulate_avalanches
    local=Path(local)
    data._DONORS=[]
    for i in range(240):
        with np.load(local/'photons'/f'{i:04}.npz') as f:data._DONORS.append({k:f[k] for k in f.files})
    data._CONFIG=load(local/'sensor-config.json')
    sipm_response.simulate_avalanches=simulate_avalanches


def physical_one(task):
    from .alpha_pair_resolution_data import joint_case
    i,case,local=task;path=Path(local)/'physical-cases'/f'{i:05}.npz'
    if path.exists():
        with np.load(path) as f:return i,(f['y'],f['q'],f['t'],f['components'],json.loads(str(f['case'])))
    row=joint_case(case)
    temporary=path.with_suffix('.partial')
    with temporary.open('wb') as stream:
        np.savez_compressed(stream,y=row[0],q=row[1],t=row[2],components=row[3],case=json.dumps(row[4]))
    temporary.replace(path)
    return i,row


def physical(workers=4):
    data=adapt_data();check_freeze();c=cache()
    (OUTPUT/'photons').mkdir(exist_ok=True);(c/'physical-cases').mkdir(exist_ok=True)
    records=data.donor_cache()
    (c/'photons').mkdir(exist_ok=True)
    for path in (OUTPUT/'photons').glob('*.npz'):
        target=c/'photons'/path.name
        if not target.exists():shutil.copyfile(path,target)
        assert digest(path)==digest(target)
    save(c/'sensor-config.json',load(PROJECT/records[0]['metadata'])['config'])
    p=load(EXPERIMENT/'protocol.json');groups={}
    for energy,name in zip(p['energies_MeV'],('low','mid','high')):
        groups[name]=[i for i,r in enumerate(records) if r['energy_MeV']==energy]
        assert len(groups[name])==80
    pairs={name+'_'+name:list(zip(group[::2],group[1::2])) for name,group in groups.items()}
    pairs['low_high']=list(zip(groups['low'][:40],groups['high'][:40]))
    pairs['high_low']=[(b,a) for a,b in pairs['low_high']]
    cases=[]
    for ci,(condition,members) in enumerate(pairs.items()):
        for gap in GAPS:
            for rep,parents in enumerate(members):
                seed=p['digitizer_seed_start']+ci*10000+rep
                cases.append(dict(kind='pair',condition=condition,gap_ns=gap,replicate=rep,
                    parents=list(parents),seed=seed,first_ns=float(np.random.default_rng(seed).uniform(96.,256.)),
                    cluster='parents-'+'-'.join(map(str,sorted(parents)))))
    for rep in range(80):
        cases.append(dict(kind='empty',condition='all',gap_ns=None,replicate=rep,parents=[],
            seed=p['digitizer_seed_start']+100000+rep,first_ns=128.,cluster=f'empty-{rep}'))
    rows=[None]*len(cases)
    # One small joint-digitization preflight before launching the remaining cases.
    init_physical(c);started=time.monotonic();i,row=physical_one((0,cases[0],c));rows[i]=row
    save(RECORD/'physical-preflight.json',dict(seconds=time.monotonic()-started,
        closure_error_mV=row[4]['component_sum_max_error_mV'],passed=True))
    with ProcessPoolExecutor(max_workers=workers,initializer=init_physical,initargs=(c,),
                             mp_context=mp.get_context('spawn')) as pool:
        futures=[pool.submit(physical_one,(i,case,c)) for i,case in enumerate(cases) if i]
        for n,future in enumerate(as_completed(futures),2):
            i,row=future.result();rows[i]=row
            if n%100==0:print('PHYSICAL',n,'/',len(rows),flush=True)
    for i,r in enumerate(records):
        with np.load(PROJECT/r['waveform']) as f:timegrid=f['time_ns'];y=f['analog_mV']
        with np.load(PROJECT/r['truth']) as f:q=float(f['charge_pe'][0]);first=float(f['first_active_ns'][0])
        pre=float(np.random.default_rng(262450000+i).uniform(96.,256.))
        start=int(np.floor((5000.-pre-timegrid[0])/2.));origin=float(timegrid[start]-.25)
        name=('low','mid','high')[p['energies_MeV'].index(r['energy_MeV'])]
        case=dict(kind='single',condition=name,gap_ns=None,replicate=r['replicate'],parents=[i],
            seed=r['seed'],cluster=f'parent-{i}',crop_first_index=start,origin_ns=origin,
            crop='Nominal trigger plus independent jitter, not true first-photon centering')
        rows.append((y[start:start+1024],np.array([q,0.]),np.array([first-origin,0.]),np.zeros((2,1024)),case))
    np.savez_compressed(c/'test.waveforms.npz',voltage_mV=np.asarray([r[0] for r in rows],np.float32))
    np.savez_compressed(c/'test.truth.npz',charge_pe=np.asarray([r[1] for r in rows]),
        times_ns=np.asarray([r[2] for r in rows]),components_mV=np.asarray([r[3] for r in rows],np.float32))
    save(c/'test.cases.json',[r[4] for r in rows])
    for path in c.glob('test.*'):shutil.copyfile(path,OUTPUT/path.name)
    save(RECORD/'physical-validation.json',dict(frames=len(rows),parents=240,workers=workers,
        max_closure_error_mV=max(r[4].get('component_sum_max_error_mV',0.) for r in rows),
        source='Joint offline digitization of fresh optical ROOTs, not native pair ROOTs',
        photon_scaling=False,input_field='analog_mV',adc_clipping=False))
    print('PHYSICAL_COMPLETE',len(rows),flush=True)


@observation_only
def test_seeds():
    """Frozen baselines receive the observation array only."""
    import torch
    from .alpha_pair_ordered import OrderedNet,predict
    check_freeze();c=cache();torch.set_num_threads(1)
    with np.load(c/'test.waveforms.npz') as f:
        assert f.files==['voltage_mV'];y=f['voltage_mV']
    seeds={};settings=load(EXPERIMENT/'settings.json')
    for j,name in enumerate(MODEL_NAMES,1):
        folder=EXPERIMENT/'models'/name
        model=OrderedNet(settings,head='physical')
        model.load_state_dict(torch.load(folder/'model.pt',map_location='cpu',weights_only=True));model.eval()
        threshold=load(folder/'selection.json')['threshold'];fits=[]
        for start in range(0,len(y),64):
            part,_=predict(model,y[start:start+64],threshold);fits.extend(part)
        seeds[f'ml_{j}']=fits
    # Store no metadata or true K/Q/t in the seed file. Fitting initialization
    # is calculated independently from voltage by the following stage.
    base=PairFitter(settings);candidates=[]
    for i,voltage in enumerate(y):
        candidates.append([{k:v for k,v in a.items() if k!='model_mV'} for a in base.candidates(voltage)])
        if (i+1)%200==0:print('BASELINE',i+1,'/',len(y),flush=True)
    seeds['candidates']=candidates
    for j in (1,2):
        seeds[f'old_hybrid_{j}']=[{k:v for k,v in base.refine(voltage,fit['times_ns'],radius=12.).items() if k!='model_mV'}
                                for voltage,fit in zip(y,seeds[f'ml_{j}'])]
    save(c/'test.seeds.json',seeds);shutil.copyfile(c/'test.seeds.json',OUTPUT/'test.seeds.json')
    save(RECORD/'test-baseline-inference.json',dict(events=len(y),truth_inputs=[],
        observation_sha256=digest(c/'test.waveforms.npz'),seeds_sha256=digest(c/'test.seeds.json')))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=('prepare','calibrate','infer','select','freeze','native','physical','seeds'))
    p.add_argument('--dataset',default='development');p.add_argument('--workers',type=int,default=4)
    p.add_argument('--limit',type=int)
    a=p.parse_args()
    if a.stage=='prepare':prepare()
    elif a.stage=='calibrate':calibration()
    elif a.stage=='infer':infer(a.dataset,a.workers,a.limit)
    elif a.stage=='select':choose_rules()
    elif a.stage=='freeze':freeze()
    elif a.stage=='native':native(a.workers)
    elif a.stage=='physical':physical(a.workers)
    elif a.stage=='seeds':test_seeds()


if __name__=='__main__':main()
