"""Count-preserving optimizer ablation with fresh independent confirmation."""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import copy
import hashlib
import importlib.metadata
import json
import multiprocessing as mp
import os
from pathlib import Path
import subprocess
import sys
import time
import unittest
import numpy as np
from .alpha_pair_blind_uncertainty_study import (
    write_new,save,save_array,load,digest,rows,now,source_closure,normalize,block_labels)
from .alpha_pair_joint_likelihood import ConservingResponse,JointNoise,waveform_candidates,select
from .alpha_pair_optimizer_polish import DEFAULTS,waveform_polish
from .alpha_pair_model_sampler import GaussianFactor,lift_innovations
from .alpha_pileup_fit import score

PROJECT=Path(os.environ.get('GAGG_POLISH_PROJECT',Path(__file__).resolve().parents[2])).resolve()
NAME='alpha_pair_optimizer_polish_20261028'
EXPERIMENT=PROJECT/'pileup/experiments'/NAME
OUTPUT=PROJECT/'jupyter/artifacts'/NAME
RECORD=PROJECT/'jupyter/run_records'/NAME
NOTEBOOK=PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
OLD_NAME='alpha_pair_blind_uncertainty_v2_20261028'
OLD=PROJECT/'jupyter/artifacts'/OLD_NAME
OLD_EXP=PROJECT/'pileup/experiments'/OLD_NAME
DEV_INDICES=[2,73,77,78,84,114,0,1,64,65,128,160,224,250,251]
METHODS=('baseline','restart','multistart','polished')

def models():
    return ConservingResponse.load(EXPERIMENT/'joint-response.npz'),JointNoise.load(EXPERIMENT/'noise.npz')

def design():
    groups=[('pair40',[18000.,800.],[200.3,240.3],64),
            ('pair200',[18000.,800.],[200.3,400.3],64),('empty',[],[],32)]
    groups += [(f'single{q}',[float(q)],[200.3],32) for q in (800,3000,8000,18800)]
    result=[]
    for band,q,t,n in groups:
        for _ in range(n):
            i=len(result);result.append(dict(index=i,band=band,count=len(q),charge_pe=q,times_ns=t,
                                             seed_words=[202610281501,3,i,997]))
    return result

def prepare():
    prior=load(PROJECT/'jupyter/run_records'/OLD_NAME/'final-validation.json')
    assert prior['passed'] and digest(NOTEBOOK)==prior['notebook_sha256']
    for folder in (EXPERIMENT,OUTPUT,RECORD):folder.mkdir(parents=True,exist_ok=False)
    for p,name in ((NOTEBOOK,'notebook.before.ipynb'),(PROJECT/'README.md','README.before.md'),
                   (PROJECT/'docs/literature-use.md','literature.before.md')):
        write_new(RECORD/name,p.read_bytes());print('BACKUP',name,flush=True)
    sources=source_closure(['alpha_pair_optimizer_polish_study','test_alpha_pair_optimizer_polish',
                           'test_alpha_pair_model_sampler'])
    for rel in ('pileup/__init__.py','pileup/analysis/__init__.py','analysis/__init__.py'):
        p=PROJECT/rel;write_new(EXPERIMENT/'source'/rel,p.read_bytes() if p.exists() else b'')
    for rel in sources:write_new(EXPERIMENT/'source'/rel,(PROJECT/rel).read_bytes())
    for name in ('joint-response.npz','noise.npz','original-gate.json'):
        write_new(EXPERIMENT/name,(OLD_EXP/name).read_bytes())
    save(EXPERIMENT/'optimizer-settings.json',DEFAULTS)
    save(EXPERIMENT/'test-design.truth.json',design())
    save(EXPERIMENT/'protocol.json',dict(created=now(),methods=list(METHODS),primary='polished',settings=DEFAULTS,
        question='Does stronger optimization of already selected K2 improve actual pair recovery?',
        count_policy='Keep original waveform_candidates/select native gate alpha0.025. Only originally selected K2 is polished; K0/K1 exact copies. No Q>=500 re-cut after polishing. Unknown true K. Cannot rescue missed K1 pairs by construction.',
        objective='Same full Gaussian quadratic+logdet,conditional mean/covariance. No new weight,prior,response or physical parameter.',
        optimization='Ordered t1/t2,Q50..50000pe,observed K1 anchor-96..+512ns enlarged to original K2 and clipped to full record. Original fit retained. Restart;6 diverse observed starts;bounded Powell. Finite search,not global proof.',
        seed_bank='Fitted times +/-1,2,4,8,16,32ns with charge identity preserved;K1 anchor,dyadic4..512ns gaps,both orders,shares.02/.05/.15/.5/.85/.95/.98. No true40/200ns starts.',
        development='15 exposed old frames,including6 lower-cost witnesses;examples not a random performance sample. Prior witnesses/truth/errors are not optimizer inputs.',
        dev_indices=DEV_INDICES,known_witness_indices=DEV_INDICES[:6],
        confirmation='288 fresh draws:64 each40/200ns,strong18000->weak800pe,t1=200.3ns;32 each empty,singles800/3000/8000/18800pe. SeedSequence[202610281501,3,index,997]. Frozen before generation;no test tuning.',
        stopping='Tests and finite,nonincreasing-cost,count-preserving development validation precede fresh generation. Configuration/primary fixed in advance;no stop based on development accuracy.',
        evaluation='Paired same-waveform comparison;all frames retained. Strict bothQ<=20%,both time<=10ns,no extra outputs. Per-band binomial intervals and paired discordance;pooled fixed strata counts only.',
        model='Purealpha fixed95/351ns,full1024x2ns float64 voltage,ample dynamic range,same conditional Gaussian model. FIR first3 residuals zero. Q equivalentpe,not depositedenergy.',
        limits='No nativeSiPM,new optical mothers,free-tau,alpha/gamma,K>=3 or rate-scan claims. Better objective need not improve truth error.',
        example_indices=[0,64],workers_limit=8,workers_used=4,blas_threads=1,new_ROOTs=0,new_training=0,default_changed=False))
    wrapper='''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
root=str(folder.parents[2])
os.environ['GAGG_POLISH_PROJECT']=root
os.environ['GAGG_UNCERTAINTY_PROJECT']=root
sys.path.insert(0,str(folder/'source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_optimizer_polish_study import main
    main()
'''
    write_new(EXPERIMENT/'run.py',wrapper.encode())
    frozen={str(p.relative_to(EXPERIMENT)):digest(p) for p in (EXPERIMENT/'source').rglob('*.py')}
    for name in ('run.py','protocol.json','optimizer-settings.json','test-design.truth.json',
                 'joint-response.npz','noise.npz','original-gate.json'):frozen[name]=digest(EXPERIMENT/name)
    save(EXPERIMENT/'freeze.json',dict(created=now(),files=frozen))
    save(RECORD/'environment.json',dict(created=now(),python=sys.version,executable=sys.executable,
        dependencies={n:importlib.metadata.version(n) for n in ('numpy','scipy','matplotlib','pandas','nbformat')},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        previous_validation_sha256=digest(PROJECT/'jupyter/run_records'/OLD_NAME/'final-validation.json'),
        existing_changes='README/docs modified;five prior uncertainty modules and experiment snapshots untracked;all preserved. New optimizer source copied exactly. Execute-bit changes ignored. Broad read-only status interrupted;scoped checks used.',
        original_inputs={str((OLD_EXP/n).relative_to(PROJECT)):digest(OLD_EXP/n) for n in
            ('joint-response.npz','noise.npz','original-gate.json')},no_commit_or_push=True))
    print('PREPARED',len(frozen),flush=True)

def check_freeze():
    for name,h in load(EXPERIMENT/'freeze.json')['files'].items():assert digest(EXPERIMENT/name)==h,name

def tests():
    check_freeze();suite=unittest.TestSuite()
    for name in ('test_alpha_pair_optimizer_polish','test_alpha_pair_model_sampler'):
        suite.addTests(unittest.defaultTestLoader.loadTestsFromName('pileup.analysis.'+name))
    result=unittest.TextTestRunner(verbosity=2).run(suite);assert result.wasSuccessful()
    save(RECORD/'tests.json',dict(passed=True,tests=result.testsRun,created=now()))

def development_inputs():
    check_freeze();assert load(RECORD/'tests.json')['passed']
    folder=OUTPUT/'development';folder.mkdir(exist_ok=False)
    by_index={r['index']:r for r in rows(OLD/'test/predictions.jsonl')}
    with np.load(OLD/'test/observations.npz') as f:y=f['voltage_mV'][DEV_INDICES]
    with np.load(OLD/'test/labels.truth.npz') as f:
        arrays={key:f[key][DEV_INDICES] for key in ('charge_pe','times_ns','count')}
    cases0=load(OLD/'test/labels.cases.json');cases=[];baselines=[]
    for i,j in enumerate(DEV_INDICES):
        cases.append(dict(cases0[j],index=i,original_index=j,known_witness=j in DEV_INDICES[:6]))
        p=by_index[j];baselines.append(dict(index=i,candidates=p['candidates'],selected=p['selected'],
            prior_exterior=p['exterior']['best_exterior'] if p['exterior'] else None))
    save_array(folder/'observations.npz',voltage_mV=y);save_array(folder/'labels.truth.npz',**arrays)
    save(folder/'labels.cases.json',cases);save(folder/'baseline.json',baselines)
    save(RECORD/'development-inputs.json',dict(passed=True,created=now(),frames=len(y),
        old_predictions_sha256=digest(OLD/'test/predictions.jsonl'),
        files={p.name:digest(p) for p in folder.iterdir()},exposed_development_only=True))

def generate():
    check_freeze();assert load(RECORD/'development-validation.json')['passed']
    folder=OUTPUT/'test';folder.mkdir(exist_ok=False)
    cases=load(EXPERIMENT/'test-design.truth.json');response,noise=models();waves=[]
    q=np.zeros((len(cases),2));t=q.copy();error=0.
    for c in cases:
        curves=response.components(c['times_ns'],c['charge_pe'])
        factor=GaussianFactor.from_noise(noise,curves,c['charge_pe'])
        z=factor.draw(np.random.default_rng(np.random.SeedSequence(c['seed_words'])))
        residual=lift_innovations(z,noise.coefficients)
        error=max(error,float(np.max(np.abs(noise.transform(residual)-z))))
        assert np.isclose(noise.score(residual,curves,c['charge_pe']),
                          np.sum(factor.whiten(z)**2)+factor.logdet,rtol=1e-10,atol=1e-8)
        waves.append(curves.sum(0)+response.baseline+noise.background_mV+residual)
        q[c['index'],:c['count']]=c['charge_pe'];t[c['index'],:c['count']]=c['times_ns']
    y=np.asarray(waves,dtype=np.float64);assert y.shape==(288,1024) and np.isfinite(y).all() and error<1e-9
    hashes=[hashlib.sha256(row.tobytes()).hexdigest() for row in y];assert len(set(hashes))==288
    old_hashes=set()
    for split in ('calibration','test'):
        with np.load(OLD/split/'observations.npz') as f:
            old_hashes.update(hashlib.sha256(row.tobytes()).hexdigest() for row in f['voltage_mV'])
    assert not old_hashes.intersection(hashes)
    save_array(folder/'observations.npz',voltage_mV=y)
    save_array(folder/'labels.truth.npz',charge_pe=q,times_ns=t,count=np.array([c['count'] for c in cases]))
    save(folder/'labels.cases.json',cases)
    save(RECORD/'test-generation.json',dict(passed=True,created=now(),frames=288,max_fir_error=error,
        source_freeze_sha256=digest(EXPERIMENT/'freeze.json'),new_seed_root=202610281501,
        unique_new_observations=True,disjoint_from_previous_432_observations=True,
        files={p.name:digest(p) for p in folder.iterdir()}))
    print('GENERATED_TEST',len(y),flush=True)

def initialize(mode):
    global RESPONSE,NOISE,OBS,GATE,CONFIG,BASELINE
    RESPONSE,NOISE=models();sys.addaudithook(block_labels)
    try:
        with (EXPERIMENT/'test-design.truth.json').open():raise AssertionError('Label guard failed')
    except PermissionError:pass
    GATE=load(EXPERIMENT/'original-gate.json');CONFIG=load(EXPERIMENT/'optimizer-settings.json')
    with np.load(OUTPUT/mode/'observations.npz') as f:OBS=f['voltage_mV']
    BASELINE=load(OUTPUT/mode/'baseline.json') if mode=='development' else None

def infer_one(index):
    started=time.monotonic();y=OBS[index];failure=None;polish_failure=None
    try:
        candidates=BASELINE[index]['candidates'] if BASELINE is not None else waveform_candidates(RESPONSE,NOISE,y)
        selected=select(candidates,GATE)
        if BASELINE is not None:assert selected==BASELINE[index]['selected']
    except (ArithmeticError,ValueError,np.linalg.LinAlgError) as exc:
        failure=type(exc).__name__+': '+str(exc)
        candidates=[];selected=dict(times_ns=[],charge_pe=[],cost=None,converged=False)
    try:
        diagnostic=waveform_polish(RESPONSE,NOISE,y,candidates,selected,CONFIG) if not failure else None
    except (ArithmeticError,ValueError,np.linalg.LinAlgError) as exc:
        diagnostic=None;polish_failure=type(exc).__name__+': '+str(exc)
    fits=dict(baseline=selected)
    fits.update(diagnostic['stages'] if diagnostic else {m:copy.deepcopy(selected) for m in METHODS[1:]})
    assert all(len(f['times_ns'])==len(selected['times_ns']) for f in fits.values())
    return dict(index=index,candidates=candidates,fits=fits,diagnostic=diagnostic,
                baseline_failure=failure,polish_failure=polish_failure,truth_inputs=[],seconds=time.monotonic()-started)

def infer(mode,workers):
    check_freeze();assert 1<=workers<=8
    assert load(RECORD/('test-generation.json' if mode=='test' else 'development-inputs.json'))['passed']
    with np.load(OUTPUT/mode/'observations.npz') as f:n=len(f['voltage_mV'])
    path=OUTPUT/mode/'predictions.jsonl';fd=os.open(path,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o644)
    offset=0;collected=[];started=time.monotonic()
    try:
        with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),
                initializer=initialize,initargs=(mode,)) as pool:
            future={pool.submit(infer_one,i):i for i in range(n)}
            for job in as_completed(future):
                r=job.result();payload=(json.dumps(normalize(r),allow_nan=False)+'\n').encode();written=0
                while written<len(payload):
                    step=os.pwrite(fd,payload[written:],offset+written)
                    if step<=0:raise OSError('Short prediction write')
                    written+=step
                offset+=len(payload);os.fsync(fd);collected.append(r)
                gain=(r['fits']['baseline']['cost']-r['fits']['polished']['cost']) if not r['baseline_failure'] else None
                print('POLISH',mode,len(collected),n,'index',r['index'],'K',len(r['fits']['baseline']['times_ns']),
                      'gain',round(gain,5) if gain is not None else None,'s',round(r['seconds'],2),flush=True)
        os.ftruncate(fd,offset);os.fsync(fd)
    finally:os.close(fd)
    check=rows(path);assert len(check)==n and {r['index'] for r in check}==set(range(n))
    assert sorted(check,key=lambda r:r['index'])==sorted(collected,key=lambda r:r['index'])
    save(RECORD/(mode+'-inference.json'),dict(passed=True,created=now(),frames=n,sha256=digest(path),
        seconds=time.monotonic()-started,workers=workers,labels_blocked=True,truth_inputs=[],
        baseline_failures=sum(bool(r['baseline_failure']) for r in check),
        polish_failures=sum(bool(r['polish_failure']) for r in check)))
    print('INFERENCE_COMPLETE',mode,n,flush=True)

def score_mode(mode):
    check_freeze();frozen=load(RECORD/(mode+'-inference.json'));assert frozen['passed']
    assert digest(OUTPUT/mode/'predictions.jsonl')==frozen['sha256']
    predictions={r['index']:r for r in rows(OUTPUT/mode/'predictions.jsonl')}
    cases=load(OUTPUT/mode/'labels.cases.json');gate=load(EXPERIMENT/'original-gate.json')
    scored=[];max_increase=0.;response,noise=models()
    with np.load(OUTPUT/mode/'observations.npz') as f:waves=f['voltage_mV']
    max_replay=0.
    for c in cases:
        r=predictions[c['index']];baseline=r['fits']['baseline'];n=c['count']
        if not r['baseline_failure']:assert select(r['candidates'],gate)==baseline
        truth=dict(first_active_ns=c['times_ns'],emission_ns=c['times_ns'],charge_pe=c['charge_pe'])
        for method in METHODS:
            fit=r['fits'][method];assert len(fit['times_ns'])==len(baseline['times_ns'])
            if len(baseline['times_ns'])!=2:assert fit==baseline
            evaluation=score(fit,truth,[0.,2048.],tolerance_ns=10.)
            strict=bool(not r['baseline_failure'] and evaluation['recovered_count']==n and evaluation['unmatched_predictions']==0)
            gain=None;distance=None;qerror=None;terror=None
            if not r['baseline_failure']:
                curves=response.components(fit['times_ns'],fit['charge_pe'])
                residual=waves[c['index']]-response.baseline-noise.background_mV-curves.sum(0)
                value=noise.score(residual,curves,fit['charge_pe'])
                max_replay=max(max_replay,abs(value-fit['cost'])/max(1.,abs(value)))
                gain=baseline['cost']-fit['cost'];max_increase=max(max_increase,-gain)
            if n==2 and len(fit['times_ns'])==2:
                qerror=np.asarray(fit['charge_pe'])/c['charge_pe']-1.
                terror=np.asarray(fit['times_ns'])-c['times_ns']
                distance=float(max(np.max(abs(qerror)/.2),np.max(abs(terror)/10.)))
                assert strict==(distance<=1.+1e-12)
            scored.append(dict(index=c['index'],original_index=c.get('original_index'),
                known_witness=c.get('known_witness',False),band=c['band'],true_count=n,method=method,
                predicted_count=len(fit['times_ns']),strict_success=strict,evaluation=evaluation,
                cost_gain=gain,accuracy_distance=distance,ordered_Q_errors=qerror,ordered_t_errors_ns=terror,
                baseline_failure=r['baseline_failure'],polish_failure=r['polish_failure'],
                converged=fit.get('converged',False),origin=fit.get('origin','baseline')))
    assert max_increase<1e-7 and max_replay<1e-9
    assert not any(r['baseline_failure'] or r['polish_failure'] for r in predictions.values())
    save(OUTPUT/mode/'scores.truth.json',scored)
    save(RECORD/(mode+'-validation.json'),dict(passed=True,created=now(),frames=len(cases),
        prediction_sha256=frozen['sha256'],scores_sha256=digest(OUTPUT/mode/'scores.truth.json'),
        count_preserved=True,non_K2_exactly_unchanged=True,maximum_cost_increase=max_increase,
        maximum_cost_replay_error=max_replay,zero_baseline_and_polish_failures=True,
        finite_search_only=True,exposed_development_only=mode=='development'))
    for band in dict.fromkeys(c['band'] for c in cases):
        print('SCORED',mode,band,{m:sum(r['strict_success'] for r in scored if r['band']==band and r['method']==m)
                               for m in METHODS},flush=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['prepare','tests','development-inputs',
        'infer-development','score-development','generate','infer-test','score-test'])
    parser.add_argument('--workers',type=int,default=4);args=parser.parse_args()
    actions={'prepare':prepare,'tests':tests,'development-inputs':development_inputs,
             'infer-development':lambda:infer('development',args.workers),'score-development':lambda:score_mode('development'),
             'generate':generate,'infer-test':lambda:infer('test',args.workers),'score-test':lambda:score_mode('test')}
    actions[args.command]()

if __name__=='__main__':main()
