"""Frozen continuous-charge training and native-transfer experiment."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
from functools import lru_cache
import importlib.metadata
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

from .alpha_pair_information_study import load,save,digest,rows
from .alpha_pair_risk_learning import (Adam,RiskNet,RISKS,TIME_GRID,INPUT_SIZE,charge_risk,
    decode,objective,single_features,regularized_refine)

PROJECT=Path(os.environ.get('GAGG_RISK_PROJECT',Path(__file__).resolve().parents[2])).resolve()
NAME='alpha_pair_risk_training_20261027'
EXPERIMENT=PROJECT/'pileup/experiments'/NAME
OUTPUT=PROJECT/'jupyter/artifacts'/NAME
RECORD=PROJECT/'jupyter/run_records'/NAME
PREVIOUS=PROJECT/'pileup/experiments/alpha_pair_charge_allocation_20261027'
OLD_OUTPUT=PROJECT/'jupyter/artifacts/alpha_pair_charge_allocation_20261027'


def models():
    from .alpha_pair_joint_likelihood import ConservingResponse,JointNoise
    return ConservingResponse.load(EXPERIMENT/'joint-response.npz'),JointNoise.load(EXPERIMENT/'joint_conserved_nll-noise.npz')


def check_freeze():
    for name,expected in load(EXPERIMENT/'freeze.json')['files'].items():assert digest(EXPERIMENT/name)==expected,name


def prepare():
    previous=load(PROJECT/'jupyter/run_records/alpha_pair_charge_allocation_20261027/final-validation.json')
    assert previous['passed']
    notebook=PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    assert digest(notebook)==previous['notebook_sha256']
    for folder in (EXPERIMENT,OUTPUT,RECORD):folder.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(notebook,RECORD/'notebook.before.ipynb')
    for name in ('joint-response.npz','joint_conserved_nll-noise.npz','baseline-gate.json','settings.json'):
        shutil.copyfile(PREVIOUS/name,EXPERIMENT/name)
    protocol=dict(created=datetime.now(timezone.utc).isoformat(),scope='Pure alpha; fixed95/351ns; complete2ns x1024 analog; local unknownK<=2; ample dynamic range',
        design='Same NumPy MLP, initialization, minibatch stream and continuous training data; absolute/logarithmic/relative charge risks compared',
        synthetic='Exact frozen conditional Gaussian likelihood sampler; Q is a conditioning parameter, not newly digitized avalanche charge or deposited energy',
        training=dict(frames=9600,validation=960,counts=[.1,.3,.6],
            Q='Single log-uniform500..25000pe; pairs half independent log-uniform500..25000, quarter weak500..3000 +strong6000..25000 in random order, quarter balanced common1000..16000 with ratio exp(uniform[-.35,.35]); reject any individual outside500..25000 or total>48000',
            gaps='Equal mixture uniform4..320ns and log-uniform4..320ns; first onset uniform96..288ns',
            prior='Declared continuous simulation-design distribution only; not a measured source spectrum or a prior fitted to exposed four-energy evaluation labels',
            seeds=dict(train=270310001,validation=270310002,evaluation=270310003),
            network_seed=270310010,minibatch_seed=270310011,steps=3000,batch=64,rate=.0007,
            hidden=[128,96],features='1021 full FIR innovations +385 observed-single-aligned innovations +log pilotQ and pilot time; asinh(z/3), training-only mean/std, standardized clip[-8,8]',
            one_pilot='Derivative start and integralQ; same single-hypothesis Gaussian NLL, radius24ns,max35 iterations; observed hypothesis, not true K',
            outputs='Two chronological candidate times on -384..384ns/2ns relative grid, locally averaged within two bins of each mode; totalQ pilot*exp(.1*tanh(z)), allocation sigmoid; count2 logit',
            loss='BCE(count2)+3*mean(pair charge risk)+.4*mean(pair time CE); charge risks ((Qhat-Q)/5000)^2,log(Qhat/Q)^2,(Qhat/Q-1)^2; targets used only during training',
            selection='Every500 steps including zero: highest validation forced-K2 strict pair success, then smallest median absolute logQ error, then lowest common count BCE; not blind test performance'),
        evaluation=dict(continuous=1000,native=760,native_gate=280,
            targeted='balanced/high_weak/weak_high x gap bands6..12,24..64,128..256ns x40 independent Gaussian draws; weak600..1500,strong12000..24000; balanced common3000..15000 times exp(+-uniform[0,.2])',
            prior_shift='240 pairs: independent uniform500..23000pe and uniform4..320ns; changes training Q/gap distributions',
            controls='80 empty and80 singles in each continuous Q bin500..1500,1500..5000,5000..12500,12500..25000',
            independence='Continuous evaluation generated only after weight/selection freeze; new parameter/noise seeds; same calibrated mean/noise model, not fresh physical optical parents. Existing native760 frames are exposed transfer diagnostics only.'),
        gate='Each output calibrates learned count2 logit on the same280 native singles, original fitted-singleQ bins1500/5000/12500; plus-one tie-preserving upper rank<=.025; candidates need eachQ>=500; preserve baselineK0/K1 presence',
        hybrids='Choose risk by the same validation key only. Compare original-NLL refinement vs NLL plus learned-logQ penalty from that candidate; widths=max(.2,validation68th percentile absolute logQ errors) per chronological slot. Same waveform enters network and fit: data-dependent regularization, NOT independent Bayesian evidence.',
        hybrid_budget='Original and learned two-component starts, each max70 iterations within24ns; preserve best initial candidate under the appropriate objective; all failures retained; count gate separately calibrated',
        selection='Do not replace default. Development promise requires continuous pair strict success at least baseline+5pp, continuous and native single controls(extra<=5%,success>=95%,medianabsQbias<=5%), and no native pooled pair loss. All comparisons retained.',
        metrics='Each trueQ>=500; eachQwithin20%; eachtwithinmin(10ns,gap/4); no extra outputs; also missing,extra,split,merge,Q bias/sigma68,t errors,forcedK2 and both-time matches',
        examples='Native high_weak replicate0 at40 and200ns; full0..2048ns; selected outputs only; absent weak output=0',
        figures=['0004-87-1','0004-88-1','0004-89-1','0004-90-1','0004-90-2','0004-91-1'],
        literature='Reuse Morad2024 detection/estimation separation and caution on trained scores vs likelihood; Kim2023(NIMA168492) supervision/weak-component and synthetic/physical distinction; project MLP, losses and regularization are not paper reproductions')
    save(EXPERIMENT/'protocol.json',protocol)
    for name in ('observations.npz','labels.truth.npz','labels.cases.json','predictions.baseline.jsonl'):
        source=OLD_OUTPUT/name
        if name.endswith('.npz'):
            with np.load(source) as f:values={k:f[k][:1040] for k in f.files}
            np.savez_compressed(OUTPUT/('native-'+name),**values)
        elif name.endswith('.json'):
            save(OUTPUT/('native-'+name),load(source)[:1040])
        else:
            with (OUTPUT/('native-'+name)).open('w') as stream:
                for r in rows(source):
                    if r['index']<1040:stream.write(json.dumps(r)+'\n')
    provenance=dict(previous_validation=previous,input_files={str((OLD_OUTPUT/n).relative_to(PROJECT)):digest(OLD_OUTPUT/n) for n in ('observations.npz','labels.truth.npz','labels.cases.json','predictions.baseline.jsonl')},
        native_files={p.name:digest(p) for p in OUTPUT.glob('native-*')},
        mother_groups=load(PROJECT/'jupyter/run_records/alpha_pair_charge_allocation_20261027/input-provenance.json')['mother_groups'],
        notebook_before_sha256=digest(notebook),new_physical_mothers=0)
    save(RECORD/'input-provenance.json',provenance)
    save(RECORD/'environment.json',dict(python=sys.version,executable=sys.executable,
        dependencies={n:importlib.metadata.version(n) for n in ('numpy','scipy','matplotlib','pandas','nbformat')},
        max_workers=4,blas_threads=1,cpu_affinity=sorted(os.sched_getaffinity(0)),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        uncommitted='Information and allocation work retained; actual new sources frozen; NumPy learner adds no external ML dependency'))
    (RECORD/'git-status.before.txt').write_text(subprocess.check_output(['git','--no-optional-locks','-c','core.fileMode=false','status','--short','--','README.md','docs/literature-use.md','pileup/analysis'],cwd=PROJECT,text=True))
    # Reuse the previously checked dependency snapshot; override only this new experiment's files.
    shutil.copytree(PREVIOUS/'inference-source',EXPERIMENT/'inference-source')
    for name in ('alpha_pair_risk_learning.py','alpha_pair_risk_study.py','test_alpha_pair_risk_learning.py'):
        shutil.copyfile(PROJECT/'pileup/analysis'/name,EXPERIMENT/'inference-source/pileup/analysis'/name)
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for name in ('RISK','INFORMATION'):
    os.environ['GAGG_'+name+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_risk_study import main
    main()
''')
    files={str(p.relative_to(EXPERIMENT)):digest(p) for p in (EXPERIMENT/'inference-source').rglob('*.py')}
    for name in ('joint-response.npz','joint_conserved_nll-noise.npz','baseline-gate.json','settings.json','protocol.json','run.py'):files[name]=digest(EXPERIMENT/name)
    save(EXPERIMENT/'freeze.json',dict(files=files,created=protocol['created'],purpose='Actual inference/training sources and old physical likelihood; no detector change'))
    print('RISK_PREPARED',len(files),'frozen files',flush=True)


@lru_cache(maxsize=1)
def design():
    return load(EXPERIMENT/'protocol.json')


def plan_case(role,index):
    seeds=design()['training']['seeds']
    seed=seeds[role]*100000+index
    rng=np.random.default_rng(seed);first=float(rng.uniform(96.,288.))
    logdraw=lambda low,high:float(np.exp(rng.uniform(np.log(low),np.log(high))))
    condition='broad';band=None
    if role=='evaluation':
        if index<360:
            group=index//40;condition=('balanced','high_weak','weak_high')[group//3]
            band=('short','middle','long')[group%3];lo,hi=((6.,12.),(24.,64.),(128.,256.))[group%3]
            gap=float(rng.uniform(lo,hi));count=2
            if condition=='balanced':
                common=logdraw(3000.,15000.);ratio=float(rng.uniform(0.,.2));q=np.array([common*np.exp(ratio),common*np.exp(-ratio)])
                if rng.random()<.5:q=q[::-1]
            else:
                q=np.array([logdraw(12000.,24000.),logdraw(600.,1500.)])
                if condition=='weak_high':q=q[::-1]
        elif index<600:
            condition='uniform_shift';band='all';count=2;gap=float(rng.uniform(4.,320.));q=rng.uniform(500.,23000.,2)
        elif index<680:condition='all';count=0;gap=None;q=np.empty(0)
        else:
            group=(index-680)//80;condition=('weak','low','mid','high')[group]
            lo,hi=((500.,1500.),(1500.,5000.),(5000.,12500.),(12500.,25000.))[group]
            count=1;gap=None;q=np.array([logdraw(lo,hi)])
    else:
        count=int(rng.choice((0,1,2),p=(.1,.3,.6)));gap=None;q=np.empty(0)
        if count==1:q=np.array([logdraw(500.,25000.)])
        elif count==2:
            gap=logdraw(4.,320.) if rng.random()<.5 else float(rng.uniform(4.,320.))
            while True:
                kind=float(rng.random())
                if kind<.5:q=np.array([logdraw(500.,25000.),logdraw(500.,25000.)])
                elif kind<.75:
                    q=np.array([logdraw(500.,3000.),logdraw(6000.,25000.)])
                    if rng.random()<.5:q=q[::-1]
                else:
                    common=logdraw(1000.,16000.);ratio=float(rng.uniform(-.35,.35));q=np.array([common*np.exp(ratio),common*np.exp(-ratio)])
                if q.min()>=500. and q.max()<=25000. and q.sum()<=48000.:break
    times=np.array([first,first+gap]) if count==2 else np.array([first]) if count==1 else np.empty(0)
    return dict(index=index,seed=seed,role=role,domain='continuous',kind='pair' if count==2 else 'single' if count==1 else 'empty',
        condition=condition,band=band,gap_ns=gap,count=count,times_ns=times.tolist(),charge_pe=q.tolist(),replicate=index%40)


def initialize():
    global RESPONSE,NOISE
    RESPONSE,NOISE=models()


def make_one(case):
    from .alpha_pair_model_sampler import GaussianFactor,lift_innovations
    rng=np.random.default_rng(np.random.SeedSequence([case['seed'],993]))
    curves=RESPONSE.components(case['times_ns'],case['charge_pe'])
    factor=GaussianFactor.from_noise(NOISE,curves,case['charge_pe'])
    residual=lift_innovations(factor.draw(rng),NOISE.coefficients)
    y=curves.sum(0)+RESPONSE.baseline+NOISE.background_mV+residual
    feature,pilot=single_features(RESPONSE,NOISE,y)
    return case['index'],y.astype(np.float32),feature,pilot


def generate(role,workers=4,pilot=False):
    check_freeze();protocol=load(EXPERIMENT/'protocol.json')
    if role=='evaluation':assert load(EXPERIMENT/'training-freeze.json')['passed']
    size=protocol['training']['frames'] if role=='train' else protocol['training']['validation'] if role=='validation' else protocol['evaluation']['continuous']
    if pilot:size=24
    if not pilot:assert load(RECORD/'pilot-validation.json')['passed']
    folder=OUTPUT/(role+('-pilot' if pilot else ''));folder.mkdir(exist_ok=False)
    cases=[plan_case(role,i) for i in range(size)]
    y=np.empty((size,1024),np.float32);x=np.empty((size,INPUT_SIZE),np.float32);pilots=[None]*size
    started=time.monotonic()
    with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize) as pool:
        futures=[pool.submit(make_one,c) for c in cases]
        for n,future in enumerate(as_completed(futures),1):
            i,wave,feature,one=future.result();y[i]=wave;x[i]=feature;pilots[i]=one
            if n%500==0:print('RISK_GENERATE',role,n,size,flush=True)
    q=np.zeros((size,2));t=np.zeros((size,2));count=np.array([c['count'] for c in cases])
    for i,c in enumerate(cases):q[i,:c['count']]=c['charge_pe'];t[i,:c['count']]=c['times_ns']
    anchor=np.array([p['times_ns'][0] for p in pilots]);oneq=np.array([p['charge_pe'][0] for p in pilots])
    outside=int(np.sum((np.abs(t-anchor[:,None])>384.)&(np.arange(2)[None,:]<count[:,None])))
    assert outside==0,'Training localization range is insufficient'
    np.savez_compressed(folder/'observations.npz',voltage_mV=y)
    np.savez_compressed(folder/'features.npz',features=x,one_Q_pe=oneq,one_time_ns=anchor)
    np.savez_compressed(folder/'labels.truth.npz',charge_pe=q,times_ns=t,count=count)
    save(folder/'labels.cases.json',cases);save(folder/'pilots.json',pilots)
    save(RECORD/(role+('-pilot' if pilot else '')+'-generation.json'),dict(passed=True,frames=size,seconds=time.monotonic()-started,
        pilot_not_converged=sum(not p['converged'] for p in pilots),outside_time_grid=outside,all_failed_pilots_retained=True,
        synthetic_not_native=True,files={str(p.relative_to(OUTPUT)):digest(p) for p in folder.iterdir()}))
    print('RISK_GENERATED',role,size,flush=True)


def validation_metrics(fits,output,truth):
    from .alpha_pileup_fit import score
    q,t,count=truth['charge_pe'],truth['times_ns'],truth['count'];success=[];errors=[]
    for i in np.flatnonzero(count==2):
        target=dict(first_active_ns=t[i],emission_ns=t[i],charge_pe=q[i])
        s=score(fits[i],target,[0.,2048.],tolerance_ns=min(10.,(t[i,1]-t[i,0])/4.))
        success.append(s['recovered_count']==2 and s['unmatched_predictions']==0)
        errors.extend(np.abs(np.log(np.asarray(fits[i]['charge_pe'])/q[i])).tolist())
    bce=float(np.mean(np.logaddexp(0.,output[:,0])-(count==2)*output[:,0]))
    return dict(forced_pair_success=float(np.mean(success)),median_abs_logQ=float(np.median(errors)),count_BCE=bce,
                key=[float(np.mean(success)),-float(np.median(errors)),-bce])


def train(kind):
    check_freeze();assert kind in RISKS
    protocol=load(EXPERIMENT/'protocol.json')['training']
    with np.load(OUTPUT/'train/features.npz') as f:training={k:f[k] for k in f.files}
    with np.load(OUTPUT/'train/labels.truth.npz') as f:truth={k:f[k] for k in f.files}
    with np.load(OUTPUT/'validation/features.npz') as f:validation={k:f[k] for k in f.files}
    with np.load(OUTPUT/'validation/labels.truth.npz') as f:valid_truth={k:f[k] for k in f.files}
    x=training['features'];mean=x.mean(0,dtype=np.float64);scale=np.maximum(x.std(0,dtype=np.float64),.05)
    model=RiskNet(mean,scale,protocol['network_seed'],hidden=tuple(protocol['hidden']))
    folder=EXPERIMENT/'models'/kind;folder.mkdir(parents=True,exist_ok=False)
    optimizer=Adam(model.parameters,rate=protocol['rate']);rng=np.random.default_rng(protocol['minibatch_seed'])
    history=[];best=None;started=time.monotonic()
    for step in range(protocol['steps']+1):
        if step:
            ids=rng.integers(0,len(x),protocol['batch']);z,saved=model.forward(x[ids],cache=True)
            loss,gradient,parts=objective(z,training['one_Q_pe'][ids],truth['charge_pe'][ids],truth['times_ns'][ids],
                                          training['one_time_ns'][ids],truth['count'][ids],kind)
            assert np.isfinite(loss) and np.isfinite(gradient).all()
            norm=optimizer.update(model.parameters,model.backward(gradient,saved))
        if step%500==0:
            output=model.forward(validation['features']);fits=decode(output,validation['one_Q_pe'],validation['one_time_ns'])
            metric=validation_metrics(fits,output,valid_truth)
            value,_,parts=objective(output,validation['one_Q_pe'],valid_truth['charge_pe'],valid_truth['times_ns'],validation['one_time_ns'],valid_truth['count'],kind)
            history.append(dict(step=step,seconds=time.monotonic()-started,validation_objective=value,parts=parts,**metric))
            if best is None or metric['key']>best['key']:
                best=dict(step=step,**metric);model.save(folder/'model.npz',dict(risk=kind,step=step,truth_inputs='Training loss only'))
                save(folder/'selection.json',best)
            save(folder/'history.json',history)
            print('RISK_TRAIN',kind,step,'forced_pair',round(metric['forced_pair_success'],4),flush=True)
    model=RiskNet.load(folder/'model.npz');z=model.forward(validation['features']);fits=decode(z,validation['one_Q_pe'],validation['one_time_ns'])
    pair=np.flatnonzero(valid_truth['count']==2)
    errors=np.abs(np.log(np.array([fits[i]['charge_pe'] for i in pair])/valid_truth['charge_pe'][pair]))
    sigma=np.maximum(.2,np.percentile(errors,68.,axis=0))
    save(folder/'validation.json',dict(metrics=validation_metrics(fits,z,valid_truth),regularization_sigma_logQ=sigma.tolist(),
        regularization='68th absolute logQ residual of every validation pair, including wrong times; floor.2; same-data penalty, not independent prior'))
    save(RECORD/('training-'+kind+'.json'),dict(passed=True,risk=kind,steps=protocol['steps'],seconds=time.monotonic()-started,
        model_sha256=digest(folder/'model.npz'),best_step=best['step'],validation_used_only=True,
        source_sha256=digest(PROJECT/'pileup/analysis/alpha_pair_risk_learning.py')))
    print('RISK_TRAINED',kind,best['step'],flush=True)


def freeze_training():
    check_freeze();selections={kind:load(EXPERIMENT/'models'/kind/'selection.json') for kind in RISKS}
    for kind in RISKS:assert load(RECORD/('training-'+kind+'.json'))['passed']
    chosen=max(RISKS,key=lambda kind:(selections[kind]['key'],-RISKS.index(kind)))
    files={str(p.relative_to(EXPERIMENT)):digest(p) for p in (EXPERIMENT/'models').rglob('*') if p.is_file()}
    save(EXPERIMENT/'training-freeze.json',dict(passed=True,chosen=chosen,selections=selections,files=files,
        frozen_at=datetime.now(timezone.utc).isoformat(),continuous_evaluation_generated=(OUTPUT/'evaluation').exists(),
        sigma=load(EXPERIMENT/'models'/chosen/'validation.json')['regularization_sigma_logQ']))
    assert not (OUTPUT/'evaluation').exists()
    print('RISK_WEIGHTS_FROZEN',chosen,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('prepare','generate','train','freeze-training'))
    parser.add_argument('--role',choices=('train','validation','evaluation'),default='train')
    parser.add_argument('--risk',choices=RISKS,default='logarithmic');parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--pilot',action='store_true');args=parser.parse_args()
    if args.command=='prepare':prepare()
    elif args.command=='generate':generate(args.role,args.workers,args.pilot)
    elif args.command=='train':train(args.risk)
    else:freeze_training()


if __name__=='__main__':main()
