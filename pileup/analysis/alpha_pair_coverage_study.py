"""Fresh, fixed-parameter Gaussian repetitions after numerical-rule freeze.

This is a small model-conditional coverage pilot, not calibrated confidence
intervals, new optical simulations, or a test of unknown-K event efficiency.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
import json
import multiprocessing as mp
import sys
import time
import unittest
import numpy as np

from .alpha_pair_contour_study import PROJECT,EXPERIMENT,OUTPUT,RECORD,models,check_freeze
from .alpha_pair_stability_study import copy_chunks
from .alpha_pair_information_study import load,save,digest,rows,serializable
from .alpha_pair_evidence_profile import profile,fit_at_charge
from .alpha_pair_profile_stability import refine_stage,midpoint_grid,point_cost,observed_starts,diverse_starts
from .alpha_pair_profile_contours import converge_contours


def freeze():
    check_freeze();gate=load(RECORD/'refinement-validation.json')
    assert gate['passed'] and gate['coverage_gate_passed'] and gate['stable']==5
    assert digest(OUTPUT/'refined-profiles.jsonl')==gate['sha256']
    folder=EXPERIMENT/'coverage-source';folder.mkdir(exist_ok=False)
    for name in ('alpha_pair_coverage_study.py','test_alpha_pair_model_sampler.py'):
        copy_chunks(PROJECT/'pileup/analysis'/name,folder/name)
    risk=PROJECT/'pileup/experiments/alpha_pair_risk_training_v2_20261027'
    for src,name in ((risk/'models/relative/model.npz','coverage-relative-model.npz'),(risk/'baseline-gate.json','coverage-baseline-gate.json')):
        copy_chunks(src,EXPERIMENT/name)
    settings=dict(initial_profile=dict(grid_count=33,starts=4,maxiter=120),dense=dict(starts=4,maxiter=150),
        contour=load(EXPERIMENT/'protocol.json')['settings'],hybrid='Old relative-MLP observed seed; unregularized physical fits from baseline/network,70 iterations each',
        controls_start=32,pilot_indices=[0,16],pair_indices=list(range(32)),all_indices=list(range(44)),
        scoring_probe=dict(starts=8,maxiter=300),workers_limit=8,blas_threads=1)
    save(EXPERIMENT/'coverage-inference-settings.json',settings)
    cases=[]
    for i in range(44):
        if i<32:
            gap=40. if i<16 else 200.;q=[18000.,800.];t=[200.3,200.3+gap];kind='pair';band=str(int(gap));condition='high_weak'
        else:
            group=(i-32)//4;gap=0.;q=[] if group==0 else [800. if group==1 else 18000.];t=[] if group==0 else [200.3]
            kind='empty' if group==0 else 'single';band=('empty','weak','strong')[group];condition=band
        cases.append(dict(index=i,seed_words=[202610271801,10000+i,994],charge_pe=q,times_ns=t,count=len(q),
            kind=kind,band=band,condition=condition,gap_ns=gap,role='independent_model_coverage_pilot',domain='fresh_gaussian'))
    save(EXPERIMENT/'coverage-design.truth.json',cases)
    protocol=dict(created=datetime.now(timezone.utc).isoformat(),frames=44,pairs=32,controls=12,pairs_per_gap=16,
        question='How often do fixed1/4/9 numerical contours include the fixed true later charge under independent repetitions of the frozen model?',
        sampling='16 independent draws at each gap40/200ns,Q=[18000,800]pe,t1=200.3ns. Four each empty,800pe single,18000pe single controls.',
        distribution='Exact frozen JointNoise Gaussian in FIR innovations; conditional prefix residual zero; full float64 voltage. No ADC clipping or optical parents.',
        seeds='NumPy default_rng SeedSequence[202610271801,10000+index,994], disjoint from previous simulation-design seeds. No reuse of a noise draw across gaps.',
        independence='New noise only after code,model,design and numeric rules freeze. Fixed previously chosen parameter points; same calibrated simulator/likelihood, not independent physical validation. Pilot first case per gap is retained unchanged.',
        profile='Initial33-log-grid plus observed-candidate nodes,4starts120iterations,forward/reverse; all midpoints4starts150; then frozen local adaptive contour routine. This complete finite algorithm, not previous curves, is what coverage tests.',
        evaluation='All32 knownK2 frames, including numerical failures, stay in denominator. Separately show conditional subset passing numerical screen, not as the primary coverage estimate.',
        inference_inputs='Observation, frozen response/noise/network, numeric settings and indices only. Truth/cases/ROOT blocked inside workers. TrueQ2 supplied only in downstream scoring probe after prediction hash freezes.',
        thresholds=[1.,4.,9.],threshold_selection='No threshold tuning from this pilot; no68/95/99% interpretation. Exact binomial95% intervals describe sampling uncertainty in observed inclusion fractions, not event charge intervals.',
        control_check='Original baseline unknownK selector and original native gate only; no recalibration.12controls are a smoke check, not a precise false-positive bound.',
        scope='Purealpha fixed95/351ns,knownK2 ordered charge,1024x2ns full voltage; local observed-trigger time domain and50-50000pe charge limits unchanged.',
        recovery='All complete-waveform profiles separately scored for Q/t accuracy,missing,extra,split,merge and sigma68; no claimed detection improvement.',
        numerical_gate_sha256=digest(RECORD/'refinement-validation.json'),primary_freeze_sha256=digest(EXPERIMENT/'freeze.json'),
        new_ROOTs=0,new_training=0,estimated_mc_precision='Only16 pairs per condition: empirical inclusion is preliminary and cannot validate nominal95% coverage.')
    save(EXPERIMENT/'coverage-protocol.json',protocol)
    (EXPERIMENT/'coverage.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for key in ('CONTOUR','STABILITY','EVIDENCE','GATE','RISK','INFORMATION'):os.environ['GAGG_'+key+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
import pileup.analysis
pileup.analysis.__path__.insert(0,str(folder/'coverage-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_coverage_study import main
    main()
''')
    names=['coverage.py','coverage-protocol.json','coverage-inference-settings.json','coverage-design.truth.json',
        'coverage-relative-model.npz','coverage-baseline-gate.json']
    files={n:digest(EXPERIMENT/n) for n in names}
    files.update({str(p.relative_to(EXPERIMENT)):digest(p) for p in folder.glob('*.py')})
    save(EXPERIMENT/'coverage-freeze.json',dict(created=protocol['created'],files=files))
    print('COVERAGE_FROZEN',len(files),flush=True)


def check_coverage():
    check_freeze()
    for n,h in load(EXPERIMENT/'coverage-freeze.json')['files'].items():assert digest(EXPERIMENT/n)==h,n


def tests():
    check_coverage();suite=unittest.defaultTestLoader.loadTestsFromName('pileup.analysis.test_alpha_pair_model_sampler')
    result=unittest.TextTestRunner(verbosity=2).run(suite);assert result.wasSuccessful()
    save(RECORD/'sampler-tests.json',dict(passed=True,tests=result.testsRun))


def generate():
    from .alpha_pair_model_sampler import GaussianFactor,lift_innovations
    check_coverage();assert load(RECORD/'sampler-tests.json')['passed']
    cases=load(EXPERIMENT/'coverage-design.truth.json');response,noise=models();folder=OUTPUT/'coverage';folder.mkdir(exist_ok=False)
    observations=[];charges=np.zeros((len(cases),2));times=charges.copy();invariants=[]
    for c in cases:
        q=c['charge_pe'];t=c['times_ns'];curves=response.components(t,q)
        factor=GaussianFactor.from_noise(noise,curves,q)
        rng=np.random.default_rng(np.random.SeedSequence(c['seed_words']))
        z=factor.draw(rng);residual=lift_innovations(z,noise.coefficients)
        y=curves.sum(0)+response.baseline+noise.background_mV+residual
        error=float(np.max(np.abs(noise.transform(residual)-z)));assert error<1e-9
        cost=noise.score(residual,curves,q);expect=float(np.sum(factor.whiten(z)**2)+factor.logdet)
        assert np.isclose(cost,expect,rtol=1e-10,atol=1e-8)
        observations.append(y);charges[c['index'],:len(q)]=q;times[c['index'],:len(t)]=t
        invariants.append(dict(index=c['index'],fir_roundtrip_max_error=error,quadratic_plus_logdet=expect,cost=float(cost)))
    y=np.asarray(observations,dtype=np.float64);assert y.shape==(44,1024) and np.isfinite(y).all()
    assert len({v.tobytes() for v in y})==44
    np.savez_compressed(folder/'observations.npz',voltage_mV=y)
    np.savez_compressed(folder/'labels.truth.npz',charge_pe=charges,times_ns=times,count=np.array([c['count'] for c in cases]))
    save(folder/'labels.cases.json',cases)
    save(RECORD/'coverage-generation.json',dict(passed=True,frames=44,source_freeze_sha256=digest(EXPERIMENT/'coverage-freeze.json'),
        created=datetime.now(timezone.utc).isoformat(),invariants=invariants,files={p.name:digest(p) for p in folder.iterdir()},new_ROOTs=0,new_noise_draws=44))
    print('COVERAGE_GENERATED',len(y),flush=True)


def initialize(probe=False):
    global RESPONSE,NOISE,NET,GATE,OBS,CONFIG,PRED
    from .alpha_pair_charge_allocation import forbid_labels
    from .alpha_pair_risk_learning import RiskNet
    RESPONSE,NOISE=models();NET=RiskNet.load(EXPERIMENT/'coverage-relative-model.npz')
    sys.addaudithook(forbid_labels)
    GATE=load(EXPERIMENT/'coverage-baseline-gate.json');CONFIG=load(EXPERIMENT/'coverage-inference-settings.json')
    with np.load(OUTPUT/'coverage/observations.npz') as f:OBS=f['voltage_mV']
    PRED={r['index']:r for r in rows(OUTPUT/'coverage/predictions.jsonl')} if probe else {}


def infer_one(index):
    from .alpha_pair_joint_likelihood import waveform_candidates,select
    from .alpha_pair_risk_learning import single_features,decode,regularized_refine
    started=time.monotonic();y=OBS[index];base=waveform_candidates(RESPONSE,NOISE,y);selected=select(base,GATE)
    if index>=CONFIG['controls_start']:
        return dict(index=index,control=True,baseline=base,baseline_selected=selected,seconds=time.monotonic()-started)
    if len(base)!=3:raise RuntimeError('Missing observed K2 seed bank; do not substitute truth')
    features,pilot=single_features(RESPONSE,NOISE,y)
    reference=decode(NET.forward(features[None,:]),[pilot['charge_pe'][0]],[pilot['times_ns'][0]])[0]
    # Only an initialization. No learned-Q penalty or truth assistance.
    initial=dict(reference,cost=point_cost(RESPONSE,NOISE,y,reference))
    candidates=[base[2],initial]+[regularized_refine(RESPONSE,NOISE,y,c,reference,np.ones(2),regularized=False,maxiter=70) for c in (base[2],reference)]
    best=min(candidates,key=lambda f:f['cost']);row=dict(sets=dict(baseline=base,mlp=base[:2]+[best]))
    first=profile(RESPONSE,NOISE,y,row,**CONFIG['initial_profile'])
    dense=refine_stage(RESPONSE,NOISE,y,row,first,midpoint_grid(first['q_grid_pe']),refit_existing=False,label='all_midpoints',**CONFIG['dense'])
    r=converge_contours(RESPONSE,NOISE,y,row,dense,CONFIG['contour'])
    final=r['stages'][r['final_stage']];error=0.
    for c in (first,*r['stages'].values()):
        for f in c['fits']:
            cost=point_cost(RESPONSE,NOISE,y,f);error=max(error,abs(cost-f['cost'])/max(1.,abs(cost)))
    assert error<1e-9
    return dict(index=index,control=False,row=row,baseline_selected=selected,coarse=first,**r,
        total_seconds=time.monotonic()-started,cost_replay_error=error,
        final_best=min(final['fits'],key=lambda f:f['cost']),primary_truth_inputs=[])


def probe_one(task):
    index,q2=task['index'],task['true_Q2_pe'];r=PRED[index];curve=r['stages'][r['final_stage']]
    seeds=observed_starts(r['row'],q2,curve['bounds_ns'],curve)
    seeds.extend(min(c['fits'],key=lambda f:f['cost']) for c in r['stages'].values())
    chosen=diverse_starts(RESPONSE,NOISE,OBS[index],q2,seeds,curve['bounds_ns'],CONFIG['scoring_probe']['starts'])
    fit=fit_at_charge(RESPONSE,NOISE,OBS[index],q2,chosen,curve['bounds_ns'],**CONFIG['scoring_probe'])
    actual=point_cost(RESPONSE,NOISE,OBS[index],fit);error=abs(actual-fit['cost'])/max(1.,abs(actual));assert error<1e-9
    return dict(**task,fit=fit,delta=fit['cost']-curve['reference_cost'],cost_replay_error=error,truth_assistance='Fixed true laterQ2 scoring only; not fed back to profiles')


def calculate(command,workers):
    check_coverage();assert 1<=workers<=8 and load(RECORD/'coverage-generation.json')['passed']
    for n,h in load(RECORD/'coverage-generation.json')['files'].items():assert digest(OUTPUT/'coverage'/n)==h
    config=load(EXPERIMENT/'coverage-inference-settings.json');prior=[];probe=command=='probe'
    tasks=config['pilot_indices'] if command=='pilot' else config['all_indices']
    if command=='infer':
        assert load(RECORD/'coverage-pilot-validation.json')['passed']
        prior=rows(OUTPUT/'coverage/pilot-predictions.jsonl');done={r['index'] for r in prior};tasks=[i for i in tasks if i not in done]
    if probe:
        validation=load(RECORD/'coverage-inference.json');assert validation['passed']
        assert digest(OUTPUT/'coverage/predictions.jsonl')==validation['sha256']
        with np.load(OUTPUT/'coverage/labels.truth.npz') as f:q=f['charge_pe']
        tasks=[dict(index=i,true_Q2_pe=float(q[i,1])) for i in config['pair_indices']]
    target=OUTPUT/'coverage'/('truth-probes.jsonl' if probe else 'pilot-predictions.jsonl' if command=='pilot' else 'predictions.jsonl')
    started=time.monotonic();result=list(prior)
    with target.open('x') as stream,ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize,initargs=(probe,)) as pool:
        for r in prior:stream.write(json.dumps(r,default=serializable,allow_nan=False)+'\n')
        futures=[pool.submit(probe_one if probe else infer_one,t) for t in tasks]
        for future in as_completed(futures):
            r=future.result();result.append(r);stream.write(json.dumps(r,default=serializable,allow_nan=False)+'\n');stream.flush()
            print('COVERAGE_'+command.upper(),len(result),len(prior)+len(tasks),r['index'],r.get('stable'),round(r.get('total_seconds',r.get('seconds',0.)),1),flush=True)
    expected=32 if probe else 2 if command=='pilot' else 44;assert len(result)==expected
    key='coverage-probe-validation' if probe else 'coverage-pilot-validation' if command=='pilot' else 'coverage-inference'
    save(RECORD/(key+'.json'),dict(passed=True,frames=len(result),workers=workers,seconds=time.monotonic()-started,sha256=digest(target),
        primary_sha256=digest(OUTPUT/'coverage/predictions.jsonl') if probe else None,primary_truth_inputs=[],label_files_blocked_in_workers=True,
        knownK2=True,stable=None if probe else sum(r.get('stable',False) for r in result),
        maximum_cost_replay_error=max(r.get('cost_replay_error',0.) for r in result),inherited_pilot=len(prior)))
    print('COVERAGE_'+command.upper()+'_COMPLETE',len(result),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('freeze','tests','generate','pilot','infer','probe'));p.add_argument('--workers',type=int,default=7)
    args=p.parse_args()
    if args.command=='freeze':freeze()
    elif args.command=='tests':tests()
    elif args.command=='generate':generate()
    else:calculate(args.command,args.workers)


if __name__=='__main__':main()
