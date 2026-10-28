"""Observation-only application of frozen continuous-charge learners."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor,as_completed
import json
import multiprocessing as mp
import sys
import time

import numpy as np

from .alpha_pair_risk_study import PROJECT,EXPERIMENT,OUTPUT,RECORD,models,check_freeze
from .alpha_pair_information_study import load,save,digest,rows
from .alpha_pair_risk_learning import RISKS,RiskNet,decode,single_features,regularized_refine

METHODS=('baseline','nn_absolute','nn_logarithmic','nn_relative','mlp_nll','mlp_regularized')


def initialize():
    global RESPONSE,NOISE
    RESPONSE,NOISE=models()
    from .alpha_pair_charge_allocation import forbid_labels
    sys.addaudithook(forbid_labels)


def feature_one(task):
    i,y=task;started=time.monotonic();feature,pilot=single_features(RESPONSE,NOISE,y)
    return i,feature,pilot,time.monotonic()-started


def native_features(workers=4):
    check_freeze();started=time.monotonic()
    with np.load(OUTPUT/'native-observations.npz') as f:y=f['voltage_mV']
    results=[None]*len(y)
    with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize) as pool:
        futures=[pool.submit(feature_one,(i,v)) for i,v in enumerate(y)]
        for n,f in enumerate(as_completed(futures),1):
            i,feature,pilot,seconds=f.result();results[i]=(feature,pilot,seconds)
            if n%200==0:print('RISK_NATIVE_FEATURES',n,len(y),flush=True)
    np.savez_compressed(OUTPUT/'native-features.npz',features=np.array([r[0] for r in results]),
        one_Q_pe=np.array([r[1]['charge_pe'][0] for r in results]),one_time_ns=np.array([r[1]['times_ns'][0] for r in results]))
    save(OUTPUT/'native-pilots.json',[r[1] for r in results])
    save(RECORD/'native-features.json',dict(passed=True,frames=len(y),seconds=time.monotonic()-started,
        frame_seconds_p16_p50_p84=np.percentile([r[2] for r in results],[16,50,84]).tolist(),
        pilot_not_converged=sum(not r[1]['converged'] for r in results),truth_inputs=[],
        feature_sha256=digest(OUTPUT/'native-features.npz'),waveform_sha256=digest(OUTPUT/'native-observations.npz')))
    print('RISK_NATIVE_FEATURES_COMPLETE',len(y),flush=True)


def baseline_one(task):
    from .alpha_pair_joint_likelihood import waveform_candidates
    i,y=task;started=time.monotonic();candidates=waveform_candidates(RESPONSE,NOISE,y)
    assert all(np.isfinite(c['cost']) for c in candidates)
    return dict(index=i,candidates=candidates,seconds=time.monotonic()-started)


def baseline(workers=4):
    check_freeze();assert load(EXPERIMENT/'training-freeze.json')['passed']
    with np.load(OUTPUT/'evaluation/observations.npz') as f:y=f['voltage_mV']
    destination=OUTPUT/'evaluation/predictions.baseline.jsonl'
    previous=rows(destination) if destination.exists() else [];done={r['index'] for r in previous}
    assert len(done)==len(previous)
    pending=[i for i in range(len(y)) if i not in done];started=time.monotonic()
    with destination.open('a') as stream:
        with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize) as pool:
            futures=[pool.submit(baseline_one,(i,y[i])) for i in pending]
            for n,f in enumerate(as_completed(futures),1):
                stream.write(json.dumps(f.result(),allow_nan=False)+'\n');stream.flush()
                if n%100==0:print('RISK_CONTINUOUS_BASELINE',n,len(pending),flush=True)
    save(RECORD/'continuous-baseline.json',dict(passed=True,frames=len(y),new_frames=len(pending),seconds=time.monotonic()-started,
        prediction_sha256=digest(destination),truth_inputs=[],source_sha256=digest(__file__)))
    print('RISK_CONTINUOUS_BASELINE_COMPLETE',len(y),flush=True)


def hybrid_one(task):
    i,y,base,reference,sigma=task;result={};started=time.monotonic()
    if len(base)<3:
        for name in ('mlp_nll','mlp_regularized'):result[name]=dict(reference,seconds=0.,skipped_no_baseline_pair=True)
        return dict(index=i,hybrids=result,seconds=time.monotonic()-started)
    observed=np.asarray(y,dtype=float)-RESPONSE.baseline-NOISE.background_mV
    target=np.log(reference['charge_pe'])
    for name,regularized in (('mlp_nll',False),('mlp_regularized',True)):
        clock=time.monotonic();initial=[]
        for candidate in (base[2],reference):
            curves=RESPONSE.components(candidate['times_ns'],candidate['charge_pe'])
            cost=NOISE.score(observed-curves.sum(0),curves,candidate['charge_pe'],details=True)
            penalty=float(np.sum(((np.log(candidate['charge_pe'])-target)/sigma)**2)) if regularized else 0.
            initial.append(dict(candidate,**cost,penalty=penalty,penalized_cost=cost['cost']+penalty,network_score=reference['network_score']))
        trials=[regularized_refine(RESPONSE,NOISE,y,c,reference,sigma,regularized=regularized) for c in initial]
        best=min(initial+trials,key=lambda c:c['penalized_cost'])
        result[name]=dict(best,seconds=time.monotonic()-clock,trials=trials,
            retained_initial=best is initial[0] or best is initial[1],original_baseline_cost=base[2]['cost'])
        if not regularized:assert result[name]['cost']<=base[2]['cost']+1e-7
    return dict(index=i,hybrids=result,seconds=time.monotonic()-started)


def predict(domain,workers=4):
    check_freeze();frozen=load(EXPERIMENT/'training-freeze.json')
    for name,expected in frozen['files'].items():assert digest(EXPERIMENT/name)==expected,name
    native=domain=='native'
    feature_path=OUTPUT/('native-features.npz' if native else 'evaluation/features.npz')
    observation_path=OUTPUT/('native-observations.npz' if native else 'evaluation/observations.npz')
    baseline_path=OUTPUT/('native-predictions.baseline.jsonl' if native else 'evaluation/predictions.baseline.jsonl')
    with np.load(feature_path) as f:features={k:f[k] for k in f.files}
    with np.load(observation_path) as f:y=f['voltage_mV']
    bases={r['index']:r['candidates'] for r in rows(baseline_path)}
    assert set(bases)==set(range(len(y)))
    # No labels have been opened in this process. Block them before model application.
    from .alpha_pair_charge_allocation import forbid_labels
    sys.addaudithook(forbid_labels)
    candidates={};outputs={};network_seconds={}
    for kind in RISKS:
        model=RiskNet.load(EXPERIMENT/'models'/kind/'model.npz');clock=time.monotonic()
        output=model.forward(features['features']);assert np.isfinite(output).all()
        fits=decode(output,features['one_Q_pe'],features['one_time_ns']);network_seconds[kind]=time.monotonic()-clock
        outputs[kind]=output;candidates['nn_'+kind]=fits
    np.savez_compressed(OUTPUT/(domain+'-network-outputs.npz'),**outputs)
    chosen=candidates['nn_'+frozen['chosen']];sigma=np.asarray(frozen['sigma'])
    destination=OUTPUT/(domain+'-predictions.jsonl');previous=rows(destination) if destination.exists() else []
    done={r['index'] for r in previous};assert len(done)==len(previous)
    pending=[i for i in range(len(y)) if i not in done];started=time.monotonic()
    with destination.open('a') as stream:
        with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize) as pool:
            futures=[pool.submit(hybrid_one,(i,y[i],bases[i],chosen[i],sigma)) for i in pending]
            for n,f in enumerate(as_completed(futures),1):
                r=f.result();i=r['index'];r['variants']={name:fits[i] for name,fits in candidates.items()};r['variants'].update(r.pop('hybrids'))
                stream.write(json.dumps(r,allow_nan=False)+'\n');stream.flush()
                if n%100==0:print('RISK_HYBRID',domain,n,len(pending),flush=True)
    save(RECORD/(domain+'-inference.json'),dict(passed=True,frames=len(y),new_frames=len(pending),hybrid_wall_seconds=time.monotonic()-started,
        network_seconds=network_seconds,feature_sha256=digest(feature_path),waveform_sha256=digest(observation_path),
        prediction_sha256=digest(destination),network_outputs_sha256=digest(OUTPUT/(domain+'-network-outputs.npz')),
        truth_inputs=[],label_access_blocked=True,chosen_from_validation=frozen['chosen'],source_sha256=digest(__file__)))
    print('RISK_INFERENCE_COMPLETE',domain,len(y),flush=True)
