"""Reproducible development run for a joint, normalized pulse likelihood."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime
import json
import multiprocessing as mp
import os
import shutil
import sys
import time

import numpy as np

from .alpha_pair_joint_data import PROJECT,EXPERIMENT,OUTPUT,RECORD
from .alpha_pair_resolution_data import load,save,digest

OLD=PROJECT/'jupyter/artifacts/alpha_pair_null_gate_20261023'
VARIANTS={'old_mean_nll':dict(mean='old',conserve=False,kappa=0.),
          'joint_conserved_nll':dict(mean='conserved',conserve=True,kappa=.95)}


def calibrate():
    from .alpha_pair_conditional import ConditionalResponse
    from .alpha_pair_joint_likelihood import ConservingResponse,JointNoise
    from .alpha_pair_joint_calibration import innovations,poles
    with np.load(OUTPUT/'calibration-ancestry.truth.npz') as f:
        q=f['charge_pe'];t=f['times_ns'];comp=f['components_mV'];fold=f['fold']
    with np.load(OUTPUT/'calibration.observations.npz') as f:y=f['voltage_mV']
    cfg=load(EXPERIMENT/'historical-sensor-config.json');coeff=poles(cfg);train=fold==0
    diagnostics={};background=y-comp.sum(1)-cfg['readout']['baseline_mV']
    zback=innovations(background,coeff)[train]
    electronic=np.correlate(coeff,coeff,'full')[len(coeff)-1:]*cfg['readout']['noise_rms_mV']**2
    dark=max(float(np.var(zback))-electronic[0],1e-5)
    noise_lags=electronic.copy();noise_lags[0]+=dark;noise_lags[1]+=.39*dark
    for name,options in VARIANTS.items():
        model=(ConservingResponse if options['conserve'] else ConditionalResponse).load(
            EXPERIMENT/('joint-response.npz' if options['conserve'] else 'response.npz'))
        mu=np.zeros_like(comp)
        for i in range(len(q)):
            ids=np.flatnonzero(q[i]>0);mu[i,ids]=model.components(t[i,ids],q[i,ids])
        r=innovations(comp-mu,coeff);m=innovations(mu,coeff)
        # Explicitly reproduce the exploratory conserved-mean array if present.
        path=OUTPUT/f'calibration-{options["mean"]}-residuals.truth.npz'
        if path.exists():
            with np.load(path) as f:
                np.testing.assert_allclose(mu,f['mean_mV'],rtol=2e-6,atol=5e-5)
                np.testing.assert_allclose(r,f['innovations'],rtol=2e-4,atol=2e-4)
        else:np.savez_compressed(path,mean_mV=mu,innovations=r.astype('f4'),innovations_mean=m.astype('f4'))
        activity=np.maximum(m[train],0.);residual=r[train];qt=q[train]
        ids=qt>0;totals=activity.sum(-1);neighbor=activity[:,::-1]
        x=np.stack((np.ones_like(qt),np.log(np.maximum(qt,50.)/3000.),
                    (activity*neighbor).sum(-1)/np.maximum(totals,1e-12)/10.),axis=-1)[ids]
        target=np.log(np.maximum((residual**2).sum(-1)[ids],1e-12)/totals[ids])
        params=np.linalg.lstsq(x,target,rcond=None)[0]
        params=np.clip(params,[-4.,-.5,-1.],[0.,.5,1.])
        v=activity*np.exp(np.clip(params[0]+params[1]*np.log(np.maximum(qt[:,:,None],50.)/3000.)+params[2]*neighbor/10.,-8.,4.))
        lag=float(np.clip(np.sum(residual[...,1:]*residual[...,:-1])/np.maximum(np.sum(np.sqrt(v[...,1:]*v[...,:-1])),1e-30),0.,.45))
        pairs=(qt>0).sum(1)==2
        cross=float(np.clip(np.sum(residual[pairs,0]*residual[pairs,1])/
                           np.maximum(np.sum(np.sqrt(v[pairs,0]*v[pairs,1])),1e-30),-.25,.25))
        noise=JointNoise(coeff,noise_lags,params,lag,cross,options['kappa'],float(background[train].mean()))
        noise.save(EXPERIMENT/f'{name}-noise.npz')
        groups=[]
        for lo,hi in [(500,1500),(1500,5000),(5000,12500),(12500,50000)]:
            keep=(qt>=lo)&(qt<hi);rr=residual[keep];vv=v[keep]
            denom=float(np.mean(vv.sum(1)+2*lag*np.sqrt(vv[:,1:]*vv[:,:-1]).sum(1)))
            groups.append(dict(q_range_pe=[lo,hi],components=int(keep.sum()),
                total_residual_median=float(np.median(rr.sum(1))),
                empirical_conditioning_fraction=1-float(np.var(rr.sum(1)))/denom))
        diagnostics[name]=dict(parameters=params.tolist(),lag1=lag,cross=cross,kappa=options['kappa'],
            noise_lags=noise_lags.tolist(),dark_innovation_variance=dark,background_mV=noise.background_mV,
            charge_groups=groups,heldout_component_rmse_mV=float(np.sqrt(np.mean((comp[fold==1]-mu[fold==1])**2))))
        print('JOINT_NOISE_CALIBRATED',name,diagnostics[name],flush=True)
    save(RECORD/'joint-noise-calibration.json',dict(models=diagnostics,training_fold=0,
        truth_used='Training means, residual moments and covariance only',source_sha256=digest(__file__)))


def prepare():
    if (EXPERIMENT/'development-freeze.json').exists():raise FileExistsError('Already frozen')
    cases=load(OUTPUT/'calibration.cases.json')
    gate_indices=[i for i,c in enumerate(cases) if c['fold']==1 and c['kind']=='single']
    assert len(gate_indices)==280
    with np.load(OUTPUT/'calibration.observations.npz') as f:gate_y=f['voltage_mV'][gate_indices]
    with np.load(OUTPUT/'calibration-ancestry.truth.npz') as f:gate_q=f['charge_pe'][gate_indices];gate_t=f['times_ns'][gate_indices]
    gate_cases=[dict(cases[i],role='calibration',pooled_index=i) for i in gate_indices]
    indices=load(EXPERIMENT/'validation-indices.json');vc=load(OLD/'test/labels.cases.json')
    with np.load(OLD/'test/observations.npz') as f:y=f['voltage_mV'][indices]
    with np.load(OLD/'test/labels.truth.npz') as f:q=f['charge_pe'][indices];t=f['times_ns'][indices]
    all_cases=gate_cases+[dict(vc[i],role='development',source_index=i) for i in indices]
    np.savez_compressed(OUTPUT/'observations.npz',voltage_mV=np.concatenate((gate_y,y)))
    np.savez_compressed(OUTPUT/'labels.truth.npz',charge_pe=np.concatenate((gate_q,q)),times_ns=np.concatenate((gate_t,t)))
    save(OUTPUT/'labels.cases.json',all_cases)
    # Identical voltage-only baseline predictions, with a newly matched gate split.
    dev=OLD/'model-development';dc=load(dev/'labels.cases.json')
    byseed={c['seed']:i for i,c in enumerate(dc) if c['kind']=='single'}
    old=[json.loads(s) for s in (dev/'predictions.model.jsonl').read_text().splitlines()]
    byindex={r['index']:r['candidates'] for r in old}
    preds=[byindex[byseed[c['seed']]] for c in gate_cases]
    with np.load(dev/'observations.npz') as f:
        np.testing.assert_array_equal(gate_y,f['voltage_mV'][[byseed[c['seed']] for c in gate_cases]])
    old={r['index']:r['candidates'] for r in map(json.loads,(OLD/'test/predictions.model.jsonl').read_text().splitlines())}
    preds += [old[i] for i in indices]
    with (OUTPUT/'predictions.old_ar.jsonl').open('w') as stream:
        for i,p in enumerate(preds):stream.write(json.dumps(dict(index=i,candidates=p))+'\n')
    plan=dict(created=datetime.now().astimezone().isoformat(),variants=VARIANTS,
        model='Gaussian surrogate after fixed electronics AR3 transform; first3 samples marginalized',
        covariance='Candidate shot activity times calibrated gain, nearest-neighbor correlation, joint cross correlation; electronics AR3 band plus dark band; optional rank2 charge conditioning',
        full_likelihood='r^T C^-1 r + logdet C; constant omitted; NLL may be negative',
        mean_area='Conserved variant: every full conditional component integrates to Q times13.618501389mVns/pe',
        conditioning='0.95 fixed conservative fraction; training moments of area-conserved components near0.995',
        gate='Raw NLL improvement, plus-one heldout native-single rank <=0.025; fitted-Q edges1500/5000/12500pe',
        selection=load(EXPERIMENT/'data-plan.json')['selection'],
        search='Raw voltage0/1/2 starts; legacy tangent residual4 starts; best4 NLL gap/share starts;20iteration screening then70iteration polish of best2',
        attribution='Mean, normalization, covariance, score and ranking change together; old_mean_nll is a partial ablation',
        training_mothers=280,gate_mothers=280,exposed_validation_mothers=320,
        inference_frames=1880,independent_test_not_started=True)
    save(EXPERIMENT/'method-plan.json',plan)
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for name in ('JOINT','CONDITIONAL','RESOLUTION','PAIR','COMPONENT'):
    os.environ['GAGG_'+name+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'development-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_joint_study import main
    main()
''')
    for namespace in ('analysis','pileup/analysis'):
        destination=EXPERIMENT/'development-source'/namespace;destination.mkdir(parents=True)
        for p in (PROJECT/namespace).glob('*.py'):shutil.copyfile(p,destination/p.name)
    if (PROJECT/'pileup/__init__.py').exists():shutil.copyfile(PROJECT/'pileup/__init__.py',EXPERIMENT/'development-source/pileup/__init__.py')
    files=[p for p in (EXPERIMENT/'development-source').rglob('*') if p.is_file()]
    files += [EXPERIMENT/n for n in ('response.npz','joint-response.npz','settings.json','method-plan.json',
        'historical-sensor-config.json','run.py',*[v+'-noise.npz' for v in VARIANTS])]
    save(EXPERIMENT/'development-freeze.json',dict(files={str(p.relative_to(EXPERIMENT)):digest(p) for p in files},
        independent_test_not_started=True,created=datetime.now().astimezone().isoformat()))
    save(RECORD/'inference-data-provenance.json',dict(gate_indices=gate_indices,validation_indices=indices,
        files={str(p.relative_to(PROJECT)):digest(p) for p in (OUTPUT/'observations.npz',OUTPUT/'labels.truth.npz',OUTPUT/'labels.cases.json',
            dev/'predictions.model.jsonl',OLD/'test/predictions.model.jsonl')}))
    print('JOINT_FROZEN',len(files),len(all_cases),flush=True)


def initialize(variant):
    global RESPONSE,NOISE
    from .alpha_pair_joint_likelihood import ConservingResponse,JointNoise
    from .alpha_pair_conditional import ConditionalResponse
    from .alpha_pair_conditional_study import forbid_labels
    options=VARIANTS[variant]
    RESPONSE=(ConservingResponse if options['conserve'] else ConditionalResponse).load(
        EXPERIMENT/('joint-response.npz' if options['conserve'] else 'response.npz'))
    NOISE=JointNoise.load(EXPERIMENT/f'{variant}-noise.npz')
    sys.addaudithook(forbid_labels)


def one(task):
    from .alpha_pair_joint_likelihood import waveform_candidates
    i,y=task;start=time.monotonic()
    return dict(index=i,candidates=waveform_candidates(RESPONSE,NOISE,y),seconds=time.monotonic()-start)


def infer(variant,workers=4,pilot=False):
    freeze=EXPERIMENT/os.environ.get('GAGG_JOINT_FREEZE','development-freeze.json')
    for name,sha in load(freeze)['files'].items():assert digest(EXPERIMENT/name)==sha,name
    with np.load(OUTPUT/'observations.npz') as f:y=f['voltage_mV']
    indices=list(range(len(y)))
    if pilot:
        cases=load(OUTPUT/'labels.cases.json')
        groups={}
        for i,c in enumerate(cases):
            key=(c['role'],c['kind'],c['condition'],c['gap_ns'])
            if key not in groups:groups[key]=i
        indices=list(groups.values())
    path=OUTPUT/f'predictions.{variant}.jsonl'
    previous=[json.loads(s) for s in path.read_text().splitlines()] if path.exists() else []
    done={r['index'] for r in previous};assert len(done)==len(previous)
    start=time.monotonic();count=0
    with path.open('a') as stream:
        with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize,initargs=(variant,)) as pool:
            futures=[pool.submit(one,(i,y[i])) for i in indices if i not in done]
            for future in as_completed(futures):
                result=future.result();stream.write(json.dumps(result,allow_nan=False)+'\n');stream.flush();count+=1
                if count%40==0:
                    os.fsync(stream.fileno());print('JOINT_INFER',variant,count,'/',len(futures),flush=True)
        os.fsync(stream.fileno())
    save(RECORD/(variant+('.pilot' if pilot else '')+'-inference.json'),dict(requested=len(indices),new_frames=count,
        seconds=time.monotonic()-start,workers=workers,truth_inputs=[],worker_labels_blocked=True,
        waveform_sha256=digest(OUTPUT/'observations.npz'),prediction_sha256=digest(path),freeze_file=freeze.name,freeze_sha256=digest(freeze)))
    print('JOINT_INFERENCE_COMPLETE',variant,len(indices),flush=True)


def select():
    from .alpha_pair_joint_likelihood import calibrate_gate,select as select_nll
    from .alpha_pair_null_gate import calibrate as old_calibrate,select as old_select
    from .alpha_pair_conditional_study import evaluate_fit
    from .alpha_pair_null_report import summarize
    cases=load(OUTPUT/'labels.cases.json')
    with np.load(OUTPUT/'labels.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
    previous=load(OLD/'test/report.json')['summary'];baselines={}
    averages=lambda rows,weak:float(np.mean([r['efficiency'] for r in rows if r['kind']=='pair' and ('weak' in r['condition'])==weak]))
    for method in ('observed_penalty','model_relative_residual'):
        rows=[r for r in previous if r['method']==method and (r['kind']!='pair' or r['gap_ns'] in (8.,40.,200.))]
        baselines[method]=dict(weak_mean=averages(rows,True),strong_mean=averages(rows,False))
    results={};calibrations={};all_summary=[];all_rows=[]
    for name in [*VARIANTS,'old_ar']:
        predictions={r['index']:r['candidates'] for r in map(json.loads,(OUTPUT/f'predictions.{name}.jsonl').read_text().splitlines())}
        assert len(predictions)==len(cases)==1880,name
        singles=[predictions[i] for i,c in enumerate(cases) if c['role']=='calibration']
        gate=(old_calibrate(singles,edges=[1500.,5000.,12500.],mode='relative_residual',alpha=.025,penalty=.001,require_birth_penalty=False)
              if name=='old_ar' else calibrate_gate(singles))
        calibrations[name]=gate;groups={}
        for i,c in enumerate(cases):
            if c['role']!='development':continue
            fit=(old_select if name=='old_ar' else select_nll)(predictions[i],gate)
            row=dict(evaluate_fit(fit,t[i],q[i],c),method=name,kind=c['kind'],condition=c['condition'],gap_ns=c['gap_ns'],
                index=i,fit=fit,true_q=q[i].tolist(),true_t=t[i].tolist())
            all_rows.append(row);groups.setdefault((c['kind'],c['condition'],c['gap_ns']),[]).append(row)
        summary=[summarize(rows) for rows in groups.values()];all_summary+=summary
        controls=[r for r in summary if r['kind']!='pair']
        control_pass=all(r['extra_fraction']<=.05 for r in controls) and all(r['efficiency']>=.95 and abs(r['q_error_p16_p50_p84'][1])<=.05 for r in controls if r['kind']=='single')
        weak=averages(summary,True);strong=averages(summary,False)
        eligible=bool(name in VARIANTS and control_pass and weak>max(b['weak_mean'] for b in baselines.values()) and strong>=baselines['observed_penalty']['strong_mean']-.02)
        results[name]=dict(controls_pass=control_pass,weak_mean=weak,strong_mean=strong,eligible=eligible,
            controls={r['condition']:dict(extra=r['extra_frames'],success=r['success'],frames=r['frames'],median_q_bias=None if r['q_error_p16_p50_p84'] is None else r['q_error_p16_p50_p84'][1]) for r in controls})
        print('JOINT_SELECTION',name,json.dumps(results[name]),flush=True)
    eligible=[name for name in VARIANTS if results[name]['eligible']]
    chosen=max(eligible,key=lambda name:results[name]['weak_mean']+results[name]['strong_mean']) if eligible else None
    save(EXPERIMENT/'development-selection.json',dict(chosen=chosen,results=results,baselines=baselines,independent_test_not_started=True))
    save(EXPERIMENT/'gate-calibrations.json',calibrations)
    save(OUTPUT/'development-report.json',dict(chosen=chosen,results=results,baselines=baselines,summary=all_summary,
        independent_performance=False,paired_parent_reuse_across_conditions=True))
    with (OUTPUT/'scored-development.jsonl').open('w') as stream:
        for row in all_rows:stream.write(json.dumps(row)+'\n')
    print('JOINT_CHOSEN',chosen,flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['calibrate','prepare','infer','select'])
    parser.add_argument('--variant',choices=list(VARIANTS),default='joint_conserved_nll');parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--pilot',action='store_true');args=parser.parse_args()
    if args.command=='infer':infer(args.variant,args.workers,args.pilot)
    else:globals()[args.command]()


if __name__=='__main__':main()
