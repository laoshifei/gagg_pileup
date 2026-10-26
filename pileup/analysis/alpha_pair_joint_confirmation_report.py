"""Independent confirmation metrics, mother-block intervals and full plots."""
from __future__ import annotations

import argparse
import csv
import json

import numpy as np

from .alpha_pair_joint_confirmation import PROJECT,EXPERIMENT,OUTPUT,RECORD,PARENT,CONDITIONS,check_freeze
from .alpha_pair_resolution_data import load,save,digest

LABELS=dict(observed_penalty='Old observed GLS',model_relative_residual='Old model GLS',
    model_matched_gate='Old model, matched gate',joint_conserved_nll='Conserved joint NLL')
COLORS=dict(observed_penalty='#cc3311',model_relative_residual='#0077bb',model_matched_gate='#888888',joint_conserved_nll='#009988')


def rows(path):return [json.loads(s) for s in path.read_text().splitlines()]


def validate():
    from .alpha_pair_joint_likelihood import select as joint_select
    from .alpha_pair_null_gate import select as old_select
    from .alpha_pair_residual import choose
    from .alpha_pair_conditional_study import evaluate_fit
    from .alpha_pair_null_report import summarize
    check_freeze();manifest=load(EXPERIMENT/'manifest.json');native=load(OUTPUT/'native-records.json')['test_single']
    assert len(manifest['runs'])==len(native)==len(manifest['planned_runs'])==320
    assert len({r['seed'] for r in native})==len({r['root'] for r in native})==320
    expected_config=load(PARENT/'historical-sensor-config.json')
    for r in manifest['runs']:
        assert r['status']=='completed' and r['returncode']==0
        assert r['validation']['runinfo_macro_verified'] and r['validation']['seed_verified']
        assert digest(PROJECT/r['root'])==r['root_sha256']
        assert digest(PROJECT/r['macro'])==r['macro_sha256']
    for r in native:assert load(PROJECT/r['metadata'])['config']==expected_config
    for r in manifest.get('failed_attempts',[]):
        assert digest(PROJECT/r['macro'])==r['macro_sha256']
        assert digest(PROJECT/r['log'])==r['log_sha256']
        assert (PROJECT/r['root']).exists()==r['root_exists']
        if r['root_exists']:assert (PROJECT/r['root']).stat().st_size==r['root_bytes']
    assert load(RECORD/'seed-independence.json')['passed']
    physical=load(RECORD/'physical-validation.json');assert physical['passed'] and physical['closure_max_mV']<1e-8
    assert digest(OUTPUT/'observations.npz')==physical['observation_sha256']
    assert load(RECORD/'physical-preflight.json')['scalar_digitizer_exact']
    dev=load(PROJECT/'jupyter/run_records'/PARENT.name/'data-provenance.json')
    old_mothers=set(dev['validation_mothers'])|set(dev['fold_mothers'][0])|set(dev['fold_mothers'][1])
    assert not old_mothers&{r['root'] for r in native}
    with np.load(OUTPUT/'labels.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
    cases=load(OUTPUT/'labels.cases.json');assert len(cases)==1600
    pred={name:{r['index']:r['candidates'] for r in rows(OUTPUT/f'predictions.{name}.jsonl')} for name in ('joint','observed','model')}
    for name,values in pred.items():
        assert len(values)==1600
        record=load(RECORD/f'{name}-inference.json');assert record['worker_labels_blocked'] and not record['truth_inputs']
        assert digest(OUTPUT/f'predictions.{name}.jsonl')==record['predictions_sha256']
    gates=load(EXPERIMENT/'gate-calibrations.json');legacy=load(EXPERIMENT/'legacy-gate-calibrations.json')['relative_residual']
    saved=rows(OUTPUT/'scores.jsonl');lookup={(r['method'],r['index']):r for r in saved}
    summary=[];known=[];convergence={}
    for method in LABELS:
        groups={};not_converged=0
        for i,c in enumerate(cases):
            if method=='observed_penalty':fit=choose(pred['observed'][i],.001)
            elif method=='model_relative_residual':fit=old_select(pred['model'][i],legacy)
            elif method=='model_matched_gate':fit=old_select(pred['model'][i],gates['old_ar'])
            else:fit=joint_select(pred['joint'][i],gates['joint_conserved_nll'])
            score=evaluate_fit(fit,t[i],q[i],c);stored=lookup[method,i]
            assert fit==stored['fit'] and all(v==stored[k] for k,v in score.items())
            groups.setdefault((c['kind'],c['condition'],c['gap_ns']),[]).append(stored)
            not_converged+=not fit.get('converged',True)
            if method=='joint_conserved_nll' and c['kind']=='pair' and c['condition']=='high_weak' and score['target']==2:
                known.append(dict(index=i,gap_ns=c['gap_ns'],blind_success=score['success'],**evaluate_fit(pred['joint'][i][2],t[i],q[i],c)))
        convergence[method]=not_converged;summary += [summarize(rr) for rr in groups.values()]
    assert summary==load(OUTPUT/'report.json')['summary']
    known_summary=[dict(gap_ns=gap,frames=sum(r['gap_ns']==gap for r in known),generated_frames=40,
        below_threshold_pairs=40-sum(r['gap_ns']==gap for r in known),
        blind_success=sum(r['blind_success'] for r in known if r['gap_ns']==gap),
        known_count_success=sum(r['success'] for r in known if r['gap_ns']==gap),
        matched_two=sum(r['matched']==2 for r in known if r['gap_ns']==gap)) for gap in (8.,40.,200.)]
    save(OUTPUT/'known-count-diagnosis.json',dict(summary=known_summary,rows=known,
        truth_used='Known K=2 only at scoring; existing blind double candidate, no new search'))
    # Common weights respect repeated use of the same optical mother block.
    block=np.full(len(cases),-1,int)
    for i,c in enumerate(cases):
        if c['kind']!='pair':continue
        blocks={native[j]['replicate']//2 for j in c['parents']};assert len(blocks)==1
        block[i]=next(iter(blocks))
    plan=load(EXPERIMENT/'report-plan.json');rng=np.random.default_rng(plan['bootstrap_seed']);B=plan['bootstrap_replicates']
    weights=np.zeros((B,40),int)
    for start in (0,20):
        draws=rng.integers(start,start+20,size=(B,20))
        for j in range(B):weights[j]+=np.bincount(draws[j],minlength=40)
    assert np.all(weights[:,:20].sum(1)==20) and np.all(weights[:,20:].sum(1)==20)
    bootstrap={};paired=[]
    for method in LABELS:
        estimates={True:[],False:[]}
        for condition in CONDITIONS:
            for gap in (8.,40.,200.):
                ids=[i for i,c in enumerate(cases) if c['kind']=='pair' and c['condition']==condition and c['gap_ns']==gap]
                denominator=np.bincount(block[ids],weights=[lookup[method,i]['target']==2 for i in ids],minlength=40)
                numerator=np.bincount(block[ids],weights=[lookup[method,i]['target']==2 and lookup[method,i]['success'] for i in ids],minlength=40)
                den=weights@denominator;assert np.all(den>0)
                estimates['weak' in condition].append((weights@numerator)/den)
                if method!='joint_conserved_nll':
                    changes=[int(lookup['joint_conserved_nll',i]['success'])-int(lookup[method,i]['success']) for i in ids if lookup[method,i]['target']==2]
                    paired.append(dict(baseline=method,condition=condition,gap_ns=gap,frames=len(changes),wins=sum(v>0 for v in changes),losses=sum(v<0 for v in changes),ties=sum(v==0 for v in changes)))
        bootstrap[method]={kind:np.mean(estimates[weak],axis=0) for weak,kind in [(True,'weak'),(False,'strong')]}
    report=load(OUTPUT/'report.json');intervals=[]
    for baseline in ('observed_penalty','model_relative_residual','model_matched_gate'):
        for kind in ('weak','strong'):
            delta=bootstrap['joint_conserved_nll'][kind]-bootstrap[baseline][kind]
            intervals.append(dict(baseline=baseline,group=kind,difference=report['results']['joint_conserved_nll'][kind+'_mean']-report['results'][baseline][kind+'_mean'],
                mother_block_bootstrap95=np.quantile(delta,[.025,.975]).tolist()))
    save(OUTPUT/'paired-confirmation.json',dict(comparisons=paired,aggregate_intervals=intervals,
        resampling=plan['paired_aggregate_uncertainty'],bootstrap_seed=plan['bootstrap_seed'],replicates=B,
        limits='Conditional on four energies, fixed response and this block design; not an uncertainty band for all detectors/rates'))
    with (OUTPUT/'summary.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(summary[0]));writer.writeheader();writer.writerows(summary)
    save(RECORD/'science-validation.json',dict(passed=True,completed_roots=320,interrupted_attempts=len(manifest.get('failed_attempts',[])),
        roots_macros_seeds_verified=True,scores_replayed=len(saved),summary_rows=len(summary),mothers_disjoint=True,
        selected_not_converged=convergence,gate_unchanged=True,freeze_unchanged=True,bootstrap_plan_sha256=digest(EXPERIMENT/'report-plan.json')))
    print('CONFIRMATION_VALIDATED',len(saved),'scores',intervals,flush=True)


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .alpha_pair_conditional import ConditionalResponse
    from .alpha_pair_joint_likelihood import ConservingResponse
    folder=OUTPUT/'figures';folder.mkdir(exist_ok=True);report=load(OUTPUT/'report.json');summary=report['summary']
    methods=['observed_penalty','model_relative_residual','joint_conserved_nll']
    fig,axes=plt.subplots(2,5,figsize=(18,8),sharex=True,sharey=True)
    for ax,condition,panel in zip(axes.flat,CONDITIONS,'abcdefghij'):
        for name in methods:
            rr=sorted([r for r in summary if r['method']==name and r['condition']==condition],key=lambda r:r['gap_ns'])
            x=[r['gap_ns'] for r in rr];y=np.array([100*r['efficiency'] for r in rr])
            ax.plot(x,y,'o-',color=COLORS[name],label=LABELS[name])
            if name=='joint_conserved_nll':
                bounds=100*np.array([r['success_wilson95'] for r in rr]).T
                ax.errorbar(x,y,yerr=np.maximum(np.vstack((y-bounds[0],bounds[1]-y)),0),fmt='none',color=COLORS[name],alpha=.4,capsize=3)
        ax.set(title=f'({panel}) {condition}',xscale='log',xticks=[8,40,200],xticklabels=['8','40','200'],ylim=(-3,103),xlabel='Separation (ns)');ax.grid(alpha=.2)
    axes[0,0].set_ylabel('Strict pair success (%)');axes[1,0].set_ylabel('Strict pair success (%)')
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3)
    fig.suptitle('Figure 0004-63-1 | Fresh optical parents; frozen methods; selected-method pointwise Wilson 95%')
    fig.tight_layout(rect=(0,.06,1,.95));fig.savefig(folder/'0004-63-1.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(16,5));conditions=['all','weak','low','mid','high']
    for j,name in enumerate(methods):
        rr=[next(r for r in summary if r['method']==name and r['kind']!='pair' and r['condition']==c) for c in conditions];x=np.arange(5)+(j-1)*.21
        axes[0].bar(x,[100*r['extra_fraction'] for r in rr],width=.2,color=COLORS[name],label=LABELS[name])
        axes[1].plot(x,[100*r['efficiency'] for r in rr],'o-',color=COLORS[name])
        axes[2].plot(x[1:],[100*r['q_error_p16_p50_p84'][1] for r in rr[1:]],'o-',color=COLORS[name])
    for ax,title,level,ylabel in zip(axes,['(a) Extra-output frames','(b) Strict success','(c) Matched charge bias'],[5,95,5],['Frames (%)','Success (%)','Median Q error (%)']):
        ax.axhline(level,color='k',ls='--',lw=.8);ax.set(title=title,xticks=np.arange(5),xticklabels=conditions,ylabel=ylabel);ax.grid(alpha=.2)
    axes[2].axhline(-5,color='k',ls='--',lw=.8)
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3)
    fig.suptitle('Figure 0004-63-2 | Independent empty/native-single controls; 80 frames per group')
    fig.tight_layout(rect=(0,.06,1,.95));fig.savefig(folder/'0004-63-2.png',dpi=150);plt.close(fig)
    old=ConditionalResponse.load(EXPERIMENT/'response.npz');new=ConservingResponse.load(EXPERIMENT/'joint-response.npz')
    with np.load(OUTPUT/'observations.npz') as f:y=f['voltage_mV']
    with np.load(OUTPUT/'labels.truth.npz') as f:components=f['components_mV']
    cases=load(OUTPUT/'labels.cases.json');scored={(r['method'],r['index']):r for r in rows(OUTPUT/'scores.jsonl')}
    fig,axes=plt.subplots(3,2,figsize=(16,11));examples=[]
    for col,gap in enumerate([40.,200.]):
        i=next(i for i,c in enumerate(cases) if c['condition']=='high_weak' and c['gap_ns']==gap and c['replicate']==0);grid=old.sample_grid
        axes[0,col].plot(grid,y[i]-old.baseline,'k',lw=.7,label='Observed total')
        axes[1,col].plot(grid,components[i,1],'k',lw=.7,label='True weak ancestry waveform')
        example=dict(index=i,case=cases[i],methods={})
        for method in methods:
            model=new if method=='joint_conserved_nll' else old;r=scored[method,i];fit=r['fit'];curves=model.components(fit['times_ns'],fit['charge_pe']);total=curves.sum(0)
            axes[0,col].plot(grid,total,color=COLORS[method],label=f'{LABELS[method]}, K={len(curves)}')
            if len(curves)==2:axes[1,col].plot(grid,curves[1],color=COLORS[method],label=LABELS[method]+'; later output')
            axes[2,col].plot(grid,y[i]-model.baseline-total,color=COLORS[method],alpha=.7,label=LABELS[method])
            example['methods'][method]=dict(fit=fit,success=r['success'],true_q=r['true_q'],true_t=r['true_t'])
        for k in range(3):
            ax=axes[k,col];ax.set(xlim=(0,2048),xlabel='Full local record time (ns)',ylabel='mV',title=f'({chr(97+2*k+col)}) high_weak, {gap:g} ns, replicate 0');ax.grid(alpha=.2);ax.legend(fontsize=8)
        examples.append(example)
    fig.suptitle('Figure 0004-64-1 | Predetermined fresh examples; blind fits, truth only as reference')
    fig.tight_layout();fig.savefig(folder/'0004-64-1.png',dpi=150);plt.close(fig)
    save(OUTPUT/'examples.json',examples)
    print('CONFIRMATION_FIGURES_COMPLETE',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['validate','figures']);args=parser.parse_args();globals()[args.command]()
