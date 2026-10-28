"""Paired event recovery, risk diagnostics and figures for continuous-Q training."""
from __future__ import annotations

from collections import defaultdict
import csv
import json

import numpy as np

from .alpha_pair_risk_study import PROJECT,EXPERIMENT,OUTPUT,RECORD,models,check_freeze
from .alpha_pair_risk_inference import METHODS
from .alpha_pair_information_study import load,save,digest,rows
from .alpha_pair_risk_learning import RISKS,RiskNet,calibrated_gate,select_candidate
from .alpha_pair_joint_likelihood import select
from .alpha_pileup_fit import score

LABELS={'baseline':'Physical NLL','nn_absolute':'MLP / absolute Q','nn_logarithmic':'MLP / log Q',
        'nn_relative':'MLP / relative Q','mlp_nll':'MLP + NLL','mlp_regularized':'MLP + regularized NLL'}
COLORS=dict(zip(METHODS,('#555555','#0077bb','#33bbee','#ee7733','#009988','#cc3311')))


def quantiles(values):return np.percentile(values,[16,50,84]).tolist() if len(values) else None


def csv_table(path,table):
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(table[0]));writer.writeheader();writer.writerows(table)


def evaluate(fit,t,q,case):
    ids=np.flatnonzero(q>=500.);truth=dict(first_active_ns=t[ids],emission_ns=t[ids],charge_pe=q[ids])
    ordinary=score(fit,truth,[0.,2048.])
    tolerance=min(10.,case['gap_ns']/4.) if case['kind']=='pair' else 10.
    strict=score(fit,truth,[0.,2048.],tolerance_ns=tolerance)
    return dict(success=strict['recovered_count']==len(ids) and ordinary['unmatched_predictions']==0,
        target=len(ids),recovered=ordinary['recovered_count'],strict_recovered=strict['recovered_count'],
        matched=ordinary['matched_count'],strict_matched=strict['matched_count'],extra=ordinary['unmatched_predictions'],
        missed=ordinary['missed_count'],split=ordinary['split_flags'],merge=ordinary['merge_flags'],
        q_errors=ordinary['relative_charge_errors'],t_errors=ordinary['time_errors_ns'],count=len(fit['times_ns']),
        converged=fit.get('converged',True))


def summarize(rr):
    first=rr[0];q=[v for r in rr for v in r['q_errors']];t=[v for r in rr for v in r['t_errors']]
    qs=quantiles(q)
    return dict(domain=first['domain'],method=first['method'],kind=first['kind'],condition=first['condition'],band=first['band'],
        frames=len(rr),success=sum(r['success'] for r in rr),efficiency=float(np.mean([r['success'] for r in rr])),
        target=sum(r['target'] for r in rr),recovered=sum(r['recovered'] for r in rr),strict_recovered=sum(r['strict_recovered'] for r in rr),
        missed=sum(r['missed'] for r in rr),extra=sum(r['extra'] for r in rr),extra_frames=sum(r['extra']>0 for r in rr),
        extra_fraction=float(np.mean([r['extra']>0 for r in rr])),split=sum(r['split'] for r in rr),merge=sum(r['merge'] for r in rr),
        q_error_p16_p50_p84=qs,q_error_sigma68=(qs[2]-qs[0])/2. if qs else None,time_error_p16_p50_p84_ns=quantiles(t),
        K2_outputs=sum(r['count']==2 for r in rr),both_times_strict=sum(r['target']==2 and r['strict_matched']==2 for r in rr),
        not_converged=sum(not r['converged'] for r in rr))


def read_domain(domain):
    native=domain=='native'
    folder=OUTPUT if native else OUTPUT/'evaluation'
    prefix='native-' if native else ''
    cases=load(folder/(prefix+'labels.cases.json'))
    with np.load(folder/(prefix+'labels.truth.npz')) as f:q=f['charge_pe'];t=f['times_ns']
    base={r['index']:r['candidates'] for r in rows(folder/(prefix+'predictions.baseline.jsonl'))}
    predicted={r['index']:r for r in rows(OUTPUT/(domain+'-predictions.jsonl'))}
    return cases,q,t,base,predicted


def analyze():
    check_freeze();frozen=load(EXPERIMENT/'training-freeze.json')
    assert frozen['passed'] and not frozen['continuous_evaluation_generated']
    for name,expected in frozen['files'].items():assert digest(EXPERIMENT/name)==expected,name
    for name in ('numerical-tests','pilot-validation','train-generation','validation-generation','evaluation-generation',
                 'native-features','continuous-baseline','native-inference','continuous-inference'):
        assert load(RECORD/(name+'.json'))['passed'],name
    native=read_domain('native');continuous=read_domain('continuous')
    assert len(native[0])==1040 and len(continuous[0])==1000
    gates={'baseline':load(EXPERIMENT/'baseline-gate.json')}
    nq=[sum(native[3][i][1]['charge_pe']) for i in range(280)]
    for method in METHODS[1:]:gates[method]=calibrated_gate([native[4][i]['variants'][method] for i in range(280)],nq)
    save(EXPERIMENT/'gate-calibrations.json',gates)
    scored=[];forced=[];cost_replay_error=0.;network_replay_error=0.
    response,noise=models()
    for domain,data in (('native',native),('continuous',continuous)):
        cases,q,t,base,predicted=data
        expected=len(cases);assert set(predicted)==set(base)==set(range(expected))
        inference=load(RECORD/(domain+'-inference.json'))
        assert digest(OUTPUT/(domain+'-predictions.jsonl'))==inference['prediction_sha256']
        assert inference['truth_inputs']==[] and inference['label_access_blocked']
        feature_path=OUTPUT/('native-features.npz' if domain=='native' else 'evaluation/features.npz')
        wave_path=OUTPUT/('native-observations.npz' if domain=='native' else 'evaluation/observations.npz')
        with np.load(feature_path) as f:features=f['features']
        with np.load(wave_path) as f:y=f['voltage_mV']
        with np.load(OUTPUT/(domain+'-network-outputs.npz')) as f:
            for kind in RISKS:
                actual=RiskNet.load(EXPERIMENT/'models'/kind/'model.npz').forward(features)
                error=float(np.max(np.abs(actual-f[kind])));network_replay_error=max(network_replay_error,error)
                assert error<1e-5
        for i,case in enumerate(cases):
            for method in ('mlp_nll','mlp_regularized'):
                fit=predicted[i]['variants'][method]
                if fit.get('skipped_no_baseline_pair'):continue
                curves=response.components(fit['times_ns'],fit['charge_pe'])
                cost=noise.score(y[i]-response.baseline-noise.background_mV-curves.sum(0),curves,fit['charge_pe'])
                error=abs(cost-fit['cost'])/max(1.,abs(cost));cost_replay_error=max(cost_replay_error,error);assert error<1e-7
                if method=='mlp_nll':assert cost<=base[i][2]['cost']+1e-6
            if domain=='native' and i<280:continue
            band=str(int(case['gap_ns'])) if domain=='native' and case['kind']=='pair' else case.get('band')
            identity=dict(index=i,domain=domain,kind=case['kind'],condition=case['condition'],band=band,
                gap_ns=case['gap_ns'],block=case.get('block',i),native_index=case.get('native_index'))
            for method in METHODS:
                candidate=base[i][2] if method=='baseline' and len(base[i])==3 else predicted[i]['variants'][method] if method!='baseline' else None
                fit=select(base[i],gates[method]) if method=='baseline' else select_candidate(base[i],candidate,gates[method])
                if len(fit['times_ns'])<2:
                    original=base[i][len(fit['times_ns'])]
                    assert fit['times_ns']==original['times_ns'] and fit['charge_pe']==original['charge_pe']
                scored.append(dict(evaluate(fit,t[i],q[i],case),**identity,method=method,fit=fit))
                if case['kind']=='pair':
                    assert candidate is not None
                    cq=np.asarray(candidate['charge_pe']);ct=np.asarray(candidate['times_ns'])
                    forced.append(dict(evaluate(candidate,t[i],q[i],case),**identity,method=method,
                        relative_Q=(cq/q[i]-1.).tolist(),time_error_ns=(ct-t[i]).tolist(),
                        weak_Q_error=float(cq[np.argmin(q[i])]/q[i].min()-1.),
                        strong_Q_error=float(cq[np.argmax(q[i])]/q[i].max()-1.),
                        below500=bool(np.any(cq<500.)),total_Q_error=float(cq.sum()/q[i].sum()-1.)))
    assert len(scored)==10560 and len(forced)==5760
    old={r['index']:r for r in rows(PROJECT/'jupyter/artifacts/alpha_pair_charge_allocation_20261027/scores.jsonl')
         if r['domain']=='native' and r['method']=='baseline' and r['gate']=='primary'}
    for r in scored:
        if r['domain']=='native' and r['method']=='baseline':
            previous=old[r['index']]
            for key in ('success','target','recovered','strict_recovered','extra','missed','split','merge','q_errors','t_errors','count'):
                assert r[key]==previous[key],(r['index'],key)
    for name,table in (('scores.jsonl',scored),('forced-K2.jsonl',forced)):
        with (OUTPUT/name).open('w') as stream:
            for row in table:stream.write(json.dumps(row,allow_nan=False)+'\n')
    groups=defaultdict(list)
    for r in scored:groups[r['domain'],r['method'],r['kind'],r['condition'],r['band']].append(r)
    summary=[summarize(rr) for rr in groups.values()]
    save(OUTPUT/'summary.json',summary);csv_table(OUTPUT/'summary.csv',summary)
    aggregates=[]
    for domain in ('native','continuous'):
        for method in METHODS:
            rr=[r for r in scored if r['domain']==domain and r['method']==method and r['kind']=='pair']
            for part,subset in [('all',rr)]+([('targeted',[r for r in rr if r['condition']!='uniform_shift']),('prior_shift',[r for r in rr if r['condition']=='uniform_shift'])] if domain=='continuous' else []):
                aggregates.append(dict(domain=domain,method=method,part=part,frames=len(subset),success=sum(r['success'] for r in subset),
                    efficiency=float(np.mean([r['success'] for r in subset])),both_times_strict=sum(r['strict_matched']==2 for r in subset)))
    diagnostics=[];groups=defaultdict(list)
    for r in forced:groups[r['domain'],r['method'],r['condition'],r['band']].append(r)
    for (domain,method,condition,band),rr in groups.items():
        diagnostics.append(dict(domain=domain,method=method,condition=condition,band=band,frames=len(rr),
            success=sum(r['success'] for r in rr),efficiency=float(np.mean([r['success'] for r in rr])),
            below500=sum(r['below500'] for r in rr),both_times_strict=sum(r['strict_matched']==2 for r in rr),
            Q_error_p16_p50_p84=quantiles([v for r in rr for v in r['relative_Q']]),
            weak_Q_error_p16_p50_p84=quantiles([r['weak_Q_error'] for r in rr]),
            strong_Q_error_p16_p50_p84=quantiles([r['strong_Q_error'] for r in rr]),
            time_error_p16_p50_p84_ns=quantiles([v for r in rr for v in r['time_error_ns']]),
            total_Q_error_p16_p50_p84=quantiles([r['total_Q_error'] for r in rr]),not_converged=sum(not r['converged'] for r in rr)))
    lookup={(r['domain'],r['method'],r['index']):r for r in scored};rng=np.random.default_rng(270310099);paired=[]
    groups=defaultdict(list)
    for r in scored:
        if r['method']=='baseline' and r['kind']=='pair':groups[r['domain'],r['condition'],r['band']].append(r)
    for (domain,condition,band),rr in groups.items():
        if domain=='native':
            blocks=np.array([r['block'] for r in rr]);unique=np.unique(blocks)
            # A group contains either 20 shared blocks twice or 40 blocks once.
            draws=rng.choice(unique,(2000,len(unique)),replace=True)
            weights=np.stack([np.count_nonzero(draws==b,axis=1) for b in blocks],axis=1)
        else:
            draws=rng.integers(0,len(rr),(2000,len(rr)));weights=np.stack([np.count_nonzero(draws==i,axis=1) for i in range(len(rr))],axis=1)
        for method in METHODS[1:]:
            difference=np.array([float(lookup[domain,method,r['index']]['success'])-float(r['success']) for r in rr])
            samples=(weights@difference)/weights.sum(1)
            paired.append(dict(domain=domain,method=method,condition=condition,band=band,difference=float(difference.mean()),
                paired95=np.percentile(samples,[2.5,97.5]).tolist(),gained=int((difference>0).sum()),lost=int((difference<0).sum())))
    save(OUTPUT/'paired-comparison.json',dict(comparisons=paired,replicates=2000,seed=270310099,
        scope='Within-model independent Gaussian draws; native exposed optical-parent block resampling per condition. Not physical/model uncertainty; degenerate intervals not certainty.'))
    decisions=[]
    for method in METHODS[1:]:
        continuous_new=next(r for r in aggregates if (r['domain'],r['method'],r['part'])==('continuous',method,'all'))
        continuous_old=next(r for r in aggregates if (r['domain'],r['method'],r['part'])==('continuous','baseline','all'))
        native_new=next(r for r in aggregates if (r['domain'],r['method'],r['part'])==('native',method,'all'))
        native_old=next(r for r in aggregates if (r['domain'],r['method'],r['part'])==('native','baseline','all'))
        checks={}
        for domain in ('native','continuous'):
            controls=[r for r in summary if r['domain']==domain and r['method']==method and r['kind']!='pair']
            singles=[r for r in controls if r['kind']=='single']
            checks[domain+'_extra']=all(r['extra_fraction']<=.05 for r in controls)
            checks[domain+'_single_success']=all(r['efficiency']>=.95 for r in singles)
            checks[domain+'_single_bias']=all(r['q_error_p16_p50_p84'] and abs(r['q_error_p16_p50_p84'][1])<=.05 for r in singles)
        checks['continuous_pair_gain']=continuous_new['efficiency']>=continuous_old['efficiency']+.05-1e-12
        checks['native_pair_no_drop']=native_new['success']>=native_old['success']
        decisions.append(dict(method=method,passed=all(checks.values()),checks=checks))
    timing=[]
    for domain,data in (('native',native),('continuous',continuous)):
        predicted=data[4];ids=range(280,1040) if domain=='native' else range(1000)
        inference=load(RECORD/(domain+'-inference.json'))
        for kind,seconds in inference['network_seconds'].items():
            timing.append(dict(domain=domain,method='nn_'+kind,mean_incremental_seconds=seconds/len(data[0]),
                incremental_seconds_p16_p50_p84=None,scope='Batch network+decode only; excludes features, baseline and gate'))
        for method in ('mlp_nll','mlp_regularized'):
            values=[predicted[i]['variants'][method]['seconds'] for i in ids]
            timing.append(dict(domain=domain,method=method,mean_incremental_seconds=float(np.mean(values)),
                incremental_seconds_p16_p50_p84=quantiles(values),scope='Additional physical search only; excludes baseline, features, network and gate'))
    report=dict(summary=summary,aggregates=aggregates,diagnostics=diagnostics,paired=paired,decisions=decisions,timing=timing,
        selected_training_risk=frozen['chosen'],regularization_sigma=frozen['sigma'],no_default_replacement=True,
        new_ROOTs=0,new_physical_mothers=0,continuous_is_conditional_Gaussian=True,training_seeds=1)
    save(OUTPUT/'report.json',report);csv_table(OUTPUT/'pair-aggregates.csv',aggregates);csv_table(OUTPUT/'forced-K2-summary.csv',diagnostics)
    hashes={p.name:digest(p) for p in (OUTPUT/'scores.jsonl',OUTPUT/'forced-K2.jsonl',OUTPUT/'summary.json',OUTPUT/'report.json',OUTPUT/'paired-comparison.json')}
    save(RECORD/'science-validation.json',dict(passed=True,evaluation_frames=1760,gate_frames=280,primary_scores=10560,forced_K2_scores=5760,
        baseline_previous_exact=True,K0_K1_outputs_preserved=True,labels_blocked_in_inference=True,all_failed_fits_retained=True,
        maximum_network_replay_error=network_replay_error,maximum_cost_replay_relative_error=cost_replay_error,
        output_hashes=hashes,source_sha256=digest(__file__)))
    print('RISK_ANALYSIS_COMPLETE',frozen['chosen'],'promising',[r['method'] for r in decisions if r['passed']],flush=True)


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.2})
    assert load(RECORD/'science-validation.json')['passed'];report=load(OUTPUT/'report.json');summary=report['summary']
    folder=OUTPUT/'figures';folder.mkdir(exist_ok=True);hashes={}
    def finish(fig,number,title,legend_axis=None):
        fig.suptitle('Figure '+number+' | '+title,fontsize=12)
        if legend_axis is not None:
            handles,labels=legend_axis.get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,fontsize=8)
        fig.tight_layout(rect=(0,.13 if legend_axis is not None else .02,1,.93))
        path=folder/(number+'.png');fig.savefig(path,dpi=145);plt.close(fig);hashes[number]=digest(path)
    fig,axes=plt.subplots(1,3,figsize=(15,4.6))
    for ax,role in zip(axes,('train','validation','evaluation')):
        with np.load(OUTPUT/role/'labels.truth.npz') as f:q=f['charge_pe'];count=f['count']
        ids=np.flatnonzero(count==2);stride=max(1,len(ids)//2000)
        ax.scatter(q[ids[::stride],0],q[ids[::stride],1],s=5,alpha=.25)
        ax.set(xscale='log',yscale='log',xlabel='First Q (pe)',ylabel='Second Q (pe)',title=role)
    finish(fig,'0004-87-1','Continuous conditioning Q; new Gaussian draws, not native avalanche rescaling')
    fig,axes=plt.subplots(1,3,figsize=(15,4.6))
    for kind in RISKS:
        history=load(EXPERIMENT/'models'/kind/'history.json')
        for ax,key,label in zip(axes,('forced_pair_success','median_abs_logQ','count_BCE'),('Forced K=2 strict success','Median absolute log Q error','Count BCE')):
            ax.plot([r['step'] for r in history],[r[key] for r in history],'-o',ms=3,label=kind)
            ax.set(xlabel='Training steps',ylabel=label)
    finish(fig,'0004-88-1','Same initialization and minibatches; common validation metrics, not test selection',axes[0])
    fig,axes=plt.subplots(2,2,figsize=(13,8))
    for ax,condition in zip(axes.flat,('balanced','high_weak','weak_high','uniform_shift')):
        bands=('all',) if condition=='uniform_shift' else ('short','middle','long')
        for j,method in enumerate(METHODS):
            rr=[next(r for r in summary if (r['domain'],r['method'],r['condition'],r['band'])==('continuous',method,condition,b)) for b in bands]
            ax.bar(np.arange(len(bands))+(j-2.5)*.13,[100*r['efficiency'] for r in rr],width=.12,label=LABELS[method],color=COLORS[method])
        ax.set(xticks=np.arange(len(bands)),xticklabels=['4-320 ns'] if condition=='uniform_shift' else ['6-12 ns','24-64 ns','128-256 ns'],
            ylim=(0,103),ylabel='Strict pair recovery (%)',title=condition)
    finish(fig,'0004-89-1','Fresh continuous Gaussian evaluation; targeted 40/group and prior-shift 240 pairs',axes[0,0])
    fig,axes=plt.subplots(1,3,figsize=(15,4.8))
    for ax,condition in zip(axes,('mid_mid','high_weak','weak_high')):
        for j,method in enumerate(METHODS):
            rr=[next(r for r in summary if (r['domain'],r['method'],r['condition'],r['band'])==('native',method,condition,b)) for b in ('8','40','200')]
            ax.bar(np.arange(3)+(j-2.5)*.13,[100*r['efficiency'] for r in rr],width=.12,label=LABELS[method],color=COLORS[method])
        ax.set(xticks=np.arange(3),xticklabels=['8','40','200'],xlabel='True separation (ns)',ylim=(0,103),ylabel='Strict pair recovery (%)',title=condition)
    finish(fig,'0004-90-1','Exposed native transfer diagnostics, 40 pairs/group; not independent physical confirmation',axes[0])
    fig,axes=plt.subplots(2,3,figsize=(15,8))
    for row,domain in enumerate(('continuous','native')):
        for j,method in enumerate(METHODS):
            rr=[next(r for r in summary if r['domain']==domain and r['method']==method and r['kind']!='pair' and r['condition']==c) for c in ('all','weak','low','mid','high')]
            x=np.arange(5)+(j-2.5)*.12
            axes[row,0].bar(x,[100*r['extra_fraction'] for r in rr],width=.11,label=LABELS[method],color=COLORS[method])
            axes[row,1].plot(x,[100*r['efficiency'] for r in rr],'-o',ms=3,color=COLORS[method])
            axes[row,2].plot(x[1:],[100*r['q_error_p16_p50_p84'][1] for r in rr[1:]],'-o',ms=3,color=COLORS[method])
        for col,label in enumerate(('Extra-output frames (%)','Strict success (%)','Matched median Q error (%)')):
            axes[row,col].set(xticks=np.arange(5),xticklabels=['empty','weak','low','mid','high'],ylabel=label,title=domain)
        axes[row,0].axhline(5.,color='black',ls='--',lw=.8);axes[row,1].axhline(95.,color='black',ls='--',lw=.8)
        for value in (-5.,5.):axes[row,2].axhline(value,color='black',ls='--',lw=.8)
    finish(fig,'0004-90-2','Empty/single controls, 80/group; same native calibration pool and separate output gates',axes[0,0])
    cases,q,t,base,predicted=read_domain('native');response,noise=models();gates=load(EXPERIMENT/'gate-calibrations.json')
    with np.load(OUTPUT/'native-observations.npz') as f:y=f['voltage_mV']
    original=PROJECT/'jupyter/artifacts/alpha_pair_joint_confirmation_20261026/labels.truth.npz'
    with np.load(original) as f:actual=f['components_mV']
    fig,axes=plt.subplots(2,3,figsize=(16,8.6));examples=[]
    displayed=('baseline','nn_'+report['selected_training_risk'],'mlp_nll','mlp_regularized')
    for row,gap in enumerate((40.,200.)):
        i=next(i for i,c in enumerate(cases) if i>=280 and c['condition']=='high_weak' and c['gap_ns']==gap and c['replicate']==0)
        mean=response.components(t[i],q[i]);weak=int(np.argmin(q[i]));strong=1-weak
        axes[row,0].plot(response.sample_grid,y[i],color='#bbbbbb',lw=.7,label='Observed voltage')
        for col,index in ((1,strong),(2,weak)):
            axes[row,col].plot(response.sample_grid,actual[cases[i]['original_index'],index],color='black',lw=.7,label='Native actual component (truth)')
            axes[row,col].plot(response.sample_grid,mean[index],color='#999999',ls='--',lw=.8,label='Conditional true-parameter mean')
        example=dict(index=i,gap_ns=gap,methods={})
        for method in displayed:
            fit=select(base[i],gates[method]) if method=='baseline' else select_candidate(base[i],predicted[i]['variants'][method],gates[method])
            curves=response.components(fit['times_ns'],fit['charge_pe'])
            axes[row,0].plot(response.sample_grid,curves.sum(0)+response.baseline+noise.background_mV,color=COLORS[method],lw=.9,label=LABELS[method])
            larger=curves[np.argmax(fit['charge_pe'])] if len(curves) else np.zeros(1024)
            smaller=curves[np.argmin(fit['charge_pe'])] if len(curves)==2 else np.zeros(1024)
            axes[row,1].plot(response.sample_grid,larger,color=COLORS[method],lw=.9,label=LABELS[method])
            axes[row,2].plot(response.sample_grid,smaller,color=COLORS[method],lw=.9,label=LABELS[method]+f' (K={len(curves)})')
            example['methods'][method]=dict(fit=fit,score=evaluate(fit,t[i],q[i],cases[i]))
        for col,title in enumerate(('Total voltage','Larger-Q selected output','Smaller-Q selected output')):
            axes[row,col].set(xlim=(0,2048),xlabel='Record time (ns)',ylabel='mV',title=f'({chr(97+row*3+col)}) {gap:g} ns: {title}')
            axes[row,col].legend(fontsize=6)
        examples.append(example)
    finish(fig,'0004-91-1','Prespecified native pair0; full record, selected output only; absent second output shown as zero')
    save(OUTPUT/'examples.json',examples)
    save(RECORD/'plot-truth-provenance.json',dict(file=str(original.relative_to(PROJECT)),sha256=digest(original),purpose='Plot only; no inference access'))
    assert list(hashes)==load(EXPERIMENT/'protocol.json')['figures']
    save(RECORD/'figure-delivery.json',dict(passed=True,figures=hashes,source_sha256=digest(__file__)))
    print('RISK_FIGURES_COMPLETE',len(hashes),flush=True)
