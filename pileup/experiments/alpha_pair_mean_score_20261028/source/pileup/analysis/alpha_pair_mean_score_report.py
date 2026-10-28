"""Known-count mean-score diagnostics; truth is opened only after inference locks."""
from __future__ import annotations

from collections import defaultdict
import json
import numpy as np

from .alpha_pair_mean_score_study import PROJECT,EXPERIMENT,OUTPUT,RECORD,INPUT,models,check_freeze
from .alpha_pair_mean_score import METHODS
from .alpha_pair_information_study import load,save,digest,rows,serializable
from .alpha_pair_risk_report import evaluate,summarize,quantiles,csv_table

LABELS=dict(baseline='Baseline NLL',wide_nll='Wider NLL search',frozen_gls='One frozen-C GLS',iterated_gls='Iterated GLS',
    known_times='True times fixed',known_charges='True charges fixed',truth_initialized_free='Truth-start free NLL',true_parameters='Truth / no fit')
COLORS=dict(zip(METHODS,('#555555','#0077bb','#ee7733','#009988')))


def advancement(scored,protocol):
    lookup={(r['method'],r['index']):r for r in scored}
    base=[r for r in scored if r['method']=='baseline'];result=[]
    pairs=[r for r in base if r['kind']=='pair'];singles=[r for r in base if r['kind']=='single']
    for method in protocol['candidate_methods']:
        subset=[lookup[method,r['index']] for r in pairs]
        changes={b:sum(int(lookup[method,r['index']]['success'])-int(r['success']) for r in pairs if r['band']==b) for b in ('pair40','pair200')}
        isolated_change=sum(int(lookup[method,r['index']]['success'])-int(r['success']) for r in singles)
        convergence=float(np.mean([r['converged'] for r in subset]))
        checks=dict(gain40=changes['pair40']>=protocol['minimum_40ns_strict_gain'],
            retain200=changes['pair200']>=-protocol['maximum_200ns_strict_loss'],
            convergence=convergence>=protocol['minimum_pair_convergence_fraction'],
            isolated=isolated_change>=-protocol['maximum_isolated_strict_loss'])
        result.append(dict(method=method,eligible=all(checks.values()),checks=checks,
            change40=changes['pair40'],change200=changes['pair200'],isolated_change=isolated_change,
            pair_convergence=convergence,pair_success=sum(r['success'] for r in subset),
            median_absolute_logQ2_error=float(np.median([abs(np.log(r['fit']['charge_pe'][1]/r['true_charge_pe'][1])) for r in subset]))))
    eligible=[r for r in result if r['eligible']]
    chosen=min(eligible,key=lambda r:(-r['pair_success'],r['median_absolute_logQ2_error'],METHODS.index(r['method'])))['method'] if eligible else None
    return dict(candidates=result,chosen=chosen,independent_confirmation_required=chosen is not None,
        no_default_change=True,scope='Predeclared development gate,not proof of improvement.')


def analyze():
    check_freeze()
    for label,name in (('inference','predictions.jsonl'),('oracle','oracle-diagnostics.jsonl')):
        v=load(RECORD/(label+'-validation.json'));assert v['passed'] and digest(OUTPUT/name)==v['sha256']
    predictions={r['index']:r for r in rows(OUTPUT/'predictions.jsonl')}
    oracle={r['index']:r for r in rows(OUTPUT/'oracle-diagnostics.jsonl')}
    assert load(RECORD/'oracle-validation.json')['primary_sha256']==digest(OUTPUT/'predictions.jsonl')
    cases={r['index']:r for r in load(INPUT/'test/labels.cases.json')}
    with np.load(INPUT/'test/labels.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
    scored=[];probed=[]
    for i,row in predictions.items():
        k=row['forced_K'];case=cases[i]
        identity=dict(index=i,domain='exposed_gaussian_known_count',kind=case['kind'],condition=case['condition'],band=case['band'],gap_ns=case['gap_ns'])
        assert row['truth_inputs']==[] and k==len(q[i,:k]) and k==case['count']
        for source,destination in ((row,scored),(oracle.get(i),probed)):
            if source is None:continue
            assert not source.get('outside_time_domain',False)
            for method,fit in source['fits'].items():
                qt=np.asarray(fit['charge_pe']);tt=np.asarray(fit['times_ns'])
                assert len(qt)==k
                destination.append(dict(evaluate(fit,t[i],q[i],case),**identity,method=method,fit=fit,
                    known_count=True,truth_assisted=source is not row,true_charge_pe=q[i,:k].tolist(),true_times_ns=t[i,:k].tolist(),
                    ordered_Q_error=(qt/q[i,:k]-1.).tolist(),ordered_t_error_ns=(tt-t[i,:k]).tolist(),
                    below500=bool(np.any(qt<500.)),charge_boundary=bool(np.any(qt<=50.0001) or np.any(qt>=49999.99))))
    assert len(scored)==176 and len(probed)==160
    groups=defaultdict(list)
    for r in scored+probed:groups[r['method'],r['band']].append(r)
    summary=[]
    for (method,band),rr in groups.items():
        r=summarize(rr);r.update(truth_assisted=rr[0]['truth_assisted'],known_count=True,
            below500=sum(v['below500'] for v in rr),charge_boundary=sum(v['charge_boundary'] for v in rr))
        r.update(Q2_error_p16_p50_p84=None,t2_error_p16_p50_p84=None,Q2_error_sigma68=None,both_Q_within20=None)
        if rr[0]['kind']=='pair':
            for name,key in (('Q2','ordered_Q_error'),('t2','ordered_t_error_ns')):
                r[name+'_error_p16_p50_p84']=quantiles([v[key][1] for v in rr])
            r['Q2_error_sigma68']=(r['Q2_error_p16_p50_p84'][2]-r['Q2_error_p16_p50_p84'][0])/2
            r['both_Q_within20']=sum(max(abs(v) for v in s['ordered_Q_error'])<=.2 for s in rr)
        summary.append(r)
    lookup={(r['method'],r['index']):r for r in scored+probed};paired=[];search=[]
    for band in ('pair40','pair200'):
        ids=[r['index'] for r in scored if r['method']=='baseline' and r['band']==band]
        for method in METHODS[1:]:
            delta=[int(lookup[method,i]['success'])-int(lookup['baseline',i]['success']) for i in ids]
            paired.append(dict(band=band,method=method,frames=len(ids),gained=delta.count(1),lost=delta.count(-1),unchanged=delta.count(0)))
        for i in ids:
            values={m:lookup[m,i] for m in ('baseline','wide_nll','truth_initialized_free','true_parameters')}
            search.append(dict(index=i,band=band,baseline_minus_wide_cost=values['baseline']['fit']['cost']-values['wide_nll']['fit']['cost'],
                wide_minus_truthstart_cost=values['wide_nll']['fit']['cost']-values['truth_initialized_free']['fit']['cost'],
                truth_minus_wide_cost=values['true_parameters']['fit']['cost']-values['wide_nll']['fit']['cost'],
                baseline_success=values['baseline']['success'],wide_success=values['wide_nll']['success'],truthstart_success=values['truth_initialized_free']['success']))
    decision=advancement(scored,load(EXPERIMENT/'protocol.json')['advancement'])
    report=dict(summary=summary,paired=paired,search=search,advancement=decision,
        primary_frames=44,pair_frames=32,oracle_frames=40,new_observations=0,
        known_count=True,exposed_development=True,no_default_change=True,no_nominal_coverage_claim=True)
    for name,table in (('scores.jsonl',scored),('oracle-scores.jsonl',probed)):
        with (OUTPUT/name).open('w') as f:
            for row in table:f.write(json.dumps(row,default=serializable,allow_nan=False)+'\n')
    save(OUTPUT/'report.json',report);csv_table(OUTPUT/'summary.csv',summary);csv_table(OUTPUT/'paired.csv',paired)
    save(EXPERIMENT/'advancement-decision.json',dict(**decision,report_sha256=digest(OUTPUT/'report.json'),primary_sha256=digest(OUTPUT/'predictions.jsonl')))
    assert all(r['baseline_minus_wide_cost']>=-1e-7 for r in search)
    save(RECORD/'science-validation.json',dict(passed=True,primary_scores=len(scored),oracle_scores=len(probed),groups=len(summary),
        all_frames_kept=True,known_count_explicit=True,oracle_separate=True,primary_sha256=digest(OUTPUT/'predictions.jsonl'),
        report_sha256=digest(OUTPUT/'report.json'),advancement=decision['chosen']))
    print('MEAN_SCORE_ANALYZED',json.dumps(decision),flush=True)


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    report=load(OUTPUT/'report.json');assert digest(OUTPUT/'report.json')==load(RECORD/'science-validation.json')['report_sha256']
    scored=rows(OUTPUT/'scores.jsonl');pred={r['index']:r for r in rows(OUTPUT/'predictions.jsonl')}
    folder=OUTPUT/'figures';folder.mkdir(exist_ok=True);saved={}
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.22})
    def finish(fig,number,title):
        for i,ax in enumerate(fig.axes):ax.set_title(f'({chr(97+i)}) '+ax.get_title(),fontsize=10)
        fig.suptitle(number+': '+title,fontsize=12);fig.tight_layout(rect=(0,0,1,.94))
        p=folder/(number+'.png');fig.savefig(p,dpi=130);plt.close(fig);saved[number]=digest(p)
    def group(method,band):return next(r for r in report['summary'] if r['method']==method and r['band']==band)
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    for ax,band in zip(axes[:2],('pair40','pair200')):
        for x,m in enumerate(METHODS):
            r=group(m,band);ax.bar(x,r['success'],color=COLORS[m]);ax.text(x,r['success']+.25,f"{r['success']}/16",ha='center')
        ax.set_xticks(range(4),[LABELS[m] for m in METHODS],rotation=25,ha='right');ax.set_ylim(0,18);ax.set_title(band[4:]+' ns / known K=2');ax.set_ylabel('Strict pair successes')
    ax=axes[2]
    for x,m in enumerate(METHODS):
        rr=[r for r in report['summary'] if r['method']==m and r['kind']=='single'];n=sum(r['success'] for r in rr)
        ax.bar(x,n,color=COLORS[m]);ax.text(x,n+.12,f'{n}/8',ha='center')
    ax.set_xticks(range(4),[LABELS[m] for m in METHODS],rotation=25,ha='right');ax.set_ylim(0,9);ax.set_title('Isolated controls / known K=1')
    finish(fig,'0004-120-1','Known-count development; same waveforms, no efficiency claim')
    fig,axes=plt.subplots(2,2,figsize=(11,8),sharex=True,sharey=True)
    for ax,m in zip(axes.flat,METHODS):
        for band,marker,color in (('pair40','o','#cc3311'),('pair200','^','#0077bb')):
            rr=[r for r in scored if r['method']==m and r['band']==band]
            ax.scatter([r['ordered_t_error_ns'][1] for r in rr],[100*r['ordered_Q_error'][1] for r in rr],label=band[4:]+' ns',marker=marker,color=color,s=36,alpha=.8)
        ax.axhspan(-20,20,color='green',alpha=.10);ax.axvspan(-10,10,color='green',alpha=.10)
        ax.set_xscale('symlog',linthresh=10);ax.set_yscale('symlog',linthresh=20)
        ax.set_xlabel('Ordered t2 error / ns');ax.set_ylabel('Ordered Q2 error / %');ax.set_title(LABELS[m]);ax.legend()
    finish(fig,'0004-121-1','All 32 pairs retained; shaded bands are marginal tolerances')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5));look={(r['method'],r['index']):r for r in scored}
    for band,color in (('pair40','#cc3311'),('pair200','#0077bb')):
        rr=[r for r in report['search'] if r['band']==band]
        axes[0].scatter([r['baseline_minus_wide_cost'] for r in rr],[100*look['wide_nll',r['index']]['ordered_Q_error'][1] for r in rr],color=color,label=band[4:]+' ns')
        axes[1].scatter([r['baseline_minus_wide_cost'] for r in rr],[r['wide_minus_truthstart_cost'] for r in rr],color=color,label=band[4:]+' ns')
    axes[0].axhspan(-20,20,color='green',alpha=.1);axes[0].set_yscale('symlog',linthresh=20);axes[0].set_ylabel('Wider NLL: ordered Q2 error / %');axes[0].set_title('Lower cost vs accurate charge')
    axes[1].axhline(0,color='black',lw=.8);axes[1].set_yscale('symlog',linthresh=1);axes[1].set_ylabel('J(wide) - J(truth-start free)');axes[1].set_title('Truth-start search diagnostic only')
    for ax in axes:ax.set_xlabel('J(baseline) - J(wider search)');ax.legend()
    finish(fig,'0004-122-1','Full likelihood cost improvements do not establish global optimality')
    methods=('baseline','wide_nll','known_times','known_charges','truth_initialized_free')
    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    for ax,band in zip(axes,('pair40','pair200')):
        for x,m in enumerate(methods):
            r=group(m,band);ax.bar(x,r['success'],color=COLORS.get(m,'#aa4499'));ax.text(x,r['success']+.25,str(r['success'])+'/16',ha='center')
        ax.set_xticks(range(len(methods)),[LABELS[m] for m in methods],rotation=25,ha='right');ax.set_ylim(0,18);ax.set_ylabel('Strict pair successes');ax.set_title(band[4:]+' ns / truth-assisted purple bars')
    finish(fig,'0004-123-1','Fixed truth probes and truth-start fit are diagnostics, not blind algorithms')
    response,noise=models()
    with np.load(INPUT/'test/observations.npz') as f:y=f['voltage_mV']
    fig,axes=plt.subplots(2,2,figsize=(13,7))
    for col,i in enumerate((320,384)):
        rr=next(r for r in scored if r['index']==i and r['method']=='baseline');t=response.sample_grid
        truth=response.components(rr['true_times_ns'],rr['true_charge_pe']);ax=axes[0,col]
        ax.plot(t,y[i]-response.baseline-noise.background_mV,color='#bbbbbb',lw=.7,label='Observed full waveform')
        for j,v in enumerate(truth):ax.plot(t,v,'--',color=('#aa4499','#cc3311')[j],label=f'True component {j+1}')
        for m in METHODS:
            fit=pred[i]['fits'][m];components=response.components(fit['times_ns'],fit['charge_pe'])
            ax.plot(t,components.sum(0),color=COLORS[m],lw=1,alpha=.8,label=LABELS[m])
            axes[1,col].plot(t,components[1],color=COLORS[m],lw=1,label=LABELS[m])
        axes[1,col].plot(t,truth[1],'--',color='#cc3311',label='True weak component')
        ax.set_title(f'Index {i}: complete {rr["gap_ns"]:g} ns pair');ax.set_ylabel('Voltage / mV');ax.legend(fontsize=7,ncol=2)
        axes[1,col].set_title('Weak component, full recorded time span');axes[1,col].set_ylabel('Voltage / mV');axes[1,col].legend(fontsize=7)
        for ax in axes[:,col]:ax.set_xlim(float(t[0]),float(t[-1]));ax.set_xlabel('Record time / ns')
    finish(fig,'0004-124-1','Fixed pilot examples, known K=2; no time-axis cropping')
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    for i,row in pred.items():
        if row['forced_K']!=2:continue
        h=row['iteration_history'];color='#cc3311' if i<384 else '#0077bb'
        for ax,key in zip(axes,('charge_log_change','time_change_ns')):
            ax.plot([r['iteration'] for r in h],[max(r[key],1e-8) for r in h],color=color,alpha=.4,lw=1,marker='.',markersize=4)
    for ax,cutoff,title in zip(axes,(.002,.02),('Undamped proposal: max |delta log Q|','Undamped proposal: max |delta t| / ns')):
        ax.axhline(cutoff,color='black',ls='--',label='Stop tolerance');ax.set_yscale('log');ax.set_xlabel('Outer GLS iteration');ax.set_title(title);ax.legend()
    finish(fig,'0004-125-1','IRGLS convergence: red 40 ns, blue 200 ns; all pairs retained')
    save(RECORD/'figure-manifest.json',dict(passed=True,figures=saved))
    print('MEAN_SCORE_FIGURES',len(saved),flush=True)


if __name__=='__main__':analyze();figures()
