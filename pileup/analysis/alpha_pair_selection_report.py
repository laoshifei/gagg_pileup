"""Selection denominators, paired gate comparisons, and full-waveform figures."""
from __future__ import annotations

from collections import defaultdict
import json
import numpy as np

from .alpha_pair_selection_study import (PROJECT, EXPERIMENT, OUTPUT, RECORD, OLD_OUT,
    METHODS, fraction, inclusion_partition, models, check_freeze)
from .alpha_pair_information_study import load, save, digest, rows, serializable
from .alpha_pair_joint_likelihood import select, evidence
from .alpha_pair_risk_report import evaluate, summarize, quantiles, csv_table

LABELS = {f'{origin}_{a:03d}':f'{"Original" if origin=="original" else "Matched"}\n{a/10:g}%'
          for origin in ('original','matched') for a in (25,50,100)}


def gates():
    v=load(RECORD/'gate-validation.json');assert v['passed']
    assert digest(EXPERIMENT/'gates.json')==v['sha256']
    return load(EXPERIMENT/'gates.json')


def selection_reason(candidates, selected):
    if len(selected['charge_pe'])==0:return 'no_presence'
    if len(candidates)<3:return 'no_K2_candidate'
    if min(candidates[2]['charge_pe'])<500.:return 'candidate_below500'
    if evidence(candidates)<=0.:return 'nonpositive_gain'
    return 'accepted' if len(selected['charge_pe'])==2 else 'rank_rejected'


def thresholds(allgates):
    result=[]
    for name,gate in allgates.items():
        for b,values in enumerate(gate['scores']):
            allowed=int(np.floor(gate['alpha']*(len(values)+1)-1+1e-12))
            cutoff=float(sorted(values,reverse=True)[allowed]) if 0<=allowed<len(values) else None
            result.append(dict(method=name,bin=b,calibration_n=len(values),alpha=gate['alpha'],
                minimum_rank=1/(len(values)+1),max_allowed_null_scores_ge_candidate=allowed,
                largest_grid_rank_at_alpha=max(0,allowed+1)/(len(values)+1),
                strict_evidence_cutoff=cutoff,also_requires='presence,bothQ>=500,gain>0'))
    return result


def old_audit(allgates):
    for name,h in load(RECORD/'input-provenance.json')['old_results'].items():assert digest(OLD_OUT/name)==h,name
    pred={r['index']:r for r in rows(OLD_OUT/'coverage/predictions.jsonl')}
    intervals=load(OLD_OUT/'coverage-intervals.json');details=[];summary=[]
    for band,indices in (('40',range(16)),('200',range(16,32))):
        for name,gate in allgates.items():
            decisions=[]
            for i in indices:
                base=pred[i]['row']['sets']['baseline'];s=select(base,gate)
                if name=='original_025':assert s==pred[i]['baseline_selected']
                decisions.append(len(s['charge_pe'])==2)
                details.append(dict(index=i,band=band,method=name,K2=len(s['charge_pe'])==2,
                    rank=s['null_rank_p'],evidence=evidence(base),raw_gain=base[1]['cost']-base[2]['cost'],
                    baseline_Q2=base[2]['charge_pe'][1],profile_Q2=pred[i]['final_best']['charge_pe'][1],
                    reason=selection_reason(base,s)))
            for level in (1.,4.,9.):
                rr=[next(r for r in intervals if r['index']==i and r['level']==level) for i in indices]
                parts=inclusion_partition([r['grid_included'] for r in rr],decisions)
                dd=[r for r in details if r['band']==band and r['method']==name]
                summary.append(dict(band=band,method=name,level=level,**parts,
                    selected_Q2_relative_error=quantiles([r['baseline_Q2']/800.-1 for r in dd if r['K2']]),
                    unselected_Q2_relative_error=quantiles([r['baseline_Q2']/800.-1 for r in dd if not r['K2']]),
                    selected_interval_span=quantiles([r['width_over_trueQ'] for r,s in zip(rr,decisions) if s]),
                    exploratory=True))
    return details,summary


def analyze():
    check_freeze();allgates=gates()
    validation=load(RECORD/'test-inference.json');assert validation['passed'] and validation['frames']==448
    folder=OUTPUT/'test';assert digest(folder/'predictions.jsonl')==validation['sha256']
    predicted={r['index']:r for r in rows(folder/'predictions.jsonl')};assert set(predicted)==set(range(448))
    cases=load(folder/'labels.cases.json')
    with np.load(folder/'labels.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
    scored=[];forced=[]
    for case in cases:
        i=case['index'];base=predicted[i]['candidates']
        identity=dict(index=i,domain='fresh_gaussian',kind=case['kind'],condition=case['condition'],band=case['band'],gap_ns=case['gap_ns'])
        for name,gate in allgates.items():
            fit=select(base,gate)
            scored.append(dict(evaluate(fit,t[i],q[i],case),**identity,method=name,fit=fit,reason=selection_reason(base,fit)))
        if case['kind']=='pair':
            assert len(base)==3
            forced.append(dict(evaluate(base[2],t[i],q[i],case),**identity,method='forced_K2_same_candidate',
                Q2_relative_error=base[2]['charge_pe'][1]/q[i,1]-1.,t2_error_ns=base[2]['times_ns'][1]-t[i,1],
                both_Q_ge500=min(base[2]['charge_pe'])>=500.,evidence=evidence(base),fit=base[2]))
    for name,rr in (('scores.jsonl',scored),('forced-K2.jsonl',forced)):
        with (OUTPUT/name).open('w') as stream:
            for r in rr:stream.write(json.dumps(r,default=serializable,allow_nan=False)+'\n')
    groups=defaultdict(list)
    for r in scored:groups[r['method'],r['band']].append(r)
    summaries=[];bias=[]
    for (method,band),rr in groups.items():
        s=summarize(rr)
        s.update(K2_fraction=fraction(s['K2_outputs'],len(rr)),success_fraction=fraction(s['success'],len(rr)),
                 extra_frame_fraction=fraction(s['extra_frames'],len(rr)),
                 reasons={k:sum(r['reason']==k for r in rr) for k in ('accepted','rank_rejected','candidate_below500','nonpositive_gain','no_presence','no_K2_candidate')})
        summaries.append(s)
        if rr[0]['kind']=='pair':
            ff={r['index']:r for r in forced if r['band']==band}
            for state in ('all','selected','unselected'):
                subset=[ff[r['index']] for r in rr if state=='all' or ((r['count']==2)==(state=='selected'))]
                bias.append(dict(method=method,band=band,state=state,frames=len(subset),
                    ordered_Q2_error_p16_p50_p84=quantiles([v['Q2_relative_error'] for v in subset]),
                    ordered_t2_error_p16_p50_p84=quantiles([v['t2_error_ns'] for v in subset]),
                    strict_success=sum(v['success'] for v in subset),interpretation='Ordered second candidate slot; includes time-mismatched candidates, not conditional matched-charge resolution.'))
    paired=[]
    for band in ('pair40','pair200'):
        a={r['index']:r for r in groups['original_025',band]};b={r['index']:r for r in groups['matched_025',band]}
        for outcome in ('selected_K2','success'):
            av={i:(r['count']==2 if outcome=='selected_K2' else r['success']) for i,r in a.items()}
            bv={i:(r['count']==2 if outcome=='selected_K2' else r['success']) for i,r in b.items()}
            paired.append(dict(band=band,outcome=outcome,both=sum(av[i] and bv[i] for i in a),
                original_only=sum(av[i] and not bv[i] for i in a),matched_only=sum(not av[i] and bv[i] for i in a),
                neither=sum(not av[i] and not bv[i] for i in a),frames=len(a)))
    caps=[dict(band=band,frames=64,success=sum(r['success'] for r in forced if r['band']==band),
        both_times=sum(r['strict_matched']==2 for r in forced if r['band']==band),
        below500=sum(not r['both_Q_ge500'] for r in forced if r['band']==band)) for band in ('pair40','pair200')]
    old_details,old_summary=old_audit(allgates)
    report=dict(summary=summaries,paired=paired,selection_bias=bias,forced_candidate_limit=caps,
        old_selection_summary=old_summary,thresholds=thresholds(allgates),calibration_frames=256,test_frames=448,
        no_new_intervals=True,no_nominal_coverage_claim=True,no_default_change=True)
    save(OUTPUT/'report.json',report);save(OUTPUT/'old-selection-frames.json',old_details)
    csv_table(OUTPUT/'summary.csv',summaries);csv_table(OUTPUT/'selection-bias.csv',bias);csv_table(OUTPUT/'thresholds.csv',report['thresholds'])
    # Selection at increasing alpha is monotone. Any new strict success must
    # already exist in the unchanged K2 candidate, never created by the gate.
    for origin in ('original','matched'):
        for i in range(448):
            ss=[next(r for r in scored if r['index']==i and r['method']==f'{origin}_{a:03d}') for a in (25,50,100)]
            assert [r['count'] for r in ss]==sorted(r['count'] for r in ss)
    for cap in caps:
        for s in summaries:
            if s['band']==cap['band']:assert s['success']<=cap['success']
    for s in old_summary:
        assert s['selected']['numerator']+s['unselected']['numerator']==s['all']['numerator']
        assert s['selected']['denominator']+s['unselected']['denominator']==16
    assert len(scored)==2688 and len(forced)==128 and len(summaries)==42
    save(RECORD/'science-validation.json',dict(passed=True,scored=len(scored),forced=len(forced),groups=len(summaries),
        old_interval_partitions=len(old_summary),predictions_sha256=validation['sha256'],report_sha256=digest(OUTPUT/'report.json'),
        original_old_decisions_reproduced=True,monotone_alpha=True,unchanged_candidate_limit=True,all_test_frames_kept=True))
    print('SELECTION_ANALYZED',len(scored),len(summaries),flush=True)


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    report=load(OUTPUT/'report.json');assert digest(OUTPUT/'report.json')==load(RECORD/'science-validation.json')['report_sha256']
    allgates=gates();folder=OUTPUT/'figures';folder.mkdir(exist_ok=True);saved={}
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.22})
    def finish(fig,number,title):
        for i,ax in enumerate(fig.axes):ax.set_title(f'({chr(97+i)}) '+ax.get_title(),fontsize=10)
        fig.suptitle(f'{number}: {title}',fontsize=12);fig.tight_layout(rect=(0,0,1,.95))
        path=folder/(number+'.png');fig.savefig(path,dpi=130);plt.close(fig);saved[number]=digest(path)
    colors=['#555555','#999999','#bbbbbb','#0077bb','#33bbee','#009988']
    old=load(OUTPUT/'old-selection-frames.json')
    fig,axes=plt.subplots(1,3,figsize=(15,4.8))
    for ax,band in zip(axes[:2],('40','200')):
        rr=[r for r in old if r['band']==band and r['method']=='original_025']
        for state,color,label in ((False,'#999999','Not selected'),(True,'#cc3311','Selected K2')):
            x=[r for r in rr if r['K2']==state]
            ax.scatter([r['evidence'] for r in x],[r['baseline_Q2']/800 for r in x],color=color,label=label)
        ax.axhspan(.8,1.2,color='#009988',alpha=.1);ax.axhline(1.,ls=':',color='black');ax.set(title=f'Old {band} ns: unchanged K2 candidates',xlabel='Guarded J1 - J2',ylabel='Ordered fitted Q2 / true Q2');ax.legend(fontsize=8)
    labels=('All intervals','Selected intervals','Selected and included')
    for j,band in enumerate(('40','200')):
        r=next(r for r in report['old_selection_summary'] if r['band']==band and r['method']=='original_025' and r['level']==4.)
        a=[r[k] for k in ('all','selected','selected_and_included')];x=np.arange(3)+(j-.5)*.30
        axes[2].bar(x,[v['fraction'] for v in a],width=.28,color=('#ee7733','#0077bb')[j],label=f'{band} ns')
        for xx,v in zip(x,a):axes[2].text(xx,v['fraction']+.025,f"{v['numerator']}/{v['denominator']}",ha='center',fontsize=8)
    axes[2].set(xticks=range(3),xticklabels=['P(C)','P(C | S)','P(C and S)'],ylim=(0,1.17),ylabel='Observed fraction',title='Old level-4 intervals: different denominators');axes[2].legend(fontsize=8,loc='upper center',ncol=2)
    finish(fig,'0004-114-1','Selection and interval inclusion in the 32 earlier diagnostic pairs')
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for b,ax in enumerate(axes.flat):
        for origin,color in (('original','#555555'),('matched','#0077bb')):
            values=np.sort(allgates[origin+'_025']['scores'][b]);unique,counts=np.unique(values,return_counts=True)
            margin=max(1.,float(np.ptp(values)))*.05
            ax.step(np.r_[unique[0]-margin,unique,unique[-1]+margin],np.r_[0.,np.cumsum(counts)/len(values),1.],where='post',color=color,label=f'{origin}, n={len(values)}')
            cutoff=next(r['strict_evidence_cutoff'] for r in report['thresholds'] if r['method']==origin+'_025' and r['bin']==b)
            if cutoff is not None:ax.axvline(cutoff,ls=':',color=color,lw=1)
        ax.set(xlabel='Guarded J1 - J2 on calibration singles',ylabel='Empirical cumulative fraction',ylim=(-.02,1.04),title=f'Observed single-fit charge bin {b}');ax.legend(fontsize=8)
    finish(fig,'0004-115-1','Single-pulse calibration; dotted lines are strict 2.5% rank cutoffs')
    fig,axes=plt.subplots(1,2,figsize=(14,5.2))
    for ax,field in zip(axes,('K2_outputs','success')):
        for j,(band,color) in enumerate((('pair40','#ee7733'),('pair200','#0077bb'))):
            rr=[next(r for r in report['summary'] if r['method']==m and r['band']==band) for m in METHODS];x=np.arange(6)+(j-.5)*.16
            vals=[fraction(r[field],64) for r in rr];means=np.array([v['fraction'] for v in vals]);bounds=np.array([v['binomial95'] for v in vals]).T
            ax.errorbar(x,means,yerr=[means-bounds[0],bounds[1]-means],fmt='o-',capsize=3,color=color,label=band.replace('pair','')+' ns')
            for xx,r in zip(x,rr):ax.annotate(str(r[field])+'/64',(xx,r[field]/64),xytext=(0,8 if r[field]<=2 or j else -14),textcoords='offset points',ha='center',fontsize=7)
        ax.set(xticks=range(6),xticklabels=[LABELS[m] for m in METHODS],ylim=(-.05,1.14),ylabel='Fraction of all 64 true pairs',title='Outputs two events' if field=='K2_outputs' else 'Strictly restores both events');ax.legend(fontsize=8)
    finish(fig,'0004-116-1','Independent confirmation: detection and strict accuracy with exact binomial bars')
    fig,axes=plt.subplots(1,2,figsize=(14,5))
    for ax,field in zip(axes,('K2_outputs','extra_frames')):
        bands=['empty','800','3000','8000','18800'];matrix=np.array([[next(r[field] for r in report['summary'] if r['method']==m and r['band']==b) for m in METHODS] for b in bands])
        ax.imshow(matrix,cmap='YlOrRd',vmin=0,vmax=max(1,int(matrix.max())),aspect='auto');ax.grid(False)
        for (i,j),v in np.ndenumerate(matrix):ax.text(j,i,str(v)+'/64',ha='center',va='center',fontsize=9,color='white' if v>.55*matrix.max() else 'black')
        ax.set(xticks=range(6),xticklabels=[LABELS[m] for m in METHODS],yticks=range(5),yticklabels=bands,ylabel='True control charge (pe); empty = 0',title='False two-event outputs' if field=='K2_outputs' else 'Frames with any time-unmatched output')
    finish(fig,'0004-116-2','Independent empty and single controls; two different false-output definitions')
    forced=rows(OUTPUT/'forced-K2.jsonl');scores=rows(OUTPUT/'scores.jsonl')
    fig,axes=plt.subplots(2,2,figsize=(13,8))
    for row,band in enumerate(('pair40','pair200')):
        ff=sorted([r for r in forced if r['band']==band],key=lambda r:r['index'])
        for col,method in enumerate(('original_025','matched_025')):
            ax=axes[row,col];selected={r['index']:r['count']==2 for r in scores if r['method']==method and r['band']==band}
            for state,color,label in ((False,'#aaaaaa','Not selected'),(True,'#0077bb','Selected K2')):
                rr=[r for r in ff if selected[r['index']]==state]
                ax.scatter([r['index'] for r in rr],[r['Q2_relative_error'] for r in rr],s=19,color=color,label=f'{label}, n={len(rr)}')
            ax.axhspan(-.2,.2,color='#009988',alpha=.1);ax.axhline(0.,color='black',ls=':',lw=1)
            ax.set_yscale('symlog',linthresh=.25)
            ax.set(xlabel='Prespecified confirmation index',ylabel='Ordered Q2 / true Q2 - 1',title=f'{band[4:]} ns, {method}; selected {sum(selected.values())}/64')
    fig.legend([Line2D([0],[0],marker='o',lw=0,color=c) for c in ('#aaaaaa','#0077bb')],['Not selected','Selected K2'],loc='upper center',bbox_to_anchor=(.5,.96),ncol=2,fontsize=9)
    finish(fig,'0004-117-1','Selection bias in all 128 fresh pairs; unmatched or failed candidates retained')
    response,noise=models();cases=load(OUTPUT/'test/labels.cases.json');pred={r['index']:r for r in rows(OUTPUT/'test/predictions.jsonl')}
    with np.load(OUTPUT/'test/observations.npz') as f:y=f['voltage_mV']
    fig,axes=plt.subplots(4,2,figsize=(14,13));examples=[]
    for row,index in enumerate(load(EXPERIMENT/'protocol.json')['example_test_indices']):
        case=cases[index];true=response.components(case['times_ns'],case['charge_pe']);x=response.sample_grid;record=dict(index=index,case=case,fits={})
        axes[row,0].plot(x,y[index],color='.65',lw=.6,label='Observed full voltage')
        axes[row,1].plot(x,true[1],color='black',lw=1,label='True conditional weak mean')
        for method,color in (('original_025','#ee7733'),('matched_025','#0077bb')):
            fit=select(pred[index]['candidates'],allgates[method]);record['fits'][method]=fit;c=response.components(fit['times_ns'],fit['charge_pe'])
            axes[row,0].plot(x,c.sum(0)+response.baseline+noise.background_mV,color=color,lw=1,label=f'{method}, K={len(c)}')
            if len(c)==2:axes[row,1].plot(x,c[1],color=color,lw=1,label=f'{method}: {fit["charge_pe"][1]:.0f} pe')
            else:axes[row,1].plot([],[],color=color,label=f'{method}: no second output')
        for col in (0,1):
            axes[row,col].set(xlim=(0,2048),xlabel='Record time (ns)',ylabel='mV',title=f'Index {index}, gap {case["gap_ns"]:g} ns: '+('total waveform' if col==0 else 'second output'));axes[row,col].legend(fontsize=7)
        examples.append(record)
    finish(fig,'0004-118-1','Four prespecified complete records; a gate changes acceptance, not fitted parameters')
    save(OUTPUT/'examples.json',examples);save(RECORD/'figure-delivery.json',dict(passed=True,figures=saved))
    print('SELECTION_FIGURES',len(saved),flush=True)


def main():
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('analyze','figures'));a=p.parse_args()
    (analyze if a.command=='analyze' else figures)()


if __name__=='__main__':main()
