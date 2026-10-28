"""Separate contour numerics, fixed-model inclusion, and event accuracy."""
from __future__ import annotations

from collections import defaultdict
import numpy as np
from scipy.stats import binomtest

from .alpha_pair_contour_study import PROJECT,EXPERIMENT,OUTPUT,RECORD,INPUT,models,paths,check_freeze
from .alpha_pair_information_study import load,save,digest,rows
from .alpha_pair_evidence_profile import compatibility_segments
from .alpha_pair_profile_stability import compare_regions
from .alpha_pair_risk_report import evaluate,summarize,quantiles,csv_table


def region(curve,level):
    return compatibility_segments(curve['q_grid_pe'],[f['cost'] for f in curve['fits']],curve['reference_cost'],level)


def analyze():
    check_freeze();check=load(RECORD/'refinement-validation.json');assert check['passed']
    assert digest(OUTPUT/'refined-profiles.jsonl')==check['sha256']
    records=rows(OUTPUT/'refined-profiles.jsonl');truth={};diagnostic=[];diagnostic_scores=[]
    for domain in ('native','continuous'):
        with np.load(paths(domain)['truth']) as f:q=f['charge_pe'];t=f['times_ns']
        truth[domain]=dict(q=q,t=t,cases=load(paths(domain)['cases']))
    for r in records:
        d,i=r['domain'],r['index'];old=r['stages']['initial'];local=r['stages']['local_search'];final=r['stages'][r['final_stage']]
        identity=dict(domain=d,index=i);data=truth[d];case=data['cases'][i]
        for label,c in (('previous',old),('refined',final)):
            best=min(c['fits'],key=lambda f:f['cost'])
            diagnostic_scores.append(dict(evaluate(best,data['t'][i],data['q'][i],case),**identity,method=label,kind=case['kind'],condition=case['condition'],band=str(case.get('band',case['gap_ns']))))
        for level in (1.,4.,9.):
            a=compare_regions(old,local,level);b=compare_regions(local,final,level);c=compare_regions(old,final,level)
            diagnostic.append(dict(**identity,level=level,stable=r['stable'],points_before=len(old['q_grid_pe']),points_after=len(final['q_grid_pe']),
                final_max_bracket_fraction=r['final_quality']['maximum_bracket_fraction'],same_grid_endpoint_fraction=float(np.expm1(a['endpoint_log_change'])) if a['endpoint_log_change'] is not None else None,
                grid_endpoint_fraction=float(np.expm1(b['endpoint_log_change'])) if b['endpoint_log_change'] is not None else None,
                full_endpoint_fraction=float(np.expm1(c['endpoint_log_change'])) if c['endpoint_log_change'] is not None else None,
                same_grid_topology=a['same_topology'],grid_topology=b['same_topology'],
                reference_gain=old['reference_cost']-final['reference_cost'],same_grid_node_gain=r['history'][0]['old_node_max_improvement'],
                segments_before=a['old']['segments'],segments_final=c['new']['segments'],
                old_strict_success=diagnostic_scores[-2]['success'],final_strict_success=diagnostic_scores[-1]['success'],
                extra_trials=sum(h['trials'] for h in r['history']),nonconverged_trials=sum(h['nonconverged_trials'] for h in r['history']),seconds=r['seconds']))
    save(OUTPUT/'refinement-summary.json',diagnostic);csv_table(OUTPUT/'refinement-summary.csv',diagnostic)
    save(OUTPUT/'refinement-scores.json',diagnostic_scores)
    report=dict(refinement=diagnostic,refinement_scores=diagnostic_scores,refinement_frames=5,refinement_stable=check['stable'],
        knownK2=True,coverage_is_model_conditional=True,formal_nominal_CI_claim=False,no_default_change=True)
    if (RECORD/'coverage-inference.json').exists():
        from .alpha_pair_coverage_study import check_coverage
        check_coverage();validation=load(RECORD/'coverage-inference.json');probecheck=load(RECORD/'coverage-probe-validation.json')
        folder=OUTPUT/'coverage';assert digest(folder/'predictions.jsonl')==validation['sha256']==probecheck['primary_sha256']
        predictions=rows(folder/'predictions.jsonl');probes={r['index']:r for r in rows(folder/'truth-probes.jsonl')}
        assert digest(folder/'truth-probes.jsonl')==probecheck['sha256']
        cases=load(folder/'labels.cases.json')
        with np.load(folder/'labels.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
        scores=[];intervals=[];frame_diagnostics=[]
        for r in predictions:
            i=r['index'];case=cases[i];identity=dict(domain='fresh_gaussian',index=i,kind=case['kind'],condition=case['condition'],band=case['band'])
            scores.append(dict(evaluate(r['baseline_selected'],t[i],q[i],case),**identity,method='original_blind_baseline'))
            if r['control']:continue
            scores.append(dict(evaluate(r['final_best'],t[i],q[i],case),**identity,method='conditional_profile_K2'))
            c=r['stages'][r['final_stage']];probe=probes[i];delta=probe['fit']['cost']-c['reference_cost']
            interp=float(np.interp(np.log(q[i,1]),np.log(c['q_grid_pe']),[f['cost']-c['reference_cost'] for f in c['fits']]))
            frame_diagnostics.append(dict(index=i,band=case['band'],stable=r['stable'],probe_delta=delta,interpolated_delta=interp,
                interpolation_minus_probe=interp-delta,probe_below_reference=delta<-.01,probe_converged=probe['fit']['converged'],
                reference_cost=c['reference_cost'],width_domain=[50.,50000.],truth_outside_time_domain=bool(np.any(t[i]<c['bounds_ns'][0]) or np.any(t[i]>c['bounds_ns'][1])),
                final_Q2_error=r['final_best']['charge_pe'][1]/q[i,1]-1.,final_t2_error_ns=r['final_best']['times_ns'][1]-t[i,1],
                strict_success=scores[-1]['success'],baseline_selected_K=len(r['baseline_selected']['times_ns']),total_seconds=r['total_seconds'],
                optimizer_trials_initial_nominal=2*4*len(r['coarse']['q_grid_pe']),
                refinement_trials=sum(h['trials'] for h in r['history'])+r['stages']['initial']['trials'],
                nonconverged_trials=r['coarse']['nonconverged_trials']+r['stages']['initial']['nonconverged_trials']+sum(h['nonconverged_trials'] for h in r['history'])))
            for level in (1.,4.,9.):
                rg=region(c,level);included=any(lo<=q[i,1]<=hi for lo,hi in rg['segments'])
                intervals.append(dict(index=i,band=case['band'],stable=r['stable'],**rg,
                    components=len(rg['segments']),width_over_trueQ=(rg['upper']-rg['lower'])/q[i,1],truth_Q2_pe=float(q[i,1]),
                    grid_included=included,probe_included=delta<=level,probe_delta=delta,interpolated_delta=interp,
                    baseline_selected_K=len(r['baseline_selected']['times_ns'])))
        groups=defaultdict(list)
        for s in scores:groups[s['method'],s['kind'],s['band']].append(s)
        metrics=[summarize(rr) for rr in groups.values()];coverage=[]
        for band in ('40','200'):
            for level in (1.,4.,9.):
                rr=[r for r in intervals if r['band']==band and r['level']==level];n=len(rr);assert n==16
                k=sum(r['grid_included'] for r in rr);ci=binomtest(k,n).proportion_ci(confidence_level=.95,method='exact')
                stable=[r for r in rr if r['stable']];selected=[r for r in rr if r['baseline_selected_K']==2]
                coverage.append(dict(band=band,level=level,frames=n,included=k,inclusion_fraction=k/n,binomial95=[ci.low,ci.high],
                    probe_included=sum(r['probe_included'] for r in rr),inclusion_disagreements=sum(r['grid_included']!=r['probe_included'] for r in rr),
                    stable_frames=len(stable),stable_included=sum(r['grid_included'] for r in stable),
                    baseline_selected_K2_frames=len(selected),baseline_selected_K2_included=sum(r['grid_included'] for r in selected),
                    width_over_trueQ=quantiles([r['width_over_trueQ'] for r in rr]),disconnected=sum(r['components']>1 for r in rr),
                    lower_boundary=sum(r['lower_boundary'] for r in rr),upper_boundary=sum(r['upper_boundary'] for r in rr)))
        assert len(predictions)==44 and len(frame_diagnostics)==32 and len(intervals)==96 and len(scores)==76
        report.update(coverage=coverage,conditional_metrics=metrics,frames=frame_diagnostics,
            all32_retained=True,stable_frames=sum(r['stable'] for r in frame_diagnostics),
            probe_below_reference=sum(r['probe_below_reference'] for r in frame_diagnostics),
            failed_truth_probes=sum(not r['probe_converged'] for r in frame_diagnostics),
            inclusion_disagreements=sum(r['grid_included']!=r['probe_included'] for r in intervals),
            new_noise_draws=44,new_ROOTs=0,new_training=0,coverage_thresholds_not_tuned=True)
        for name,data in (('coverage-intervals.json',intervals),('coverage-scores.json',scores),('coverage-frames.json',frame_diagnostics)):
            save(OUTPUT/name,data)
        csv_table(OUTPUT/'coverage-summary.csv',coverage);csv_table(OUTPUT/'conditional-metrics.csv',metrics)
    save(OUTPUT/'report.json',report)
    names=['refinement-summary.json','refinement-summary.csv','refinement-scores.json','report.json']
    if 'coverage' in report:names+=['coverage-intervals.json','coverage-scores.json','coverage-frames.json','coverage-summary.csv','conditional-metrics.csv']
    save(RECORD/'science-validation.json',dict(passed=True,refinement_frames=5,coverage_frames=32 if 'coverage' in report else 0,
        coverage_controls=12 if 'coverage' in report else 0,primary_inclusion='Membership in computed interval set; actual fixed-trueQ probe reported separately',
        output_hashes={n:digest(OUTPUT/n) for n in names}))
    print('CONTOUR_COVERAGE_SCORED',report['refinement_stable'],report.get('stable_frames'),flush=True)


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.2})
    report=load(OUTPUT/'report.json');rr=rows(OUTPUT/'refined-profiles.jsonl');folder=OUTPUT/'figures';folder.mkdir(exist_ok=True);saved={}
    def finish(fig,number,title):
        for i,ax in enumerate(fig.axes):ax.set_title(f'({chr(97+i)}) '+ax.get_title())
        fig.suptitle('Figure '+number+' | '+title,fontsize=12);fig.tight_layout(rect=(0,0,1,.95))
        p=folder/(number+'.png');fig.savefig(p,dpi=145);plt.close(fig);saved[number]=digest(p)
    fig,axes=plt.subplots(2,2,figsize=(13,8))
    for row,index in enumerate((160,190)):
        r=next(r for r in rr if r['domain']=='continuous' and r['index']==index);final=r['stages'][r['final_stage']];ref=final['reference_cost']
        oldregion=region(r['stages']['initial'],9.);newregion=region(final,9.)
        oldends=np.asarray(oldregion['segments']).reshape(-1);newends=np.asarray(newregion['segments']).reshape(-1)
        center=float(newends[np.argmax(np.abs(np.log(newends/oldends)))]) if len(newends)==len(oldends) else float(newends[0])
        for s,color,label in (('initial','#777777','Previous grid'),('local_search','#ee7733','Same-grid refit'),(r['final_stage'],'#0077bb','Adaptive + verified')):
            c=r['stages'][s];q=c['q_grid_pe'];cost=np.asarray([f['cost'] for f in c['fits']])-ref
            for ax in axes[row]:ax.plot(q,cost,'.-',ms=3,lw=1,color=color,label=label)
        for ax in axes[row]:ax.axhline(9.,color='black',ls=':',lw=.8);ax.set(xscale='log',xlabel='Later Q2 (pe)',ylabel='J - final sampled minimum');ax.legend(fontsize=8)
        axes[row,0].set(xlim=(50,50000),ylim=(-.1,14.),title=f'Gaussian index {index}: full Q domain')
        axes[row,1].set(xlim=(center/1.3,center*1.3),ylim=(7.,11.),title=f'Index {index}: most shifted level-9 endpoint')
    finish(fig,'0004-108-1','Two previously changing contours; evaluated nodes versus interpolation')
    fig,axes=plt.subplots(1,2,figsize=(13,4.6))
    for r in sorted(rr,key=lambda r:(r['domain'],r['index'])):
        label=f"{r['domain']} {r['index']}";hist=r['history'];x=np.arange(len(hist))
        axes[0].plot(x,[100*h['quality']['maximum_bracket_fraction'] for h in hist],'o-',ms=3,label=label)
        axes[1].plot(x,[r['stages']['initial']['reference_cost']-r['stages'][h['stage']]['reference_cost'] for h in hist],'o-',ms=3,label=label)
    axes[0].axhline(1.,ls=':',color='black');axes[0].set(yscale='log',ylabel='Largest crossing bracket (%)',xlabel='Refinement stage (0 = same-grid refit)',title='All five selected diagnostics')
    axes[1].set(ylabel='Decrease in minimum J',xlabel='Refinement stage',title='Reference minimum changes separately');axes[0].legend(fontsize=8)
    finish(fig,'0004-109-1','Boundary resolution and optimum cost are distinct numerical checks')
    if 'coverage' in report:
        fig,axes=plt.subplots(1,3,figsize=(15,4.7));levels=(1.,4.,9.)
        for k,(band,color) in enumerate((('40','#ee7733'),('200','#0077bb'))):
            a=[next(v for v in report['coverage'] if v['band']==band and v['level']==l) for l in levels];pos=np.arange(3)+(k-.5)*.12
            mean=np.array([v['inclusion_fraction'] for v in a]);bounds=np.array([v['binomial95'] for v in a]).T
            axes[0].errorbar(pos,mean,yerr=np.vstack((mean-bounds[0],bounds[1]-mean)),fmt='o',capsize=4,color=color,label=f'{band} ns, n=16')
            for x,v in zip(pos,a):axes[0].annotate(f"{v['included']}/16",(x,v['inclusion_fraction']),xytext=(0,-15 if k==0 else 9),textcoords='offset points',ha='center',fontsize=8,color=color)
            w=np.array([v['width_over_trueQ'] for v in a]).T
            axes[1].errorbar(pos,w[1],yerr=np.vstack((w[1]-w[0],w[2]-w[1])),fmt='o-',capsize=3,color=color,label=f'{band} ns')
            fs=[f for f in report['frames'] if f['band']==band]
            for j,method in enumerate(('original_blind_baseline','conditional_profile_K2')):
                m=next(m for m in report['conditional_metrics'] if m['band']==band and m['method']==method and m['kind']=='pair')
                axes[2].bar(k+(j-.5)*.28,m['success'],width=.26,color=color,alpha=.45 if j==0 else 1.,label=method if k==0 else None)
                axes[2].text(k+(j-.5)*.28,m['success']+.3,str(m['success']),ha='center',fontsize=9)
        axes[0].set(xticks=range(3),xticklabels=['1','4','9'],xlabel='Delta J threshold',ylabel='Computed interval contains true Q2',ylim=(-.03,1.15),title='Exact binomial 95% sampling bars');axes[0].legend(fontsize=8,loc='lower right')
        axes[1].set(xticks=range(3),xticklabels=['1','4','9'],xlabel='Delta J threshold',ylabel='Total interval span / true Q2',title='Median and 16-84% spread');axes[1].legend(fontsize=8)
        axes[2].set(xticks=[0,1],xticklabels=['40 ns','200 ns'],ylabel='Strictly recovered pairs / 16',ylim=(0,17.5),title='Accuracy differs from inclusion');axes[2].legend(fontsize=7)
        finish(fig,'0004-110-1','Fresh fixed-model repetitions; not nominal 68/95/99% charge intervals')
        intervals=load(OUTPUT/'coverage-intervals.json');fig,axes=plt.subplots(1,2,figsize=(13,7.5))
        for ax,band in zip(axes,('40','200')):
            ids=sorted({v['index'] for v in intervals if v['band']==band})
            for j,index in enumerate(ids):
                for k,(level,color) in enumerate(zip(levels,('#009988','#0077bb','#ee7733'))):
                    v=next(v for v in intervals if v['index']==index and v['level']==level)
                    for lo,hi in v['segments']:ax.plot([lo/v['truth_Q2_pe'],hi/v['truth_Q2_pe']],[j+(k-1)*.19]*2,color=color,lw=1.8)
            ax.axvline(1.,color='black',ls=':',lw=1);ax.set(xscale='log',xlabel='Later Q2 / 800 pe (truth scale only)',ylabel='Prespecified frame index',yticks=range(len(ids)),yticklabels=ids,title=f'{band} ns: all 16 repetitions')
        from matplotlib.lines import Line2D
        axes[0].legend([Line2D([0],[0],color=c,lw=2) for c in ('#009988','#0077bb','#ee7733')],['Delta <= 1','Delta <= 4','Delta <= 9'],fontsize=8)
        finish(fig,'0004-111-1','Every disconnected interval retained; no exclusion based on numerical or fit outcome')
        pred={r['index']:r for r in rows(OUTPUT/'coverage/predictions.jsonl')};cases=load(OUTPUT/'coverage/labels.cases.json');response,noise=models()
        with np.load(OUTPUT/'coverage/observations.npz') as f:y=f['voltage_mV']
        fig,axes=plt.subplots(4,2,figsize=(14,13));example=[]
        for row,i in enumerate((0,1,16,17)):
            r=pred[i];c=cases[i];fit=r['final_best'];true=response.components(c['times_ns'],c['charge_pe']);fitted=response.components(fit['times_ns'],fit['charge_pe']);grid=response.sample_grid
            axes[row,0].plot(grid,y[i],color='.7',lw=.7,label='Observed');axes[row,0].plot(grid,fitted.sum(0)+response.baseline+noise.background_mV,color='#0077bb',lw=1,label='Conditional K2 fit')
            axes[row,1].plot(grid,true[1],color='black',lw=1,label='True-parameter conditional component');axes[row,1].plot(grid,fitted[1],color='#0077bb',lw=1,label=f"Fitted later Q2 = {fit['charge_pe'][1]:.0f} pe")
            for col in (0,1):axes[row,col].set(xlim=(0,2048),xlabel='Record time (ns)',ylabel='mV',title=f"Index {i}; gap {c['gap_ns']:g} ns; "+('full total voltage' if col==0 else 'later component, known K2'));axes[row,col].legend(fontsize=7)
            example.append(dict(index=i,true_times_ns=c['times_ns'],true_charge_pe=c['charge_pe'],fit=fit))
        finish(fig,'0004-112-1','Four prespecified fresh examples, complete record; black curves are conditional means')
        save(OUTPUT/'examples.json',example)
    save(RECORD/'figure-delivery.json',dict(passed=True,figures=saved,source_sha256=digest(__file__)))
    print('CONTOUR_COVERAGE_FIGURES',len(saved),flush=True)
