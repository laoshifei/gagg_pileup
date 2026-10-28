"""Threshold sweep and conditional profile scoring; truth is used only here."""
from __future__ import annotations

from collections import defaultdict
import json

import numpy as np

from .alpha_pair_evidence_study import PROJECT,EXPERIMENT,OUTPUT,RECORD,INPUT,paths,models,check_freeze
from .alpha_pair_evidence_profile import decision_endpoints,choose_endpoint,compatibility_segments
from .alpha_pair_physical_gate import decide
from .alpha_pair_information_study import load,save,digest,rows
from .alpha_pair_risk_report import evaluate,summarize,quantiles,csv_table

METHODS=('baseline','learned_gate','local_physical','mlp_physical')
LABELS=dict(baseline='Original physical NLL',learned_gate='ML seed / learned gate',local_physical='Local / physical gate',mlp_physical='ML seed / physical gate')
COLORS=dict(zip(METHODS,('#555555','#ee7733','#0077bb','#009988')))


def load_data():
    data={}
    for domain in ('native','continuous'):
        p=paths(domain)
        with np.load(p['truth']) as f:q=f['charge_pe'];t=f['times_ns']
        data[domain]=dict(cases=load(p['cases']),q=q,t=t,cached={r['index']:r for r in rows(p['candidates'])})
    return data


def sweep():
    check_freeze();gates=load(EXPERIMENT/'gate-calibrations.json');protocol=load(EXPERIMENT/'protocol.json')
    for name,value in load(RECORD/'input-provenance.json')['inputs'].items():assert digest(PROJECT/name)==value,name
    data=load_data();records=[];old={(r['domain'],r['index'],r['method']):r for r in rows(INPUT/'scores.jsonl')}
    for domain,d in data.items():
        for i,case in enumerate(d['cases']):
            if domain=='native' and i<280:continue
            band=str(int(case['gap_ns'])) if domain=='native' and case['kind']=='pair' else case.get('band')
            for method in METHODS:
                end=decision_endpoints(d['cached'][i],method,gates[method])
                scored={side:evaluate(end[side],d['t'][i],d['q'][i],case) for side in ('low','high')}
                chosen=scored[choose_endpoint(end,.025)]
                for key in ('success','target','recovered','strict_recovered','matched','strict_matched','extra','missed','split','merge','q_errors','t_errors','count'):
                    assert chosen[key]==old[domain,i,method][key],(domain,i,method,key)
                records.append(dict(domain=domain,index=i,method=method,kind=case['kind'],condition=case['condition'],band=band,
                    gap_ns=case['gap_ns'],block=case.get('block',i),endpoints=end,metrics=scored))
    assert len(records)==7040
    with (OUTPUT/'decision-records.jsonl').open('w') as f:
        for row in records:f.write(json.dumps(row,allow_nan=False)+'\n')
    alphas=protocol['threshold']['alpha_grid'];all_summaries=[];aggregate=[];at_default=[]
    grouped=defaultdict(list)
    for r in records:grouped[r['domain'],r['method'],r['kind'],r['condition'],r['band']].append(r)
    for alpha in alphas:
        summary=[]
        for key,rr in grouped.items():
            selected=[dict(r['metrics'][choose_endpoint(r['endpoints'],alpha)],domain=r['domain'],method=r['method'],kind=r['kind'],condition=r['condition'],band=r['band']) for r in rr]
            summary.append(dict(summarize(selected),alpha=alpha))
        all_summaries.extend(summary)
        if alpha==.025:at_default=summary
        for domain in ('native','continuous'):
            for method in METHODS:
                pp=[r for r in summary if r['domain']==domain and r['method']==method and r['kind']=='pair']
                cc=[r for r in summary if r['domain']==domain and r['method']==method and r['kind']!='pair']
                target='40' if domain=='native' else 'middle';long='200' if domain=='native' else 'long'
                short_result=next(r for r in pp if r['condition']=='high_weak' and r['band']==target)
                long_result=next(r for r in pp if r['condition']=='high_weak' and r['band']==long)
                aggregate.append(dict(domain=domain,method=method,alpha=alpha,frames=sum(r['frames'] for r in pp),
                    strict_success=sum(r['success'] for r in pp),target_success=short_result['success'],long_success=long_result['success'],
                    target_K2=short_result['K2_outputs'],target_both_times=short_result['both_times_strict'],
                    long_K2=long_result['K2_outputs'],long_both_times=long_result['both_times_strict'],
                    worst_control_extra=max(r['extra_fraction'] for r in cc),control_extra={r['condition']:r['extra_frames'] for r in cc},
                    min_single_success=min(r['efficiency'] for r in cc if r['kind']=='single'),
                    missed=sum(r['missed'] for r in pp),extra=sum(r['extra'] for r in pp),split=sum(r['split'] for r in pp),merge=sum(r['merge'] for r in pp)))
    csv_table(OUTPUT/'threshold-summary.csv',all_summaries);csv_table(OUTPUT/'threshold-aggregate.csv',aggregate)
    lookup={(r['domain'],r['method'],r['alpha']):r for r in aggregate};envelopes=[]
    for method in METHODS:
        allowed=[a for a in alphas if all(lookup[d,method,a]['worst_control_extra']<=.05 for d in ('native','continuous'))]
        alpha=max(allowed)
        envelopes.append(dict(method=method,alpha=alpha,development_only=True,
            native=lookup['native',method,alpha],continuous=lookup['continuous',method,alpha],
            default_native=lookup['native',method,.025],default_continuous=lookup['continuous',method,.025],
            all_native=lookup['native',method,1.],all_continuous=lookup['continuous',method,1.]))
    candidates=[]
    for r in records:
        if r['domain']=='native' and r['kind']=='pair' and r['condition']=='high_weak' and r['band']=='40' and r['metrics']['high']['success']:
            alpha=r['endpoints']['rank'];native=lookup['native',r['method'],alpha];continuous=lookup['continuous',r['method'],alpha]
            candidates.append(dict(method=r['method'],index=r['index'],needed_alpha=alpha,
                native_controls=native['control_extra'],continuous_controls=continuous['control_extra'],
                worst_extra=max(native['worst_control_extra'],continuous['worst_control_extra'])))
    for domain in ('native','continuous'):
        for method in METHODS:
            rr=[lookup[domain,method,a] for a in alphas]
            assert np.all(np.diff([r['strict_success'] for r in rr])>=0)
    result=dict(aggregate=aggregate,default_summary=at_default,envelopes=envelopes,short_correct_candidates=candidates,
        alpha_grid=alphas,all_data_exposed=True,no_default_replacement=True)
    save(OUTPUT/'sweep.json',result)
    save(RECORD/'sweep-validation.json',dict(passed=True,evaluation_frames=1760,decision_records=7040,alpha_count=len(alphas),
        equivalent_selected_scores=7040*len(alphas),summary_rows=len(all_summaries),default_previous_exact=True,
        strict_pair_success_monotone=True,parameter_estimates_unchanged=True,
        output_hashes={n:digest(OUTPUT/n) for n in ('decision-records.jsonl','threshold-summary.csv','threshold-aggregate.csv','sweep.json')}))
    print('EVIDENCE_SWEEP_COMPLETE',len(alphas),len(all_summaries),[(r['method'],r['alpha'],r['native']['target_success']) for r in envelopes],flush=True)


def analyze():
    check_freeze();assert load(RECORD/'sweep-validation.json')['passed']
    inference=load(RECORD/'inference.json');probe_record=load(RECORD/'probe-validation.json')
    assert inference['passed'] and probe_record['passed']
    assert digest(OUTPUT/'profiles.jsonl')==inference['sha256']==probe_record['primary_profiles_sha256']
    assert digest(OUTPUT/'truth-probes.jsonl')==probe_record['sha256']
    curves=rows(OUTPUT/'profiles.jsonl');probes={(r['domain'],r['index']):r for r in rows(OUTPUT/'truth-probes.jsonl')}
    data=load_data();gates=load(EXPERIMENT/'gate-calibrations.json');values=[]
    for curve in curves:
        domain=curve['domain'];i=curve['index'];d=data[domain];case=d['cases'][i];q=d['q'][i];t=d['t'][i]
        band=str(int(case['gap_ns'])) if domain=='native' else case['band'];probe=probes[domain,i]
        costs=[f['cost'] for f in curve['fits']];best=min(curve['fits'],key=lambda r:r['cost'])
        selected=decide(d['cached'][i],'mlp_physical',gates['mlp_physical'])
        for level in (1.,4.,9.):
            region=compatibility_segments(curve['q_grid_pe'],costs,curve['reference_cost'],level)
            slice_region=compatibility_segments(curve['q_grid_pe'],curve['slice_cost'],min(curve['slice_cost']),level)
            assert region['segments']
            values.append(dict(domain=domain,index=i,band=band,**region,
                width_over_trueQ=(region['upper']-region['lower'])/q[1],span_ratio=region['upper']/region['lower'],
                slice_width_over_trueQ=(slice_region['upper']-slice_region['lower'])/q[1],
                components=len(region['segments']),contains500=any(a<=500.<=b for a,b in region['segments']),
                truth_Q2_pe=float(q[1]),truth_delta=probe['delta_cost'],truth_included=probe['delta_cost']<=level,
                truth_probe_converged=probe['fit']['converged'],probe_improved_reference=probe['delta_cost']<-.01,
                true_times_in_domain=bool(min(t)>=curve['bounds_ns'][0] and max(t)<=curve['bounds_ns'][1]),
                selected_K2=len(selected['times_ns'])==2,selected_strict_success=evaluate(selected,t,q,case)['success'],
                conditional_best_Q2_error=best['charge_pe'][1]/q[1]-1.,conditional_best_success=evaluate(best,t,q,case)['success'],
                reverse_gain=curve['reference_cost']-min(curve['forward_cost']),
                outside_cache=curve['cached_outside_local_domain'],nonconverged_trials=curve['nonconverged_trials'],seconds=curve['seconds']))
    # A true-Q probe must never alter the primary curve or the default decision.
    assert digest(OUTPUT/'profiles.jsonl')==inference['sha256']
    groups=defaultdict(list)
    for row in values:groups[row['domain'],row['band'],row['level']].append(row)
    summary=[]
    for (domain,band,level),rr in sorted(groups.items()):
        selected=[r for r in rr if r['selected_K2']]
        summary.append(dict(domain=domain,band=band,level=level,frames=len(rr),truth_included=sum(r['truth_included'] for r in rr),
            truth_inclusion_fraction=float(np.mean([r['truth_included'] for r in rr])),
            width_over_trueQ=quantiles([r['width_over_trueQ'] for r in rr]),span_ratio=quantiles([r['span_ratio'] for r in rr]),
            slice_width_over_trueQ=quantiles([r['slice_width_over_trueQ'] for r in rr]),
            lower_boundary=sum(r['lower_boundary'] for r in rr),upper_boundary=sum(r['upper_boundary'] for r in rr),
            disconnected=sum(r['components']>1 for r in rr),contains500=sum(r['contains500'] for r in rr),
            selected_frames=len(selected),selected_truth_included=sum(r['truth_included'] for r in selected),
            selected_strict_success=sum(r['selected_strict_success'] for r in rr),
            conditional_best_success=sum(r['conditional_best_success'] for r in rr),
            conditional_best_Q2_error=quantiles([r['conditional_best_Q2_error'] for r in rr]),
            truth_delta=quantiles([r['truth_delta'] for r in rr]),probe_lower_than_reference=sum(r['probe_improved_reference'] for r in rr),
            failed_trueQ_probes=sum(not r['truth_probe_converged'] for r in rr),
            true_times_outside=sum(not r['true_times_in_domain'] for r in rr),cached_outside=sum(r['outside_cache'] for r in rr),
            nonconverged_trials=sum(r['nonconverged_trials'] for r in rr),seconds_median=float(np.median([r['seconds'] for r in rr]))))
    assert len(curves)==160 and len(values)==480 and len(summary)==12
    save(OUTPUT/'profile-scores.json',values);save(OUTPUT/'profile-summary.json',summary);csv_table(OUTPUT/'profile-summary.csv',summary)
    result=dict(sweep=load(OUTPUT/'sweep.json'),profiles=summary,all_data_exposed=True,formal_confidence_intervals=False,
        no_default_replacement=True,new_ROOTs=0,new_training=0,known_K2_profiles=True,true_Q2_probes_for_scoring_only=True)
    save(OUTPUT/'report.json',result)
    save(RECORD/'science-validation.json',dict(passed=True,primary_profiles=160,truth_Q2_probes=160,profile_summaries=12,
        profile_level_scores=480,default_decisions_unchanged=True,primary_profile_hash_unchanged_after_truth_probe=True,
        output_hashes={n:digest(OUTPUT/n) for n in ('profile-scores.json','profile-summary.json','report.json')},
        maximum_cost_replay_error=inference['maximum_cost_replay_error']))
    print('EVIDENCE_PROFILES_SCORED',[(r['domain'],r['band'],r['truth_included'],r['width_over_trueQ'][1]) for r in summary if r['level']==4.],flush=True)


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.2})
    report=load(OUTPUT/'report.json');sweep=report['sweep'];figures={};folder=OUTPUT/'figures';folder.mkdir(exist_ok=True)
    def finish(fig,number,title,legend=None):
        for i,ax in enumerate(fig.axes):ax.set_title(f'({chr(97+i)}) '+ax.get_title())
        fig.suptitle('Figure '+number+' | '+title,fontsize=12)
        if legend is not None:
            handles,labels=legend.get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=4,fontsize=8)
        fig.tight_layout(rect=(0,.09 if legend is not None else .02,1,.94))
        p=folder/(number+'.png');fig.savefig(p,dpi=145);plt.close(fig);figures[number]=digest(p)
    fig,axes=plt.subplots(2,2,figsize=(13,8))
    for row,domain in enumerate(('native','continuous')):
        for col,(field,title) in enumerate((('target_success','40 ns' if domain=='native' else '24-64 ns'),('long_success','200 ns' if domain=='native' else '128-256 ns'))):
            ax=axes[row,col]
            for method in METHODS:
                rr=[r for r in sweep['aggregate'] if r['domain']==domain and r['method']==method]
                ax.plot([100*r['worst_control_extra'] for r in rr],[100*r[field]/40 for r in rr],color=COLORS[method],label=LABELS[method])
                r=next(r for r in rr if r['alpha']==.025);ax.scatter(100*r['worst_control_extra'],100*r[field]/40,color=COLORS[method],marker='o',s=35,zorder=4)
            ax.axvline(5.,color='black',ls='--',lw=.8);ax.set(xlabel='Worst empty/single extra-output frames (%)',ylabel='Strict pair recovery (%)',ylim=(-1,102),title=domain+' high then weak: '+title)
    finish(fig,'0004-98-1','Fixed candidates; all attainable rank thresholds; circles mark alpha=0.025',axes[0,0])
    fig,axes=plt.subplots(2,3,figsize=(15,8))
    for row,(prefix,title) in enumerate((('target','40 ns'),('long','200 ns'))):
        for col,field in enumerate((prefix+'_K2',prefix+'_both_times',prefix+'_success')):
            for method in METHODS:
                rr=[r for r in sweep['aggregate'] if r['domain']=='native' and r['method']==method]
                axes[row,col].step([r['alpha'] for r in rr],[r[field] for r in rr],where='post',label=LABELS[method],color=COLORS[method])
            axes[row,col].axvline(.025,color='black',ls='--',lw=.8)
            axes[row,col].set(xlabel='Allowed upper-tail rank alpha',ylabel='Pairs out of 40',ylim=(-1,42),title=title+': '+('selected K2','both times','Q and times correct')[col])
    finish(fig,'0004-98-2','Native high then weak: relaxed decisions do not repair incorrect parameters',axes[0,0])
    curves={(r['domain'],r['index']):r for r in rows(OUTPUT/'profiles.jsonl')};data=load_data();values=load(OUTPUT/'profile-scores.json')
    groups=[('native','40'),('native','200'),('continuous','middle'),('continuous','long')]
    fig,axes=plt.subplots(2,2,figsize=(13,8));ratios=np.geomspace(.1,20.,180)
    for ax,(domain,band) in zip(axes.flat,groups):
        rr=[r for r in values if (r['domain'],r['band'],r['level'])==(domain,band,4.)];p=[];s=[]
        for r in rr:
            c=curves[domain,r['index']];x=np.array(c['q_grid_pe'])/r['truth_Q2_pe']
            p.append(np.interp(np.log(ratios),np.log(x),np.array([f['cost'] for f in c['fits']])-c['reference_cost']))
            s.append(np.interp(np.log(ratios),np.log(x),np.array(c['slice_cost'])-min(c['slice_cost'])))
        low,median,high=np.percentile(p,[16,50,84],axis=0)
        ax.fill_between(ratios,low,high,color='#009988',alpha=.2,label='Profile: 16-84% across frames')
        ax.plot(ratios,median,color='#009988',label='Profile: median')
        ax.plot(ratios,np.median(s,axis=0),color='#ee7733',ls='--',label='Frozen nuisance slice: median')
        ax.axvline(1.,color='black',lw=.8);ax.axhline(4.,color='black',ls=':',lw=.8)
        ax.set(xscale='log',xlim=(.1,20.),ylim=(0,20),xlabel='Later Q2 / true Q2 (plot normalization only)',ylabel='Delta cost',title=domain+' / '+band)
    finish(fig,'0004-99-1','Forced K2 profiles; cost=quadratic+logdet, equal to -2 log L up to a constant',axes[0,0])
    fig,axes=plt.subplots(2,2,figsize=(13,8));x=np.arange(4);labels=['native 40','native 200','Gaussian middle','Gaussian long']
    for j,(level,color) in enumerate(zip((1.,4.,9.),('#0077bb','#009988','#ee7733'))):
        rr=[next(r for r in report['profiles'] if (r['domain'],r['band'],r['level'])==(d,b,level)) for d,b in groups]
        pos=x+(j-1)*.22
        axes[0,0].bar(pos,[100*r['truth_inclusion_fraction'] for r in rr],width=.21,color=color,label=f'Delta cost <= {level:g}')
        axes[0,1].plot(pos,[r['width_over_trueQ'][1] for r in rr],'-o',color=color,ms=4)
        axes[1,0].bar(pos,[100*r['lower_boundary']/40 for r in rr],width=.21,color=color)
    rr=[next(r for r in report['profiles'] if (r['domain'],r['band'],r['level'])==(d,b,4.)) for d,b in groups]
    axes[1,1].bar(x-.16,[100*r['truth_inclusion_fraction'] for r in rr],width=.3,label='All 40 pairs',color='#009988')
    axes[1,1].bar(x+.16,[100*r['selected_truth_included']/max(r['selected_frames'],1) for r in rr],width=.3,label='Only K2 accepted at alpha=.025',color='#555555')
    for i,r in enumerate(rr):axes[1,1].text(i+.16,100*r['selected_truth_included']/max(r['selected_frames'],1)+1,f'n={r["selected_frames"]}',ha='center',fontsize=8)
    for ax,title,ylabel in zip(axes.flat,('True-Q probe inclusion','Profile range width','Range touches lower Q=50 bound','Selection changes the evaluated subset'),('Frames (%)','Median (upper-lower) / true Q2','Frames (%)','Delta<=4 inclusion (%)')):
        ax.set(xticks=x,xticklabels=labels,title=title,ylabel=ylabel);ax.tick_params(axis='x',labelrotation=15)
    axes[0,1].set_yscale('log');axes[1,1].set_ylim(0,110);axes[1,1].legend(fontsize=8)
    finish(fig,'0004-100-1','Empirical compatibility only; levels 1/4/9 are not named confidence levels',axes[0,0])
    response,noise=models();fig,axes=plt.subplots(4,2,figsize=(14,13));examples=[]
    observations={}
    for domain in ('native','continuous'):
        with np.load(paths(domain)['observation']) as f:observations[domain]=f['voltage_mV']
    for j,task in enumerate(load(EXPERIMENT/'protocol.json')['examples']):
        domain=task['domain'];i=task['index'];c=curves[domain,i];d=data[domain];q=d['q'][i];t=d['t'][i]
        accepted=[f for f in c['fits'] if f['cost']<=c['reference_cost']+4.]
        low=min(accepted,key=lambda f:f['charge_pe'][1]);high=max(accepted,key=lambda f:f['charge_pe'][1]);best=min(c['fits'],key=lambda f:f['cost'])
        truth=response.components(t,q);axes[j,0].plot(response.sample_grid,observations[domain][i],color='#bbbbbb',lw=.7,label='Observed full voltage')
        axes[j,1].plot(response.sample_grid,truth[1],color='black',lw=1.,label='True-parameter conditional component')
        chosen={}
        for label,fit,color,style in (('Best sampled K2',best,'#009988','-'),('Low Q2, delta<=4',low,'#0077bb','--'),('High Q2, delta<=4',high,'#ee7733',':')):
            m=response.components(fit['times_ns'],fit['charge_pe']);axes[j,0].plot(response.sample_grid,m.sum(0)+response.baseline+noise.background_mV,color=color,ls=style,lw=1.,label=label)
            axes[j,1].plot(response.sample_grid,m[1],color=color,ls=style,lw=1.,label=label+f' ({fit["charge_pe"][1]:.0f} pe)');chosen[label]=fit
        for col,ax in enumerate(axes[j]):
            ax.set(xlim=(0,2048),xlabel='Record time (ns)',ylabel='mV',title=f'{domain} index{i}: '+('total waveform' if col==0 else 'later component (forced K2)'));ax.legend(fontsize=7)
        examples.append(dict(task,models=chosen,true_Q2_pe=float(q[1]),true_gap_ns=float(t[1]-t[0])))
    finish(fig,'0004-101-1','Prespecified full records: different conditional K2 explanations of the same voltage')
    save(OUTPUT/'examples.json',examples)
    assert list(figures)==load(EXPERIMENT/'protocol.json')['figures']
    save(RECORD/'figure-delivery.json',dict(passed=True,figures=figures,source_sha256=digest(__file__)))
    print('EVIDENCE_FIGURES_COMPLETE',len(figures),flush=True)
