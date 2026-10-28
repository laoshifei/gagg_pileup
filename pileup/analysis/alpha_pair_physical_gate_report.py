"""Paired scoring and diagnostics for physical count decision with ML seeds."""
from __future__ import annotations

from collections import defaultdict
import json

import numpy as np

from .alpha_pair_physical_gate_study import PROJECT,EXPERIMENT,OUTPUT,RECORD,INPUT,input_paths,models,check_freeze
from .alpha_pair_physical_gate import METHODS,PRIMARY,bank_for,decide,gates_from_controls,rejection
from .alpha_pair_information_study import load,save,digest,rows
from .alpha_pair_risk_report import evaluate,summarize,quantiles,csv_table

LABELS=dict(baseline='Original physical NLL',learned_gate='ML seed / learned gate',local_physical='Local refit / physical gate',
            mlp_physical='ML seed / physical gate',mlp_fixed_gate='ML seed / original threshold')
COLORS=dict(zip(METHODS,('#555555','#ee7733','#0077bb','#009988','#bb4499')))


def read_data(domain):
    paths=input_paths(domain);cases=load(paths['cases'])
    with np.load(paths['truth']) as f:q=f['charge_pe'];t=f['times_ns']
    cached={r['index']:r for r in rows(OUTPUT/(domain+'-candidates.jsonl'))}
    assert set(cached)==set(range(len(cases)))
    return cases,q,t,cached


def analyze():
    check_freeze();inference=load(RECORD/'inference.json');assert inference['passed']
    for name,expected in inference['outputs'].items():assert digest(OUTPUT/name)==expected,name
    for name,expected in load(RECORD/'input-provenance.json')['inputs'].items():assert digest(PROJECT/name)==expected,name
    data={domain:read_data(domain) for domain in ('native','continuous')}
    controls=[data['native'][3][i] for i in range(280)]
    gates=gates_from_controls(controls,load(EXPERIMENT/'baseline-gate.json'),load(EXPERIMENT/'previous-gates.json')['mlp_nll'])
    save(EXPERIMENT/'gate-calibrations.json',gates)
    old_scores={(r['domain'],r['index'],r['method']):r for r in rows(INPUT/'scores.jsonl') if r['method'] in ('baseline','mlp_nll')}
    scored=[];forced=[];search=[]
    for domain,(cases,q,t,cached) in data.items():
        for i,case in enumerate(cases):
            if domain=='native' and i<280:continue
            row=cached[i];band=str(int(case['gap_ns'])) if domain=='native' and case['kind']=='pair' else case.get('band')
            identity=dict(index=i,domain=domain,kind=case['kind'],condition=case['condition'],band=band,
                gap_ns=case['gap_ns'],block=case.get('block',i))
            diagnostic={}
            if case['kind']=='pair':
                for name,bank in row['sets'].items():
                    assert len(bank)==3
                    value=evaluate(bank[2],t[i],q[i],case);diagnostic[name]=value
                    forced.append(dict(value,**identity,bank=name,below500=min(bank[2]['charge_pe'])<500.,
                        weak_Q_error=float(np.asarray(bank[2]['charge_pe'])[np.argmin(q[i])]/q[i].min()-1.),
                        total_Q_error=float(sum(bank[2]['charge_pe'])/q[i].sum()-1.)))
            if len(row['sets']['baseline'])==3:
                a,b,c=(row['sets'][k][2] for k in ('baseline','local','mlp'))
                search.append(dict(**identity,local_cost_gain=a['cost']-b['cost'],network_cost_gain=b['cost']-c['cost'],
                    network_lower_cost=b['cost']-c['cost']>1e-6,network_time_change_ns=float(np.max(np.abs(np.array(b['times_ns'])-c['times_ns']))),
                    network_Q_change_pe=float(np.max(np.abs(np.array(b['charge_pe'])-c['charge_pe'])))))
            for method in METHODS:
                fit=decide(row,method,gates[method]);value=evaluate(fit,t[i],q[i],case)
                reason=rejection(row,method,fit)
                assert (reason=='accepted')==(len(fit['times_ns'])==2)
                value.update(identity,method=method,sensitivity=method not in PRIMARY,fit=fit,rejection=reason)
                if case['kind']=='pair':
                    d=diagnostic[bank_for(method)];value['forced_success']=d['success']
                    value['correct_candidate_rejected']=bool(d['success'] and reason!='accepted')
                    assert value['success']==(d['success'] and reason=='accepted')
                scored.append(value)
                if method in ('baseline','learned_gate'):
                    old=old_scores[domain,i,'baseline' if method=='baseline' else 'mlp_nll']
                    for key in ('success','target','recovered','strict_recovered','matched','extra','missed','split','merge','q_errors','t_errors','count'):
                        assert value[key]==old[key],(domain,i,method,key)
    assert len(scored)==8800 and len(forced)==2880
    for filename,values in (('scores.jsonl',scored),('forced-K2.jsonl',forced),('search-diagnostics.jsonl',search)):
        with (OUTPUT/filename).open('w') as stream:
            for row in values:stream.write(json.dumps(row,allow_nan=False)+'\n')
    groups=defaultdict(list)
    for r in scored:groups[r['domain'],r['method'],r['kind'],r['condition'],r['band']].append(r)
    summary=[dict(summarize(v),sensitivity=v[0]['sensitivity']) for v in groups.values()]
    aggregates=[];stages=[];forced_summary=[]
    for domain in ('native','continuous'):
        for method in METHODS:
            rr=[r for r in scored if r['domain']==domain and r['method']==method and r['kind']=='pair']
            parts=[('all',rr)]
            if domain=='continuous':parts += [('targeted',[r for r in rr if r['condition']!='uniform_shift']),('prior_shift',[r for r in rr if r['condition']=='uniform_shift'])]
            for part,subset in parts:
                aggregates.append(dict(domain=domain,method=method,part=part,frames=len(subset),success=sum(r['success'] for r in subset),
                    efficiency=float(np.mean([r['success'] for r in subset])),both_times_strict=sum(r['strict_matched']==2 for r in subset),
                    correct_candidate_rejected=sum(r['correct_candidate_rejected'] for r in subset)))
    for key,rr in groups.items():
        if key[2]!='pair':continue
        correct=[r for r in rr if r['forced_success']]
        counts={reason:sum(r['rejection']==reason for r in correct) for reason in ('accepted','subthreshold_Q','no_presence','no_pair','nonpositive_evidence','rank_rejected')}
        assert sum(counts.values())==len(correct)
        stages.append(dict(domain=key[0],method=key[1],condition=key[3],band=key[4],frames=len(rr),forced_correct=len(correct),
            selected_K2=sum(r['count']==2 for r in rr),both_times_strict=sum(r['strict_matched']==2 for r in rr),
            strict_success=sum(r['success'] for r in rr),correct_candidate_outcomes=counts))
    groups=defaultdict(list)
    for r in forced:groups[r['domain'],r['bank'],r['condition'],r['band']].append(r)
    for (domain,bank,condition,band),rr in groups.items():
        forced_summary.append(dict(domain=domain,bank=bank,condition=condition,band=band,frames=len(rr),success=sum(r['success'] for r in rr),
            both_times_strict=sum(r['strict_matched']==2 for r in rr),below500=sum(r['below500'] for r in rr),
            weak_Q_error_p16_p50_p84=quantiles([r['weak_Q_error'] for r in rr]),total_Q_error_p16_p50_p84=quantiles([r['total_Q_error'] for r in rr]),
            not_converged=sum(not r['converged'] for r in rr)))
    lookup={(r['domain'],r['method'],r['index']):r for r in scored};paired=[];rng=np.random.default_rng(270327062)
    groups=defaultdict(list)
    for r in scored:
        if r['method']=='baseline' and r['kind']=='pair':groups[r['domain'],r['condition'],r['band']].append(r)
    comparisons=(('local_physical','baseline'),('mlp_physical','baseline'),('mlp_physical','learned_gate'),('mlp_physical','local_physical'),('mlp_fixed_gate','baseline'))
    for (domain,condition,band),rr in groups.items():
        blocks=np.array([r['block'] if domain=='native' else j for j,r in enumerate(rr)]);unique=np.unique(blocks)
        draws=rng.choice(unique,(2000,len(unique)),replace=True)
        weights=np.stack([np.count_nonzero(draws==b,axis=1) for b in blocks],axis=1)
        for method,reference in comparisons:
            difference=np.array([float(lookup[domain,method,r['index']]['success'])-float(lookup[domain,reference,r['index']]['success']) for r in rr])
            samples=weights@difference/weights.sum(1)
            paired.append(dict(domain=domain,condition=condition,band=band,method=method,reference=reference,
                difference=float(difference.mean()),paired95=np.percentile(samples,[2.5,97.5]).tolist(),
                gained=int((difference>0).sum()),lost=int((difference<0).sum())))
    aggregate={(r['domain'],r['method'],r['part']):r for r in aggregates}
    decisions=[]
    for method in ('local_physical','mlp_physical','mlp_fixed_gate'):
        checks={}
        for domain in ('native','continuous'):
            controls=[r for r in summary if r['domain']==domain and r['method']==method and r['kind']!='pair']
            singles=[r for r in controls if r['kind']=='single']
            checks[domain+'_extra']=all(r['extra_fraction']<=.05 for r in controls)
            checks[domain+'_single_success']=all(r['efficiency']>=.95 for r in singles)
            checks[domain+'_single_bias']=all(r['q_error_p16_p50_p84'] and abs(r['q_error_p16_p50_p84'][1])<=.05 for r in singles)
        checks['native200_recovered']=next(r for r in summary if (r['domain'],r['method'],r['condition'],r['band'])==('native',method,'high_weak','200'))['success']>=20
        checks['native_total_preserved']=aggregate['native',method,'all']['success']>=241
        checks['continuous_total_preserved']=aggregate['continuous',method,'all']['success']>=500
        decisions.append(dict(method=method,gate_repair_passed=all(checks.values()),checks=checks,sensitivity=method not in PRIMARY))
    value_checks={d:aggregate[d,'mlp_physical','all']['efficiency']>=aggregate[d,'local_physical','all']['efficiency']+.05-1e-12 for d in ('native','continuous')}
    value_checks['controls']=next(r for r in decisions if r['method']=='mlp_physical')['gate_repair_passed']
    inference_counts={}
    for domain in ('native','continuous'):
        rr=[r for r in search if r['domain']==domain]
        inference_counts[domain]=dict(frames=len(rr),network_lower_cost=sum(r['network_lower_cost'] for r in rr),
            network_cost_gain_p16_p50_p84=quantiles([r['network_cost_gain'] for r in rr]))
    report=dict(summary=summary,aggregates=aggregates,stages=stages,forced_summary=forced_summary,paired=paired,
        decisions=decisions,network_practical_value=dict(passed=all(value_checks.values()),checks=value_checks),search=inference_counts,
        all_data_exposed=True,new_ROOTs=0,new_training=0,new_optimization=0,no_default_replacement=True,
        count_gate_rule='Physical NLL difference; same calibration/control rules; original numerical threshold sensitivity separately labelled')
    save(OUTPUT/'report.json',report);save(OUTPUT/'summary.json',summary);csv_table(OUTPUT/'summary.csv',summary)
    csv_table(OUTPUT/'pair-aggregates.csv',aggregates);csv_table(OUTPUT/'stage-diagnostics.csv',stages);csv_table(OUTPUT/'forced-K2-summary.csv',forced_summary)
    save(OUTPUT/'paired-comparison.json',dict(comparisons=paired,replicates=2000,seed=270327062,
        scope='Exposed-data descriptive paired intervals per condition; native optical-parent blocks, continuous independent draws; degenerate intervals are not certainty'))
    save(RECORD/'science-validation.json',dict(passed=True,frames=2040,evaluation_frames=1760,calibration_frames=280,
        primary_scores=7040,all_scores=8800,forced_K2_scores=2880,old_baseline_and_learned_scores_exact=True,K0_K1_preserved=True,
        all_correct_candidate_rejections_partitioned=True,physical_gate_uses_no_network_score=True,all_failed_fits_retained=True,
        input_cost_replay_relative_error=inference['maximum_cost_replay_relative_error'],
        output_hashes={n:digest(OUTPUT/n) for n in ('scores.jsonl','forced-K2.jsonl','search-diagnostics.jsonl','report.json','summary.json','paired-comparison.json')},source_sha256=digest(__file__)))
    print('PHYSICAL_GATE_ANALYSIS_COMPLETE','repair',[(r['method'],r['gate_repair_passed']) for r in decisions],'network_value',report['network_practical_value']['passed'],flush=True)


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.2})
    assert load(RECORD/'science-validation.json')['passed']
    report=load(OUTPUT/'report.json');summary=report['summary'];hashes={}
    folder=OUTPUT/'figures';folder.mkdir(exist_ok=True)
    def finish(fig,number,title,legend_axis=None):
        fig.suptitle('Figure '+number+' | '+title,fontsize=12)
        for i,ax in enumerate(fig.axes):ax.set_title(f'({chr(97+i)}) '+ax.get_title())
        if legend_axis is not None:
            handles,labels=legend_axis.get_legend_handles_labels()
            fig.legend(handles,labels,loc='lower center',ncol=3,fontsize=8)
        fig.tight_layout(rect=(0,.12 if legend_axis is not None else .02,1,.92))
        path=folder/(number+'.png');fig.savefig(path,dpi=145);plt.close(fig);hashes[number]=digest(path)
    def grouped(ax,domain,condition,bands):
        for j,method in enumerate(METHODS):
            values=[next(r for r in summary if (r['domain'],r['method'],r['condition'],r['band'])==(domain,method,condition,b)) for b in bands]
            ax.bar(np.arange(len(bands))+(j-2)*.15,[100*r['efficiency'] for r in values],width=.14,
                label=LABELS[method],color=COLORS[method],hatch='//' if method=='mlp_fixed_gate' else None)
        ax.set(xticks=np.arange(len(bands)),ylim=(0,104),ylabel='Strict pair recovery (%)',title=condition)
    fig,axes=plt.subplots(1,3,figsize=(15,4.8))
    for ax,condition in zip(axes,('mid_mid','high_weak','weak_high')):
        grouped(ax,'native',condition,('8','40','200'));ax.set(xticklabels=['8','40','200'],xlabel='True separation (ns)')
    finish(fig,'0004-93-1','Exposed native development: same physical count rule, 40 pairs/group',axes[0])
    fig,axes=plt.subplots(2,2,figsize=(13,8))
    for ax,condition in zip(axes.flat,('balanced','high_weak','weak_high','uniform_shift')):
        bands=('all',) if condition=='uniform_shift' else ('short','middle','long')
        grouped(ax,'continuous',condition,bands)
        ax.set(xticklabels=['4-320 ns'] if condition=='uniform_shift' else ['6-12 ns','24-64 ns','128-256 ns'])
    finish(fig,'0004-93-2','Previously exposed continuous Gaussian draws; original threshold is sensitivity only',axes[0,0])
    fig,axes=plt.subplots(2,2,figsize=(13,8));outcomes=('accepted','rank_rejected','subthreshold_Q','other')
    outcome_colors=('#009988','#999999','#ee7733','#bb4499')
    for row,band in enumerate(('40','200')):
        stages=[next(r for r in report['stages'] if (r['domain'],r['method'],r['condition'],r['band'])==('native',m,'high_weak',band)) for m in PRIMARY]
        for j,method in enumerate(PRIMARY):
            axes[row,0].bar(np.arange(3)+(j-1.5)*.20,[stages[j][k] for k in ('selected_K2','both_times_strict','strict_success')],
                width=.19,color=COLORS[method],label=LABELS[method])
        axes[row,0].set(xticks=np.arange(3),xticklabels=['Selected K=2','Both times','Both Q and times'],ylim=(0,42),ylabel='Pairs out of 40',title=band+' ns: selected-output stages')
        bottom=np.zeros(len(PRIMARY))
        for outcome,color in zip(outcomes,outcome_colors):
            values=[s['correct_candidate_outcomes'].get(outcome,0) if outcome!='other' else sum(s['correct_candidate_outcomes'][k] for k in ('no_presence','no_pair','nonpositive_evidence')) for s in stages]
            axes[row,1].bar(np.arange(len(PRIMARY)),values,bottom=bottom,color=color,label=outcome);bottom+=values
        for j,n in enumerate(bottom):axes[row,1].text(j,n+.6,str(int(n)),ha='center',fontsize=9)
        axes[row,1].set(xticks=np.arange(4),xticklabels=['baseline','learned','local','ML + physical'],ylim=(0,42),ylabel='Correct K=2 candidates',title=band+' ns: fate of correct candidates')
        axes[row,1].legend(fontsize=7,loc='upper left')
    finish(fig,'0004-94-1','Native high then weak: candidate correctness and blind decision are separate',axes[0,0])
    fig,axes=plt.subplots(2,3,figsize=(15,8))
    for row,domain in enumerate(('native','continuous')):
        for j,method in enumerate(PRIMARY):
            values=[next(r for r in summary if (r['domain'],r['method'],r['condition'],r['kind'])==(domain,method,c,'empty' if c=='all' else 'single')) for c in ('all','weak','low','mid','high')]
            x=np.arange(5)+(j-1.5)*.14
            axes[row,0].bar(x,[100*r['extra_fraction'] for r in values],width=.13,color=COLORS[method],label=LABELS[method])
            axes[row,1].plot(x,[100*r['efficiency'] for r in values],'-o',ms=3,color=COLORS[method])
            axes[row,2].plot(x[1:],[100*r['q_error_p16_p50_p84'][1] for r in values[1:]],'-o',ms=3,color=COLORS[method])
        for col,label in enumerate(('Extra-output frames (%)','Strict success (%)','Matched median Q error (%)')):
            axes[row,col].set(xticks=np.arange(5),xticklabels=['empty','weak','low','mid','high'],title=domain,ylabel=label)
        axes[row,0].axhline(5.,color='black',ls='--',lw=.8);axes[row,1].axhline(95.,color='black',ls='--',lw=.8)
        for value in (-5.,5.):axes[row,2].axhline(value,color='black',ls='--',lw=.8)
    finish(fig,'0004-95-1','Empty and single controls, 80/group; same 280 native calibration singles',axes[0,0])
    cases,q,t,cached=read_data('native');response,noise=models();gates=load(EXPERIMENT/'gate-calibrations.json')
    with np.load(input_paths('native')['observation']) as f:y=f['voltage_mV']
    original=PROJECT/'jupyter/artifacts/alpha_pair_joint_confirmation_20261026/labels.truth.npz'
    with np.load(original) as f:actual=f['components_mV']
    fig,axes=plt.subplots(2,3,figsize=(16,8.6));examples=[]
    for row,gap in enumerate((40.,200.)):
        i=next(i for i,c in enumerate(cases) if i>=280 and c['condition']=='high_weak' and c['gap_ns']==gap and c['replicate']==0)
        mean=response.components(t[i],q[i]);weak=int(np.argmin(q[i]));strong=1-weak
        axes[row,0].plot(response.sample_grid,y[i],color='#bbbbbb',lw=.7,label='Observed voltage')
        for col,index in ((1,strong),(2,weak)):
            axes[row,col].plot(response.sample_grid,actual[cases[i]['original_index'],index],color='black',lw=.7,label='Native actual component (truth)')
            axes[row,col].plot(response.sample_grid,mean[index],color='#999999',ls='--',lw=.8,label='Conditional true-parameter mean')
        example=dict(index=i,gap_ns=gap,methods={})
        for method in PRIMARY:
            fit=decide(cached[i],method,gates[method]);curves=response.components(fit['times_ns'],fit['charge_pe'])
            axes[row,0].plot(response.sample_grid,curves.sum(0)+response.baseline+noise.background_mV,color=COLORS[method],lw=.9,label=LABELS[method])
            larger=curves[np.argmax(fit['charge_pe'])] if len(curves) else np.zeros(1024)
            smaller=curves[np.argmin(fit['charge_pe'])] if len(curves)==2 else np.zeros(1024)
            axes[row,1].plot(response.sample_grid,larger,color=COLORS[method],lw=.9,label=LABELS[method])
            axes[row,2].plot(response.sample_grid,smaller,color=COLORS[method],lw=.9,label=LABELS[method]+f' (K={len(curves)})')
            example['methods'][method]=dict(fit=fit,score=evaluate(fit,t[i],q[i],cases[i]))
        for col,title in enumerate(('Total voltage','Larger-Q selected output','Smaller-Q selected output')):
            axes[row,col].set(xlim=(0,2048),xlabel='Record time (ns)',ylabel='mV',title=f'{gap:g} ns: {title}');axes[row,col].legend(fontsize=6)
        examples.append(example)
    finish(fig,'0004-96-1','Prespecified native pair0; full record and final selected components; missing second output=0')
    save(OUTPUT/'examples.json',examples)
    save(RECORD/'plot-truth-provenance.json',dict(file=str(original.relative_to(PROJECT)),sha256=digest(original),purpose='Plot only; no selection/inference access'))
    assert list(hashes)==load(EXPERIMENT/'protocol.json')['figures']
    save(RECORD/'figure-delivery.json',dict(passed=True,figures=hashes,source_sha256=digest(__file__)))
    print('PHYSICAL_GATE_FIGURES_COMPLETE',len(hashes),flush=True)
