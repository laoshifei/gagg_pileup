"""Paired scoring and failure analysis for the charge-allocation study."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import json

import numpy as np

from .alpha_pair_allocation_study import PROJECT,EXPERIMENT,OUTPUT,RECORD,CLOSURE,models,check_freeze
from .alpha_pair_information_study import load,save,digest,rows
from .alpha_pair_charge_allocation import VARIANTS
from .alpha_pair_joint_likelihood import calibrate_gate,select
from .alpha_pileup_fit import score

METHODS=('baseline',*VARIANTS)
LABELS={'baseline':'Baseline','component_local':'Component / local','allocation_local':'Total-share / local',
        'component_profile':'Component / profile','allocation_profile':'Total-share / profile'}
COLORS=dict(zip(METHODS,('#555555','#0077bb','#ee7733','#33bbee','#cc3311')))
CONDITIONS=('mid_mid','high_weak','weak_high')
GAPS=(8.,40.,200.)


def quantiles(values):return np.percentile(values,[16,50,84]).tolist() if len(values) else None


def table_csv(path,records):
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)


def evaluate(fit,t,q,case):
    ids=np.flatnonzero(q>=500.)
    truth=dict(first_active_ns=t[ids],emission_ns=t[ids],charge_pe=q[ids])
    normal=score(fit,truth,[0.,2048.])
    tolerance=min(10.,case['gap_ns']/4.) if case['kind']=='pair' else 10.
    strict=score(fit,truth,[0.,2048.],tolerance_ns=tolerance)
    return dict(success=strict['recovered_count']==len(ids) and normal['unmatched_predictions']==0,
        target=len(ids),recovered=normal['recovered_count'],strict_recovered=strict['recovered_count'],
        matched=normal['matched_count'],strict_matched=strict['matched_count'],extra=normal['unmatched_predictions'],missed=normal['missed_count'],
        split=normal['split_flags'],merge=normal['merge_flags'],q_errors=normal['relative_charge_errors'],
        t_errors=normal['time_errors_ns'],count=len(fit['times_ns']),converged=fit.get('converged',True))


def summarize(rr):
    first=rr[0];q=[v for r in rr for v in r['q_errors']];t=[v for r in rr for v in r['t_errors']]
    qs=quantiles(q);ts=quantiles(t)
    return dict(domain=first['domain'],method=first['method'],gate=first['gate'],kind=first['kind'],
        condition=first['condition'],gap_ns=first['gap_ns'],frames=len(rr),success=sum(r['success'] for r in rr),
        efficiency=float(np.mean([r['success'] for r in rr])),target=sum(r['target'] for r in rr),
        recovered=sum(r['recovered'] for r in rr),strict_recovered=sum(r['strict_recovered'] for r in rr),
        missed=sum(r['missed'] for r in rr),extra=sum(r['extra'] for r in rr),
        extra_frames=sum(r['extra']>0 for r in rr),extra_fraction=float(np.mean([r['extra']>0 for r in rr])),
        split=sum(r['split'] for r in rr),merge=sum(r['merge'] for r in rr),
        q_error_p16_p50_p84=qs,q_error_sigma68=(qs[2]-qs[0])/2. if qs else None,
        time_error_p16_p50_p84_ns=ts,not_converged=sum(not r['converged'] for r in rr),
        K2_outputs=sum(r['count']==2 for r in rr),both_times10=sum(r['target']==2 and r['matched']==2 for r in rr),
        both_times_strict=sum(r['target']==2 and r['strict_matched']==2 for r in rr))


def read_data():
    cases=load(OUTPUT/'labels.cases.json')
    with np.load(OUTPUT/'labels.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
    base={r['index']:r['candidates'] for r in rows(OUTPUT/'predictions.baseline.jsonl')}
    predicted={r['index']:r for r in rows(OUTPUT/'predictions.allocation.jsonl')}
    return cases,q,t,base,predicted


def analyze():
    check_freeze();cases,q,t,base,predicted=read_data()
    for record in ('numerical-tests','pilot-validation','inference'):
        assert load(RECORD/(record+'.json'))['passed'],record
    assert digest(OUTPUT/'predictions.allocation.jsonl')==load(RECORD/'inference.json')['prediction_sha256']
    assert set(predicted)==set(range(1800)) and len(cases)==len(q)==len(t)==1800
    gates={'baseline':load(EXPERIMENT/'baseline-gate.json')}
    assert calibrate_gate([base[i] for i in range(280)])==gates['baseline']
    for method in VARIANTS:gates[method]=calibrate_gate([predicted[i]['variants'][method]['candidates'] for i in range(280)])
    save(EXPERIMENT/'gate-calibrations.json',gates)
    scored=[];candidate_rows=[];cost_error=0.;profile_count=profile_not_converged=0
    response,noise=models()
    with np.load(OUTPUT/'observations.npz') as f:voltage=f['voltage_mV']
    for i in range(1800):
        c=cases[i];p=predicted[i]
        cost_error=max(cost_error,p['maximum_cost_replay_error'])
        for row in p['profile']['profiles']:
            assert row['cost']<=row['unprofiled_cost']+1e-9
            profile_count+=1;profile_not_converged+=not row['converged']
        for method in METHODS:
            candidates=base[i] if method=='baseline' else p['variants'][method]['candidates']
            assert candidates[:2]==base[i][:2]
            if len(base[i])==3:
                assert candidates[2]['cost']<=base[i][2]['cost']+1e-10
                if method!='baseline':
                    best=candidates[2];curves=response.components(best['times_ns'],best['charge_pe'])
                    replay=noise.score(voltage[i]-response.baseline-noise.background_mV-curves.sum(0),curves,best['charge_pe'])
                    err=abs(replay-best['cost'])/max(1.,abs(replay));cost_error=max(cost_error,err);assert err<1e-7
            if c['role']=='gate':continue
            for gate_name,gate in [('primary',gates[method])]+([] if method=='baseline' else [('original',gates['baseline'])]):
                fit=select(candidates,gate)
                scored.append(dict(evaluate(fit,t[i],q[i],c),index=i,native_index=c['native_index'],domain=c['domain'],
                    method=method,gate=gate_name,kind=c['kind'],condition=c['condition'],gap_ns=c['gap_ns'],block=c['block'],
                    fit=fit))
            if c['kind']=='pair':
                assert len(candidates)==3 and np.all(q[i]>=500.)
                fit=candidates[2];order=np.argsort(fit['times_ns']);fq=np.asarray(fit['charge_pe'])[order];ft=np.asarray(fit['times_ns'])[order]
                weak=int(np.argmin(q[i]));strong=1-weak
                candidate_rows.append(dict(evaluate(fit,t[i],q[i],c),index=i,native_index=c['native_index'],domain=c['domain'],
                    method=method,condition=c['condition'],gap_ns=c['gap_ns'],block=c['block'],
                    cost_gain=base[i][2]['cost']-fit['cost'],below500=bool(np.any(fq<500.)),fitted_gap_ns=float(ft[1]-ft[0]),
                    total_q_error=float(fq.sum()/q[i].sum()-1.),weak_q_error=float(fq[weak]/q[i,weak]-1.),
                    weak_delta_Q_pe=float(fq[weak]-q[i,weak]),strong_delta_Q_pe=float(fq[strong]-q[i,strong]),
                    weak_time_error_ns=float(ft[weak]-t[i,weak]),strong_time_error_ns=float(ft[strong]-t[i,strong])))
    assert len(scored)==13680 and len(candidate_rows)==3600
    original={(r['domain'],r['index']):r for r in rows(CLOSURE/'scores.jsonl') if r['method'] in ('native/native_gate','matched/native_gate')}
    for r in scored:
        if r['method']!='baseline':continue
        old_index=r['native_index']+(280 if r['domain']=='matched' else 0)
        old=original[r['domain'],old_index]
        for key in ('success','target','recovered','strict_recovered','matched','extra','missed','split','merge','q_errors','t_errors','count'):
            assert r[key]==old[key],(r['index'],key)
    with (OUTPUT/'scores.jsonl').open('w') as stream:
        for r in scored:stream.write(json.dumps(r,allow_nan=False)+'\n')
    grouped=defaultdict(list)
    for r in scored:grouped[r['domain'],r['method'],r['gate'],r['kind'],r['condition'],r['gap_ns']].append(r)
    summary=[summarize(rr) for rr in grouped.values()]
    assert len(summary)==252
    save(OUTPUT/'summary.json',summary);table_csv(OUTPUT/'summary.csv',summary)
    primary=[r for r in scored if r['gate']=='primary']
    assert len(primary)==7600
    lookup={(r['domain'],r['method'],r['condition'],r['gap_ns']):r for r in summary if r['gate']=='primary'}

    decisions=[]
    for method in VARIANTS:
        s=lambda condition,gap:lookup['native',method,condition,gap]
        b=lambda condition,gap:lookup['native','baseline',condition,gap]
        controls=[s(c,None) for c in ('all','weak','low','mid','high')]
        singles=controls[1:]
        weak_mean=float(np.mean([s(c,g)['efficiency'] for c in ('high_weak','weak_high') for g in GAPS]))
        old_mean=float(np.mean([b(c,g)['efficiency'] for c in ('high_weak','weak_high') for g in GAPS]))
        checks=dict(controls_extra=all(r['extra_fraction']<=.05 for r in controls),
            single_success=all(r['efficiency']>=.95 for r in singles),
            single_bias=all(r['q_error_p16_p50_p84'] is not None and abs(r['q_error_p16_p50_p84'][1])<=.05 for r in singles),
            high_weak40_gain=s('high_weak',40.)['efficiency']>=b('high_weak',40.)['efficiency']+.10-1e-12,
            high_weak200_no_drop=s('high_weak',200.)['success']>=b('high_weak',200.)['success'],
            weak_high8_drop=s('weak_high',8.)['efficiency']>=b('weak_high',8.)['efficiency']-.05-1e-12,
            mid_mid=all(s('mid_mid',g)['efficiency']>=.95 for g in GAPS),weak_mean=weak_mean>=old_mean-1e-12)
        decisions.append(dict(method=method,passed=all(checks.values()),checks=checks,weak_mean=weak_mean,
            baseline_weak_mean=old_mean,high_weak40_plus200=s('high_weak',40.)['success']+s('high_weak',200.)['success']))
    eligible=[r for r in decisions if r['passed']]
    chosen=max(eligible,key=lambda r:(r['high_weak40_plus200'],r['weak_mean'],-VARIANTS.index(r['method'])))['method'] if eligible else None

    diagnostic_summary=[]
    grouped=defaultdict(list)
    primary_lookup={(r['method'],r['index']):r for r in primary}
    for r in candidate_rows:grouped[r['domain'],r['method'],r['condition'],r['gap_ns']].append(r)
    for (domain,method,condition,gap),rr in grouped.items():
        diagnostic_summary.append(dict(domain=domain,method=method,condition=condition,gap_ns=gap,frames=len(rr),
            known_K_success=sum(r['success'] for r in rr),below500=sum(r['below500'] for r in rr),
            positive_cost_gain=sum(r['cost_gain']>1e-6 for r in rr),cost_gain_p16_p50_p84=quantiles([r['cost_gain'] for r in rr]),
            weak_q_error_p16_p50_p84=quantiles([r['weak_q_error'] for r in rr]),
            total_q_error_p16_p50_p84=quantiles([r['total_q_error'] for r in rr]),
            weak_time_error_p16_p50_p84_ns=quantiles([r['weak_time_error_ns'] for r in rr]),
            not_converged=sum(not r['converged'] for r in rr),fitted_gap_p16_p50_p84_ns=quantiles([r['fitted_gap_ns'] for r in rr]),
            gap_under2ns=sum(r['fitted_gap_ns']<2. for r in rr),
            gained_correct=sum(primary_lookup[method,r['index']]['success'] and not primary_lookup['baseline',r['index']]['success'] for r in rr),
            lost_correct=sum(not primary_lookup[method,r['index']]['success'] and primary_lookup['baseline',r['index']]['success'] for r in rr)))
    save(OUTPUT/'candidate-diagnostics.json',dict(rows=candidate_rows,summary=diagnostic_summary,
        warning='Forced K=2 candidates, including subthreshold/failed fits; not blind selected performance; chronological truth assignment'))
    table_csv(OUTPUT/'candidate-summary.csv',diagnostic_summary)

    rng=np.random.default_rng(2026102703)
    weights=np.array([np.bincount(np.r_[rng.choice(np.arange(20),20,replace=True),rng.choice(np.arange(20,40),20,replace=True)],minlength=40) for _ in range(2000)])
    paired=[]
    for domain in ('native','matched'):
        for condition in CONDITIONS:
            for gap in GAPS:
                rr=sorted([r for r in primary if r['domain']==domain and r['method']=='baseline' and (r['condition'],r['gap_ns'])==(condition,gap)],key=lambda r:r['index'])
                blocks=np.array([r['block'] for r in rr]);observed=np.array([r['success'] for r in rr],float)
                denominator=weights[:,blocks].sum(1)
                assert np.all(denominator>0)
                for method in VARIANTS:
                    difference=np.array([primary_lookup[method,r['index']]['success'] for r in rr],float)-observed
                    draws=(weights[:,blocks]@difference)/denominator
                    paired.append(dict(domain=domain,condition=condition,gap_ns=gap,method=method,
                        difference=float(difference.mean()),paired_block95=np.percentile(draws,[2.5,97.5]).tolist(),
                        gained=int(np.sum(difference>0)),lost=int(np.sum(difference<0))))
    save(OUTPUT/'paired-comparison.json',dict(comparisons=paired,replicates=2000,seed=2026102703,
        scope='Descriptive exposed-parent development intervals; first20/last20 stratified blocks; common weights for both domains and all methods; degenerate intervals are not certainty'))

    timing=[]
    for domain in ('native','matched'):
        ids=[i for i,c in enumerate(cases) if c['role']=='evaluation' and c['domain']==domain]
        for method in VARIANTS:
            shared='_profile' in method
            costs=[predicted[i]['variants'][method]['seconds']+(predicted[i]['profile']['seconds'] if shared else 0.) for i in ids]
            evaluations=[predicted[i]['variants'][method]['evaluations']+(predicted[i]['profile']['evaluations'] if shared else 0) for i in ids]
            timing.append(dict(domain=domain,method=method,frames=len(ids),incremental_seconds_p16_p50_p84=quantiles(costs),
                incremental_seconds_mean=float(np.mean(costs)),objective_evaluations_mean=float(np.mean(evaluations)),
                shared_profile_seconds_mean=float(np.mean([predicted[i]['profile']['seconds'] for i in ids])) if shared else 0.))
    save(OUTPUT/'timing.json',dict(summary=timing,warning='Incremental to cached baseline; profile preprocessing charged to each standalone profile variant; shared work actually computed once; concurrent wall times are not real-time throughput'))
    save(OUTPUT/'report.json',dict(summary=[r for r in summary if r['gate']=='primary'],original_gate_summary=[r for r in summary if r['gate']=='original'],
        decisions=decisions,chosen=chosen,diagnostics=diagnostic_summary,timing=timing,paired=paired,
        no_default_replacement=True,new_roots=0,independent_new_mothers=0,matched_draws_used=[0]))
    hashes={name:digest(OUTPUT/name) for name in ('scores.jsonl','summary.json','candidate-diagnostics.json','paired-comparison.json','timing.json','report.json')}
    save(RECORD/'science-validation.json',dict(passed=True,frames=1800,variant_candidate_sets=7200,primary_scores=7600,
        all_scores=13680,summary_rows=252,forced_K2_scores=3600,baseline_matches_previous=True,K0_K1_unchanged=True,
        original_K2_cost_not_worsened=True,maximum_cost_replay_error=cost_error,scalar_profile_points=profile_count,
        scalar_profile_not_converged=profile_not_converged,no_failed_fit_exclusion=True,active_mothers_disjoint=True,
        gate_calibration_uses_only280_native_singles=True,inputs_unchanged=True,freeze_unchanged=True,output_hashes=hashes,
        source_sha256=digest(__file__)))
    print('ALLOCATION_ANALYSIS_COMPLETE','chosen',chosen,'scores',len(scored),flush=True)


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.2,'figure.dpi':120})
    assert load(RECORD/'science-validation.json')['passed']
    report=load(OUTPUT/'report.json');summary=report['summary']
    folder=OUTPUT/'figures';folder.mkdir(exist_ok=True);hashes={}

    def finish(fig,number,title,legend_axis=None):
        fig.suptitle('Figure '+number+' | '+title,fontsize=13)
        bottom=.03
        if legend_axis is not None:
            h,l=legend_axis.get_legend_handles_labels();fig.legend(h,l,loc='lower center',ncol=3,fontsize=9);bottom=.10
        fig.tight_layout(rect=(0,bottom,1,.94))
        p=folder/(number+'.png');fig.savefig(p,dpi=155);plt.close(fig);hashes[number]=digest(p)

    fig,axes=plt.subplots(2,3,figsize=(14,8),sharey=True)
    for row,domain in enumerate(('native','matched')):
        for col,condition in enumerate(CONDITIONS):
            ax=axes[row,col]
            for j,method in enumerate(METHODS):
                rr=[next(r for r in summary if (r['domain'],r['method'],r['condition'],r['gap_ns'])==(domain,method,condition,g)) for g in GAPS]
                ax.bar(np.arange(3)+(j-2)*.16,[100*r['efficiency'] for r in rr],width=.15,color=COLORS[method],label=LABELS[method])
            ax.set(xticks=np.arange(3),xticklabels=['8','40','200'],ylim=(0,103),xlabel='True separation (ns)',
                   ylabel='Strict pair recovery (%)',title=f'({chr(97+row*3+col)}) {domain}: {condition}')
    finish(fig,'0004-82-1','Same gate calibration rule, 40 exposed pairs per point; no truth input',axes[0,0])

    fig,axes=plt.subplots(2,3,figsize=(14,8))
    controls=('all','weak','low','mid','high')
    for row,domain in enumerate(('native','matched')):
        for j,method in enumerate(METHODS):
            rr=[next(r for r in summary if r['domain']==domain and r['method']==method and r['kind']!='pair' and r['condition']==c) for c in controls]
            x=np.arange(5)+(j-2)*.13
            axes[row,0].bar(x,[100*r['extra_fraction'] for r in rr],width=.12,color=COLORS[method],label=LABELS[method])
            axes[row,1].plot(x,[100*r['efficiency'] for r in rr],'-o',color=COLORS[method],label=LABELS[method],ms=3)
            axes[row,2].plot(x[1:],[100*r['q_error_p16_p50_p84'][1] for r in rr[1:]],'-o',color=COLORS[method],ms=3)
        for col,label in enumerate(('Extra-output frames (%)','Strict success (%)','Matched median Q error (%)')):
            axes[row,col].set(xticks=np.arange(5),xticklabels=['empty','weak','low','mid','high'],ylabel=label,title=f'({chr(97+row*3+col)}) {domain}')
        axes[row,0].axhline(5.,ls='--',color='black',lw=.8)
        axes[row,1].axhline(95.,ls='--',color='black',lw=.8)
        for v in (-5.,5.):axes[row,2].axhline(v,ls='--',color='black',lw=.8)
    finish(fig,'0004-83-1','Empty and single controls: 80 frames per input/group; preserve false-split accounting',axes[0,0])

    fig,axes=plt.subplots(2,3,figsize=(14,8))
    diagnostics=report['diagnostics']
    for row,domain in enumerate(('native','matched')):
        for j,method in enumerate(VARIANTS):
            rr=[next(r for r in diagnostics if (r['domain'],r['method'],r['condition'],r['gap_ns'])==(domain,method,'high_weak',g)) for g in GAPS]
            x=np.arange(3)+(j-1.5)*.19
            axes[row,0].bar(x,[100*r['positive_cost_gain']/r['frames'] for r in rr],width=.18,color=COLORS[method],label=LABELS[method])
            axes[row,1].bar(x,[100*r['known_K_success']/r['frames'] for r in rr],width=.18,color=COLORS[method])
            timing=next(r for r in report['timing'] if (r['domain'],r['method'])==(domain,method))
            axes[row,2].bar(j,timing['incremental_seconds_mean'],color=COLORS[method])
        for col,label in enumerate(('Lower K=2 cost (%)','Forced K=2 strict success (%)','Mean incremental wall time (s)')):
            axes[row,col].set(ylabel=label,title=f'({chr(97+row*3+col)}) {domain}')
            if col<2:axes[row,col].set(xticks=np.arange(3),xticklabels=['8','40','200'],xlabel='True high_weak separation (ns)',ylim=(0,103))
            else:axes[row,col].set(xticks=np.arange(4),xticklabels=['C local','S local','C profile','S profile'])
    finish(fig,'0004-83-2','Cost improvement is not event recovery; C: component, S: total-share; cached baseline cost excluded',axes[0,0])

    candidate=load(OUTPUT/'candidate-diagnostics.json')['rows']
    fig,axes=plt.subplots(2,3,figsize=(14,8))
    for row,domain in enumerate(('native','matched')):
        for col,method in enumerate(('baseline','component_profile','allocation_profile')):
            ax=axes[row,col]
            rr=[r for r in candidate if (r['domain'],r['method'],r['condition'],r['gap_ns'])==(domain,method,'high_weak',40.)]
            for success,color,label in [(False,'#aa3377','Strict failure'),(True,'#009988','Strict success')]:
                selected=[r for r in rr if r['success']==success]
                ax.scatter([r['weak_delta_Q_pe']/1000 for r in selected],[r['strong_delta_Q_pe']/1000 for r in selected],s=22,color=color,alpha=.8,label=label)
            limit=max(1.,max(abs(r[k])/1000 for r in rr for k in ('weak_delta_Q_pe','strong_delta_Q_pe')))*1.1
            ax.plot([-limit,limit],[limit,-limit],'k:',lw=1,label='Unchanged total Q')
            ax.axhline(0.,color='gray',lw=.6);ax.axvline(0.,color='gray',lw=.6)
            ax.set(xlim=(-limit,limit),ylim=(-limit,limit),xlabel='Weak Q error (1000 pe)',ylabel='Strong Q error (1000 pe)',
                   title=f'({chr(97+row*3+col)}) {domain}: {LABELS[method]}')
    finish(fig,'0004-84-1','Forced K=2 at true 40 ns: complete charge-error ranges, including rejected and wrong-time candidates',axes[0,0])

    cases,q,t,base,predicted=read_data();response,noise=models();sample=response.sample_grid
    with np.load(OUTPUT/'observations.npz') as f:y=f['voltage_mV']
    gates=load(EXPERIMENT/'gate-calibrations.json')
    original=PROJECT/'jupyter/artifacts/alpha_pair_joint_confirmation_20261026/labels.truth.npz'
    with np.load(original) as f:actual_components=f['components_mV']
    examples=[]
    fig,axes=plt.subplots(2,3,figsize=(16,8.6))
    for row,gap in enumerate((40.,200.)):
        i=next(i for i,c in enumerate(cases) if c['role']=='evaluation' and c['domain']=='native' and c['condition']=='high_weak' and c['gap_ns']==gap and c['replicate']==0)
        c=cases[i];truth_curves=response.components(t[i],q[i]);weak=int(np.argmin(q[i]));strong=1-weak
        axes[row,0].plot(sample,y[i],color='#aaaaaa',lw=.7,label='Observed piled-up voltage')
        axes[row,0].plot(sample,truth_curves.sum(0)+response.baseline+noise.background_mV,color='black',ls='--',lw=.8,label='Conditional true-parameter mean')
        for col,j in ((1,strong),(2,weak)):
            axes[row,col].plot(sample,actual_components[c['original_index'],j],color='black',lw=.7,label='Native actual component (truth)')
            axes[row,col].plot(sample,truth_curves[j],color='#009988',ls='--',lw=.8,label='Conditional true-parameter mean')
        example=dict(index=i,native_index=c['native_index'],gap_ns=gap,methods={})
        for method in ('baseline','allocation_profile'):
            candidates=base[i] if method=='baseline' else predicted[i]['variants'][method]['candidates']
            fit=select(candidates,gates[method]);curves=response.components(fit['times_ns'],fit['charge_pe'])
            axes[row,0].plot(sample,curves.sum(0)+response.baseline+noise.background_mV,color=COLORS[method],lw=.9,label=LABELS[method])
            # Plot actually selected outputs. A selected single has no weak output.
            fitted_strong=curves[int(np.argmax(fit['charge_pe']))] if len(curves) else np.zeros_like(sample)
            fitted_weak=curves[int(np.argmin(fit['charge_pe']))] if len(curves)==2 else np.zeros_like(sample)
            axes[row,1].plot(sample,fitted_strong,color=COLORS[method],lw=.9,label=LABELS[method])
            axes[row,2].plot(sample,fitted_weak,color=COLORS[method],lw=.9,label=LABELS[method]+f' (K={len(curves)})')
            example['methods'][method]=dict(fit=fit,score=evaluate(fit,t[i],q[i],c))
        for col,title in enumerate(('Total voltage','Larger-Q selected output','Smaller-Q selected output')):
            ax=axes[row,col];ax.set(xlim=(0,2048),xlabel='Record time (ns)',ylabel='mV',title=f'({chr(97+row*3+col)}) {gap:g} ns: {title}')
            for onset in t[i]:ax.axvline(onset,color='black',ls=':',lw=.6,alpha=.3)
            ax.legend(fontsize=7)
        examples.append(example)
    finish(fig,'0004-85-1','Prespecified native pair 0, complete 0-2048 ns; no weak output when selected K=1')
    save(OUTPUT/'examples.json',examples)
    save(RECORD/'plot-truth-provenance.json',dict(file=str(original.relative_to(PROJECT)),sha256=digest(original),
        purpose='Display actual native avalanche components only; never opened by inference workers',examples=[r['index'] for r in examples]))
    assert list(hashes)==load(EXPERIMENT/'protocol.json')['figures']
    save(RECORD/'figure-delivery.json',dict(passed=True,figures=hashes,source_sha256=digest(__file__)))
    print('ALLOCATION_FIGURES_COMPLETE',len(hashes),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('analyze','figures'));args=parser.parse_args()
    if args.command=='analyze':analyze()
    else:figures()


if __name__=='__main__':main()
