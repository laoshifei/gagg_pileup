"""Auditable plots and checks for the joint likelihood development study."""
from __future__ import annotations

import argparse
import csv
import json

import numpy as np

from .alpha_pair_joint_data import PROJECT,EXPERIMENT,OUTPUT,RECORD
from .alpha_pair_joint_study import OLD,VARIANTS
from .alpha_pair_resolution_data import load,save,digest

LABELS=dict(observed_penalty='Old observed GLS',model_relative_residual='Old model GLS',old_ar='Old GLS, matched gate',
    old_mean_nll='Old mean + new NLL',joint_conserved_nll='Conserved joint mean + NLL')
COLORS=dict(observed_penalty='#cc3311',model_relative_residual='#888888',old_ar='#777777',
    old_mean_nll='#0077bb',joint_conserved_nll='#009988')


def rows(path):return [json.loads(s) for s in path.read_text().splitlines()]


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .alpha_pair_conditional import ConditionalResponse
    from .alpha_pair_joint_likelihood import ConservingResponse
    folder=OUTPUT/'figures';folder.mkdir(exist_ok=True)
    diagnosis=load(OUTPUT/'joint-diagnosis.json');calibration=load(RECORD/'joint-noise-calibration.json')['models']
    fig,axes=plt.subplots(1,2,figsize=(13,4.8))
    selected=[r for r in diagnosis['summary'] if r['fold']==1 and r['model']=='joint' and r['condition']=='high_weak' and r['gap_ns'] in (8.,40.,200.)]
    selected.sort(key=lambda r:r['gap_ns'])
    axes[0].bar([str(int(r['gap_ns'])) for r in selected],[100*r['cross_variance_fraction'] for r in selected],color='#4477aa')
    axes[0].axhline(0,color='k',lw=.7)
    axes[0].set(title='(a) Strong -> weak: cross contribution',xlabel='Separation (ns)',ylabel='2 tr(C12) / [tr(C1)+tr(C2)] (%)')
    for condition,kind in [('high','single'),('weak','single'),('high_weak','pair')]:
        r=next(r for r in diagnosis['summary'] if r['fold']==1 and r['model']=='joint' and r['kind']==kind and r['condition']==condition and (kind=='single' or r['gap_ns']==40.))
        axes[1].plot(np.arange(9)*2,r['innovation_lag_correlation'],'o-',label=condition)
    axes[1].set(title='(b) Residual correlation after electronics filter',xlabel='Sample lag (ns)',ylabel='Correlation');axes[1].legend()
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('Figure 0004-58-1 | Held-out component diagnostics; truth used only here')
    fig.tight_layout();fig.savefig(folder/'0004-58-1.png',dpi=150);plt.close(fig)

    old=ConditionalResponse.load(EXPERIMENT/'response.npz');new=ConservingResponse.load(EXPERIMENT/'joint-response.npz')
    grid=np.arange(-16.,6000.01,.25);qgrid=np.geomspace(550.,22000.,50)
    fig,axes=plt.subplots(1,2,figsize=(13,4.8));area_checks=[]
    for label,other,delay,style in [('isolated',0.,0.,'-'),('after 17000 pe, 40 ns',17000.,-40.,'--')]:
        for model,name,color in [(old,'Old mean','#cc3311'),(new,'Conserved mean','#009988')]:
            ratio=np.asarray([np.trapezoid(model.unit(grid,q,other,delay),grid)/model.base.area_per_pe for q in qgrid])
            axes[0].plot(qgrid,100*(ratio-1),style,color=color,label=f'{name}, {label}')
            area_checks.append(dict(model=name,condition=label,q_pe=qgrid.tolist(),relative_area=ratio.tolist()))
    axes[0].set(xscale='log',title='(a) Full-template area versus realized Q',xlabel='Candidate Q (pe)',ylabel='Area / expected area - 1 (%)')
    axes[0].legend(fontsize=8);axes[0].axhline(0,color='k',lw=.6)
    labels=['weak','low','mid','high'];x=np.arange(4)
    for name,offset in [('old_mean_nll',-.18),('joint_conserved_nll',.18)]:
        fractions=[1-r['empirical_conditioning_fraction'] for r in calibration[name]['charge_groups']]
        axes[1].bar(x+offset,fractions,width=.36,color=COLORS[name],label=LABELS[name])
    axes[1].set(xticks=x,xticklabels=labels,yscale='log',title='(b) Variance of integrated component residual',ylabel='Measured / unconditioned short-range variance')
    axes[1].legend(fontsize=8)
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('Figure 0004-59-1 | Charge conservation: a mean-model correction, not free lifetimes')
    fig.tight_layout();fig.savefig(folder/'0004-59-1.png',dpi=150);plt.close(fig)
    save(OUTPUT/'area-diagnosis.json',dict(rows=area_checks,numerical_grid_ns=.25,
        maximum_conserved_area_error=max(abs(v-1) for r in area_checks if r['model']=='Conserved mean' for v in r['relative_area'])))

    report=load(OUTPUT/'development-report.json');summary=report['summary']
    legacy=[r for r in load(OLD/'test/report.json')['summary'] if r['method'] in ('observed_penalty','model_relative_residual') and (r['kind']!='pair' or r['gap_ns'] in (8.,40.,200.))]
    plot_rows=legacy+summary
    methods=['observed_penalty','model_relative_residual','old_mean_nll','joint_conserved_nll']
    conditions=['low_low','mid_mid','high_high','low_high','high_low','low_weak','weak_low','high_weak','weak_high','weak_weak']
    fig,axes=plt.subplots(2,5,figsize=(18,8),sharex=True,sharey=True)
    for ax,condition,panel in zip(axes.flat,conditions,'abcdefghij'):
        for name in methods:
            rr=sorted([r for r in plot_rows if r['method']==name and r['condition']==condition],key=lambda r:r['gap_ns'])
            ax.plot([r['gap_ns'] for r in rr],[100*r['efficiency'] for r in rr],marker='o',lw=1.4,
                color=COLORS[name],label=LABELS[name])
        ax.set(title=f'({panel}) {condition}',xscale='log',xticks=[8,40,200],xticklabels=['8','40','200'],ylim=(-3,103))
        ax.set_xlabel('Separation (ns)');ax.grid(alpha=.2)
    axes[0,0].set_ylabel('Strict pair success (%)');axes[1,0].set_ylabel('Strict pair success (%)')
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=4,fontsize=9)
    fig.suptitle('Figure 0004-60-1 | Exposed development; 40 pairs per point, reused mothers across points')
    fig.tight_layout(rect=(0,.06,1,.95));fig.savefig(folder/'0004-60-1.png',dpi=150);plt.close(fig)

    fig,axes=plt.subplots(1,3,figsize=(17,5));controls=['all','weak','low','mid','high']
    for j,name in enumerate(methods):
        rr=[next(r for r in plot_rows if r['method']==name and r['kind']!='pair' and r['condition']==c) for c in controls]
        x=np.arange(5)+(j-1.5)*.17
        axes[0].bar(x,[100*r['extra_fraction'] for r in rr],width=.16,label=LABELS[name],color=COLORS[name])
        axes[1].plot(x,[100*r['efficiency'] for r in rr],'o-',color=COLORS[name],label=LABELS[name])
        axes[2].plot(x[1:],[100*r['q_error_p16_p50_p84'][1] for r in rr[1:]],'o-',color=COLORS[name],label=LABELS[name])
    for ax,title,level,ylabel in zip(axes,['(a) Extra-output frames','(b) Strict success','(c) Matched charge bias'],[5,95,5],['Frames (%)','Success (%)','Median Q error (%)']):
        ax.axhline(level,color='k',ls='--',lw=.8);ax.set(title=title,xticks=np.arange(5),xticklabels=controls,ylabel=ylabel);ax.grid(alpha=.2)
    axes[2].axhline(-5,color='k',ls='--',lw=.8)
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=4,fontsize=9)
    fig.suptitle('Figure 0004-60-2 | Empty and native-single controls; 80 frames in each group')
    fig.tight_layout(rect=(0,.06,1,.95));fig.savefig(folder/'0004-60-2.png',dpi=150);plt.close(fig)

    cases=load(OUTPUT/'labels.cases.json');scored=rows(OUTPUT/'scored-development.jsonl')
    lookup={(r['method'],r['index']):r for r in scored}
    with np.load(OUTPUT/'observations.npz') as f:y=f['voltage_mV']
    with np.load(OLD/'test/labels.truth.npz') as f:truth=f['components_mV']
    fig,axes=plt.subplots(3,2,figsize=(16,11));examples=[]
    for col,gap in enumerate([40.,200.]):
        i=next(i for i,c in enumerate(cases) if c['role']=='development' and c['condition']=='high_weak' and c['gap_ns']==gap and c['replicate']==0)
        actual=truth[cases[i]['source_index']];grid=old.sample_grid
        axes[0,col].plot(grid,y[i]-old.baseline,'k',lw=.7,label='Observed total')
        axes[1,col].plot(grid,actual[1],'k',lw=.7,label='True weak ancestry component')
        example=dict(index=i,case=cases[i],methods={})
        for name,model in [('old_ar',old),('old_mean_nll',old),('joint_conserved_nll',new)]:
            row=lookup[name,i];fit=row['fit'];curves=model.components(fit['times_ns'],fit['charge_pe']);total=curves.sum(0)
            axes[0,col].plot(grid,total,color=COLORS[name],label=f"{LABELS[name]}, K={len(curves)}")
            if len(curves)==2:axes[1,col].plot(grid,curves[1],color=COLORS[name],label=LABELS[name]+'; later output')
            axes[2,col].plot(grid,y[i]-model.baseline-total,color=COLORS[name],alpha=.7,label=LABELS[name])
            example['methods'][name]=dict(fit=fit,success=row['success'],true_q=row['true_q'],true_t=row['true_t'])
        for k in range(3):
            axes[k,col].set(xlim=(0,2048),xlabel='Full local record time (ns)',ylabel='mV',title=f'({chr(97+2*k+col)}) high_weak, {gap:g} ns, replicate 0')
            axes[k,col].grid(alpha=.2);axes[k,col].legend(fontsize=8)
        examples.append(example)
    fig.suptitle('Figure 0004-61-1 | Predetermined complete records; fits are blind, truth is a plotting reference')
    fig.tight_layout();fig.savefig(folder/'0004-61-1.png',dpi=150);plt.close(fig)
    save(OUTPUT/'examples.json',examples)
    with (OUTPUT/'development-summary.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(summary[0]));writer.writeheader();writer.writerows(summary)
    print('JOINT_FIGURES_COMPLETE',flush=True)


def validate():
    from .alpha_pair_joint_likelihood import select as new_select,calibrate_gate
    from .alpha_pair_null_gate import select as old_select,calibrate as old_calibrate
    from .alpha_pair_conditional_study import evaluate_fit
    from .alpha_pair_null_report import summarize
    for freeze in ('development-freeze.json','development-freeze-v2.json'):
        for p,sha in load(EXPERIMENT/freeze)['files'].items():assert digest(EXPERIMENT/p)==sha,(freeze,p)
    ancestry=load(RECORD/'ancestry-correction.json')
    assert digest(PROJECT/ancestry['corrected_file'])==ancestry['corrected_sha256']
    assert len(ancestry['checks'])==560 and max(r['closure_error_mV'] for r in ancestry['checks'])<1e-8
    provenance=load(RECORD/'data-provenance.json')
    mother0,mother1=map(set,provenance['fold_mothers']);valid=set(provenance['validation_mothers'])
    assert len(mother0)==len(mother1)==280 and len(valid)==320
    assert not mother0&mother1 and not (mother0|mother1)&valid
    for p,sha in provenance['input_files'].items():assert digest(PROJECT/p)==sha,p
    for name,sha in provenance['array_hashes'].items():assert digest(OUTPUT/name)==sha,name
    for p,sha in load(RECORD/'inference-data-provenance.json')['files'].items():assert digest(PROJECT/p)==sha,p
    cases=load(OUTPUT/'labels.cases.json');gates=load(EXPERIMENT/'gate-calibrations.json');stored=load(OUTPUT/'development-report.json')
    with np.load(OUTPUT/'labels.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
    all_summary=[];count=0;diagnostics={};comparisons=[];all_scored={};known_rows=[];known_groups={}
    for name in [*VARIANTS,'old_ar']:
        pred=rows(OUTPUT/f'predictions.{name}.jsonl');byindex={r['index']:r for r in pred}
        assert len(pred)==len(byindex)==len(cases)==1880
        singles=[byindex[i]['candidates'] for i,c in enumerate(cases) if c['role']=='calibration']
        gate=(old_calibrate(singles,edges=[1500.,5000.,12500.],mode='relative_residual',alpha=.025,penalty=.001,require_birth_penalty=False)
              if name=='old_ar' else calibrate_gate(singles))
        assert gate==gates[name]
        groups={};scores={};not_converged=0;calls=[]
        for i,c in enumerate(cases):
            p=byindex[i]['candidates'];assert len(p) in (1,3)
            for k,candidate in enumerate(p):
                assert len(candidate['times_ns'])==len(candidate['charge_pe'])==k
                assert np.isfinite(candidate['cost'])
                if name in VARIANTS:
                    assert np.isclose(candidate['cost'],candidate['quadratic']+candidate['logdet'])
            if c['role']!='development':continue
            fit=(old_select if name=='old_ar' else new_select)(p,gate)
            row=dict(evaluate_fit(fit,t[i],q[i],c),method=name,kind=c['kind'],condition=c['condition'],gap_ns=c['gap_ns'],true_q=q[i].tolist(),true_t=t[i].tolist())
            scores[i]=row;groups.setdefault((c['kind'],c['condition'],c['gap_ns']),[]).append(row);count+=1
            if c['kind']=='pair' and c['condition']=='high_weak':
                # Scoring-side oracle count only: exactly reuse the blind K=2 candidate.
                assert len(p)==3
                known=evaluate_fit(p[2],t[i],q[i],c)
                result=dict(index=i,method=name,gap_ns=c['gap_ns'],blind_success=row['success'],
                    known_count_success=known['success'],matched_two=known['matched']==2,
                    recovered_two=known['recovered']==2,strict_recovered_two=known['strict_recovered']==2,
                    truth_q=q[i].tolist(),truth_t=t[i].tolist(),fit=p[2])
                known_rows.append(result);known_groups.setdefault((name,c['gap_ns']),[]).append(result)
            not_converged+=not fit.get('converged',True)
            if 'seconds' in byindex[i]:calls.append(byindex[i]['seconds'])
        all_scored[name]=scores;all_summary += [summarize(rr) for rr in groups.values()]
        diagnostics[name]=dict(selected_not_converged=not_converged,development_frames=1600,
            seconds_median=None if not calls else float(np.median(calls)),seconds_p95=None if not calls else float(np.quantile(calls,.95)))
    assert all_summary==stored['summary']
    for name in VARIANTS:
        groups={}
        for i,r in all_scored[name].items():
            if r['kind']=='pair':groups.setdefault((r['condition'],r['gap_ns']),[]).append(int(r['success'])-int(all_scored['old_ar'][i]['success']))
        comparisons += [dict(method=name,condition=k[0],gap_ns=k[1],wins=sum(x>0 for x in v),losses=sum(x<0 for x in v),ties=sum(x==0 for x in v)) for k,v in groups.items()]
    save(OUTPUT/'paired-development.json',dict(comparisons=comparisons,independent_performance=False))
    save(OUTPUT/'known-count-diagnosis.json',dict(truth_used='True K=2 at scoring only; no new fitting, no true times/charges as starts',
        summary=[dict(method=k[0],gap_ns=k[1],frames=len(v),blind_success=sum(r['blind_success'] for r in v),
            known_count_success=sum(r['known_count_success'] for r in v),matched_two=sum(r['matched_two'] for r in v),
            strict_recovered_two=sum(r['strict_recovered_two'] for r in v)) for k,v in known_groups.items()],rows=known_rows))
    save(RECORD/'science-validation.json',dict(passed=True,scores_replayed=count,summaries=len(all_summary),
        diagnostics=diagnostics,source_freezes_preserved=True,input_hashes_preserved=True,mothers_disjoint=True,
        native_ancestry_closure_max_mV=max(r['closure_error_mV'] for r in ancestry['checks']),
        conserved_area_error=load(OUTPUT/'area-diagnosis.json')['maximum_conserved_area_error']))
    print('JOINT_VALIDATED',count,diagnostics,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['figures','validate']);args=parser.parse_args()
    globals()[args.command]()
