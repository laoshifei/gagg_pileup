"""Post-inference evaluation, paired comparisons and numbered local figures."""
from __future__ import annotations

from collections import defaultdict
import csv
import gzip
import json
from pathlib import Path

import numpy as np

from . import alpha_pair_conditional_study as study
from .alpha_pair_conditional import ConditionalResponse,select_count,selective_refit
from .alpha_pair_fitting import PairFitter
from .alpha_pair_resolution_data import load,save,digest
from .alpha_pair_resolution_report import wilson,quantiles


def assemble(dataset):
    c=study.cache();seed=load(c/(dataset+'.seeds.json'))
    pred={r['index']:r for r in map(json.loads,(c/(dataset+'.predictions.jsonl')).read_text().splitlines())}
    cases=load(c/(dataset+'.cases.json'))
    assert set(pred)==set(range(len(cases)))
    penalty_path=study.EXPERIMENT/'fitting-selection.json'
    oldpenalty=load(penalty_path if penalty_path.exists() else study.PARENT/'fitting-selection.json')['penalty']
    selection=load(study.EXPERIMENT/'selection.json')
    methods={k:seed[k] for k in ('ml_1','ml_2','old_hybrid_1','old_hybrid_2')}
    methods['blind_fit']=[PairFitter.choose(pred[i]['original_candidates'],oldpenalty) for i in range(len(cases))]
    methods['conditional_fit']=[select_count(pred[i]['candidates'],selection['penalty']) for i in range(len(cases))]
    methods['conditional_fit_only']=[select_count(pred[i]['fitting_only_candidates'],selection['fitting_only_penalty']) for i in range(len(cases))]
    for j in (1,2):
        for mode,short in (('mean_only_hybrids','mean_hybrid'),('correlated_hybrids','correlated_hybrid')):
            methods[f'{short}_{j}']=[pred[i][mode][str(j)] for i in range(len(cases))]
        rule=selection['gates'][str(j)]['chosen']
        methods[f'selective_{j}']=[selective_refit(seed[f'ml_{j}'][i],pred[i][rule['mode']][str(j)],
            gap_threshold_ns=rule['gap_threshold_ns'],max_log_charge_change=rule['max_log_charge_change'])
            for i in range(len(cases))]
    return methods,cases


def report(dataset):
    c=study.cache();methods,cases=assemble(dataset)
    with np.load(c/(dataset+'.truth.npz')) as f:
        q=f['charge_pe'];t=f['times_ns'];valid=~f['train'] if dataset=='development' else np.ones(len(q),bool)
    rows=[];groups=defaultdict(list)
    for method,fits in methods.items():
        for i,case in enumerate(cases):
            if not valid[i]:continue
            row=dict(study.evaluate_fit(fits[i],t[i],q[i],case),index=i,method=method,
                kind=case['kind'],condition=case['condition'],gap_ns=case['gap_ns'],cluster=case['cluster'],
                true_q=q[i].tolist(),true_t=t[i].tolist(),
                returned_below500=sum(v<500. for v in fits[i]['charge_pe']))
            groups[(method,case['kind'],case['condition'],case['gap_ns'])].append(row);rows.append(row)
    summary=[]
    for (method,kind,condition,gap),members in groups.items():
        n=len(members);s=sum(r['success'] for r in members);extra=sum(r['extra']>0 for r in members)
        qe=[x for r in members for x in r['q_errors']];te=[x for r in members for x in r['t_errors']]
        summary.append(dict(method=method,kind=kind,condition=condition,gap_ns=gap,frames=n,
            success=s,efficiency=s/n,wilson95=wilson(s,n),extra_frames=extra,extra_fraction=extra/n,
            extra_wilson95=wilson(extra,n),budget_pass=extra/n<=.05,
            target=sum(r['target'] for r in members),recovered=sum(r['recovered'] for r in members),
            strict_recovered=sum(r['strict_recovered'] for r in members),
            missed=sum(r['missed'] for r in members),extra=sum(r['extra'] for r in members),
            split=sum(r['split'] for r in members),merge=sum(r['merge'] for r in members),
            returned_below500=sum(r['returned_below500'] for r in members),
            q_error_p16_p50_p84=quantiles(qe),q_error_std=float(np.std(qe,ddof=1)) if len(qe)>1 else None,
            time_error_p16_p50_p84_ns=quantiles(te),
            accepted_refits=sum(r['accepted'] is True for r in members),
            q1_p16_p50_p84_pe=quantiles([r['true_q'][0] for r in members]) if kind!='empty' else None,
            q2_p16_p50_p84_pe=quantiles([r['true_q'][1] for r in members]) if kind=='pair' else None))
    comparisons=[];rng=np.random.default_rng(262490000)
    pairs=(('ml_1','old_hybrid_1'),('ml_1','selective_1'),('old_hybrid_1','selective_1'),
           ('blind_fit','conditional_fit'),('ml_2','selective_2'))
    for left,right in pairs:
        for condition,gap in sorted({(r['condition'],r['gap_ns']) for r in rows if r['kind']=='pair'}):
            a=groups[(left,'pair',condition,gap)];b=groups[(right,'pair',condition,gap)]
            assert [r['index'] for r in a]==[r['index'] for r in b]
            d=np.array([int(y['success'])-int(x['success']) for x,y in zip(a,b)])
            boot=d[rng.integers(0,len(d),(2000,len(d)))].mean(1)
            comparisons.append(dict(left=left,right=right,condition=condition,gap_ns=gap,
                n=len(d),wins=int(sum(d>0)),losses=int(sum(d<0)),difference=float(d.mean()),
                paired_bootstrap95=np.percentile(boot,[2.5,97.5]).tolist(),
                unit='Independent parent pairs within this one condition/gap only; percentile intervals'))
    result=dict(dataset=dataset,summary=summary,comparisons=comparisons,
        protocol=load(study.EXPERIMENT/'protocol.json') if (study.EXPERIMENT/'protocol.json').exists() else None,
        selection=load(study.EXPERIMENT/'selection.json'),
        truth_sha256=digest(c/(dataset+'.truth.npz')),waveform_sha256=digest(c/(dataset+'.waveforms.npz')),
        event_rows=len(rows),independence='Same parents reused across gaps and some conditions. No pooled independent-exposure claim.')
    save(study.OUTPUT/(dataset+'-report.json'),result)
    with gzip.open(study.OUTPUT/(dataset+'-scores.jsonl.gz'),'wt') as f:
        for row in rows:f.write(json.dumps(row,allow_nan=False)+'\n')
    with (study.OUTPUT/(dataset+'-summary.csv')).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(summary[0]));writer.writeheader();writer.writerows(summary)
    return result


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    result=load(study.OUTPUT/'test-report.json');summary=result['summary']
    folder=study.OUTPUT/'figures';folder.mkdir(exist_ok=True)
    selected=('blind_fit','ml_1','old_hybrid_1','conditional_fit','selective_1','selective_2')
    labels=dict(blind_fit='Original blind fit',ml_1='Raw ML seed1',old_hybrid_1='Old ML + fixed fit',
                conditional_fit='Conditional GLS + ML starts',selective_1='Selective refit seed1',selective_2='Selective refit seed2')
    colors=dict(zip(selected,('#777777','#228833','#cc3311','#0077bb','#aa3377','#ee7733')))
    model=ConditionalResponse.load(study.EXPERIMENT/'response.npz')
    s=np.arange(0.,600.,2.);fig,axes=plt.subplots(1,3,figsize=(15,4))
    for q in (2700.,9500.,16000.):axes[0].plot(s,model.unit(s,q),label=f'{q:g} pe')
    axes[0].plot(s,np.interp(s,model.base.grid,model.base.kernel),'k--',label='Original template')
    for q,delay in ((0.,0.),(9500.,-8.),(9500.,-40.),(9500.,8.),(9500.,40.)):
        axes[1].plot(s,model.unit(s,9500.,q,delay),label='Single' if q==0 else f'Neighbor delay {delay:g} ns')
    axes[0].set(title='(a) Single-event charge dependence',xlabel='Time from first photon (ns)',ylabel='mV / pe')
    axes[1].set(title='(b) Inferred neighbor-state surrogate',xlabel='Time from first photon (ns)',ylabel='mV / pe')
    dev=load(study.OUTPUT/'development-report.json')['summary']
    for method in ('ml_1','old_hybrid_1','conditional_fit','selective_1'):
        data=sorted([r for r in dev if r['method']==method and r['condition']=='mid_mid'],key=lambda r:r['gap_ns'])
        axes[2].plot([r['gap_ns'] for r in data],[r['efficiency'] for r in data],'.-',label=labels[method])
    axes[2].set(title='(c) Development validation only',xlabel='Separation (ns)',ylabel='Strict pair recovery',xscale='log',ylim=(-.02,1.04))
    for ax in axes:ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Figure 0004-41-1 | Response calibration and development validation')
    fig.tight_layout();fig.savefig(folder/'0004-41-1.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(15,9))
    for ax,condition,letter in zip(axes.flat,study.load(study.EXPERIMENT/'protocol.json')['conditions'],'abcde'):
        for method in selected:
            data=sorted([r for r in summary if r['method']==method and r['condition']==condition],key=lambda r:r['gap_ns'])
            ax.plot([r['gap_ns'] for r in data],[r['efficiency'] for r in data],'.-',label=labels[method],color=colors[method])
        ax.axhline(.9,color='black',ls=':',lw=1);ax.set(title=f'({letter}) {condition}',xlabel='Separation (ns)',ylabel='Strict pair recovery',xscale='log',ylim=(-.03,1.04));ax.grid(alpha=.2)
    axes.flat[-1].axis('off');handles,legend=axes.flat[0].get_legend_handles_labels();axes.flat[-1].legend(handles,legend,loc='center')
    fig.suptitle('Figure 0004-42-1 | Independent optical parents, 40 pairs per point')
    fig.tight_layout();fig.savefig(folder/'0004-42-1.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(13,5));categories=[('empty','all'),('single','low'),('single','mid'),('single','high')]
    for mi,method in enumerate(selected):
        data=[next(r for r in summary if (r['method'],r['kind'],r['condition'])==(method,k,c)) for k,c in categories]
        x=np.arange(4)+(mi-2.5)*.09
        p=np.array([r['extra_fraction'] for r in data]);ci=np.array([r['extra_wilson95'] for r in data])
        axes[0].errorbar(x,p,yerr=np.array([p-ci[:,0],ci[:,1]-p]),fmt='.',label=labels[method],color=colors[method],capsize=2)
        axes[1].plot(x,[r['efficiency'] for r in data],'.-',label=labels[method],color=colors[method])
    axes[0].axhline(.05,color='black',ls=':');axes[0].set(title='(a) Extra-output frames (Wilson95)',ylabel='Fraction')
    axes[1].set(title='(b) Strict control success, no extras',ylabel='Fraction',ylim=(-.03,1.04))
    for ax in axes:ax.set_xticks(range(4),['Empty','Low single','Mid single','High single']);ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle('Figure 0004-42-2 | Empty and native single-event controls, n=80 each')
    fig.tight_layout();fig.savefig(folder/'0004-42-2.png',dpi=150);plt.close(fig)
    methods,cases=assemble('test');c=study.cache()
    with np.load(c/'test.waveforms.npz') as f:y=f['voltage_mV']
    with np.load(c/'test.truth.npz') as f:true=f['components_mV'];t=f['times_ns'];q=f['charge_pe']
    # Predefined indices, never selected by success/failure.
    examples=[next(i for i,r in enumerate(cases) if r['kind']=='pair' and r['condition']=='mid_mid' and r['gap_ns']==gap and r['replicate']==0) for gap in (8.,40.)]
    fig,axes=plt.subplots(3,2,figsize=(14,10));grid=np.arange(1024)*2.+.25
    detail=[]
    for col,i in enumerate(examples):
        fit=methods['selective_1'][i];curves=model.components(fit['times_ns'],fit['charge_pe'])
        old=methods['old_hybrid_1'][i];oldcurve=model.base.columns(grid,old['times_ns'])@np.asarray(old['charge_pe'])
        axes[0,col].plot(grid,y[i]-model.baseline,color='black',label='Observed total')
        axes[0,col].plot(grid,curves.sum(0),label='Selective output rendered by calibrated response')
        axes[0,col].plot(grid,oldcurve,ls='--',label='Old fixed-template refit')
        for j in range(2):
            axes[1,col].plot(grid,true[i,j],label=f'True component {j+1}')
            if j<len(curves):axes[1,col].plot(grid,curves[j],ls='--',label=f'Recovered component {j+1}')
        axes[2,col].plot(grid,y[i]-model.baseline-curves.sum(0),label='Calibrated-response residual')
        for row in range(3):
            axes[row,col].set(title=f'({"abcdef"[row*2+col]}) gap={cases[i]["gap_ns"]:g} ns',xlabel='Local record time (ns)',ylabel='mV',xlim=(0,2048));axes[row,col].legend(fontsize=8);axes[row,col].grid(alpha=.2)
        detail.append(dict(index=i,case=cases[i],true_t=t[i].tolist(),true_q=q[i].tolist(),fit=fit,
                           score=study.evaluate_fit(fit,t[i],q[i],cases[i])))
    fig.suptitle('Figure 0004-43-1 | Full records and individual components; fixed replicate 0')
    fig.tight_layout();fig.savefig(folder/'0004-43-1.png',dpi=150);plt.close(fig)
    save(study.OUTPUT/'examples.json',detail)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('stage',choices=('development','test','figures'))
    args=p.parse_args()
    if args.stage=='figures':figures()
    else:report(args.stage)
