"""Evaluate frozen residual-search methods, including unsuccessful conditions."""
from __future__ import annotations

import json
import gzip
from collections import Counter

import numpy as np

from . import alpha_pair_residual_study as study
from .alpha_pair_conditional import ConditionalResponse,select_count,selective_refit
from .alpha_pair_fitting import PairFitter
from .alpha_pair_residual import choose


LABELS=dict(blind_fit='Original blind fit',old_hybrid_1='Network + fixed fit',
    conditional_fit='Previous conditional count fit',selective_1='Previous GLS, seed1',
    regularized_1='Regularized GLS, same network K',regularized_count='Regularized existing candidates',
    residual_original='Residual search, original precision',residual_regularized='Residual search + regularization',
    residual_waveform='Residual search, waveform starts only')


def same_count_output(raw,refined):
    # The zero-only waveform gate skips all searches, including neural refits.
    # Keep the raw network result in this comparison if no same-K refit exists.
    # This does not change the new count-selection method's empty output.
    if len(raw['times_ns'])!=len(refined['times_ns']):
        return dict(raw,accepted_refit=False,refit_unavailable='Waveform zero-only gate')
    return selective_refit(raw,refined,gap_threshold_ns=0.,max_log_charge_change=10.)


def assemble(dataset):
    c=study.cache();seeds=study.load(c/(dataset+'.seeds.json'));cases=study.load(c/(dataset+'.cases.json'))
    pred={r['index']:r for r in map(json.loads,(c/(dataset+'.predictions.jsonl')).read_text().splitlines())}
    assert set(pred)==set(range(len(cases)))
    if dataset=='validation':
        old={r['index']:r for r in map(json.loads,(c/'validation.previous.jsonl').read_text().splitlines())}
    else:old={i:r['previous'] for i,r in pred.items()}
    original=study.load(study.EXPERIMENT/'selection.json');selection=study.load(study.EXPERIMENT/'residual-selection.json')
    chosen=selection['chosen_policy'];penalties=selection['penalties']
    methods={key:seeds[key] for key in ('ml_1','ml_2','old_hybrid_1','old_hybrid_2')}
    penalty=study.load(study.EXPERIMENT/'fitting-selection.json')['penalty']
    methods['blind_fit']=[PairFitter.choose(seeds['candidates'][i],penalty) for i in range(len(cases))]
    methods['conditional_fit']=[select_count(old[i]['candidates'],original['penalty']) for i in range(len(cases))]
    for j in (1,2):
        gate=original['gates'][str(j)]['chosen']
        methods[f'selective_{j}']=[selective_refit(seeds[f'ml_{j}'][i],old[i][gate['mode']][str(j)],
            gap_threshold_ns=gate['gap_threshold_ns'],max_log_charge_change=gate['max_log_charge_change']) for i in range(len(cases))]
        methods[f'regularized_{j}']=[same_count_output(seeds[f'ml_{j}'][i],
            pred[i]['policies'][chosen]['hybrids'][str(j)]) for i in range(len(cases))]
    for method,policy,variant in (('regularized_count',chosen,'existing'),('residual_original','original','expanded'),
            ('residual_regularized',chosen,'expanded'),('residual_waveform',chosen,'waveform_only')):
        methods[method]=[choose(pred[i]['policies'][policy][variant],penalties[policy][variant]['penalty']) for i in range(len(cases))]
    return methods,cases


def report(dataset):
    # Reuse the established scoring implementation and uncertainty convention.
    from . import alpha_pair_conditional_report as reference
    reference.study=study.adapt();reference.assemble=assemble
    result=reference.report(dataset)
    result['residual_selection']=study.load(study.EXPERIMENT/'residual-selection.json')
    result['methods_note']='14 workflows; two inherited neural weight seeds, no new neural training'
    with gzip.open(study.OUTPUT/(dataset+'-scores.jsonl.gz'),'rt') as stream:
        scores=[json.loads(line) for line in stream]
    # Within a fixed condition/gap the optical parent pairs are independent.
    # Do not pool reused parents across gaps into a spurious exposure count.
    groups={}
    for row in scores:
        if row['kind']=='pair':
            groups.setdefault((row['method'],row['condition'],row['gap_ns']),[]).append(row)
    contrasts=[('selective_1','regularized_1'),('regularized_1','residual_regularized'),
        ('regularized_count','residual_regularized'),('residual_original','residual_regularized'),
        ('residual_waveform','residual_regularized'),('conditional_fit','residual_regularized')]
    rng=np.random.default_rng(263990001);comparisons=[]
    points=sorted({(r['condition'],r['gap_ns']) for r in scores if r['kind']=='pair'})
    for left,right in contrasts:
        for condition,gap in points:
            a=groups[(left,condition,gap)];b=groups[(right,condition,gap)]
            assert [r['index'] for r in a]==[r['index'] for r in b]
            delta=np.array([int(y['success'])-int(x['success']) for x,y in zip(a,b)])
            boot=delta[rng.integers(0,len(delta),(2000,len(delta)))].mean(1)
            comparisons.append(dict(left=left,right=right,condition=condition,gap_ns=gap,
                n=len(delta),wins=int(sum(delta>0)),losses=int(sum(delta<0)),
                difference=float(delta.mean()),paired_bootstrap95=np.percentile(boot,[2.5,97.5]).tolist(),
                unit='One fixed condition/gap; no cross-gap pooling or simultaneous interval'))
    result['residual_comparisons']=comparisons
    study.save(study.OUTPUT/(dataset+'-report.json'),result)
    counts=[]
    for method in ('ml_1','ml_2','selective_1','selective_2','regularized_1','regularized_2',
                   'conditional_fit','regularized_count','residual_original','residual_regularized','residual_waveform'):
        rows=groups[(method,'high_low',200.)]
        histogram=Counter(r['count'] for r in rows)
        counts.append(dict(method=method,frames=len(rows),output_K0=histogram[0],output_K1=histogram[1],
            output_K2=histogram[2],strict_pair_success=sum(r['success'] for r in rows)))
    study.save(study.OUTPUT/(dataset+'-diagnostics.json'),dict(high_low_200ns_counts=counts,
        report_sha256=study.digest(study.OUTPUT/(dataset+'-report.json')),
        purpose='Frozen-output explanation only; no threshold optimization'))
    return result


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output=study.OUTPUT;folder=output/'figures';folder.mkdir(exist_ok=True)
    selected=study.load(study.EXPERIMENT/'residual-selection.json')
    summary=study.load(output/'test-report.json')['summary']
    fig,axes=plt.subplots(1,2,figsize=(13,4.5))
    for policy,params in study.POLICIES.items():
        rows=[r for r in selected['grid'] if r['policy']==policy and r['variant']=='expanded']
        rho=params['rho'] or study.load(study.EXPERIMENT/'calibration.json')['residual_calibration']['rho']
        reference=next(r for r in rows if r['penalty']==.003)
        axes[0].scatter([rho],[100*reference['controls']['low']['q_median_relative_error']],label=policy)
        axes[1].plot([max(g['extra_fraction'] for g in r['controls'].values()) for r in rows],
            [r['pair_success'] for r in rows],'.-',label=policy)
    axes[0].axhline(0,color='gray',lw=1);axes[0].axhspan(-5,5,color='gray',alpha=.12)
    axes[0].set(title='(a) Low-single charge bias, penalty=0.003',xlabel='Assumed AR correlation',ylabel='Median relative Q error (%)')
    chosen=selected['chosen'];x=max(g['extra_fraction'] for g in chosen['controls'].values())
    axes[1].plot([x],[chosen['pair_success']],'k*',ms=14,label='Frozen choice')
    axes[1].axvline(.05,color='gray',ls=':');axes[1].set(title='(b) Development trade-off',
        xlabel='Largest extra-output fraction among controls',ylabel='Strict pair success')
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Figure 0004-45-1 | Development only: precision and count selection')
    fig.tight_layout();fig.savefig(folder/'0004-45-1.png',dpi=150);plt.close(fig)
    methods=['old_hybrid_1','selective_1','conditional_fit','regularized_1','residual_regularized','residual_waveform']
    colors=dict(zip(methods,['#888888','#cc3311','#ee7733','#228833','#0077bb','#aa3377']))
    fig,axes=plt.subplots(2,3,figsize=(15,9))
    for ax,condition,letter in zip(axes.flat,['low_low','mid_mid','high_high','low_high','high_low'],'abcde'):
        for method in methods:
            rows=sorted([r for r in summary if r['condition']==condition and r['method']==method],key=lambda r:r['gap_ns'])
            ax.plot([r['gap_ns'] for r in rows],[r['efficiency'] for r in rows],'.-',label=LABELS[method],color=colors[method])
        ax.axhline(.9,color='black',ls=':',lw=1)
        ax.set(title=f'({letter}) {condition}',xlabel='Separation (ns)',ylabel='Strict pair success',xscale='log',ylim=(-.03,1.04));ax.grid(alpha=.2)
    axes.flat[-1].axis('off');h,l=axes.flat[0].get_legend_handles_labels();axes.flat[-1].legend(h,l,loc='center',fontsize=9)
    fig.suptitle('Figure 0004-46-1 | New independent optical parents, 40 pairs per point')
    fig.tight_layout();fig.savefig(folder/'0004-46-1.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(16,4.8));controls=['all','low','mid','high']
    displayed=['old_hybrid_1','selective_1','regularized_1','residual_regularized','residual_waveform']
    for mi,method in enumerate(displayed):
        rows=[next(r for r in summary if r['method']==method and r['kind']!='pair' and r['condition']==condition) for condition in controls]
        x=np.arange(3)+(mi-2)*.08
        med=np.array([r['q_error_p16_p50_p84'][1] for r in rows[1:]])*100
        q=np.array([r['q_error_p16_p50_p84'] for r in rows[1:]])*100
        axes[0].errorbar(x,med,yerr=[med-q[:,0],q[:,2]-med],fmt='.',capsize=2,label=LABELS[method],color=colors[method])
        x=np.arange(4)+(mi-2)*.08;p=np.array([r['extra_fraction'] for r in rows]);ci=np.array([r['extra_wilson95'] for r in rows])
        axes[1].errorbar(x,p,yerr=[p-ci[:,0],ci[:,1]-p],fmt='.',capsize=2,color=colors[method])
        axes[2].plot(x,[r['efficiency'] for r in rows],'.-',color=colors[method])
    axes[0].set_xticks(range(3),['Low','Mid','High']);axes[0].axhline(0,color='gray',lw=1)
    axes[0].set(title='(a) Matched-charge error: median, p16-p84',ylabel='Relative Q error (%)');axes[0].legend(fontsize=7)
    for ax in axes[1:]:ax.set_xticks(range(4),['Empty','Low','Mid','High'])
    axes[1].axhline(.05,color='gray',ls=':');axes[1].set(title='(b) Extra-output frames, Wilson95',ylabel='Fraction')
    axes[2].set(title='(c) Strict control success',ylabel='Fraction',ylim=(-.03,1.04))
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('Figure 0004-46-2 | Native single-event and empty controls, n=80 each')
    fig.tight_layout();fig.savefig(folder/'0004-46-2.png',dpi=150);plt.close(fig)
    fits,cases=assemble('test');c=study.cache();response=ConditionalResponse.load(study.EXPERIMENT/'response.npz')
    with np.load(c/'test.waveforms.npz') as f:y=f['voltage_mV']
    with np.load(c/'test.truth.npz') as f:q=f['charge_pe'];t=f['times_ns'];components=f['components_mV']
    examples=[next(i for i,r in enumerate(cases) if r['condition']=='high_low' and r['gap_ns']==gap and r['replicate']==0) for gap in (40.,200.)]
    fig,axes=plt.subplots(3,2,figsize=(14,10));details=[];grid=response.sample_grid
    for col,i in enumerate(examples):
        new=fits['residual_regularized'][i];old=fits['selective_1'][i]
        curves=response.components(new['times_ns'],new['charge_pe']);oldcurve=response.components(old['times_ns'],old['charge_pe']).sum(0)
        axes[0,col].plot(grid,y[i]-response.baseline,'k',label='Observed total')
        axes[0,col].plot(grid,oldcurve,'--',label='Previous same-K GLS')
        axes[0,col].plot(grid,curves.sum(0),label='Residual search + regularization')
        for j in range(2):
            axes[1,col].plot(grid,components[i,j],label=f'True component {j+1}')
            if j<len(curves):axes[1,col].plot(grid,curves[j],'--',label=f'Recovered component {j+1}')
        axes[2,col].plot(grid,y[i]-response.baseline-oldcurve,label='Previous residual')
        axes[2,col].plot(grid,y[i]-response.baseline-curves.sum(0),label='New residual')
        for row in range(3):
            axes[row,col].set(title=f'({"abcdef"[2*row+col]}) Strong then weak, gap={cases[i]["gap_ns"]:g} ns',
                xlabel='Local record time (ns)',ylabel='mV',xlim=(0,2048))
            axes[row,col].grid(alpha=.2);axes[row,col].legend(fontsize=8)
        details.append(dict(index=i,case=cases[i],true_q=q[i].tolist(),true_t=t[i].tolist(),previous=old,fit=new,
            score=study.previous.evaluate_fit(new,t[i],q[i],cases[i]),selection='Predetermined high_low replicate0; not chosen after viewing test results'))
    fig.suptitle('Figure 0004-47-1 | Full records and recovered components, fixed examples')
    fig.tight_layout();fig.savefig(folder/'0004-47-1.png',dpi=150);plt.close(fig)
    study.save(output/'examples.json',details)


if __name__=='__main__':
    import sys
    if sys.argv[1]=='figures':figures()
    else:report(sys.argv[1])
