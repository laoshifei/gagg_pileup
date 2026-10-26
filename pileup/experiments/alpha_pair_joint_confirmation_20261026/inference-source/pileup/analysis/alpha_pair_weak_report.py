"""Plots and auditable development summaries for the weak-pulse diagnosis."""
from __future__ import annotations

import argparse
import csv
import gzip
import json

import numpy as np

from . import alpha_pair_weak_diagnostic_study as diagnostic
from . import alpha_pair_empirical_study as empirical
from .alpha_pair_resolution_data import load, save, digest

LABELS = dict(observed_previous='Previous observed GLS', blind_previous='Previous model GLS',
    known_count_previous='True K', fixed_true_times='True K and t',
    fixed_true_times_ols='True K and t, OLS', truth_start_joint='True K/Q/t start',
    truth_start_joint_ols='True K/Q/t start, OLS', oracle_subtract_other='Exact other subtraction',
    diagonal='Empirical diagonal',rank16_half='Rank 16, fraction 0.5',
    rank16_full='Rank 16, fraction 0.9',rank32_full='Rank 32, fraction 0.9',old_ar='Old AR, matched gate')
COLORS=dict(observed_previous='#cc3311',blind_previous='#0077bb',known_count_previous='#888888',
    fixed_true_times='#ee7733',fixed_true_times_ols='#ddaa33',truth_start_joint='#009988',
    truth_start_joint_ols='#44aa99',oracle_subtract_other='#aa3377',
    diagonal='#ee7733',rank16_half='#009988',rank16_full='#aa3377',rank32_full='#33bbee',old_ar='#777777')


def read_rows(path):
    return [json.loads(s) for s in path.read_text().splitlines()]


def figures_diagnostic():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from .alpha_pair_conditional import ConditionalResponse
    folder=diagnostic.OUTPUT/'figures';folder.mkdir(exist_ok=True)
    summary=load(diagnostic.OUTPUT/'report.json')['summary']
    rows=read_rows(diagnostic.OUTPUT/'scores.jsonl')
    methods=['blind_previous','known_count_previous','fixed_true_times',
             'truth_start_joint','truth_start_joint_ols','oracle_subtract_other']
    short=['Blind','True K','True K,t','True start\nGLS','True start\nOLS','Subtract\ntrue other']
    fig,axes=plt.subplots(1,3,figsize=(16,5),sharey=True)
    for ax,gap,panel in zip(axes,(8.,40.,200.),'abc'):
        for j,method in enumerate(methods):
            r=next(r for r in summary if r['condition']=='high_weak' and r['gap_ns']==gap and r['method']==method)
            ax.bar(j,100*r['efficiency'],color=COLORS[method],alpha=.8)
            lo,hi=100*np.asarray(r['success_wilson95']);value=100*r['efficiency']
            ax.errorbar(j,value,yerr=[[value-lo],[hi-value]],fmt='none',color='k',capsize=3)
            ax.text(j,max(10.,value)+3,f"{r['success']}/{r['eligible_frames']}",ha='center',fontsize=9)
        ax.set(xticks=range(6),xticklabels=short,title=f'({panel}) Strong -> weak, {gap:g} ns',ylim=(0,116))
        ax.tick_params(axis='x',labelsize=8);ax.grid(axis='y',alpha=.2)
    axes[0].set_ylabel('Strict pair success (%)')
    fig.suptitle('Figure 0004-53-1 | Exposed development; only the first bar is blind')
    fig.tight_layout();fig.savefig(folder/'0004-53-1.png',dpi=150);plt.close(fig)

    cases=load(diagnostic.INPUT/'labels.cases.json')
    diagnostics={r['index']:r['diagnostics'] for r in read_rows(diagnostic.OUTPUT/'diagnostics.jsonl')}
    response=ConditionalResponse.load(diagnostic.EXPERIMENT/'response.npz')
    with np.load(diagnostic.INPUT/'observations.npz') as f:y=f['voltage_mV']
    with np.load(diagnostic.INPUT/'labels.truth.npz') as f:
        q=f['charge_pe'];t=f['times_ns'];components=f['components_mV']
    fig,axes=plt.subplots(3,2,figsize=(16,11));examples=[]
    for col,spec in enumerate(load(diagnostic.EXPERIMENT/'protocol.json')['examples']):
        i=next(i for i,c in enumerate(cases) if all(c[k]==v for k,v in spec.items()))
        d=diagnostics[i];grid=response.sample_grid;observed=y[i]-response.baseline
        known=d['known_count_previous'];joint=d['truth_start_joint'];oracle=d['oracle_subtract_other']
        known_model=response.components(known['times_ns'],known['charge_pe'])
        joint_model=response.components(joint['times_ns'],joint['charge_pe'])
        true_model=response.components(t[i],q[i]);oracle_model=response.components(oracle['times_ns'],oracle['charge_pe'])
        axes[0,col].plot(grid,observed,'k',lw=.8,label='Observed total')
        axes[0,col].plot(grid,known_model.sum(0),'--',label='Known-K fit')
        axes[0,col].plot(grid,true_model.sum(0),alpha=.7,label='Model at true Q/t')
        axes[1,col].plot(grid,components[i,1],color='k',lw=.8,label='True weak ancestry waveform')
        axes[1,col].plot(grid,known_model[1],color=COLORS['known_count_previous'],label='Known-K weak fit')
        axes[1,col].plot(grid,joint_model[1],color=COLORS['truth_start_joint'],label='True-start weak fit')
        axes[1,col].plot(grid,oracle_model[1],color=COLORS['oracle_subtract_other'],ls='--',label='After exact other subtraction')
        axes[2,col].plot(grid,observed-true_model.sum(0),label='Residual at true Q/t',alpha=.7)
        axes[2,col].plot(grid,observed-known_model.sum(0),label='Known-K fit residual',alpha=.7)
        for k in range(3):
            ax=axes[k,col];ax.set(xlim=(0,2048),xlabel='Full local record time (ns)',ylabel='mV',
                title=f'({chr(97+2*k+col)}) high_weak, {spec["gap_ns"]:g} ns; replicate 0')
            ax.grid(alpha=.2);ax.legend(fontsize=8)
        examples.append(dict(index=i,case=cases[i],true_q=q[i].tolist(),true_t=t[i].tolist(),diagnostics=d))
    fig.suptitle('Figure 0004-54-1 | Predetermined full records; all fitted curves here use truth assistance')
    fig.tight_layout();fig.savefig(folder/'0004-54-1.png',dpi=150);plt.close(fig)
    save(diagnostic.OUTPUT/'examples.json',examples)

    fig,axes=plt.subplots(1,3,figsize=(15,4.8));cost_counts=[]
    all_ratios=[d['truth_start_joint']['cost']/d['exact_truth_cost']
        for i,d in diagnostics.items() if cases[i]['condition']=='high_weak']
    for ax,gap,panel in zip(axes,(8.,40.,200.),'abc'):
        selected=[r for r in rows if r['condition']=='high_weak' and r['gap_ns']==gap and r['method']=='truth_start_joint']
        values=[]
        for r in selected:
            i=r['index'];d=diagnostics[i];fit=d['truth_start_joint']
            values.append([fit['times_ns'][1]-t[i,1],100*(fit['charge_pe'][1]/q[i,1]-1.),
                           fit['cost']/d['exact_truth_cost']])
        v=np.asarray(values);tol=min(10.,gap/4.)
        ax.add_patch(Rectangle((-tol,-20),2*tol,40,facecolor='#ccddcc',alpha=.45,label='Weak Q/t tolerance'))
        points=ax.scatter(v[:,0],v[:,1],c=v[:,2],vmin=min(all_ratios),vmax=1,cmap='viridis',edgecolor='k',lw=.3)
        ax.axhline(0,color='k',lw=.5);ax.axvline(0,color='k',lw=.5)
        ax.set(title=f'({panel}) {gap:g} ns, true-start GLS',xlabel='Weak timing error (ns)',ylabel='Weak Q error (%)')
        if np.max(np.abs(v[:,1]))>200:
            ax.set_yscale('symlog',linthresh=50)
        ax.grid(alpha=.2)
        failed=[r for r in selected if not r['success']]
        cost_counts.append(dict(gap_ns=gap,frames=len(selected),failed=len(failed),
            failed_cost_below_truth=sum(diagnostics[r['index']]['truth_start_joint']['cost']<diagnostics[r['index']]['exact_truth_cost']-1e-6 for r in failed)))
    fig.colorbar(points,ax=axes.tolist(),label='Fitted cost / cost at true parameters',fraction=.02,pad=.02)
    fig.suptitle('Figure 0004-54-2 | Moving away from truth can lower the fixed GLS objective')
    fig.subplots_adjust(top=.84,bottom=.15,right=.88,wspace=.3)
    fig.savefig(folder/'0004-54-2.png',dpi=150);plt.close(fig)
    save(diagnostic.OUTPUT/'cost-diagnosis.json',dict(counts=cost_counts,
        interpretation='A lower sample residual does not imply correct Q/t; no theorem of intrinsic non-identifiability'))
    print('DIAGNOSTIC_FIGURES_COMPLETE',flush=True)


def assemble_empirical():
    from .alpha_pair_conditional_study import evaluate_fit
    from .alpha_pair_null_gate import select
    from .alpha_pair_null_report import summarize
    cases=load(empirical.OUTPUT/'labels.cases.json')
    with np.load(empirical.OUTPUT/'labels.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
    calibrations=load(empirical.EXPERIMENT/'gate-calibrations.json')
    stored=load(empirical.OUTPUT/'development-report.json')
    scores=[];summaries=[];groups={}
    for method in [*empirical.VARIANTS,'old_ar']:
        rows=read_rows(empirical.OUTPUT/f'predictions.{method}.jsonl')
        predictions={r['index']:r['candidates'] for r in rows}
        assert len(rows)==len(predictions)==len(cases)
        for i,c in enumerate(cases):
            if c['role']!='development':continue
            fit=select(predictions[i],calibrations[method])
            score=dict(evaluate_fit(fit,t[i],q[i],c),index=i,source_index=c['source_index'],method=method,
                kind=c['kind'],condition=c['condition'],gap_ns=c['gap_ns'],true_q=q[i].tolist(),true_t=t[i].tolist())
            scores.append(score);groups.setdefault((method,c['kind'],c['condition'],c['gap_ns']),[]).append(score)
    summaries=[summarize(rows) for rows in groups.values()]
    assert summaries==stored['summary'],'Independent score replay disagrees with development selection'
    with gzip.open(empirical.OUTPUT/'development-scores.jsonl.gz','wt') as stream:
        for r in scores:stream.write(json.dumps(r,allow_nan=False)+'\n')
    with (empirical.OUTPUT/'development-summary.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(summaries[0]));writer.writeheader();writer.writerows(summaries)
    comparisons=[]
    for method in empirical.VARIANTS:
        for (m,kind,condition,gap),new in groups.items():
            if m!=method or kind!='pair':continue
            old=groups[('old_ar',kind,condition,gap)]
            assert [r['index'] for r in new]==[r['index'] for r in old]
            change=[int(b['success'])-int(a['success']) for a,b in zip(old,new) if b['target']==2]
            comparisons.append(dict(method=method,baseline='old_ar',condition=condition,gap_ns=gap,
                frames=len(change),wins=sum(v>0 for v in change),losses=sum(v<0 for v in change),
                mean_change=float(np.mean(change))))
    save(empirical.OUTPUT/'paired-development.json',dict(comparisons=comparisons,
        scope='Descriptive exposed-development paired comparisons; no pooling across reused parents'))
    return stored


def figures_empirical():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .alpha_pair_empirical_noise import EmpiricalNoise
    folder=empirical.OUTPUT/'figures';folder.mkdir(exist_ok=True)
    noise=EmpiricalNoise.load(empirical.EXPERIMENT/'empirical-noise.npz')
    fig,axes=plt.subplots(1,2,figsize=(13,5))
    for j,name in enumerate(('weak','low','mid','high')):
        axes[0].plot(noise.grid,np.sqrt(noise.variance[j]),label=f'{name}, N={noise.counts[j]}')
    axes[0].set(xlim=(-64,800),yscale='log',xlabel='Time relative to true first photon (ns)',
                ylabel='Residual standard deviation (mV)',title='(a) Isolated training pulses')
    axes[0].legend();axes[0].grid(alpha=.2)
    ids=np.flatnonzero((noise.grid>=0)&(noise.grid<=400));modes=noise.modes[-1,ids,:16]*np.sqrt(.9)
    diagonal=np.maximum(noise.variance[-1,ids]-(modes*modes).sum(1),noise.noise_floor)
    cov=np.diag(diagonal)+modes@modes.T
    correlation=cov/np.sqrt(np.diag(cov)[:,None]*np.diag(cov)[None,:])
    img=axes[1].imshow(correlation,origin='lower',extent=[0,400,0,400],vmin=-1,vmax=1,cmap='RdBu_r')
    axes[1].set(title='(b) High-Q bin: rank 16, fraction 0.9',xlabel='Relative time (ns)',ylabel='Relative time (ns)')
    fig.colorbar(img,ax=axes[1],label='Regularized correlation')
    fig.suptitle('Figure 0004-55-1 | Empirical covariance calibration; mean response unchanged')
    fig.tight_layout();fig.savefig(folder/'0004-55-1.png',dpi=150);plt.close(fig)

    report=load(empirical.OUTPUT/'development-report.json');summary=report['summary']
    old=load(diagnostic.OUTPUT/'report.json')['summary']
    methods=['observed_previous','blind_previous','old_ar',*empirical.VARIANTS]
    combined=old+summary
    # Preserve source condition order; do not invent a missing combination.
    conditions=list(dict.fromkeys(r['condition'] for r in summary if r['kind']=='pair'))
    assert len(conditions)==10
    fig,axes=plt.subplots(2,5,figsize=(19,8),sharex=True,sharey=True)
    for ax,condition,panel in zip(axes.flat,conditions,'abcdefghij'):
        for method in methods:
            points=sorted([r for r in combined if r['method']==method and r['kind']=='pair' and r['condition']==condition],key=lambda r:r['gap_ns'])
            ax.plot([r['gap_ns'] for r in points],[100*r['efficiency'] for r in points],
                'o-' if method in empirical.VARIANTS else 'x--',color=COLORS[method],label=LABELS[method],ms=4)
        ax.set(title=f'({panel}) {condition}',xscale='log',ylim=(-3,103),xticks=[8,40,200],xticklabels=['8','40','200'])
        ax.grid(alpha=.2);ax.set_xlabel('First-photon separation (ns)');ax.set_ylabel('Strict pair success (%)')
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=4,fontsize=9)
    fig.suptitle('Figure 0004-56-1 | All ten combinations, 40 generated pairs per point; development only')
    fig.tight_layout(rect=(0,.09,1,.95));fig.savefig(folder/'0004-56-1.png',dpi=150);plt.close(fig)

    fig,axes=plt.subplots(1,3,figsize=(16,5))
    categories=['all','weak','low','mid','high']
    for j,method in enumerate(methods):
        selected=[next(r for r in combined if r['method']==method and r['kind']!='pair' and r['condition']==c) for c in categories]
        x=np.arange(5)+(j-3)*.10
        axes[0].plot(x,[100*r['extra_fraction'] for r in selected],'o',color=COLORS[method],label=LABELS[method],ms=4)
        axes[1].plot(x,[100*r['efficiency'] for r in selected],'o',color=COLORS[method],ms=4)
        axes[2].plot(x[1:],[100*r['q_error_p16_p50_p84'][1] for r in selected[1:]],'o',color=COLORS[method],ms=4)
    axes[0].axhline(5,ls='--',color='k',lw=.8);axes[1].axhline(95,ls='--',color='k',lw=.8)
    for bound in (-5,5):axes[2].axhline(bound,ls='--',color='k',lw=.8)
    for ax,title,ylabel in zip(axes,['(a) Extra-output frames','(b) Strict single/empty success','(c) Matched Q bias'],
                              ['Fraction (%)','Success (%)','Median relative error (%)']):
        ax.set(title=title,ylabel=ylabel,xticks=range(5),xticklabels=categories);ax.grid(alpha=.2)
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=4,fontsize=9)
    fig.suptitle('Figure 0004-56-2 | Control tradeoffs; N=80 per group; dashed lines are selection limits')
    fig.tight_layout(rect=(0,.12,1,.94));fig.savefig(folder/'0004-56-2.png',dpi=150);plt.close(fig)
    print('EMPIRICAL_FIGURES_COMPLETE',flush=True)


def validate():
    from .alpha_pair_conditional_study import evaluate_fit
    diagnostic.check_freeze()
    freeze=load(empirical.EXPERIMENT/'development-freeze.json')
    for path,sha in freeze['files'].items():assert digest(empirical.EXPERIMENT/path)==sha,path
    expected=load(diagnostic.EXPERIMENT/'protocol.json')['indices']
    rows=read_rows(diagnostic.OUTPUT/'diagnostics.jsonl')
    assert len(rows)==1600 and sorted(r['index'] for r in rows)==expected
    cases=load(diagnostic.INPUT/'labels.cases.json')
    with np.load(diagnostic.INPUT/'labels.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
    scored={(r['index'],r['method']):r for r in read_rows(diagnostic.OUTPUT/'scores.jsonl')}
    checked=0
    for row in rows:
        i=row['index'];d=row['diagnostics']
        assert d['metric_replay_relative_error']<1e-7
        assert d['truth_start_joint']['cost']<=d['exact_truth_cost']*(1+1e-10)+1e-7
        for method in ('known_count_previous','fixed_true_times','fixed_true_times_ols',
                       'truth_start_joint','truth_start_joint_ols','oracle_subtract_other'):
            if method not in d:continue
            actual=evaluate_fit(d[method],t[i],q[i],cases[i])
            assert all(scored[(i,method)][k]==v for k,v in actual.items()),(i,method)
            checked+=1
    provenance=load(empirical.RECORD/'covariance-calibration.json')
    assert len(provenance['training_seeds'])==len(provenance['gate_seeds'])==280
    assert not set(provenance['training_seeds'])&set(provenance['gate_seeds'])
    validation=load(empirical.OLD/'test/native-records.json')['test_single']
    assert not (set(provenance['training_seeds'])|set(provenance['gate_seeds']))&{r['seed'] for r in validation}
    waveform_hash=digest(empirical.OUTPUT/'observations.npz')
    for method in empirical.VARIANTS:
        run=load(empirical.RECORD/(method+'-inference.json'))
        assert run['completed'] and run['requested']==1880 and run['truth_inputs']==[]
        assert run['worker_labels_blocked'] and run['waveform_sha256']==waveform_hash
        assert run['prediction_sha256']==digest(empirical.OUTPUT/f'predictions.{method}.jsonl')
    for path,sha in load(empirical.RECORD/'data-provenance.json')['input_files'].items():
        assert digest(empirical.PROJECT/path)==sha,path
    assemble_empirical()
    inputs=load(diagnostic.RECORD/'input-provenance.json')['files']
    for path,sha in inputs.items():assert digest(diagnostic.PROJECT/path)==sha,path
    save(empirical.RECORD/'science-validation.json',dict(passed=True,diagnostic_frames=1600,
        diagnostic_scores_rechecked=checked,empirical_blind_frames=4*1880,empirical_score_rows=8000,
        calibration_gate_validation_parents_disjoint=True,new_native_roots=0,
        development_freeze_files=len(freeze['files']),old_inputs_unchanged=True,
        scope='Development only; no independent test generated or performance certified',
        result_hashes={str(p.relative_to(empirical.PROJECT)):digest(p) for p in
            (diagnostic.OUTPUT/'report.json',empirical.OUTPUT/'development-report.json',empirical.EXPERIMENT/'development-selection.json')}))
    print('WEAK_DEVELOPMENT_VALIDATED',checked,'diagnostic and 8000 blind score rows',flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['diagnostic-figures','empirical-figures','validate'])
    args=parser.parse_args()
    if args.command=='diagnostic-figures':figures_diagnostic()
    elif args.command=='empirical-figures':figures_empirical()
    else:validate()


if __name__=='__main__':main()
