"""Paired recovery, truth-conditioned residuals and figures for model closure."""
from __future__ import annotations

import argparse
import csv

import numpy as np

from .alpha_pair_model_closure_study import PROJECT, EXPERIMENT, OUTPUT, RECORD, POOL, CONDITIONS, models, rows, check_freeze
from .alpha_pair_resolution_data import load, save, digest

METHODS = ('native/native_gate', 'matched/native_gate', 'matched/matched_gate', 'native/matched_gate')
LABELS = {'native/native_gate':'Native, original gate', 'matched/native_gate':'Gaussian, original gate',
          'matched/matched_gate':'Gaussian, Gaussian gate', 'native/matched_gate':'Native, Gaussian gate'}
COLORS = dict(zip(METHODS, ('#555555','#009988','#aa4499','#0077bb')))


def read_data():
    data={}
    for domain in ('native','matched'):
        prefix='native-' if domain=='native' else ''
        with np.load(OUTPUT/(prefix+'observations.npz')) as f:y=f['voltage_mV']
        with np.load(OUTPUT/(prefix+'labels.truth.npz')) as f:q=f['charge_pe'];t=f['times_ns']
        cases=load(OUTPUT/('native.cases.json' if domain=='native' else 'labels.cases.json'))
        predictions={r['index']:r for r in rows(OUTPUT/(domain+'.predictions.jsonl'))}
        data[domain]=dict(y=y,q=q,t=t,cases=cases,predictions=predictions)
    return data


def analyze():
    from .alpha_pair_joint_likelihood import calibrate_gate, select
    from .alpha_pair_conditional_study import evaluate_fit
    from .alpha_pair_null_report import summarize
    from .alpha_pair_model_sampler import GaussianFactor, conditional_mean
    from .alpha_pair_profile_diagnostics import cost_at
    check_freeze();response,noise=models();data=read_data();protocol=load(EXPERIMENT/'protocol.json')
    provenance=load(RECORD/'input-provenance.json')
    for path,expected in provenance['files'].items():assert digest(PROJECT/path)==expected,path
    for filename,expected in provenance['prepared_files'].items():assert digest(OUTPUT/filename)==expected,filename
    groups=[set(v) for v in provenance['mother_groups'].values()]
    assert all(not groups[i]&groups[j] for i in range(len(groups)) for j in range(i))
    for record in ('numerical-tests','sampler-preflight','pilot-validation','matched-inference','diagnostic-run'):
        assert load(RECORD/(record+'.json'))['passed'],record
    generation=load(RECORD/'generation.json')
    assert generation['waveform_sha256']==digest(OUTPUT/'observations.npz')
    assert load(RECORD/'matched-inference.json')['prediction_sha256']==digest(OUTPUT/'matched.predictions.jsonl')
    for domain,expected in [('native',760),('matched',1800)]:
        assert set(data[domain]['predictions'])==set(range(expected))
    original_gate=load(EXPERIMENT/'baseline-gate.json')
    original_pool={r['index']:r['candidates'] for r in rows(POOL/'predictions.baseline.jsonl')}
    assert calibrate_gate([original_pool[i] for i in range(280)])==original_gate
    gates=dict(native_gate=original_gate,matched_gate=calibrate_gate([data['matched']['predictions'][i]['candidates'] for i in range(280)]))
    save(EXPERIMENT/'gate-comparison.json',gates)
    scores=[];summary=[];cost_error=0.;selected_nonconverged={};gate_excluded=0
    for method in METHODS:
        domain,gate_name=method.split('/');d=data[domain];grouped={};nonconverged=0
        for i,c in enumerate(d['cases']):
            if c['role']=='gate':gate_excluded+=1;continue
            candidates=d['predictions'][i]['candidates'];fit=select(candidates,gates[gate_name])
            value=dict(evaluate_fit(fit,d['t'][i],d['q'][i],c),method=method,domain=domain,gate=gate_name,
                index=i,native_index=c['native_index'],block=c['block'],draw=c['draw'],kind=c['kind'],condition=c['condition'],
                gap_ns=c['gap_ns'],fit=fit,true_q=d['q'][i].tolist(),true_t=d['t'][i].tolist())
            scores.append(value);grouped.setdefault((c['kind'],c['condition'],c['gap_ns']),[]).append(value)
            nonconverged+=not fit.get('converged',True)
        summary += [summarize(rr) for rr in grouped.values()]
        selected_nonconverged[method]=int(nonconverged)
    assert len(scores)==4560 and len(summary)==56 and gate_excluded==560
    for domain,d in data.items():
        for i,r in d['predictions'].items():
            for fit in r['candidates']:
                replay=cost_at(response,noise,d['y'][i],fit['times_ns'],fit['charge_pe'])['cost']
                error=abs(replay-fit['cost'])/max(1.,abs(replay));cost_error=max(cost_error,error)
                assert error<1e-7,(domain,i)
    old=load(POOL/'report.json')['summary']
    for s in summary:
        if s['method']!='native/native_gate':continue
        old_row=next(r for r in old if r['method']=='baseline' and r['dataset_role']=='latest'
                     and (r['kind'],r['condition'],r['gap_ns'])==(s['kind'],s['condition'],s['gap_ns']))
        for key in ('success','eligible_frames','extra_frames','q_error_p16_p50_p84','time_error_p16_p50_p84_ns'):
            assert s[key]==old_row[key],(key,s)
    with (OUTPUT/'scores.jsonl').open('w') as stream:
        import json
        for s in scores:stream.write(json.dumps(s,allow_nan=False)+'\n')
    with (OUTPUT/'summary.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(summary[0]));writer.writeheader();writer.writerows(summary)
    save(OUTPUT/'summary.json',summary)

    scored={(r['method'],r['index']):r for r in scores};diagnostic_rows=[];cost_diagnostics=[];diagnostic_cost_error=0.
    for domain in ('native','matched'):
        d=data[domain];diagnostics=rows(OUTPUT/(domain+'-diagnostics.truth.jsonl'))
        assert len(diagnostics)==(240 if domain=='native' else 480)
        assert len({r['index'] for r in diagnostics})==len(diagnostics)
        for diagnostic in diagnostics:
            i=diagnostic['index'];c=d['cases'][i];baseline=d['predictions'][i]['candidates'][2]
            assert baseline==diagnostic['baseline_double']
            fits={'blind_native_gate':scored[domain+'/native_gate',i]['fit'],'known_K':baseline,
                  'true_start':diagnostic['truth_start'],'fixed_times':diagnostic['fixed_times']}
            for name,fit in fits.items():
                replay=cost_at(response,noise,d['y'][i],fit['times_ns'],fit['charge_pe'])['cost']
                error=abs(replay-fit['cost'])/max(1.,abs(replay));diagnostic_cost_error=max(diagnostic_cost_error,error)
                assert error<1e-7,(domain,i,name)
                value=evaluate_fit(fit,d['t'][i],d['q'][i],c)
                diagnostic_rows.append(dict(value,domain=domain,method=name,index=i,condition=c['condition'],gap_ns=c['gap_ns'],
                    native_index=c['native_index'],block=c['block'],fit=fit,truth_used=name!='blind_native_gate'))
            a=evaluate_fit(baseline,d['t'][i],d['q'][i],c)['success']
            b=evaluate_fit(diagnostic['truth_start'],d['t'][i],d['q'][i],c)['success']
            cost_diagnostics.append(dict(domain=domain,index=i,condition=c['condition'],gap_ns=c['gap_ns'],
                original_success=a,true_start_success=b,
                original_minus_true_start=baseline['cost']-diagnostic['truth_start']['cost'],
                true_start_lower_cost=diagnostic['truth_start']['cost']<baseline['cost']-1e-6,
                new_correct_but_higher_cost=b and not a and diagnostic['truth_start']['cost']>=baseline['cost']-1e-6))
    diagnostic_summary=[]
    for domain in ('native','matched'):
        for condition in ('high_weak','weak_high'):
            for gap in (8.,40.,200.):
                for method in ('blind_native_gate','known_K','true_start','fixed_times'):
                    rr=[r for r in diagnostic_rows if (r['domain'],r['condition'],r['gap_ns'],r['method'])==(domain,condition,gap,method)]
                    diagnostic_summary.append(dict(domain=domain,condition=condition,gap_ns=gap,method=method,frames=len(rr),
                        success=sum(r['success'] for r in rr),efficiency=float(np.mean([r['success'] for r in rr])),
                        not_converged=sum(not r['fit'].get('converged',True) for r in rr)))
    assert len(diagnostic_rows)==2880 and len(diagnostic_summary)==48
    save(OUTPUT/'diagnostic-summary.json',dict(summary=diagnostic_summary,rows=diagnostic_rows,cost_comparisons=cost_diagnostics,
        warning='Known K/true initialization/fixed true times are diagnostic and not achievable blind bounds'))

    # Truth-conditioned Gaussian coordinates. These never enter blind fitting.
    residual_rows=[];white_arrays={};white_indices={};residual_cost_error=0.
    for domain,d in data.items():
        white=[];indices=[]
        for i,c in enumerate(d['cases']):
            if c['role']=='gate':continue
            valid=d['q'][i]>0.;q=d['q'][i,valid];t=d['t'][i,valid];curves=response.components(t,q)
            r=noise.transform(d['y'][i]-conditional_mean(response,noise,t,q));factor=GaussianFactor.from_noise(noise,curves,q)
            w=factor.whiten(r);value=noise.score(r,noise.transform(curves),q,transformed=True,details=True)
            error=abs(float(w@w)-value['quadratic']);residual_cost_error=max(residual_cost_error,error);assert error<1e-7
            centered=w-w.mean();variance=float(np.mean(centered**2))
            lag={str(k):float(centered[k:]@centered[:-k]/np.sqrt((centered[k:]@centered[k:])*(centered[:-k]@centered[:-k])))
                 for k in (1,2,4,8,16)}
            residual_rows.append(dict(domain=domain,index=i,kind=c['kind'],condition=c['condition'],gap_ns=c['gap_ns'],
                native_index=c['native_index'],draw=c['draw'],block=c['block'],quadratic_per_sample=value['quadratic']/len(w),
                white_mean=float(w.mean()),white_variance=variance,
                white_skew=float(np.mean(centered**3)/variance**1.5),white_excess_kurtosis=float(np.mean(centered**4)/variance**2-3.),
                outside6_fraction=float(np.mean(abs(w)>6.)),lag_correlations=lag))
            white.append(w);indices.append(i)
        white_arrays[domain]=np.asarray(white);white_indices[domain]=np.asarray(indices)
    np.savez_compressed(OUTPUT/'whitened-residuals.truth.npz',**white_arrays,
                        **{k+'_indices':v for k,v in white_indices.items()})
    residual_summary=[]
    for domain in data:
        groups={}
        for r in residual_rows:
            if r['domain']==domain:groups.setdefault((r['kind'],r['condition'],r['gap_ns']),[]).append(r)
        for key,rr in groups.items():
            residual_summary.append(dict(domain=domain,kind=key[0],condition=key[1],gap_ns=key[2],frames=len(rr),
                quadratic_per_sample_p16_p50_p84=np.quantile([r['quadratic_per_sample'] for r in rr],[.16,.5,.84]).tolist(),
                median_white_skew=float(np.median([r['white_skew'] for r in rr])),
                median_white_excess_kurtosis=float(np.median([r['white_excess_kurtosis'] for r in rr])),
                outside6_fraction=float(np.mean([r['outside6_fraction'] for r in rr]))))
    save(OUTPUT/'residual-diagnostics.json',dict(rows=residual_rows,summary=residual_summary,
        reference='At true parameters, exact Gaussian squared whitened norm has chi-square1021 distribution; not fitted chi-square/dof',
        scope='Native inputs are exposed and the covariance is a fitted surrogate; residual mismatch is diagnostic, not a detector resolution bound'))

    paired=paired_summary(data,scores,protocol)
    save(OUTPUT/'paired-comparison.json',paired)
    boundary=boundary_summary(data,gates,scored,response,noise)
    save(OUTPUT/'boundary-comparison.json',boundary)
    report=dict(summary=summary,diagnostic_summary=diagnostic_summary,paired=paired['comparisons'],
        interpretation='Frozen algorithm; exposed native760 and matched2x760 conditional Gaussian draws, plus280 matched gate points; no independent detector-performance claim',
        no_default_method_change=True,new_roots=0,selected_nonconverged=selected_nonconverged,
        runtime_seconds_quantiles=np.quantile([r['seconds'] for r in data['matched']['predictions'].values()],[.16,.5,.84,1.]).tolist(),
        boundary_changed_selected_count=boundary['changed_selected_count'],boundary_changed_strict_success=boundary['changed_strict_success'])
    save(OUTPUT/'report.json',report)
    save(RECORD/'science-validation.json',dict(passed=True,scores=len(scores),summary_rows=len(summary),diagnostic_scores=len(diagnostic_rows),
        diagnostic_summary_rows=len(diagnostic_summary),residual_frames=len(residual_rows),original_inputs_unchanged=True,
        native_scores_match_previous=True,mother_groups_disjoint=True,freeze_unchanged=True,
        maximum_candidate_relative_cost_error=cost_error,maximum_residual_quadratic_error=residual_cost_error,
        maximum_diagnostic_relative_cost_error=diagnostic_cost_error,
        source_sha256=digest(__file__),figures_pending=True,notebook_pending=True))
    print('MODEL_CLOSURE_ANALYZED',len(scores),len(summary),len(diagnostic_rows),flush=True)


def paired_summary(data,scores,protocol):
    bootstrap=protocol['bootstrap'];rng=np.random.default_rng(bootstrap['seed']);B=bootstrap['replicates']
    weights=np.zeros((B,40),int)
    for start in (0,20):
        draws=rng.integers(start,start+20,size=(B,20))
        for j in range(B):weights[j]+=np.bincount(draws[j],minlength=40)
    groups={}
    for r in scores:
        if r['kind']=='pair':groups.setdefault((r['method'],r['condition'],r['gap_ns']),[]).append(r)
    means={};estimates={}
    for key,rr in groups.items():
        blocks=np.array([r['block'] for r in rr]);den=np.bincount(blocks,weights=[r['target']==2 for r in rr],minlength=40)
        num=np.bincount(blocks,weights=[r['success'] and r['target']==2 for r in rr],minlength=40)
        assert np.all(weights@den>0)
        estimates[key]=(weights@num)/(weights@den);means[key]=float(num.sum()/den.sum())
    comparisons=[]
    for method,reference in [('matched/native_gate','native/native_gate'),('matched/matched_gate','native/native_gate'),
                             ('matched/matched_gate','matched/native_gate'),('native/matched_gate','native/native_gate')]:
        for condition in CONDITIONS:
            for gap in (8.,40.,200.):
                key=(method,condition,gap);base=(reference,condition,gap);delta=estimates[key]-estimates[base]
                comparisons.append(dict(method=method,reference=reference,condition=condition,gap_ns=gap,
                    difference=means[key]-means[base],paired_block95=np.quantile(delta,[.025,.975]).tolist()))
    return dict(comparisons=comparisons,**bootstrap,
        warning='Noise draws are independent conditional on repeated parameters; mother blocks preserve native reuse. Conditional diagnostic intervals, not independent detector validation.')


def boundary_summary(data,gates,scored,response,noise):
    from .alpha_pair_joint_likelihood import select
    from .alpha_pair_conditional_study import evaluate_fit
    with np.load(OUTPUT/'boundary-observations.npz') as f:y=f['voltage_mV']
    metadata=load(OUTPUT/'boundary.cases.json');predictions=rows(OUTPUT/'boundary.predictions.jsonl');result=[]
    assert len(predictions)==16
    d=data['matched']
    for r in predictions:
        i=r['index'];m=metadata[i];j=m['main_index'];c=d['cases'][j]
        original=scored['matched/native_gate',j];fit=select(r['candidates'],gates['native_gate'])
        score=evaluate_fit(fit,d['t'][j],d['q'][j],c)
        error=float(np.max(abs(noise.transform(y[i]-d['y'][j]))));assert error<1e-9
        result.append(dict(index=i,main_index=j,condition=c['condition'],gap_ns=c['gap_ns'],prefix_mV=m['prefix_mV'],
            same_count=len(fit['charge_pe'])==original['count'],same_strict_success=score['success']==original['success'],
            original_fit=original['fit'],perturbed_fit=fit,transform_max_difference=error))
    return dict(rows=result,changed_selected_count=sum(not r['same_count'] for r in result),
        changed_strict_success=sum(not r['same_strict_success'] for r in result),
        interpretation='Finite search sensitivity to a likelihood-null boundary component;16 prespecified perturbations, not a full hardware boundary model')


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    assert load(RECORD/'science-validation.json')['passed']
    folder=OUTPUT/'figures';folder.mkdir(exist_ok=True)
    report=load(OUTPUT/'report.json');summary=report['summary'];data=read_data();response,noise=models()
    scores={(r['method'],r['index']):r for r in rows(OUTPUT/'scores.jsonl')}

    fig,axes=plt.subplots(1,3,figsize=(15,5.2),sharey=True)
    for ax,condition,panel in zip(axes,CONDITIONS,'abc'):
        for method in METHODS:
            rr=sorted([r for r in summary if r['method']==method and r['condition']==condition],key=lambda r:r['gap_ns'])
            ax.plot([r['gap_ns'] for r in rr],[100*r['efficiency'] for r in rr],
                    'o--' if method=='native/matched_gate' else 'o-',color=COLORS[method],label=LABELS[method],lw=1.5)
        ax.set(title=f'({panel}) {condition}',xscale='log',xticks=[8,40,200],xticklabels=['8','40','200'],
               xlabel='Separation (ns)',ylim=(-3,103));ax.grid(alpha=.2)
    axes[0].set_ylabel('Strict pair recovery (%)')
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=2,fontsize=9)
    fig.suptitle('Figure 0004-72-1 | Frozen blind inference: 40 native pairs / 80 Gaussian draws per point; repeated parameters')
    fig.tight_layout(rect=(0,.14,1,.93));fig.savefig(folder/'0004-72-1.png',dpi=150);plt.close(fig)

    fig,axes=plt.subplots(1,3,figsize=(16,5.5));controls=['all','weak','low','mid','high']
    for j,method in enumerate(METHODS):
        rr=[next(r for r in summary if r['method']==method and r['condition']==c and r['kind']!='pair') for c in controls]
        x=np.arange(5)+(j-1.5)*.16
        axes[0].bar(x,[100*r['extra_fraction'] for r in rr],width=.15,color=COLORS[method],label=LABELS[method])
        axes[1].plot(x,[100*r['efficiency'] for r in rr],'o-',color=COLORS[method])
        axes[2].plot(x[1:],[100*r['q_error_p16_p50_p84'][1] for r in rr[1:]],'o-',color=COLORS[method])
    for ax,title,panel in zip(axes,['Extra-output frames (%)','Strict success (%)','Matched median Q error (%)'],'abc'):
        ax.set(title=f'({panel}) {title}',xticks=range(5),xticklabels=['empty','weak','low','mid','high']);ax.grid(alpha=.2)
    axes[0].axhline(5.,color='black',ls='--',lw=.8);axes[1].axhline(95.,color='black',ls='--',lw=.8)
    for value in (-5.,5.):axes[2].axhline(value,color='black',ls='--',lw=.8)
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=2,fontsize=9)
    fig.suptitle('Figure 0004-72-2 | Empty / single controls: 80 native frames and 160 conditional Gaussian draws per group')
    fig.tight_layout(rect=(0,.15,1,.92));fig.savefig(folder/'0004-72-2.png',dpi=150);plt.close(fig)

    ds=load(OUTPUT/'diagnostic-summary.json')['summary']
    methods=['blind_native_gate','known_K','true_start','fixed_times']
    labels=['Blind, original gate','Existing K=2 candidate','True Q/t initialization','Fixed true times; fit Q']
    colors=['#555555','#0077bb','#aa4499','#ee7733']
    fig,axes=plt.subplots(2,2,figsize=(13,9),sharey=True)
    for row,domain in enumerate(('native','matched')):
        for col,condition in enumerate(('high_weak','weak_high')):
            ax=axes[row,col]
            for j,(method,label,color) in enumerate(zip(methods,labels,colors)):
                rr=sorted([r for r in ds if r['domain']==domain and r['condition']==condition and r['method']==method],key=lambda r:r['gap_ns'])
                ax.bar(np.arange(3)+(j-1.5)*.19,[100*r['efficiency'] for r in rr],width=.18,label=label,color=color)
            panel='abcd'[row*2+col]
            ax.set(title=f'({panel}) {domain}: {condition}',xticks=range(3),xticklabels=['8','40','200'],
                   xlabel='Separation (ns)',ylabel='Strict pair recovery (%)',ylim=(0,103));ax.grid(axis='y',alpha=.2)
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=2,fontsize=10)
    fig.suptitle('Figure 0004-73-1 | Truth-assisted diagnosis: 40 native / 80 Gaussian per point; assisted fits are not blind limits')
    fig.tight_layout(rect=(0,.09,1,.95));fig.savefig(folder/'0004-73-1.png',dpi=150);plt.close(fig)

    residual=load(OUTPUT/'residual-diagnostics.json')
    with np.load(OUTPUT/'whitened-residuals.truth.npz') as f:
        whitened={domain:f[domain] for domain in ('native','matched')}
        residual_indices={domain:f[domain+'_indices'] for domain in ('native','matched')}
    fig,axes=plt.subplots(2,2,figsize=(16,10));domain_colors={'native':'#555555','matched':'#009988'}
    ordered=[(condition,gap) for condition in CONDITIONS for gap in (8.,40.,200.)]
    for offset,domain in [(-.09,'native'),(.09,'matched')]:
        rr=[next(r for r in residual['summary'] if r['domain']==domain and (r['condition'],r['gap_ns'])==key) for key in ordered]
        quantile=np.array([r['quadratic_per_sample_p16_p50_p84'] for r in rr])
        axes[0,0].errorbar(np.arange(9)+offset,quantile[:,1],yerr=[quantile[:,1]-quantile[:,0],quantile[:,2]-quantile[:,1]],
                           fmt='o',color=domain_colors[domain],label=domain,capsize=3)
        values=whitened[domain].ravel();bins=np.linspace(-6.,6.,121);hist,edges=np.histogram(values,bins=bins)
        axes[0,1].plot((edges[1:]+edges[:-1])/2,hist/len(values)/np.diff(edges),color=domain_colors[domain],
                       label=f'{domain}; outside +/-6: {100*np.mean(abs(values)>6.):.3f}%')
        rr=[r for r in residual['rows'] if r['domain']==domain and r['kind']=='pair']
        lags=(1,2,4,8,16)
        axes[1,0].plot(lags,[np.mean([r['lag_correlations'][str(k)] for r in rr]) for k in lags],
                       'o-',color=domain_colors[domain],label=domain)
        aligned=[];relative=np.arange(-50.,501.,2.)
        d=data[domain]
        for row,i in enumerate(residual_indices[domain]):
            c=d['cases'][i]
            if c['kind']=='pair' and c['condition']=='high_weak' and c['gap_ns']==40.:
                aligned.append(np.interp(relative+d['t'][i,0],response.sample_grid[3:],whitened[domain][row]))
        axes[1,1].plot(relative,np.mean(aligned,axis=0),color=domain_colors[domain],label=domain)
    axes[0,0].axhline(1.,ls='--',color='black',lw=.8)
    axes[0,0].set(title='(a) Truth quadratic / 1021; median and 16-84% range',xticks=range(9),
                  xticklabels=[f'{c}\n{g:g} ns' for c,g in ordered],ylabel='Quadratic / sample count')
    axes[0,0].tick_params(axis='x',labelsize=8,rotation=35)
    x=np.linspace(-6.,6.,400);axes[0,1].plot(x,np.exp(-x*x/2)/np.sqrt(2*np.pi),'--',color='#ee7733',label='N(0,1)')
    axes[0,1].set(title='(b) Truth-whitened samples; pooled evaluation frames',xlabel='Whitened residual',ylabel='Density (full-sample normalization)',yscale='log',ylim=(1e-5,1.))
    axes[1,0].axhline(0.,color='black',lw=.8);axes[1,0].set(title='(c) Mean within-frame lag correlation; all pairs',
        xlabel='Lag (2 ns samples)',ylabel='Correlation',xscale='log',xticks=[1,2,4,8,16],xticklabels=['1','2','4','8','16'])
    axes[1,1].axhline(0.,color='black',lw=.8);axes[1,1].set(title='(d) high_weak 40 ns: aligned mean whitened residual',
        xlabel='Time relative to true strong onset (ns)',ylabel='Mean whitened residual')
    for ax in axes.flat:ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Figure 0004-74-1 | Residuals at true parameters; Gaussian reference is conditional on this fitted surrogate')
    fig.tight_layout(rect=(0,0,1,.95));fig.savefig(folder/'0004-74-1.png',dpi=150);plt.close(fig)

    original=PROJECT/'jupyter/artifacts/alpha_pair_joint_confirmation_20261026/labels.truth.npz'
    with np.load(original) as f:native_components=f['components_mV']
    fig,axes=plt.subplots(4,3,figsize=(18,15));examples=[];sample=response.sample_grid
    for pair,gap in enumerate((40.,200.)):
        native_index=next(i for i,c in enumerate(data['native']['cases']) if c['kind']=='pair' and c['condition']=='high_weak'
                          and c['gap_ns']==gap and c['replicate']==0)
        matched_index=next(i for i,c in enumerate(data['matched']['cases']) if c['role']=='evaluation' and c['draw']==0 and c['native_index']==native_index)
        for offset,(domain,index) in enumerate([('native',native_index),('matched',matched_index)]):
            row=pair*2+offset;d=data[domain];c=d['cases'][index];score=scores[domain+'/native_gate',index];fit=score['fit']
            valid=d['q'][index]>0.;mean_curves=response.components(d['t'][index,valid],d['q'][index,valid])
            mean=mean_curves.sum(0)+response.baseline+noise.background_mV
            curves=response.components(fit['times_ns'],fit['charge_pe']);predicted=curves.sum(0)+response.baseline+noise.background_mV
            weak=curves[int(np.argmin(fit['charge_pe']))] if len(curves)==2 else np.zeros_like(sample)
            true_weak=int(np.argmin(d['q'][index]))
            axes[row,0].plot(sample,d['y'][index],color='#aaaaaa',lw=.8,label='Observed voltage')
            axes[row,0].plot(sample,mean,color='#ee7733',ls='--',lw=1.,label='Conditional true-parameter mean')
            axes[row,0].plot(sample,predicted,color='#0077bb',lw=1.,label='Blind selected fit')
            if domain=='native':axes[row,1].plot(sample,native_components[c['original_index'],true_weak],color='black',lw=.8,label='Native weak avalanche component')
            axes[row,1].plot(sample,mean_curves[true_weak],color='#ee7733',ls='--',lw=1.,label='Conditional weak mean')
            axes[row,1].plot(sample,weak,color='#0077bb',lw=1.,label='Blind smaller-Q K=2 output')
            axes[row,2].plot(sample,d['y'][index]-predicted,color='#0077bb',lw=.7,label='Observed - selected fit')
            for col,title in enumerate(('Total voltage','Weak component / conditional mean','Fit residual')):
                ax=axes[row,col];panel=chr(ord('a')+row*3+col)
                ax.set(title=f'({panel}) {domain}, {gap:g} ns: {title}',xlim=(0,2048),ylabel='mV');ax.grid(alpha=.2)
                for onset in d['t'][index]:ax.axvline(onset,ls=':',color='black',alpha=.3,lw=.8)
                if row==3:ax.set_xlabel('Record time (ns)')
                if row<2:ax.legend(fontsize=7)
            examples.append(dict(domain=domain,index=index,native_index=native_index,gap_ns=gap,case=c,score=score))
    fig.suptitle('Figure 0004-75-1 | Prespecified high_weak pair 0; full 0-2048 ns; Gaussian draws have no unique avalanche component truth')
    fig.tight_layout(rect=(0,0,1,.96));fig.savefig(folder/'0004-75-1.png',dpi=150);plt.close(fig)
    save(OUTPUT/'examples.json',examples)
    save(RECORD/'figure-delivery.json',dict(passed=True,figures={str(p.relative_to(PROJECT)):digest(p) for p in sorted(folder.glob('*.png'))},
        source_sha256=digest(__file__)))
    print('MODEL_CLOSURE_FIGURES',len(list(folder.glob('*.png'))),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=('analyze','figures'))
    args=parser.parse_args();globals()[args.command]()
