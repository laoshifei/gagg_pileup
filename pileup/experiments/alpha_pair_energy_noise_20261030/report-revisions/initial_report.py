"""Read-only scoring and standalone figures for conditional energy/noise checks."""
from __future__ import annotations

from collections import defaultdict
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import solve_banded

from .alpha_pair_energy_noise_study import PROJECT, EXPERIMENT, OUTPUT, RECORD, load, save, digest
from .alpha_pair_recovery_energy import RecoveryEnergyResponse

CONDITIONS = [('weak', None), ('high', None), ('high_weak', 8.), ('high_weak', 40.),
              ('high_weak', 200.), ('weak_high', 40.)]
LABELS = ['Weak single','Strong single','Strong → weak\n8 ns','Strong → weak\n40 ns',
          'Strong → weak\n200 ns','Weak → strong\n40 ns']
DOMAINS = ['physical','mean_replaced','gaussian']
COLORS = ['#2878b5','#e59d25','#408950']
NAMES = ['Repeated digitization','Mean replaced (oracle)','Matched Gaussian']


def write_csv(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def stats(values):
    x = np.asarray(values)
    return dict(bias_mean=float(x.mean()), bias_median=float(np.median(x)),
                sigma68=float((np.percentile(x,84)-np.percentile(x,16))/2),
                rmse=float(np.sqrt(np.mean(x*x))), sd=float(x.std(ddof=1)) if len(x)>1 else None)


def interaction_decomposition(response, results):
    """Exact waveform identity; local truth-point energy projection is diagnostic."""
    single = {}
    for r in results:
        t=r['task']
        if len(t['parents']) != 1:continue
        with np.load(OUTPUT/'cases'/f"{t['id']:02}.npz") as data:
            single[t['parents'][0]]=(t['times_ns'][0],data['components'][:r['n_reference'],0].mean(0))
    rows=[]; errors=[]
    for r in results:
        t=r['task']
        if len(t['parents']) != 2:continue
        with np.load(OUTPUT/'cases'/f"{t['id']:02}.npz") as data:
            physical=data['components'][:r['n_reference']].mean(0).sum(0)
            factor=data['factor'].copy()
        separated=[]
        for parent,arrival in zip(t['parents'],t['times_ns']):
            origin,curve=single[parent]
            separated.append(np.interp(response.sample_grid-(arrival-origin),response.sample_grid,curve,left=0.,right=0.))
        empirical_isolated=np.sum(separated,axis=0)
        model_isolated=sum(response.components([arrival],[energy])[0][0]
                           for arrival,energy in zip(t['times_ns'],t['energy_MeV']))
        model_pair=response.components(t['times_ns'],t['energy_MeV'])[0].sum(0)
        isolated=empirical_isolated-model_isolated
        interaction=(physical-empirical_isolated)-(model_pair-model_isolated)
        total=physical-model_pair
        err=float(np.max(abs(total-isolated-interaction)));errors.append(err);assert err<1e-10
        def whiten(value):
            return solve_banded((len(factor)-1,0),factor,response.transform(value),check_finite=False)
        columns=[]
        for j in range(2):
            up=np.asarray(t['energy_MeV']).copy();down=up.copy();up[j]*=np.exp(1e-4);down[j]*=np.exp(-1e-4)
            derivative=(response.components(t['times_ns'],up)[0].sum(0)-response.components(t['times_ns'],down)[0].sum(0))/(2e-4)
            columns.append(whiten(derivative))
        inverse=np.linalg.pinv(np.column_stack(columns),rcond=1e-10)
        projected=[inverse@whiten(x) for x in [isolated,interaction,total]]
        np.testing.assert_allclose(projected[0]+projected[1],projected[2],atol=1e-11)
        for j,label in enumerate(t['labels']):
            rows.append(dict(case=t['id'],block=t['block'],condition=t['condition'],gap_ns=t['gap_ns'],component=j,label=label,
                isolated_mean_logE=float(projected[0][j]),interaction_mean_logE=float(projected[1][j]),
                total_component_mean_logE=float(projected[2][j]),waveform_identity_max_error_mV=err))
    return rows,max(errors)


def make_report():
    for path, sha in load(RECORD/'result-freeze.json')['files'].items(): assert digest(OUTPUT/path) == sha, path
    response = RecoveryEnergyResponse.load(EXPERIMENT)
    results = [load(p) for p in sorted((OUTPUT/'cases').glob('*.json'))]
    assert len(results) == 24
    frames, components, means, projections = [], [], [], []
    minimum_q = float('inf'); qualified = injected = 0; noise_closure = 0.
    for result in results:
        task = result['task']; truth = np.asarray(task['energy_MeV'])
        base = dict(case=task['id'], block=task['block'], condition=task['condition'], gap_ns=task['gap_ns'])
        with np.load(OUTPUT/'cases'/f"{task['id']:02}.npz") as arrays:
            q = arrays['charges'][result['n_reference']:]
            minimum_q = min(minimum_q,float(q.min())); qualified += int((q>=500).sum()); injected += q.size
        noise_closure = max(noise_closure, max(r['closure_mV'] for r in result['simulation_records']))
        for row in result['fits']:
            fit = row['fit']; error = np.asarray(fit['energy_MeV'])/truth-1
            np.testing.assert_array_equal(fit['times_ns'],task['times_ns'])
            common = dict(base,domain=row['domain'],weight=row['weight'],repeat=row['repeat'])
            frames.append(dict(common,success5=int(np.all(abs(error)<=.05)), success10=int(np.all(abs(error)<=.1)),
                success20=int(np.all(abs(error)<=.2)), not_converged=int(not fit['converged']),boundary=int(fit['boundary']),cost=fit['cost']))
            for j, value in enumerate(error):
                components.append(dict(common,component=j,label=task['labels'][j],energy_MeV=truth[j],
                    predicted_energy_MeV=fit['energy_MeV'][j],error=float(value),success10=int(abs(value)<=.1)))
        for row in result['mean_fits']:
            err = np.asarray(row['fit']['energy_MeV'])/truth-1
            for j,value in enumerate(err):
                means.append(dict(base,domain=row['domain'],weight=row['weight'],component=j,label=task['labels'][j],
                                  relative_energy_error=float(value),not_converged=int(not row['fit']['converged'])))
        for j,label in enumerate(task['labels']):
            p=result['projection']
            projections.append(dict(base,component=j,label=label,condition_number=p['condition_number'],
                expected_gaussian_logE_sd=p['gaussian_local_logE_sd'][j],
                measured_repeat_logE_sd=p['evaluation_projected_logE_sd'][j],
                mean_projected_logE=p['reference_mean_projected_logE'][j],
                reference_mean_se=p['reference_mean_projected_se'][j],
                whitened_repeat_variance=p['whitened_evaluation_variance_mean']))
    groups = defaultdict(list)
    for f in frames:groups[f['condition'],f['gap_ns'],f['domain'],f['weight']].append(f)
    summaries = []
    for (condition,gap,domain,weight), rr in groups.items():
        selected = [x for x in components if (x['condition'],x['gap_ns'],x['domain'],x['weight'])==(condition,gap,domain,weight)]
        weak = [r['error'] for r in selected if r['label']=='weak']
        summaries.append(dict(condition=condition,gap_ns=gap,domain=domain,weight=weight,draws=len(rr),mother_blocks=len({r['block'] for r in rr}),
            **{key:sum(r[key] for r in rr) for key in ['success5','success10','success20','not_converged','boundary']},
            component_bias=stats([r['error'] for r in selected])['bias_mean'],
            weak_bias=stats(weak)['bias_mean'] if weak else None,weak_sigma68=stats(weak)['sigma68'] if weak else None,
            weak_rmse=stats(weak)['rmse'] if weak else None))
    parent_rows = []
    for summary in summaries:
        for block in [0,10,20,30]:
            rr = [f for f in frames if all(f[k]==summary[k] for k in ['condition','gap_ns','domain','weight']) and f['block']==block]
            ee = [r['error'] for r in components if all(r[k]==summary[k] for k in ['condition','gap_ns','domain','weight']) and r['block']==block and r['label']=='weak']
            parent_rows.append(dict(condition=summary['condition'],gap_ns=summary['gap_ns'],domain=summary['domain'],weight=summary['weight'],
                block=block,draws=len(rr),success10=sum(f['success10'] for f in rr),weak_bias=float(np.mean(ee)) if ee else None))
    table_dir = EXPERIMENT/'tables'; table_dir.mkdir(exist_ok=True)
    for name, rows in [('frames',frames),('components',components),('summary',summaries),('parent-summary',parent_rows),
                       ('mean-fits',means),('local-projection',projections)]: write_csv(table_dir/(name+'.csv'),rows)
    interactions,interaction_closure=interaction_decomposition(response,results)
    write_csv(table_dir/'interaction-decomposition.csv',interactions)
    assert len(frames)==1152 and all(r['draws']==32 and r['mother_blocks']==4 for r in summaries)
    artifacts = []
    def finish(fig,name):
        fig.tight_layout()
        for suffix in ['png','svg']:
            path = OUTPUT/(name+'.'+suffix);fig.savefig(path,dpi=155,bbox_inches='tight');artifacts.append(path)
        plt.close(fig)

    lookup = {(r['condition'],r['gap_ns'],r['domain'],r['weight']):r for r in summaries}
    fig,ax=plt.subplots(figsize=(12,4.8)); x=np.arange(6)
    for k,(domain,color,name) in enumerate(zip(DOMAINS,COLORS,NAMES)):
        y=[lookup[c,g,domain,'gls']['success10']/32*100 for c,g in CONDITIONS]
        yd=[lookup[c,g,domain,'diagonal']['success10']/32*100 for c,g in CONDITIONS]
        ax.bar(x+(k-1)*.25,y,.24,color=color,label=name+' / GLS',alpha=.85)
        ax.scatter(x+(k-1)*.25,yd,marker='_',color='black',s=110,zorder=3,label='Diagonal weights' if k==0 else None)
    ax.set(xticks=x,xticklabels=LABELS,ylim=(0,109),ylabel='All injected energies within ±10% (%)',
           title='Figure 0007-40-1: Known K/t conditional energy recovery')
    ax.legend(ncol=2,fontsize=8,loc='upper right');ax.grid(axis='y',alpha=.2)
    fig.text(.5,-.02,'4 optical pairs; 8 evaluation draws per condition/pair. Repeats are not independent optical events.',ha='center',fontsize=9)
    finish(fig,'0007-40-1_conditional_recovery')

    fig,axes=plt.subplots(2,2,figsize=(12,8),sharey=True)
    for ax,(condition,gap) in zip(axes.flat,CONDITIONS[2:]):
        groups=[]
        for domain in DOMAINS:
            groups.append([r['error']*100 for r in components if r['condition']==condition and r['gap_ns']==gap and
                           r['domain']==domain and r['weight']=='gls' and r['label']=='weak'])
        bp=ax.boxplot(groups,tick_labels=['Physical','Mean replaced','Gaussian'],patch_artist=True,showfliers=True)
        for patch,color in zip(bp['boxes'],COLORS):patch.set_facecolor(color);patch.set_alpha(.5)
        for k,domain in enumerate(DOMAINS):
            pp=[r for r in parent_rows if (r['condition'],r['gap_ns'],r['domain'],r['weight'])==(condition,gap,domain,'gls')]
            ax.scatter(np.full(4,k+1)+np.linspace(-.07,.07,4),[100*r['weak_bias'] for r in pp],color='black',s=16,zorder=4)
        ax.axhspan(-10,10,color='gray',alpha=.13);ax.axhline(0,color='gray',lw=.8)
        ax.set(title=f'{condition}, {gap:g} ns',ylabel='Weak deposited-energy error (%)');ax.grid(axis='y',alpha=.2)
    fig.suptitle('Figure 0007-41-1: Weak-energy errors; black dots = four mother-block means',y=1.01)
    finish(fig,'0007-41-1_weak_energy')

    selected=results[3]; task=selected['task']; assert task['id']==3 and task['condition']=='high_weak' and task['gap_ns']==40
    with np.load(OUTPUT/'cases/03.npz') as f: data={k:f[k].copy() for k in f.files}
    nr=selected['n_reference']; physical=data['physical']; reference=physical[:nr]; example=physical[nr]
    fit=next(r['fit'] for r in selected['fits'] if r['domain']=='physical' and r['weight']=='gls' and r['repeat']==0)
    pred=response.components(task['times_ns'],fit['energy_MeV'])[0]
    true_model=response.components(task['times_ns'],task['energy_MeV'])[0]
    truth_reference=data['components'][:nr].mean(0)
    time=response.sample_grid
    fig,axes=plt.subplots(4,1,figsize=(12,11),sharex=True)
    axes[0].plot(time,example,lw=.7,alpha=.5,label='Evaluation waveform 0')
    axes[0].plot(time,reference.mean(0),label='Independent reference mean')
    axes[0].plot(time,data['model_mean'],'--',label='Model at true energies')
    axes[0].plot(time,pred.sum(0)+response.baseline+response.background,':',label='GLS fit, true times')
    axes[0].set(title='(a) Complete observed waveform',ylabel='Voltage (mV)');axes[0].legend(fontsize=8)
    residual=reference.mean(0)-data['model_mean'];se=reference.std(0,ddof=1)/np.sqrt(nr)
    axes[1].plot(time,residual,label='Reference mean − model at truth');axes[1].fill_between(time,residual-se,residual+se,alpha=.25,label='±1 conditional reference SE')
    axes[1].axhline(0,color='gray',lw=.8);axes[1].set(title='(b) Conditional mean discrepancy',ylabel='Residual (mV)');axes[1].legend(fontsize=8)
    for j,ax in enumerate(axes[2:]):
        ax.plot(time,truth_reference[j],label='Mean true avalanche component')
        ax.plot(time,true_model[j],'--',label='Model component at true energy')
        ax.plot(time,pred[j],':',label='Fitted component, evaluation waveform 0')
        ax.set(title=f'({chr(99+j)}) {task["labels"][j]} component; true Edep={task["energy_MeV"][j]:.4f} MeV',ylabel='Component (mV)');ax.legend(fontsize=8)
    for ax in axes:ax.set_xlim(0,2048);ax.grid(alpha=.2)
    axes[-1].set_xlabel('Time from window start (ns)')
    fig.suptitle('Figure 0007-42-1: Fixed example, mother block 0, strong → weak at 40 ns',y=1.005)
    finish(fig,'0007-42-1_full_waveform')

    fig,ax=plt.subplots(figsize=(11,4.8));x=np.arange(4)
    for k,(key,color,label) in enumerate(zip(['isolated_mean_logE','interaction_mean_logE','total_component_mean_logE'],
                         ['#2878b5','#e59d25','#444444'],['Isolated mother-shape mismatch','Additional interaction mismatch','Total component-mean mismatch'])):
        values=[[r[key]*100 for r in interactions if (r['condition'],r['gap_ns'])==(c,g) and r['label']=='weak'] for c,g in CONDITIONS[2:]]
        ax.bar(x+(k-1)*.25,[np.mean(v) for v in values],.24,color=color,alpha=.65,label=label)
        for j,v in enumerate(values):ax.scatter(j+(k-1)*.25+np.linspace(-.06,.06,4),v,color=color,s=16,zorder=4)
    ax.axhline(0,color='black',lw=.7);ax.set(xticks=x,xticklabels=LABELS[2:],ylabel='100 × local projected delta ln(Eweak)',
        title='Figure 0007-43-1: Conditional mean discrepancy, projected at true parameters')
    ax.legend(fontsize=8);ax.grid(axis='y',alpha=.2)
    fig.text(.5,-.02,'Truth-assisted mean components; finite reference samples. Local projection is not fitted energy error or a bound.',ha='center',fontsize=9)
    finish(fig,'0007-43-1_mean_decomposition')

    result=dict(scope=load(EXPERIMENT/'protocol.json')['scope'],conditions=24,optical_mothers=8,optical_pairs=4,
        reference_draws_per_condition=16,evaluation_draws_per_condition=8,simulations=576,pilot_simulations=8,
        physical_evaluation_frames=192,gaussian_evaluation_frames=192,mean_replaced_frames=192,
        fitted_rows=1152,noiseless_and_reference_mean_fits=96,minimum_physical_evaluation_Q_pe=minimum_q,
        physical_evaluation_components=injected,physical_evaluation_components_Q_ge_500=qualified,
        max_component_sum_error_mV=noise_closure,
        mean_interaction_identity_max_error_mV=interaction_closure,
        model_noiseless_max_relative_energy_error=max(abs(r['relative_energy_error']) for r in means if r['domain']=='model_noiseless'),
        summary=summaries,parent_summary=parent_rows,
        means=means,local_projections=projections,interaction_decomposition=interactions,figures=[str(p.relative_to(PROJECT)) for p in artifacts],
        caveats=load(EXPERIMENT/'protocol.json')['limitations'],count_errors='not assessed: supplied K/t',
        no_rate_claim=True,no_default_update=True)
    save(EXPERIMENT/'result-summary.json',result)
    save(RECORD/'report-validation.json',dict(passed=True,known_times_unchanged=True,scored_rows=len(frames),
        summary_groups=len(summaries),independent_reference_indices=[0,15],evaluation_indices=[16,23],
        source_results_verified=True,figures={p.name:digest(p) for p in artifacts},
        no_optical_independence_claim=True,no_blind_count_efficiency_claim=True))
    print(json.dumps({k:result[k] for k in ['conditions','fitted_rows','minimum_physical_evaluation_Q_pe','model_noiseless_max_relative_energy_error']},indent=2))
    for row in summaries:
        if row['weight']=='gls':print(json.dumps(row))


if __name__=='__main__':make_report()
