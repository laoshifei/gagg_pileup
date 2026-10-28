"""Truth-separated scoring and graphics for the finite profile stability audit."""
from __future__ import annotations

from collections import defaultdict
import json
import numpy as np

from .alpha_pair_stability_study import PROJECT,EXPERIMENT,OUTPUT,RECORD,INPUT,paths,models,check_freeze
from .alpha_pair_evidence_profile import compatibility_segments
from .alpha_pair_profile_stability import compare_regions,point_cost
from .alpha_pair_information_study import load,save,digest,rows
from .alpha_pair_risk_report import evaluate,summarize,quantiles,csv_table

STAGES=('original','search','dense','contour')
COLORS=dict(zip(STAGES,('#777777','#0077bb','#ee7733','#009988')))
LABELS=dict(original='Original',search='More starts',dense='All midpoints',contour='Contour refinement')
GROUPS=(('native','40'),('native','200'),('continuous','middle'),('continuous','long'))
GROUP_LABELS=('Native 40 ns','Native 200 ns','Gaussian 24-64 ns','Gaussian 128-256 ns')


def data():
    out={}
    for domain in ('native','continuous'):
        with np.load(paths(domain)['truth']) as f:q=f['charge_pe'];t=f['times_ns']
        out[domain]=dict(q=q,t=t,cases=load(paths(domain)['cases']))
    return out


def analyze():
    check_freeze();run=load(RECORD/'inference.json');pr=load(RECORD/'probe-validation.json')
    assert run['passed'] and pr['passed']
    assert digest(OUTPUT/'profiles.jsonl')==run['sha256']==pr['primary_profiles_sha256']
    assert digest(OUTPUT/'truth-probes.jsonl')==pr['sha256']
    curves=rows(OUTPUT/'profiles.jsonl');probes={(r['domain'],r['index']):r for r in rows(OUTPUT/'truth-probes.jsonl')}
    old_probes={(r['domain'],r['index']):r for r in rows(INPUT/'truth-probes.jsonl')}
    truth=data();values=[];scored=[];changes=[];frames=[]
    response,noise=models();observed={};probe_replay_error=0.
    for domain in ('native','continuous'):
        with np.load(paths(domain)['observation']) as f:observed[domain]=f['voltage_mV']
    for record in curves:
        domain,i=record['domain'],record['index'];d=truth[domain];q,t=d['q'][i],d['t'][i];case=d['cases'][i]
        band=str(int(case['gap_ns'])) if domain=='native' else case['band']
        identity=dict(domain=domain,index=i,band=band)
        probe=probes[domain,i];old_probe=old_probes[domain,i]
        # Both probes are scoring-only. Retaining the better actual cost avoids
        # changing inclusion merely because one of the finite probe fits fails.
        score_probe=min((probe,old_probe),key=lambda p:p['fit']['cost'])
        for p in (probe,old_probe):
            actual=point_cost(response,noise,observed[domain][i],p['fit'])
            probe_replay_error=max(probe_replay_error,abs(actual-p['fit']['cost'])/max(1.,abs(actual)))
            assert probe_replay_error<1e-9
        final=record['stages']['contour']
        frames.append(dict(**identity,seconds=record['seconds'],
            final_stable_all_levels=all(c['stable'] for c in record['comparisons']['dense_contour']),
            original_reference=record['stages']['original']['reference_cost'],final_reference=final['reference_cost'],
            total_reference_gain=record['stages']['original']['reference_cost']-final['reference_cost'],
            probe_reference_gain=final['reference_cost']-score_probe['fit']['cost'],
            original_probe_cost=old_probe['fit']['cost'],new_probe_cost=probe['fit']['cost'],
            probe_converged=score_probe['fit']['converged'],new_probe_converged=probe['fit']['converged'],
            cost_replay_error=record['cost_replay_error']))
        for stage,curve in record['stages'].items():
            costs=[f['cost'] for f in curve['fits']];best=min(curve['fits'],key=lambda f:f['cost'])
            scored.append(dict(evaluate(best,t,q,case),**identity,method=stage,kind=case['kind'],condition=case['condition']))
            for level in (1.,4.,9.):
                region=compatibility_segments(curve['q_grid_pe'],costs,curve['reference_cost'],level)
                assert region['segments']
                values.append(dict(**identity,stage=stage,**region,
                    width_over_trueQ=(region['upper']-region['lower'])/q[1],components=len(region['segments']),
                    truth_Q2_pe=float(q[1]),truth_delta=score_probe['fit']['cost']-curve['reference_cost'],
                    truth_included=score_probe['fit']['cost']-curve['reference_cost']<=level,
                    points=len(curve['q_grid_pe']),reference_cost=curve['reference_cost'],
                    nonconverged_trials=curve.get('nonconverged_trials',0),trials=curve.get('trials'),
                    conditional_best_Q2_error=best['charge_pe'][1]/q[1]-1.,
                    conditional_best_success=evaluate(best,t,q,case)['success']))
        for pair,rr in record['comparisons'].items():
            a,b=pair.split('_')
            for change in rr:
                actual=compare_regions(record['stages'][a],record['stages'][b],change['level'])
                assert change==actual
                changes.append(dict(**identity,comparison=pair,**actual))
    assert len(curves)==32 and len(values)==384 and len(changes)==288 and len(scored)==128
    grouped=defaultdict(list)
    for r in values:grouped[r['domain'],r['band'],r['stage'],r['level']].append(r)
    profile_summary=[]
    for (domain,band,stage,level),rr in sorted(grouped.items()):
        assert len(rr)==8
        profile_summary.append(dict(domain=domain,band=band,stage=stage,level=level,frames=len(rr),
            width_over_trueQ=quantiles([r['width_over_trueQ'] for r in rr]),truth_included=sum(r['truth_included'] for r in rr),
            lower_boundary=sum(r['lower_boundary'] for r in rr),upper_boundary=sum(r['upper_boundary'] for r in rr),
            disconnected=sum(r['components']>1 for r in rr),points=quantiles([r['points'] for r in rr]),
            conditional_best_Q2_error=quantiles([r['conditional_best_Q2_error'] for r in rr]),
            conditional_best_success=sum(r['conditional_best_success'] for r in rr),
            nonconverged_trials=sum(r['nonconverged_trials'] for r in rr)))
    grouped=defaultdict(list)
    for r in changes:grouped[r['domain'],r['band'],r['comparison'],r['level']].append(r)
    stability=[]
    for (domain,band,pair,level),rr in sorted(grouped.items()):
        endpoints=[np.expm1(r['endpoint_log_change']) for r in rr if r['endpoint_log_change'] is not None]
        stability.append(dict(domain=domain,band=band,comparison=pair,level=level,frames=len(rr),stable=sum(r['stable'] for r in rr),
            topology_changed=sum(not r['same_topology'] for r in rr),boundary_changed=sum(not r['same_boundaries'] for r in rr),
            reference_gain=quantiles([r['reference_gain'] for r in rr]),reference_gain_max=max(r['reference_gain'] for r in rr),
            endpoint_fraction_change=quantiles(endpoints),endpoint_fraction_change_max=max(endpoints) if endpoints else None))
    grouped=defaultdict(list)
    for r in scored:grouped[r['domain'],r['band'],r['method']].append(r)
    metrics=[summarize(rr) for rr in grouped.values()]
    report=dict(profile_summary=profile_summary,stability_summary=stability,conditional_metrics=metrics,frames=frames,
        stable_all_levels=sum(f['final_stable_all_levels'] for f in frames),
        probe_improves_reference=sum(f['probe_reference_gain']>.01 for f in frames),
        no_default_replacement=True,known_K2=True,all_data_exposed=True,formal_confidence_intervals=False,
        new_ROOTs=0,new_training=0,new_random_draws=0,original160_profiles_unchanged=True)
    for name,content in (('profile-scores.json',values),('changes.json',changes),('conditional-scores.json',scored),('report.json',report)):
        save(OUTPUT/name,content)
    for name,content in (('profile-summary.csv',profile_summary),('stability-summary.csv',stability),('conditional-metrics.csv',metrics)):
        csv_table(OUTPUT/name,content)
    files=('profile-scores.json','changes.json','conditional-scores.json','report.json','profile-summary.csv','stability-summary.csv','conditional-metrics.csv')
    save(RECORD/'science-validation.json',dict(passed=True,frames=32,stages=4,level_scores=384,comparisons=288,
        known_K2_scored_fits=128,old_curves_unchanged=True,scoring_probes_not_fed_back=True,
        maximum_cost_replay_error=run['maximum_cost_replay_error'],maximum_probe_cost_replay_error=probe_replay_error,
        output_hashes={n:digest(OUTPUT/n) for n in files}))
    print('STABILITY_SCORED',report['stable_all_levels'],'of32 stable; trueQ lower-reference probes',report['probe_improves_reference'],flush=True)


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.2})
    report=load(OUTPUT/'report.json');values=load(OUTPUT/'profile-scores.json')
    records=rows(OUTPUT/'profiles.jsonl');lookup={(r['domain'],r['index']):r for r in records}
    truth=data();folder=OUTPUT/'figures';folder.mkdir(exist_ok=True);saved={}
    examples=load(EXPERIMENT/'protocol.json')['examples']
    def group_of(r):
        case=truth[r['domain']]['cases'][r['index']]
        return r['domain'],str(int(case['gap_ns'])) if r['domain']=='native' else case['band']
    def finish(fig,number,title,legend=True):
        for i,ax in enumerate(fig.axes):ax.set_title(f'({chr(97+i)}) '+ax.get_title())
        fig.suptitle('Figure '+number+' | '+title,fontsize=12)
        if legend:
            fig.legend([Line2D([0],[0],color=COLORS[s],lw=2) for s in STAGES],
                       [LABELS[s] for s in STAGES],loc='lower center',ncol=4,fontsize=9)
        fig.tight_layout(rect=(0,.07 if legend else .02,1,.94))
        path=folder/(number+'.png');fig.savefig(path,dpi=145);plt.close(fig);saved[number]=digest(path)
    fig,axes=plt.subplots(2,2,figsize=(13,8))
    for ax,group,label in zip(axes.flat,GROUPS,GROUP_LABELS):
        rr=sorted((r for r in records if group_of(r)==group),key=lambda r:r['index'])
        for r in rr:
            gain=[r['stages']['original']['reference_cost']-r['stages'][s]['reference_cost'] for s in STAGES]
            ax.plot(range(4),gain,'o-',alpha=.55,ms=3,lw=1,label=str(r['index']))
        ax.set(xticks=range(4),xticklabels=[LABELS[s] for s in STAGES],ylabel='Decrease in minimum cost J',title=label)
        ax.tick_params(axis='x',labelrotation=12);ax.legend(title='Frame index',fontsize=7,ncol=2)
    finish(fig,'0004-103-1','All eight prespecified frames per group; numerical search changes',False)
    fig,axes=plt.subplots(2,2,figsize=(13,8))
    for ax,task,label in zip(axes.flat,examples,GROUP_LABELS):
        r=lookup[task['domain'],task['index']];ref=r['stages']['contour']['reference_cost']
        for s in STAGES:
            c=r['stages'][s]
            ax.plot(c['q_grid_pe'],np.array([f['cost'] for f in c['fits']])-ref,color=COLORS[s],
                    lw=1.5,ls='--' if s=='original' else '-',label=LABELS[s])
        tq=truth[task['domain']]['q'][task['index'],1]
        ax.axvline(tq,color='black',lw=.8,ls=':',label='True Q2 (plot only)')
        ax.axhline(4.,color='black',lw=.7,ls=':')
        ax.set(xscale='log',xlim=(50.,50000.),ylim=(-.1,20.),xlabel='Later Q2 (pe)',ylabel='J - final sampled minimum',
               title=label+f'; index {task["index"]}')
    finish(fig,'0004-104-1','Prespecified profiles; all curves use the same final reference; vertical dotted line is truth')
    fig,axes=plt.subplots(2,2,figsize=(13,10))
    for ax,group,label in zip(axes.flat,GROUPS,GROUP_LABELS):
        ids=sorted(r['index'] for r in records if group_of(r)==group)
        for j,i in enumerate(ids):
            for k,s in enumerate(STAGES):
                v=next(v for v in values if (v['domain'],v['index'],v['stage'],v['level'])==(group[0],i,s,4.))
                for lo,hi in v['segments']:
                    ax.plot([lo/v['truth_Q2_pe'],hi/v['truth_Q2_pe']],[j+(k-1.5)*.15]*2,color=COLORS[s],lw=2)
        ax.axvline(1.,color='black',lw=.8,ls=':')
        ax.set(xscale='log',xlabel='Later Q2 / true Q2 (scoring scale only)',ylabel='Prespecified frame index',
               yticks=range(len(ids)),yticklabels=ids,title=label,ylim=(-.6,7.6))
    finish(fig,'0004-104-2','All separate intervals at delta J <= 4; each stage uses its own minimum')
    fig,axes=plt.subplots(2,2,figsize=(13,8));x=np.arange(4);short=['Native 40','Native 200','Gaussian middle','Gaussian long']
    for j,(level,color) in enumerate(zip((1.,4.,9.),('#0077bb','#009988','#ee7733'))):
        rr=[next(v for v in report['stability_summary'] if (v['domain'],v['band'],v['comparison'],v['level'])==(d,b,'dense_contour',level)) for d,b in GROUPS]
        axes[0,0].bar(x+(j-1)*.22,[v['stable'] for v in rr],width=.2,color=color,label=f'Delta J <= {level:g}')
    axes[0,0].legend(fontsize=8,loc='upper left');axes[0,0].set_ylim(0,11);axes[0,0].set_yticks([0,2,4,6,8])
    for j,s in enumerate(STAGES):
        rr=[next(v for v in report['profile_summary'] if (v['domain'],v['band'],v['stage'],v['level'])==(d,b,s,4.)) for d,b in GROUPS]
        pos=x+(j-1.5)*.18
        axes[0,1].bar(pos,[v['width_over_trueQ'][1] for v in rr],width=.17,color=COLORS[s])
        axes[1,0].bar(pos,[v['disconnected'] for v in rr],width=.17,color=COLORS[s])
        axes[1,1].bar(pos,[v['conditional_best_success'] for v in rr],width=.17,color=COLORS[s])
    for ax,title,ylabel in zip(axes.flat,('Last grid step meets numerical rule','Median contour width / true Q2','Frames with disconnected support','Strict recovery at conditional best K2'),
                             ('Frames out of 8','Width / true Q2, delta J <= 4','Frames out of 8, delta J <= 4','Frames out of 8; known K2')):
        ax.set(xticks=x,xticklabels=short,title=title,ylabel=ylabel);ax.tick_params(axis='x',labelrotation=12)
    axes[0,1].set_yscale('log')
    for ax in (axes[1,0],axes[1,1]):ax.set_ylim(0,8.6);ax.set_yticks([0,2,4,6,8])
    finish(fig,'0004-105-1','Numerical stability, compatibility width and conditional accuracy are separate')
    response,noise=models();fig,axes=plt.subplots(4,2,figsize=(14,13));sample=[];observed={}
    for d in ('native','continuous'):
        with np.load(paths(d)['observation']) as f:observed[d]=f['voltage_mV']
    for row,(task,label) in enumerate(zip(examples,GROUP_LABELS)):
        d,i=task['domain'],task['index'];r=lookup[d,i];tt=truth[d]['t'][i];qq=truth[d]['q'][i]
        original=min(r['stages']['original']['fits'],key=lambda f:f['cost'])
        final=min(r['stages']['contour']['fits'],key=lambda f:f['cost'])
        actual=response.components(tt,qq);grid=response.sample_grid
        axes[row,0].plot(grid,observed[d][i],color='.75',lw=.7,label='Observed voltage')
        axes[row,1].plot(grid,actual[1],color='black',lw=1,label='True-parameter conditional component')
        for s,fit in (('original',original),('contour',final)):
            component=response.components(fit['times_ns'],fit['charge_pe'])
            axes[row,0].plot(grid,component.sum(0)+response.baseline+noise.background_mV,color=COLORS[s],lw=1,label=LABELS[s])
            axes[row,1].plot(grid,component[1],color=COLORS[s],lw=1,label=f'{LABELS[s]} ({fit["charge_pe"][1]:.0f} pe)')
        for col in (0,1):
            axes[row,col].set(xlim=(0,2048),xlabel='Record time (ns)',ylabel='mV',
                title=f'{label}; index {i}: '+('full total voltage' if col==0 else 'later component, known K2'))
            axes[row,col].legend(fontsize=7)
        sample.append(dict(**task,truth_times_ns=tt.tolist(),truth_charge_pe=qq.tolist(),original=original,contour=final))
    finish(fig,'0004-106-1','Full-record conditional decompositions before and after numerical refinement',False)
    save(OUTPUT/'examples.json',sample)
    save(RECORD/'figure-delivery.json',dict(passed=True,figures=saved,source_sha256=digest(__file__)))
    print('STABILITY_FIGURES_COMPLETE',len(saved),flush=True)
