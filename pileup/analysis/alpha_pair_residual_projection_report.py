"""Tables and fixed figures for the frozen truth-assisted residual study."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from pileup.analysis import alpha_pair_residual_projection_study as s
from pileup.analysis.alpha_pair_residual_projection import TERMS, WEIGHTS, paired_bootstrap

CONDITIONS = (('high_weak',8.),('high_weak',40.),('high_weak',200.),
              ('high_weak',1200.),('weak_high',40.))
LABELS = ('H→W 8','H→W 40','H→W 200','H→W 1200','W→H 40')
TERM_LABELS = ('Strong isolated','Weak isolated','Shared correction','Dark/electronics')
COLORS = ('#4267a9','#dc8f28','#389075','#9c5a9e')


def csv_write(name,rows):
    path=s.EXP/'tables'/name;path.parent.mkdir(exist_ok=True)
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def load():
    cases=s.read(s.EXP/'cases.json')
    projections={c['id']:s.read(s.ART/'projections'/f"{c['id']:04d}.json") for c in cases}
    fits=[s.read(s.ART/'fits'/f"{t['id']:04d}.json") for t in s.read(s.EXP/'fit-plan.json')]
    previous=[]
    for c in cases:
        for j,method in enumerate(s.METHODS):
            row=s.read(s.PARENTART/'pairs'/f"{288+(c['id']-16)*6+j:04d}.json")
            e=np.asarray(c['energy_dep_MeV']);bias=100*(np.asarray(row['energy_MeV'])/e-1)
            weak=int(np.argmin(e))
            row.update(case=c['id'],condition=c['condition'],gap_ns=c['gap_ns'],block=c['block'],
                method=method,weight='full',start='voltage',delete_term=None,variant='original',
                weak_index=weak,bias_percent=bias.tolist(),weak_bias_percent=float(bias[weak]),
                strict=bool(row['converged'] and np.max(abs(bias))<=10),true_energy_MeV=e.tolist())
            previous.append(row)
    return cases,projections,fits,previous


def variant(row):
    if row.get('variant')=='original':return 'original_3starts'
    if row['delete_term'] is not None:return 'delete_'+row['delete_term']
    return 'true_E_start' if row['start']=='truth' else row['weight']+'_3starts'


def summarize(rows):
    weak=np.array([r['weak_bias_percent'] for r in rows])
    strong=np.array([r['bias_percent'][1-r['weak_index']] for r in rows])
    return dict(n=len(rows),strict=sum(r['strict'] for r in rows),
        weak_bias_percent=float(weak.mean()),weak_rmse_percent=float(np.sqrt(np.mean(weak**2))),
        weak_error_sd_percent=float(np.std(weak,ddof=1)),strong_bias_percent=float(strong.mean()),
        strong_rmse_percent=float(np.sqrt(np.mean(strong**2))),
        nonconverged=sum(not r['converged'] for r in rows),boundary=sum(r['boundary'] for r in rows))


def tables():
    s.verify();cases,projections,fits,previous=load()
    projection_rows=[];condition_rows=[];intervals=[];split_rows=[]
    for c in cases:
        for row in projections[c['id']]['rows']:
            weak=row['weak_index'];contributions=np.asarray(row['contributions_log']).mean(0)[:,weak]
            out={k:row[k] for k in ('case','condition','gap_ns','block','method','weight',
                'column_cosine','normalized_column_condition','profiled_weak_norm','outside_energy_bounds')}
            out.update(weak_delta_log=row['delta_log_energy'][weak],
                weak_projected_bias_percent=row['projected_bias_percent'][weak],
                weak_bounded_bias_percent=row['bounded_bias_percent'][weak],
                weak_old_fit_bias_percent=row['old_fit_bias_percent'][weak],
                **{term+'_log':float(v) for term,v in zip(TERMS,contributions)})
            projection_rows.append(out)
        for weight,change in projections[c['id']]['symmetric_changes'].items():
            weak=int(np.argmin(c['energy_dep_MeV']))
            split_rows.append(dict(case=c['id'],condition=c['condition'],gap_ns=c['gap_ns'],
                block=c['block'],weight=weight,**{key:val[weak] for key,val in change.items()}))
    for i,condition in enumerate(CONDITIONS):
        group=[c for c in cases if (c['condition'],c['gap_ns'])==condition]
        values=[]
        for c in group:
            rows=projections[c['id']]['rows'];weak=rows[0]['weak_index']
            values.append(np.stack([np.asarray(r['contributions_log'])[:,:,weak] for r in rows],axis=1))
        values=np.asarray(values)  # mother,repeat,6 model/weight arms,4 terms
        samples,ci=paired_bootstrap(values,seed=731100+i)
        total_ci=np.percentile(samples.sum(-1),[2.5,97.5],axis=0)
        for j,(method,weight) in enumerate((m,w) for m in s.METHODS for w in WEIGHTS):
            selected=[r for r in projection_rows if (r['condition'],r['gap_ns'])==condition
                      and r['method']==method and r['weight']==weight]
            out=dict(condition=condition[0],gap_ns=condition[1],method=method,weight=weight,mothers=4,
                mean_weak_delta_log=float(values[:,:,j].sum(-1).mean()),
                total_log_ci_low=float(total_ci[0,j]),total_log_ci_high=float(total_ci[1,j]),
                mean_weak_projected_bias_percent=float(np.mean([r['weak_projected_bias_percent'] for r in selected])),
                mean_weak_bounded_bias_percent=float(np.mean([r['weak_bounded_bias_percent'] for r in selected])),
                mean_weak_old_fit_bias_percent=float(np.mean([r['weak_old_fit_bias_percent'] for r in selected])),
                mean_column_cosine=float(np.mean([r['column_cosine'] for r in selected])),
                mean_normalized_condition=float(np.mean([r['normalized_column_condition'] for r in selected])))
            for k,term in enumerate(TERMS):
                out[term+'_log']=float(values[:,:,j,k].mean())
                out[term+'_ci_low']=float(ci[0,j,k]);out[term+'_ci_high']=float(ci[1,j,k])
            condition_rows.append(out)
        # Identical resampling across the two model arms and all four terms.
        intervals.append(dict(condition=condition[0],gap_ns=condition[1],
            full_spatial_minus_uniform_total_log=float((values[:,:,3]-values[:,:,0]).sum(-1).mean()),
            paired_log_ci95=np.percentile((samples[:,3]-samples[:,0]).sum(-1),[2.5,97.5]).tolist()))
    nonlinear_rows=[]
    all_rows=previous+fits
    for condition in CONDITIONS:
        for method in s.METHODS:
            selected=[r for r in all_rows if (r['condition'],r['gap_ns'])==condition and r['method']==method]
            for label in dict.fromkeys(variant(r) for r in selected):
                subset=[r for r in selected if variant(r)==label]
                nonlinear_rows.append(dict(condition=condition[0],gap_ns=condition[1],method=method,
                    variant=label,**summarize(subset)))
    old_lookup={(r['case'],r['method']):r for r in previous};search_rows=[]
    for r in fits:
        if r['start']!='truth' or r['delete_term'] is not None:continue
        old=old_lookup[(r['case'],r['method'])]
        tolerance=max(1e-5,1e-7*abs(old['cost']))
        search_rows.append(dict(case=r['case'],condition=r['condition'],gap_ns=r['gap_ns'],method=r['method'],
            old_cost=old['cost'],truth_start_cost=r['cost'],cost_difference=r['cost']-old['cost'],
            numerical_tolerance=tolerance,lower_cost=r['cost']<old['cost']-tolerance,
            higher_cost=r['cost']>old['cost']+tolerance,
            maximum_energy_change_percent=float(100*np.max(abs(np.asarray(r['energy_MeV'])/old['energy_MeV']-1))),
            old_weak_bias_percent=old['weak_bias_percent'],truth_start_weak_bias_percent=r['weak_bias_percent']))
    csv_write('projection-by-case.csv',projection_rows)
    csv_write('projection-by-condition.csv',condition_rows)
    csv_write('symmetric-model-change.csv',split_rows)
    csv_write('nonlinear-summary.csv',nonlinear_rows)
    csv_write('search-comparison.csv',search_rows)
    csv_write('fit-by-case.csv',[dict(id=r['id'],case=r['case'],condition=r['condition'],gap_ns=r['gap_ns'],
        block=r['block'],method=r['method'],variant=variant(r),cost=r['cost'],strict=r['strict'],
        weak_bias_percent=r['weak_bias_percent'],strong_bias_percent=r['bias_percent'][1-r['weak_index']],
        converged=r['converged'],boundary=r['boundary']) for r in fits])
    result=dict(projections=condition_rows,paired_intervals=intervals,nonlinear=nonlinear_rows,
        search=dict(n=len(search_rows),lower_cost=sum(r['lower_cost'] for r in search_rows),
                    higher_cost=sum(r['higher_cost'] for r in search_rows),
                    maximum_energy_change_percent=max(r['maximum_energy_change_percent'] for r in search_rows)),
        formal_fits=len(fits),nonconverged=sum(not r['converged'] for r in fits),
        pilot_fits=4,replay_fits=4,known_K_t=True,new_ROOT=0,new_digitizations=0,default_changed=False,
        scope='Four exposed mother pairs,8-repeat means. All projections and deletions truth-assisted; weight contrasts known K/t only. No blind efficiency/merge/split/selection metrics measured.',
        recommendation='No automatic parameter update, default replacement or expansion. Interpret residual/weight/search contrasts jointly; changes require a new frozen development study.')
    s.save(s.EXP/'result-summary.json',result)
    return result


def savefig(fig,name):
    folder=s.ART/'figures';folder.mkdir(exist_ok=True)
    fig.savefig(folder/(name+'.png'),dpi=160)
    fig.savefig(folder/(name+'.svg'))
    plt.close(fig)


def plots():
    cases,projections,fits,previous=load()
    result=s.read(s.EXP/'result-summary.json')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(14,5.2),layout='constrained')
    for ax,method,letter in zip(axes,s.METHODS,('a','b')):
        rows=[next(r for r in result['projections'] if (r['condition'],r['gap_ns'])==condition
                   and r['method']==method and r['weight']=='full') for condition in CONDITIONS]
        positive=np.zeros(5);negative=np.zeros(5);x=np.arange(5)
        for term,color,label in zip(TERMS,COLORS,TERM_LABELS):
            val=100*np.array([r[term+'_log'] for r in rows]);bottom=np.where(val>=0,positive,negative)
            ax.bar(x,val,bottom=bottom,color=color,label=label,width=.65)
            positive+=np.maximum(val,0);negative+=np.minimum(val,0)
        actual=[]
        for condition in CONDITIONS:
            group=[r for r in previous if r['method']==method and (r['condition'],r['gap_ns'])==condition]
            actual.append(100*np.mean([np.log1p(r['weak_bias_percent']/100) for r in group]))
        ax.plot(x,actual,'ko',label='Original nonlinear fit')
        ax.plot(x,[100*r['mean_weak_delta_log'] for r in rows],'k_',ms=16,label='Linear total')
        ax.axhline(0,color='0.4',lw=.7);ax.set_xticks(x,LABELS,rotation=20)
        ax.set(title=f'({letter}) {method}',ylabel='100 × weak Δln(E) [log-energy units]')
    axes[0].legend(fontsize=8,loc='best')
    fig.suptitle('Figure 0007-94-1: Correlated residual contributions at true energies (full weights)')
    savefig(fig,'fig_0007_94_1_residual_contributions')

    case=next(c for c in cases if c['id']==17)
    with np.load(s.ART/'projections/0017.npz') as f:arrays={k:f[k].copy() for k in f.files}
    s.initialize();background=s.MODELS['uniform'].background;t=s.ANCHOR.sample_grid
    fig,axes=plt.subplots(2,2,figsize=(14,8),layout='constrained')
    ax=axes[0,0];ax.plot(t,arrays['voltage_mean_mV'],color='0.55',lw=.7,label='Observed mean (8 repeats)')
    for method,color in zip(s.METHODS,('#315b9a','#d77422')):
        ax.plot(t,arrays[method+'_mean_mV']+background,lw=1,label=method+' at true E',color=color)
    ax.set(title='(a) Complete observed window',ylabel='Voltage [mV]');ax.legend(fontsize=8)
    terms=arrays['spatial12_terms_mV'].mean(0)
    for values,color,label in zip(terms,COLORS,TERM_LABELS):
        axes[0,1].plot(t,values,color=color,label=label,lw=.8)
    axes[0,1].set(title='(b) Spatial model residual terms',ylabel='Residual [mV]');axes[0,1].legend(fontsize=8)
    for method,color in zip(s.METHODS,('#315b9a','#d77422')):
        axes[1,0].plot(t,arrays[method+'_full_influence'][1],color=color,label=method,lw=.8)
    axes[1,0].set(title='(c) Raw-voltage weak-energy influence',ylabel='Δln(E) / mV');axes[1,0].legend(fontsize=8)
    influence=arrays['spatial12_full_influence'][1]
    for values,color,label in zip(terms,COLORS,TERM_LABELS):
        axes[1,1].plot(t,100*np.cumsum(influence*values),color=color,label=label)
    axes[1,1].set(title='(d) Cumulative spatial contributions',ylabel='100 × weak Δln(E)')
    for ax in axes.ravel():
        ax.set(xlim=(0,2048),xlabel='Time [ns]');ax.axhline(0,color='0.6',lw=.5)
        for arrival in case['times_ns']:ax.axvline(arrival,color='0.7',lw=.7,ls=':')
    fig.suptitle('Figure 0007-95-1: Fixed first mother, H→W 40 ns; true-energy diagnostic')
    savefig(fig,'fig_0007_95_1_full_window_projection')

    fig,axes=plt.subplots(1,2,figsize=(13,5.4),layout='constrained')
    for ax,gap,letter in zip(axes,(8.,40.),('a','b')):
        for block,color in zip((5,15,25,35),COLORS):
            values=[]
            for method,label in [('uniform','original_3starts'),('spatial12','original_3starts'),
                                 ('spatial12','diagonal_3starts'),('spatial12','unit_fir_3starts')]:
                row=next(r for r in previous+fits if r['block']==block and r['method']==method
                         and r['condition']=='high_weak' and r['gap_ns']==gap and variant(r)==label)
                values.append(row['weak_bias_percent'])
            ax.plot(range(4),values,'o-',color=color,label=f'Mother block {block}',lw=1)
        ax.axhspan(-10,10,color='0.9',zorder=0);ax.axhline(0,color='0.3',lw=.7)
        ax.set_xticks(range(4),('U/full','S/full','S/diagonal','S/unit FIR'),rotation=15)
        ax.set(title=f'({letter}) H→W {gap:g} ns',ylabel='Weak energy error [%]')
    axes[0].legend(fontsize=8)
    fig.suptitle('Figure 0007-96-1: Fixed weighting contrasts, same means and three voltage-derived starts')
    savefig(fig,'fig_0007_96_1_weight_contrasts')

    fig,axes=plt.subplots(1,2,figsize=(14,5.6),layout='constrained')
    variants=['true_E_start']+['delete_'+term for term in TERMS]
    for ax,gap,letter in zip(axes,(8.,40.),('a','b')):
        for block,color in zip((5,15,25,35),COLORS):
            values=[next(r['weak_bias_percent'] for r in fits if r['block']==block
                and r['method']=='spatial12' and r['condition']=='high_weak' and r['gap_ns']==gap
                and variant(r)==label) for label in variants]
            ax.plot(range(5),values,'o-',color=color,label=f'Mother block {block}',lw=1)
        ax.axhspan(-10,10,color='0.9',zorder=0);ax.axhline(0,color='0.3',lw=.7)
        ax.set_xticks(range(5),('Raw mean','−Strong iso.','−Weak iso.','−Shared','−Dark/elec.'),rotation=20)
        ax.set(title=f'({letter}) H→W {gap:g} ns',ylabel='Weak energy error [%]')
    axes[0].legend(fontsize=8)
    fig.suptitle('Figure 0007-97-1: Oracle removal of one residual term; fixed full weights, true-E start')
    savefig(fig,'fig_0007_97_1_oracle_deletions')


def main():
    result=tables();plots()
    print(json.dumps(dict(formal_fits=result['formal_fits'],nonconverged=result['nonconverged'],
                          search=result['search']),ensure_ascii=False))


if __name__=='__main__':main()
