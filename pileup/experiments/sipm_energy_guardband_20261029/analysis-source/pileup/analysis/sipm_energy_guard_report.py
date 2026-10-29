"""Figures for guard-band calibration, independent singles and historical replay."""
from __future__ import annotations

import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import beta

from .sipm_energy_guard_data import PROJECT, EXPERIMENT, OUTPUT, OLD, load, save, digest
from .sipm_energy_guard_study import lines, FIT_METHODS, DIAGNOSTIC
from .sipm_energy_calibration import invert_charge, error_summary

LABELS={'linear_energy':'Linear E','yield_scale':'Scaled light yield',
        'yield_saturation':'Saturating light yield','empirical_pchip':'Monotone PCHIP',
        'model_relative_residual':'Model GLS','joint_conserved_nll':'Joint NLL',DIAGNOSTIC:'True Q diagnostic'}
CONDITIONS=['weak_weak','low_weak','weak_low','high_weak','weak_high','low_low','mid_mid','high_high','low_high','high_low']
GAPS=[8.,40.,200.]


def csv_file(name,rows):
    with (OUTPUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def sensitivity_table():
    """Local error propagation for the frozen curve, not a wide-error bound."""
    model=load(EXPERIMENT/'energy-model.json');rows=load(OUTPUT/'calibration.records.json')
    e=np.asarray(model['energy_MeV']);q=np.asarray(model['charge_pe']);slope=np.gradient(q,e)
    assert np.all(slope>0)
    result=[]
    levels=sorted({r['source_energy_MeV'] for r in rows if r['stage']=='select'})
    for level in levels:
        er=float(np.mean([r['edep_MeV'] for r in rows if r['stage']=='select' and r['source_energy_MeV']==level]))
        qr=float(np.interp(er,e,q));gain=float(qr/er/np.interp(er,e,slope))
        result.append(dict(source_energy_MeV=level,reference_edep_MeV=er,curve_charge_pe=qr,
            relative_inverse_gain=gain,local_Q_tolerance_for_10percent_E=.1/gain))
    csv_file('local-energy-sensitivity.csv',result)
    return result


def report_figures():
    report=load(OUTPUT/'report.json');cal=load(OUTPUT/'calibration.records.json');test=load(OUTPUT/'confirmation.labels.json')
    pilot=load(OUTPUT/'pilot.records.json');models=load(OUTPUT/'candidate-models.json');chosen=report['chosen']
    # Same previously exposed selection events, to avoid conflating sample mix
    # with the effect of wider calibration data. This is not a new test.
    common=[r for r in cal if r.get('fold')==1];assert len(common)==280
    boundary=[]
    for version,candidates in [('previous',load(OLD/'candidate-models.json')),('guard_band',models)]:
        for method,model in candidates.items():
            for level in [1.625,2.5,5.486,8.]:
                rr=[r for r in common if r['source_energy_MeV']==level]
                stats=error_summary([r['edep_MeV'] for r in rr],invert_charge(model,[r['charge_pe'] for r in rr]))
                boundary.append(dict(calibration=version,method=method,source_energy_MeV=level,**stats))
    csv_file('common-selection-comparison.csv',boundary)
    sensitivity_table()
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'figure.dpi':120,'savefig.dpi':160,'svg.fonttype':'none'})
    figures=[]
    def finish(fig,name):
        for suffix in ['.png','.svg']:fig.savefig(OUTPUT/(name+suffix),bbox_inches='tight')
        plt.close(fig);figures.append(name)
    stages=[('fit','Fit: 376','#527fa3'),('select','Select: 344','#e5a046'),('confirmation','New confirmation: 320','#338665')]
    allrows=cal+test;levels=sorted({r['source_energy_MeV'] for r in allrows})
    fig,ax=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
    for j,(stage,label,color) in enumerate(stages):
        rr=[r for r in allrows if r['stage']==stage]
        ax[0].scatter([r['source_energy_MeV'] for r in rr],[r['edep_MeV'] for r in rr],s=12,alpha=.45,color=color,label=label)
        count=[sum(r['source_energy_MeV']==e for r in rr) for e in levels]
        ax[1].bar(np.arange(len(levels))+(j-1)*.23,count,width=.23,color=color,label=label)
    ax[0].scatter([r['source_energy_MeV'] for r in pilot],[r['edep_MeV'] for r in pilot],s=28,marker='x',color='black',label='Pilot: 16 (not used for fitting)')
    ax[0].set(xlabel='Source energy / MeV',ylabel='Deposited energy / MeV',title='(a) Native transport through the front foil')
    ax[0].legend(fontsize=8)
    ax[1].set_xticks(range(len(levels)),[f'{e:g}' for e in levels],rotation=45,ha='right')
    ax[1].set(xlabel='Source energy / MeV',ylabel='Independent optical parents',title='(b) Guard points 1.5/9; unseen energies 1.8/7')
    fig.suptitle('Fig. 0007-8-1 | Wider calibration support and independent sampling')
    finish(fig,'0007-8-1_coverage')

    fig,ax=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
    for stage,label,color in stages[:2]:
        rr=[r for r in cal if r['stage']==stage]
        ax[0].scatter([r['edep_MeV'] for r in rr],[r['charge_pe'] for r in rr],s=12,alpha=.35,label=label,color=color)
    for name,m in models.items():ax[0].plot(m['energy_MeV'],m['charge_pe'],lw=2 if name==chosen else 1,label=LABELS[name]+(' [selected]' if name==chosen else ''))
    ax[0].set(xlabel='Deposited energy / MeV',ylabel='Primary realized charge / equivalent pe',title='(a) Frozen response mappings')
    ax[0].legend(fontsize=8)
    choices=[r for r in report['mapping_selection']['candidates'] if r['valid']]
    losses=[r['selection_loss'] for r in choices]
    ax[1].scatter(losses,np.arange(len(choices)),s=65,
                  color=['#c65038' if r['method']==chosen else '#527fa3' for r in choices],zorder=3)
    ax[1].set_yticks(range(len(choices)),[LABELS[r['method']] for r in choices])
    for i,value in enumerate(losses):ax[1].annotate(f'{value:.6g}',(value,i),xytext=(8,0),textcoords='offset points',va='center',fontsize=9)
    ax[1].set_xscale('log');ax[1].set_xlim(min(losses)*.7,max(losses)*3)
    ax[1].set_ylim(len(choices)-.5,-.5);ax[1].grid(axis='x',alpha=.15)
    ax[1].set(xlabel='Internal-target relative-MSE loss, including rejection',title='(b) Selection on 344 separate parents')
    fig.suptitle('Fig. 0007-9-1 | Calibration chosen before new confirmation')
    finish(fig,'0007-9-1_selection')

    controls=[]
    for r in report['single_confirmation']:
        n,k=r['n'],r['valid']
        controls.append(dict(r,valid_fraction=k/n,valid_fraction_low95=float(beta.ppf(.025,k,n-k+1)) if k else 0.,
            valid_fraction_high95=float(beta.ppf(.975,k+1,n-k)) if k<n else 1.))
    csv_file('independent-single-controls.csv',controls)
    fig,ax=plt.subplots(1,2,figsize=(12.5,4.8),layout='constrained')
    for j,(method,label,color) in enumerate([('old_four_band_calibration','Previous calibration','#929292'),
                                           ('new_guard_calibration','Guard-band calibration','#c65038')]):
        rr=[r for r in controls if r['method']==method];x=np.arange(len(rr))+(j-.5)*.2
        median=np.array([r['median_relative'] if r['median_relative'] is not None else np.nan for r in rr])*100
        ci=np.array([[v if v is not None else np.nan for v in r['median_bootstrap95']] for r in rr])*100
        ax[0].errorbar(x,median,yerr=[median-ci[:,0],ci[:,1]-median],marker='o',ls='none',capsize=3,color=color,label=label)
        ax[1].bar(x,[r['within_10']/r['n'] for r in rr],width=.2,color=color,label=label)
    energies=[r['source_energy_MeV'] for r in controls if r['method']=='new_guard_calibration']
    ticks=[f'{e:g}'+('*' if e in [1.8,7.] else '') for e in energies]
    for a in ax:a.set_xticks(range(len(energies)),ticks,rotation=45,ha='right');a.set_xlabel('Source energy / MeV; * = unseen energy in calibration')
    ax[0].axhline(0,color='gray',lw=.8);ax[0].axhspan(-5,5,color='gray',alpha=.12)
    ax[0].set(ylabel='Median relative E error / %',title='(a) Median and conditional bootstrap 95% interval')
    ax[0].legend(fontsize=8)
    ax[1].set(ylabel='Fraction within 10% E, denominator 32 per group',ylim=(0,1.08),title='(b) Known isolated primary Q; not blind detection')
    fig.suptitle('Fig. 0007-10-1 | Fresh independent energy-calibration confirmation')
    finish(fig,'0007-10-1_independent_confirmation')

    summaries=report['historical_summary'];lookup={(r['method'],r['kind'],r['condition'],r['gap_ns']):r for r in summaries}
    fig,axes=plt.subplots(1,3,figsize=(12,6.4),layout='constrained')
    for panel,(ax,method) in enumerate(zip(axes,FIT_METHODS+[DIAGNOSTIC])):
        values=np.array([[lookup[method,'pair',c,g]['frame_success_10'] for g in GAPS] for c in CONDITIONS])
        im=ax.imshow(values/40,vmin=0,vmax=1,cmap='viridis',aspect='auto')
        for i in range(10):
            for j in range(3):ax.text(j,i,f'{values[i,j]}/40',ha='center',va='center',color='white' if values[i,j]<23 else 'black',fontsize=9)
        ax.set_xticks(range(3),[str(int(g)) for g in GAPS]);ax.set_xlabel('Separation / ns')
        ax.set_yticks(range(10),[c.replace('_',' -> ') for c in CONDITIONS]);ax.set_title(f'({chr(97+panel)}) {LABELS[method]}')
    fig.colorbar(im,ax=axes,label='Both E within 10%, timing tolerance satisfied, no extras',shrink=.8)
    fig.suptitle('Fig. 0007-11-1 | Historical fixed fits with the new energy mapping')
    finish(fig,'0007-11-1_historical_map')
    csv_file('historical-condition-summary.csv',summaries)
    pair_comparison=[]
    for method in FIT_METHODS+[DIAGNOSTIC]:
        for group in ['weak','strong']:
            before=report['previous_aggregate'][method][group];after=report['historical_aggregate'][method][group]
            pair_comparison.append(dict(method=method,group=group,old_calibration=before,new_calibration=after,difference=after-before))
    csv_file('historical-calibration-comparison.csv',pair_comparison)

    cases=load(OLD/'evaluation.cases.json')
    with np.load(OLD/'evaluation.labels.npz') as f:energy=f['energy_MeV'];charge=f['charge_pe']
    newpred=lines(OUTPUT/'historical-energy.predictions.jsonl');oldpred=lines(OLD/'energy.predictions.jsonl')
    scores=lines(OUTPUT/'historical-energy.scores.jsonl')
    fig,axes=plt.subplots(3,1,figsize=(11,9),layout='constrained',sharex=True);spectra=[];bins=np.linspace(0,9,91)
    for panel,(ax,gap) in enumerate(zip(axes,GAPS)):
        ids={i for i,c in enumerate(cases) if c['kind']=='pair' and c['gap_ns']==gap};ix=sorted(ids)
        truth=energy[ix][charge[ix]>=500]
        ax.hist(truth,bins=bins,histtype='step',lw=2,color='black',label=f'True E: {len(truth)} events')
        for method,color in zip(FIT_METHODS,['#3d70a8','#c65038']):
            pp=[r for r in newpred if r['method']==method and r['index'] in ids]
            vals=[e for r in pp for e in r['energy_MeV'] if e is not None];invalid=sum(e is None for r in pp for e in r['energy_MeV'])
            ss=[s for s in scores if s['method']==method and s['index'] in ids]
            spectra.append(dict(method=method,gap_ns=gap,truth_events=len(truth),outputs=sum(s['output_count'] for s in ss),
                finite_energy=len(vals),out_of_domain=invalid,missed=sum(s['missed'] for s in ss),extra=sum(s['extra'] for s in ss),
                strict_energy10=sum(s['strict_10'] for s in ss)))
            ax.hist(vals,bins=bins,histtype='step',lw=1.4,color=color,label=f'{LABELS[method]} new E: {len(vals)} finite, {invalid} invalid')
        vals=[e for p in oldpred if p['method']=='joint_conserved_nll' and p['index'] in ids for e in p['energy_MeV'] if e is not None]
        ax.hist(vals,bins=bins,histtype='step',lw=1,color='gray',linestyle='--',label='Same NLL fits, previous energy mapping')
        ax.set(ylabel='Counts / 0.1 MeV',title=f'({chr(97+panel)}) Separation {gap:g} ns; controlled mixture, not a source-rate scan')
        ax.legend(fontsize=8);ax.grid(alpha=.13)
    axes[-1].set_xlabel('Deposited energy / MeV')
    fig.suptitle('Fig. 0007-12-1 | Absolute counts after calibration update')
    finish(fig,'0007-12-1_energy_counts')
    csv_file('historical-spectrum-counts.csv',spectra)
    save(OUTPUT/'delivery-summary.json',dict(figures=figures,single_controls=controls,historical_comparison=pair_comparison,
        spectra=spectra,figure_sha256={f:digest(OUTPUT/(f+'.png')) for f in figures},
        confirmation='Fresh isolated primary Q, not blind waveform inference',historical='Same Q/t/count fits, new energy mapping'))
    print('GUARD_FIGURES_COMPLETE',figures,flush=True)


if __name__=='__main__':report_figures()
