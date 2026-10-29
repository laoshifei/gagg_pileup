"""Figures and diagnostic tables for the frozen SiPM energy-transfer study.

Run from the project root. This module reads frozen outputs; it never fits or
chooses a calibration, event-count gate, or waveform model.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .sipm_energy_calibration import invert_charge, error_summary
from .sipm_energy_study import PROJECT, NAME, OUTPUT, EXPERIMENT, RECORD, FIT_METHODS, load, save, digest, lines


LABELS = {'linear_energy':'Linear E', 'yield_scale':'Scaled light yield',
          'yield_saturation':'Saturating light yield', 'empirical_pchip':'Monotone PCHIP',
          'model_relative_residual':'Model GLS', 'joint_conserved_nll':'Joint NLL',
          'true_charge_diagnostic':'True Q diagnostic'}
METHODS = FIT_METHODS + ['true_charge_diagnostic']
COLORS = ['#3d70a8', '#c65038', '#348c70']
CONDITIONS = ['weak_weak','low_weak','weak_low','high_weak','weak_high',
              'low_low','mid_mid','high_high','low_high','high_low']
GAPS = [8.,40.,200.]
LEVELS = [1.625,2.5,5.486,8.]


def write_csv(name,rows):
    with (OUTPUT/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def make_report():
    report=load(OUTPUT/'report.json'); model=load(EXPERIMENT/'energy-model.json')
    calibration=load(OUTPUT/'calibration-records.json')
    parents=load(OUTPUT/'evaluation-parent-labels.json')
    cases=load(OUTPUT/'evaluation.cases.json'); scores=lines(OUTPUT/'energy.scores.jsonl')
    predictions=lines(OUTPUT/'energy.predictions.jsonl')
    models=load(OUTPUT/'candidate-models.json')
    summary=report['summary'];lookup={(r['method'],r['kind'],r['condition'],r['gap_ns']):r for r in summary}
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'figure.dpi':120,'savefig.dpi':160,'svg.fonttype':'none'})
    figures=[]
    def finish(fig,name):
        fig.savefig(OUTPUT/(name+'.png'),bbox_inches='tight')
        fig.savefig(OUTPUT/(name+'.svg'),bbox_inches='tight')
        plt.close(fig);figures.append(name)

    # Calibration selection: both axes retain independent selection failures.
    fig,ax=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
    for fold,marker in [(0,'o'),(1,'x')]:
        rows=[r for r in calibration if r['fold']==fold]
        ax[0].scatter([r['edep_MeV'] for r in rows],[r['charge_pe'] for r in rows],s=12,
                      alpha=.35,marker=marker,color='gray',label=['Fit: 280','Select: 280'][fold])
    for name,m in models.items():
        ax[0].plot(m['energy_MeV'],m['charge_pe'],lw=2 if name==model['method'] else 1,
                   label=LABELS[name]+(' [selected]' if name==model['method'] else ''))
    ax[0].set(xlabel='Deposited energy / MeV',ylabel='Primary realized charge / equivalent pe',title='(a) Independent single-event calibration')
    ax[0].legend(fontsize=8)
    valid=[r for r in report['mapping_selection'] if r['valid']]
    names=[LABELS[r['method']] for r in valid]
    ax[1].barh(names,[r['selection_loss'] for r in valid],color=['#c65038' if r['method']==model['method'] else '#9eb6ce' for r in valid])
    ax[1].set(xlabel='Selection loss: relative MSE + domain rejection penalty',title='(b) Frozen selection rule; lower is better')
    ax[1].invert_yaxis()
    fig.suptitle('Fig. 0007-2-1 | Charge-to-energy calibration and selection')
    finish(fig,'0007-2-1_calibration')

    # Single-event control: primary Q conversion versus blind estimated Q.
    control=[]
    for method in METHODS:
        for level in LEVELS:
            rr=[r for r in scores if r['method']==method and r['kind']=='single'
                and parents[cases[r['index']]['parents'][0]]['source_energy_MeV']==level]
            errors=[m['relative_energy_error'] for r in rr for m in r['matches'] if m['relative_energy_error'] is not None]
            control.append(dict(method=method,source_MeV=level,frames=len(rr),
                target=sum(r['target'] for r in rr),matched=sum(r['matched'] for r in rr),
                recovered10=sum(r['strict_10'] for r in rr),frame_success10=sum(r['frame_success_10'] for r in rr),
                missed=sum(r['missed'] for r in rr),extra=sum(r['extra'] for r in rr),
                out_of_domain=sum(r['energy_out_of_domain'] for r in rr),
                median_percent=100*float(np.median(errors)) if errors else None,
                sigma68_percent=50*float(np.quantile(errors,.84)-np.quantile(errors,.16)) if errors else None))
    write_csv('single-controls.csv',control)
    fig,ax=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for j,(method,color) in enumerate(zip(METHODS,COLORS)):
        rr=[r for r in control if r['method']==method];x=np.arange(4)+(j-1)*.17
        ax[0].errorbar(x,[r['median_percent'] for r in rr],yerr=[r['sigma68_percent'] for r in rr],
                       marker='o',ls='none',capsize=3,color=color,label=LABELS[method])
        ax[1].bar(x,[r['recovered10']/r['target'] for r in rr],width=.17,color=color,label=LABELS[method])
    ax[0].axhline(0,color='gray',lw=.8);ax[0].axhspan(-5,5,color='gray',alpha=.1)
    ax[0].set(ylabel='Conditional median relative E error +/- sigma68 / %',title='(a) Matched, finite estimates only')
    ax[1].set(ylabel='Strict event recovery within 10% E',ylim=(0,1.08),title='(b) All 80 true events per source group')
    for a in ax:a.set_xticks(range(4),[str(x) for x in LEVELS]);a.set_xlabel('Source energy / MeV (not deposited energy)')
    ax[1].legend(fontsize=8)
    fig.suptitle('Fig. 0007-3-1 | Independent single-event energy controls')
    finish(fig,'0007-3-1_single_controls')

    # Same 1200 pairs, same fitting gates; each cell has 40 pairs.
    fig,axes=plt.subplots(1,3,figsize=(12,6.3),layout='constrained')
    for ax,method in zip(axes,METHODS):
        counts=np.array([[lookup[method,'pair',c,g]['frame_success_10'] for g in GAPS] for c in CONDITIONS])
        im=ax.imshow(counts/40,vmin=0,vmax=1,cmap='viridis',aspect='auto')
        for i in range(10):
            for j in range(3):ax.text(j,i,f'{counts[i,j]}/40',ha='center',va='center',fontsize=9,color='white' if counts[i,j]<23 else 'black')
        ax.set_xticks(range(3),[str(int(g)) for g in GAPS]);ax.set_xlabel('Separation / ns')
        ax.set_yticks(range(10),[c.replace('_',' -> ') for c in CONDITIONS]);ax.set_title(LABELS[method])
    fig.colorbar(im,ax=axes,label='Both events within 10% E and timing tolerance; no extras',shrink=.8)
    fig.suptitle('Fig. 0007-4-1 | Energy recovery over charge order and separation')
    finish(fig,'0007-4-1_condition_map')

    fig,ax=plt.subplots(1,2,figsize=(11.5,4.4),layout='constrained')
    groups=['weak','strong']
    for j,(method,color) in enumerate(zip(METHODS,COLORS)):
        ax[0].bar(np.arange(2)+(j-1)*.22,[report['aggregate_strict_frame10'][method][g] for g in groups],
                  width=.22,color=color,label=LABELS[method])
    ax[0].set_xticks(range(2),['Pairs containing weak alpha','Other pairs'])
    ax[0].set(ylabel='Mean strict pair recovery within 10% E',ylim=(0,1.06),title='(a) Equal mean over 15 conditions per group')
    ax[0].legend(fontsize=8)
    for i,g in enumerate(groups):
        r=report['paired_energy10'][g];v=100*r['difference'];lo,hi=np.array(r['mother_block95'])*100
        ax[1].errorbar(v,i,xerr=[[v-lo],[hi-v]],fmt='o',color='#c65038',capsize=5)
    ax[1].axvline(0,color='gray',lw=.8);ax[1].set_yticks(range(2),['Contains weak alpha','Other pairs'])
    ax[1].set(xlabel='Joint NLL - Model GLS / percentage points',title='(b) Paired mother-block 95% bootstrap interval')
    fig.suptitle('Fig. 0007-4-2 | Paired energy comparison on identical waveforms')
    finish(fig,'0007-4-2_paired_comparison')

    # Deposited spectra use absolute counts; invalid energy outputs remain in table.
    with np.load(OUTPUT/'evaluation.labels.npz') as f:energy=f['energy_MeV'];charge=f['charge_pe']
    fig,axes=plt.subplots(3,1,figsize=(11,9),layout='constrained',sharex=True)
    spectrum=[];bins=np.linspace(0,8.5,86)
    for ax,g in zip(axes,GAPS):
        ids=[i for i,c in enumerate(cases) if c['kind']=='pair' and c['gap_ns']==g]
        truth=energy[ids][charge[ids]>=500]
        ax.hist(truth,bins=bins,histtype='step',color='black',lw=2,label=f'True E: {len(truth)} events')
        for method,color in zip(FIT_METHODS,COLORS):
            pp=[p for p in predictions if p['method']==method and p['index'] in ids]
            values=[e for p in pp for e in p['energy_MeV'] if e is not None]
            invalid=sum(e is None for p in pp for e in p['energy_MeV'])
            ss=[s for s in scores if s['method']==method and s['index'] in ids]
            spectrum.append(dict(method=method,gap_ns=g,truth_events=len(truth),outputs=sum(s['output_count'] for s in ss),
                finite_energy=len(values),out_of_domain=invalid,missed=sum(s['missed'] for s in ss),
                extra=sum(s['extra'] for s in ss),energy10_recovered=sum(s['strict_10'] for s in ss)))
            ax.hist(values,bins=bins,histtype='step',color=color,lw=1.5,
                    label=f'{LABELS[method]}: {len(values)} finite, {invalid} out of domain')
        ax.set(ylabel='Counts / 0.1 MeV',title=f'Separation {g:g} ns; controlled energy mixture, not a rate scan')
        ax.legend(fontsize=8);ax.grid(alpha=.15)
    axes[-1].set_xlabel('Deposited energy / MeV')
    fig.suptitle('Fig. 0007-5-1 | Absolute reconstructed energy counts')
    finish(fig,'0007-5-1_energy_counts')
    write_csv('spectrum-counts.csv',spectrum)
    write_csv('condition-summary.csv',summary)

    # Diagnose edge-domain selection from independent selection data alone.
    boundary=[]
    for name,m in models.items():
        qmin,qmax=m['charge_pe'][0],m['charge_pe'][-1]
        for level in LEVELS:
            rr=[r for r in calibration if r['fold']==1 and r['source_energy_MeV']==level]
            q=np.array([r['charge_pe'] for r in rr]);e=np.array([r['edep_MeV'] for r in rr])
            s=error_summary(e,invert_charge(m,q))
            boundary.append(dict(method=name,source_MeV=level,n=len(rr),below_Q_domain=int((q<qmin).sum()),
                above_Q_domain=int((q>qmax).sum()),**{k:v for k,v in s.items() if k!='n'}))
    write_csv('selection-boundary-diagnostic.csv',boundary)
    save(OUTPUT/'delivery-summary.json',dict(figures=figures,single_controls=control,spectrum_counts=spectrum,
        selection_boundary=boundary,warning='True-Q diagnostic uses known K/t and is not a blind result or a rigorous upper bound',
        figure_sha256={name:digest(OUTPUT/(name+'.png')) for name in figures}))
    print('ENERGY_REPORT_FIGURES',figures,flush=True)


if __name__=='__main__':make_report()
