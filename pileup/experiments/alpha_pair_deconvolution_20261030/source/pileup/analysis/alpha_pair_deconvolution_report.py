"""Score the frozen pilot separately from inference; retain every component."""
from __future__ import annotations

from collections import defaultdict
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .alpha_pair_deconvolution import AgeRecovery, EnergyProfiles, FiniteConvolution, spe_kernel
from .alpha_pair_deconvolution_study import ARMS, ART, EXP, NAME, OLDART, REC, ROOT, hashes, read_json, sha, write_json

BASELINE = 'new_mean_new_noise'
METHODS = (BASELINE, *ARMS)
LABELS = dict(new_mean_new_noise='Calibrated voltage GLS', rl128='RL 128', recovery128='RL 128 + recovery',
              rl512='RL 512', recovery512='RL 512 + recovery')
COLORS = ['#20252b', '#3389b8', '#db7046', '#79bdb5', '#b989ad']
CONDITIONS = ['weak', 'low', 'mid', 'high', 'high_weak_8', 'high_weak_40', 'high_weak_200', 'weak_high_40']
DISPLAY = ['Single weak', 'Single low', 'Single mid', 'Single high', 'High→weak 8 ns', 'High→weak 40 ns',
           'High→weak 200 ns', 'Weak→high 40 ns']


def condition(row):
    return row['condition'] if row['kind'] == 'single' else row['condition']+'_'+str(int(row['gap_ns']))


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n'); w.writeheader(); w.writerows(rows)


def save(fig, number, name):
    fig.savefig(ART/f'0007-{number}-1_{name}.png', dpi=160, bbox_inches='tight')
    fig.savefig(ART/f'0007-{number}-1_{name}.svg', bbox_inches='tight')
    plt.close(fig)


def main():
    labels = read_json(ART/'labels.json')
    rows, components, outcomes = [], [], {}
    for label in labels:
        i = label['id']; result = read_json(ART/'cases'/f'{i:04d}.json')
        old = read_json(OLDART/'cases'/f'{i:04d}.json')['fits'][BASELINE]
        fits = {BASELINE: old, **result['fits']}
        assert result['worker_truth_read_guard'] and result['known_K_t'] and not result['truth_E_Q_read']
        for method, fit in fits.items():
            error = np.asarray(fit['energy_MeV'])/label['energy_MeV']-1.
            key = condition(label); success = bool(np.all(np.abs(error) <= .1))
            outcomes[i, method] = success
            rows.append(dict(id=i, condition=key, method=method, count=len(error), success5=bool(np.all(np.abs(error)<=.05)),
                             success10=success, success20=bool(np.all(np.abs(error)<=.2)), converged=fit['converged'], boundary=fit['boundary']))
            for j, value in enumerate(error):
                components.append(dict(id=i, condition=key, method=method, component=j, true_energy_MeV=label['energy_MeV'][j],
                    fitted_energy_MeV=fit['energy_MeV'][j], relative_error=float(value), weak=(j==int(np.argmin(label['energy_MeV'])))))
    summary = []
    for key in CONDITIONS:
        for method in METHODS:
            frames = [r for r in rows if r['condition']==key and r['method']==method]
            errors = np.array([r['relative_error'] for r in components if r['condition']==key and r['method']==method and r['weak']])*100
            ids = [r['id'] for r in frames]
            summary.append(dict(condition=key, method=method, frames=len(frames),
                success5=sum(r['success5'] for r in frames), success10=sum(r['success10'] for r in frames), success20=sum(r['success20'] for r in frames),
                gained=sum(outcomes[i, method] and not outcomes[i, BASELINE] for i in ids),
                lost=sum(not outcomes[i, method] and outcomes[i, BASELINE] for i in ids),
                bias_pct=float(np.mean(errors)), median_pct=float(np.median(errors)), sigma68_pct=float((np.percentile(errors,84)-np.percentile(errors,16))/2),
                rmse_pct=float(np.sqrt(np.mean(errors**2))), boundary=sum(r['boundary'] for r in frames), nonconverged=sum(not r['converged'] for r in frames)))
    lookup = {(r['condition'], r['method']): r for r in summary}
    oldpairs = sum(lookup[k, BASELINE]['success10'] for k in CONDITIONS[4:])
    newpairs = sum(lookup[k, 'recovery128']['success10'] for k in CONDITIONS[4:])
    singles_ok = all(lookup[k, 'recovery128']['success10'] >= lookup[k, BASELINE]['success10']-1 for k in CONDITIONS[:4])
    pair_ok = newpairs >= oldpairs
    elapsed = read_json(REC/'evaluation-run.json')['seconds']
    gate = dict(passed=bool(pair_ok and singles_ok and elapsed*480/64*1.5 < 1800),
                pairs_ok=pair_ok, singles_ok=singles_ok, baseline_pair_success=oldpairs, primary_pair_success=newpairs,
                estimated_480_seconds=elapsed*480/64*1.5,
                note='Gate is a pilot decision, not a significance test. No parameters tuned on these outcomes.')
    diagnostics = [read_json(ART/'cases'/f"{r['id']:04d}.json") for r in labels]
    report = dict(scope='Known TRUE K/t, pure alpha; pilot64 with historical exposure and mother reuse', summary=summary,
        expansion_gate=gate, evaluated_frames=len(labels), new_fits=4*len(labels), reference_fits=len(labels),
        inverse_failures=sum(not d['converged'] for r in diagnostics for d in r['recovery'].values()),
        new_fit_nonconvergence=sum(not r['converged'] for r in rows if r['method']!=BASELINE),
        inverse_correction_range=[min(d['correction'] for r in diagnostics for d in r['recovery'].values()),max(d['correction'] for r in diagnostics for d in r['recovery'].values())],
        components=len(components), unknown_K_t_metrics='Not measured; no claim of zero misses/extras/split/merge', primary='recovery128')
    write_json(ART/'report.json', report)
    tables = EXP/'tables'; tables.mkdir(exist_ok=True)
    write_csv(tables/'summary.csv', summary); write_csv(tables/'frames.csv', rows); write_csv(tables/'components.csv', components)
    write_json(EXP/'pilot-decision.json', gate)

    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)
    for ax, keys, subtitle in zip(axes, (CONDITIONS[:4], CONDITIONS[4:]), ('(a) Isolated controls', '(b) Known-count/time pairs')):
        for j, method in enumerate(METHODS):
            values = [lookup[k, method]['success10'] for k in keys]
            x = np.arange(len(keys))+(j-2)*.15
            ax.bar(x, values, width=.145, color=COLORS[j], label=LABELS[method])
            for p, v in zip(x, values): ax.text(p, v+.08, str(v), ha='center', fontsize=8)
        ax.set(xticks=np.arange(len(keys)), xticklabels=[DISPLAY[CONDITIONS.index(k)] for k in keys], ylim=(0, 9), ylabel='All E within ±10% / 8', title=subtitle)
        ax.tick_params(axis='x', rotation=18); ax.grid(axis='y', alpha=.2)
    fig.legend(*axes[0].get_legend_handles_labels(), loc='outside lower center', ncol=5, fontsize=9)
    fig.suptitle('Figure 0007-58-1: Energy recovery pilot — fixed K/t, all cases retained', y=1.02)
    save(fig, 58, 'success')

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for panel, (ax, key) in enumerate(zip(axes.flat, CONDITIONS[4:])):
        values = [[r['relative_error']*100 for r in components if r['condition']==key and r['method']==m and r['weak']] for m in METHODS]
        boxes = ax.boxplot(values, patch_artist=True, tick_labels=['GLS','RL128','+ recovery','RL512','+ recovery'], showfliers=True)
        for patch, color in zip(boxes['boxes'], COLORS): patch.set_facecolor(color); patch.set_alpha(.5)
        for j, val in enumerate(values): ax.scatter(np.full(len(val), j+1), val, color=COLORS[j], s=10, alpha=.7)
        ax.axhspan(-10, 10, color='green', alpha=.1); ax.axhline(0, color='grey', lw=.7)
        ax.set(title=f'({chr(97+panel)}) '+DISPLAY[CONDITIONS.index(key)], ylabel='Weak E error [%]'); ax.tick_params(axis='x', rotation=15)
    fig.suptitle('Figure 0007-59-1: Weak-component energy errors; independent y scales')
    save(fig, 59, 'weak_errors')

    models = {}
    for arm in ARMS:
        with np.load(EXP/(arm+'-model.npz')) as f: models[arm] = EnergyProfiles(**{k:f[k] for k in f.files})
    with np.load(ART/'transformed.npz') as f: transformed = {k:f[k] for k in f.files}
    with np.load(ART/'inputs/observations.npz') as f: voltage=f['voltage_mV']; ids=f['ids']
    latent_ns=(np.arange(1024)+.5)*2.; sampled_ns=np.arange(1024)*2.+.25
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    for col, identifier in enumerate((40, 120)):
        idx=list(ids).index(identifier); label=labels[idx]
        for method in ('rl128','recovery128'):
            axes[0,col].plot(latent_ns, transformed[method][idx], label=LABELS[method], lw=.9)
        axes[0,col].set(title=f'({chr(97+col)}) id {identifier}, {condition(label)}', ylabel='Equivalent counts per 2 ns bin', xlim=(0,2048))
        axes[0,col].legend(fontsize=9)
        for iteration in (128,512):
            ratio=np.divide(transformed['recovery'+str(iteration)][idx], transformed['rl'+str(iteration)][idx],
                            out=np.ones(1024), where=transformed['rl'+str(iteration)][idx]>1e-10)
            axes[1,col].plot(latent_ns, ratio, label=f'{iteration} RL iterations')
        axes[1,col].set(title=f'({chr(99+col)}) Correction factor', xlabel='Time from window start [ns]', ylabel='Recovery-corrected / released', xlim=(0,2048)); axes[1,col].legend()
    fig.suptitle('Figure 0007-60-1: Released-charge vs PDE-equivalent-light coordinates (different meanings)')
    save(fig, 60, 'transformed')

    configuration=read_json(EXP/'sensor-config.json'); op=FiniteConvolution(spe_kernel(configuration),1024)
    recovery=AgeRecovery(configuration['sensor']); background=read_json(EXP/'protocol.json')['background_mV']
    with np.load(OLDART/'labels.components.npz') as f: truth_components=f['components_mV'][ids]
    fig, axes=plt.subplots(3,2,figsize=(14,10),constrained_layout=True)
    for col, identifier in enumerate((40,120)):
        idx=list(ids).index(identifier); label=labels[idx]; result=diagnostics[idx]; weak=int(np.argmin(label['energy_MeV']))
        axes[0,col].plot(sampled_ns,voltage[idx]-background,color='k',lw=.8,label='Observed - background')
        axes[1,col].plot(sampled_ns,truth_components[idx,weak],color='k',lw=.8,label='True weak contribution (score only)')
        for method,color in zip(('rl128','recovery128'),COLORS[1:3]):
            model=models[method]; fit=result['fits'][method]
            curves=model.components(fit['times_ns'],fit['energy_MeV']); total=curves.sum(0)+model.background
            released=recovery.forward(total) if method.startswith('recovery') else total
            signal=op.forward(released)
            weak_released=released*np.divide(curves[weak],total,out=np.zeros_like(total),where=total>0)
            weak_voltage=op.forward(weak_released)
            axes[0,col].plot(sampled_ns,signal,label=LABELS[method],color=color,lw=1.)
            axes[1,col].plot(sampled_ns,weak_voltage,label=LABELS[method],color=color,lw=1.)
            axes[2,col].plot(sampled_ns,voltage[idx]-background-signal,label=LABELS[method],color=color,lw=.7,alpha=.8)
        for row in range(3):
            axes[row,col].set(xlim=(0,2048),ylabel=['Signal [mV]','Weak contribution [mV]','Residual [mV]'][row]); axes[row,col].legend(fontsize=8)
        axes[0,col].set_title(f'({chr(97+col)}) id {identifier}: {condition(label)}; true E={label["energy_MeV"][0]:.3f}, {label["energy_MeV"][1]:.3f} MeV')
        axes[1,col].set_title(f'({chr(99+col)}) Weak component')
        axes[2,col].set_title(f'({chr(101+col)}) Observed minus refolded model')
        axes[2,col].set_xlabel('Time from window start [ns]')
    fig.suptitle('Figure 0007-61-1: Full windows and refolded temporal decompositions')
    save(fig,61,'full_waveforms')
    sources=[ROOT/'pileup/analysis/alpha_pair_deconvolution_report.py']
    outputs=list(tables.glob('*.csv'))+[ART/'report.json',EXP/'pilot-decision.json']+list(ART.glob('0007-*.png'))+list(ART.glob('0007-*.svg'))
    write_json(EXP/'report-freeze.json',hashes(sources+outputs))
    print(json.dumps(dict(gate=gate, inverse_failures=report['inverse_failures'], fit_nonconvergence=report['new_fit_nonconvergence']),indent=2))
    for k in CONDITIONS: print(k,[(m,lookup[k,m]['success10'],round(lookup[k,m]['bias_pct'],2),round(lookup[k,m]['rmse_pct'],2)) for m in METHODS])


if __name__=='__main__': main()
