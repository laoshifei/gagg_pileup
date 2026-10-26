"""Report weak-event recovery and false births without test-time selection."""
from __future__ import annotations

from collections import Counter
import csv
import gzip
import json
import math

import numpy as np

from . import alpha_pair_null_study as study
from .alpha_pair_conditional_study import evaluate_fit
from .alpha_pair_null_gate import select
from .alpha_pair_residual import choose

LABELS = dict(observed_penalty='Observed weights + penalty', model_penalty='Model weights + penalty',
    model_relative_residual='Model weights + residual-ratio gate',
    model_relative_signal='Model weights + signal-ratio gate',
    known_count_model='Known-count diagnostic')
COLORS = dict(observed_penalty='#cc3311', model_penalty='#228833',
    model_relative_residual='#0077bb', model_relative_signal='#aa3377', known_count_model='#888888')


def wilson(success, count):
    if count == 0: return None
    z = 1.959963984540054
    p = success/count; scale = 1.+z*z/count
    center = (p+z*z/(2.*count))/scale
    half = z*math.sqrt(p*(1.-p)/count+z*z/(4.*count*count))/scale
    return [max(0., center-half), min(1., center+half)]


def assemble(stage):
    observed = study.new_rows(stage)
    model = study.new_rows(stage, 'model')
    assert len(observed) == len(model)
    calibration = study.load(study.EXPERIMENT/'model-null-calibration.json')
    fits = dict(observed_penalty=[choose(r['candidates'], .001) for r in observed],
                model_penalty=[choose(r['candidates'], .001) for r in model])
    for mode, cal in calibration.items():
        fits['model_'+mode] = [select(r['candidates'], cal) for r in model]
    # Explicit diagnostic: the evaluator supplies true above-threshold K.
    # It never enters a blind output, calibration, or method selection.
    fits['known_count_model'] = [r['candidates'][min(int(np.count_nonzero(r['charge'] >= 500.)), len(r['candidates'])-1)]
                                for r in model]
    return model, fits


def summarize(rows):
    first = rows[0]
    pair = first['kind'] == 'pair'
    eligible = [r for r in rows if not pair or r['target'] == 2]
    successes = sum(r['success'] for r in eligible)
    extra = sum(r['extra'] > 0 for r in rows)
    errors = [e for r in rows for e in r['q_errors']]
    timing = [e for r in rows for e in r['t_errors']]
    total = sum(r['target'] for r in rows)
    counts = Counter(r['count'] for r in rows)
    def quantile(values):
        return np.percentile(values, [16., 50., 84.]).tolist() if len(values) else None
    return dict(method=first['method'], kind=first['kind'], condition=first['condition'], gap_ns=first['gap_ns'],
        frames=len(rows), eligible_frames=len(eligible), below_threshold_pairs=len(rows)-len(eligible),
        success=successes, efficiency=successes/len(eligible) if eligible else None,
        efficiency_all_frames=successes/len(rows), success_wilson95=wilson(successes, len(eligible)),
        extra_frames=extra, extra_fraction=extra/len(rows), extra_wilson95=wilson(extra, len(rows)),
        targets=total, recovered=sum(r['recovered'] for r in rows),
        event_efficiency=sum(r['recovered'] for r in rows)/total if total else None,
        missed=sum(r['missed'] for r in rows), not_recovered=total-sum(r['recovered'] for r in rows),
        extras=sum(r['extra'] for r in rows), split=sum(r['split'] for r in rows), merge=sum(r['merge'] for r in rows),
        q_error_p16_p50_p84=quantile(errors), q_error_std=float(np.std(errors)) if errors else None,
        time_error_p16_p50_p84_ns=quantile(timing),
        q1_p16_p50_p84_pe=quantile([r['true_q'][0] for r in rows if r['true_q'][0] > 0.]),
        q2_p16_p50_p84_pe=quantile([r['true_q'][1] for r in rows if r['true_q'][1] > 0.]),
        output_K0=counts[0], output_K1=counts[1], output_K2=counts[2],
        budget_pass=extra/len(rows) <= .05)


def report(stage='test'):
    if stage == 'test': study.check_freeze()
    rows, fits = assemble(stage)
    all_scores = []; groups = {}
    for method, outputs in fits.items():
        for row, fit in zip(rows, outputs):
            case = row['case']
            score = dict(evaluate_fit(fit, row['times'], row['charge'], case), index=row['index'],
                method=method, kind=case['kind'], condition=case['condition'], gap_ns=case['gap_ns'],
                true_q=row['charge'].tolist(), true_t=row['times'].tolist(),
                null_rank_p=fit.get('null_rank_p'), known_count=method == 'known_count_model')
            all_scores.append(score)
            groups.setdefault((method, case['kind'], case['condition'], case['gap_ns']), []).append(score)
    summary = [summarize(v) for v in groups.values()]
    selection = study.load(study.EXPERIMENT/'model-null-selection.json')
    chosen = selection['chosen']
    primary = 'model_'+chosen if chosen else None
    controls = [r for r in summary if r['method'] == primary and r['kind'] != 'pair']
    accepted = bool(controls) and all(r['budget_pass'] for r in controls) and all(
        r['efficiency'] >= .95 and abs(r['q_error_p16_p50_p84'][1]) <= .05
        for r in controls if r['kind'] == 'single')
    comparisons = []
    rng = np.random.default_rng(265990001)
    if primary:
        for (method, kind, condition, gap), new in groups.items():
            if method != primary or kind != 'pair': continue
            for baseline in ('observed_penalty', 'model_penalty'):
                old = groups[(baseline, kind, condition, gap)]
                assert [r['index'] for r in old] == [r['index'] for r in new]
                delta = np.array([int(b['success'])-int(a['success']) for a, b in zip(old,new) if b['target'] == 2])
                boot = delta[rng.integers(0,len(delta),(2000,len(delta)))].mean(1) if len(delta) else []
                comparisons.append(dict(baseline=baseline, primary=primary, condition=condition, gap_ns=gap,
                    n=len(delta), wins=int(np.sum(delta > 0)), losses=int(np.sum(delta < 0)),
                    mean_difference=float(np.mean(delta)) if len(delta) else None,
                    paired_bootstrap95=np.percentile(boot,[2.5,97.5]).tolist() if len(boot) else None))
    output = study.OUTPUT/stage
    with gzip.open(output/'scores.jsonl.gz','wt') as stream:
        for row in all_scores: stream.write(json.dumps(row, allow_nan=False)+'\n')
    with (output/'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0])); writer.writeheader(); writer.writerows(summary)
    result = dict(stage=stage, frames=len(rows), score_rows=len(all_scores), summary=summary,
        primary=primary, all_control_requirements_passed=accepted, selection=selection,
        paired_comparisons=comparisons, candidate_count=5, models_retrained=False,
        uncertainty='Within a condition/gap only; same optical parents recur across gaps and conditions',
        charge='Realized avalanche-equivalent pe; conditional matched-charge error, not energy resolution',
        pair_denominator='Both true Q>=500pe; all generated frames and below-threshold counts also shown',
        known_count_diagnostic='Evaluator chooses candidate K using true above-threshold count; not a blind method')
    study.save(output/'report.json', result)
    print('REPORT', stage, 'PRIMARY', primary, 'CONTROLS_PASS', accepted, flush=True)
    for row in controls:
        print('CONTROL', row['condition'], row['success'], '/', row['eligible_frames'],
              'extra', row['extra_frames'], 'charge', row['q_error_p16_p50_p84'], flush=True)
    return result


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .alpha_pair_conditional import ConditionalResponse
    calibration = study.load(study.EXPERIMENT/'model-null-calibration.json')
    result = study.load(study.OUTPUT/'test/report.json')
    primary = result['primary']; summary = result['summary']
    folder = study.OUTPUT/'figures'; folder.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    for ax, (mode, cal), letter in zip(axes, calibration.items(), 'ab'):
        for label, values in zip(('Weak-Q bin', 'Low-Q bin', 'Mid-Q bin', 'High-Q bin'), cal['scores']):
            values = np.asarray(values)
            survival = (len(values)-np.arange(len(values)))/(len(values)+1.)
            line, = ax.step(values, survival, where='post', label=f'{label}, n={len(values)}')
            rank = int(np.floor(cal['alpha']*(len(values)+1)))
            if rank:
                ax.axvline(values[-rank], color=line.get_color(), lw=.8, ls=':', alpha=.6)
        ax.axhline(.025, color='black', ls='--', lw=1.)
        ax.set(title=f'({letter}) {mode}', xlabel='Post-search evidence score', ylabel='Calibration upper-tail rank',
               yscale='log', ylim=(.003,1.1))
        ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle('Figure 0004-49-1 | Isolated-event calibration; dashed line: rank 0.025')
    fig.tight_layout(); fig.savefig(folder/'0004-49-1.png', dpi=150); plt.close(fig)
    methods = ['observed_penalty', 'model_penalty', primary, 'known_count_model']
    for number, conditions, title in [('0004-50-1', study.CONDITIONS[:5], 'Original charge combinations'),
                                       ('0004-50-2', study.CONDITIONS[5:], 'Real weak-alpha combinations')]:
        fig, axes = plt.subplots(2, 3, figsize=(15, 9))
        for ax, condition, letter in zip(axes.flat, conditions, 'abcde'):
            for method in methods:
                rr = sorted([r for r in summary if r['method']==method and r['kind']=='pair' and r['condition']==condition],
                            key=lambda r:r['gap_ns'])
                ax.plot([r['gap_ns'] for r in rr], [r['efficiency'] for r in rr],
                        '.--' if method=='known_count_model' else '.-', color=COLORS[method], label=LABELS[method])
            ax.set(title=f'({letter}) {condition}', xlabel='First-active-photon separation (ns)',
                   ylabel='Strict pair recovery, both true Q >= 500 pe', xscale='log', ylim=(-.03,1.04))
            ax.axhline(.9, color='black', ls=':', lw=.8); ax.grid(alpha=.2)
        axes.flat[-1].axis('off')
        handles, labels = axes.flat[0].get_legend_handles_labels()
        axes.flat[-1].legend(handles, labels, loc='center', fontsize=9)
        axes.flat[-1].text(.05,.1,'Known-count curve uses truth K for diagnosis.\nIt is not a blind recovery result.', fontsize=9)
        fig.suptitle(f'Figure {number} | {title}; independent test, 40 generated pairs per point')
        fig.tight_layout(); fig.savefig(folder/(number+'.png'), dpi=150); plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    conditions = ['all','weak','low','mid','high']
    for j, method in enumerate(methods[:3]):
        rr = [next(r for r in summary if r['method']==method and r['kind']!='pair' and r['condition']==c) for c in conditions]
        x = np.arange(5)+(j-1)*.12
        p = np.array([r['extra_fraction'] for r in rr]); ci = np.array([r['extra_wilson95'] for r in rr])
        axes[0].errorbar(x,p,yerr=np.maximum([p-ci[:,0],ci[:,1]-p],0.),fmt='.',capsize=3,color=COLORS[method],label=LABELS[method])
        axes[1].plot(x,[r['efficiency'] for r in rr],'.-',color=COLORS[method])
        q = np.array([r['q_error_p16_p50_p84'] for r in rr[1:]])*100
        axes[2].errorbar(x[1:],q[:,1],yerr=[q[:,1]-q[:,0],q[:,2]-q[:,1]],fmt='.',capsize=3,color=COLORS[method])
    for ax in axes: ax.grid(alpha=.2)
    axes[0].set(title='(a) Extra-output frames, Wilson95',ylabel='Fraction'); axes[0].axhline(.05,color='gray',ls=':')
    axes[1].set(title='(b) Strict whole-frame recovery',ylabel='Fraction',ylim=(-.03,1.04))
    axes[2].set(title='(c) Matched-charge errors, median and p16-p84',ylabel='Relative error (%)')
    axes[2].axhspan(-5,5,color='gray',alpha=.12); axes[2].axhline(0,color='gray',lw=.8)
    for ax in axes[:2]: ax.set_xticks(range(5),['Empty','Weak','Low','Mid','High'])
    axes[2].set_xticks(range(1,5),['Weak','Low','Mid','High'])
    axes[0].legend(fontsize=7)
    fig.suptitle('Figure 0004-50-3 | Independent empty/native-single controls, 80 per group')
    fig.tight_layout(); fig.savefig(folder/'0004-50-3.png',dpi=150); plt.close(fig)
    rows, fits = assemble('test')
    response = ConditionalResponse.load(study.EXPERIMENT/'response.npz')
    with np.load(study.OUTPUT/'test/observations.npz') as f: voltage=f['voltage_mV']
    with np.load(study.OUTPUT/'test/labels.truth.npz') as f: components=f['components_mV']
    protocol = study.load(study.EXPERIMENT/'test-plan.json')
    examples = []
    specifications = protocol['examples_predeclared']
    fig, axes = plt.subplots(4,len(specifications),figsize=(7*len(specifications),13),squeeze=False)
    for col, spec in enumerate(specifications):
        i = next(i for i,r in enumerate(rows) if all(r['case'][k]==v for k,v in spec.items()))
        row=rows[i]; new=fits[primary][i]; old=fits['observed_penalty'][i]
        new_components=response.components(new['times_ns'],new['charge_pe'])
        old_model=response.components(old['times_ns'],old['charge_pe']).sum(0)
        grid=response.sample_grid; y=voltage[i]-response.baseline
        axes[0,col].plot(grid,y,'k',lw=1,label='Observed total')
        axes[0,col].plot(grid,old_model,'--',label='Previous observed-weight fit')
        axes[0,col].plot(grid,new_components.sum(0),label='Model weights + null gate')
        for k in range(2):
            axes[k+1,col].plot(grid,components[i,k],label=f'True component {k+1}')
            if k < len(new_components): axes[k+1,col].plot(grid,new_components[k],'--',label=f'Recovered {k+1}')
        axes[3,col].plot(grid,y-old_model,label='Previous residual')
        axes[3,col].plot(grid,y-new_components.sum(0),label='New residual')
        for level in range(4):
            ax=axes[level,col]
            ax.set(title=f'({chr(97+level*len(specifications)+col)}) {spec["condition"]}, {spec["gap_ns"]:g} ns',
                   xlabel='Full local record time (ns)',ylabel='mV',xlim=(0,2048))
            ax.grid(alpha=.2); ax.legend(fontsize=8)
        examples.append(dict(index=i, case=row['case'], true_q=row['charge'].tolist(), true_t=row['times'].tolist(),
            old=old, new=new, old_score=evaluate_fit(old,row['times'],row['charge'],row['case']),
            new_score=evaluate_fit(new,row['times'],row['charge'],row['case'])))
    fig.suptitle('Figure 0004-51-1 | Predetermined full-record examples; no success-based selection')
    fig.tight_layout(); fig.savefig(folder/'0004-51-1.png',dpi=150); plt.close(fig)
    study.save(study.OUTPUT/'test/examples.json',examples)


if __name__ == '__main__':
    import sys
    if sys.argv[1]=='figures': figures()
    else: report(sys.argv[1])
