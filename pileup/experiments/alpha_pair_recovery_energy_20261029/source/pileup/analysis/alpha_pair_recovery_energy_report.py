"""Post-freeze diagnostics for the recovery-aware energy development study.

No model selection or refitting occurs here. All confidence intervals describe
paired retrospective replay and exclude calibration/model-selection uncertainty.
"""
from __future__ import annotations

from collections import defaultdict
import csv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .alpha_pair_recovery_energy import RecoveryEnergyResponse
from .alpha_pair_recovery_energy_study import (
    PROJECT, EXPERIMENT, OUTPUT, RECORD, GAPS, VARIANTS,
    load, lines, save, write_lines, digest, verify,
)
from .sipm_energy_calibration import time_matches, error_summary

CONDITIONS = ['weak_weak', 'low_weak', 'weak_low', 'high_weak', 'weak_high',
              'low_low', 'mid_mid', 'high_high', 'low_high', 'high_low']
LABELS = ['Original NLL + static map', 'Coupled charge conversion',
          'Isolated E/t refit', 'Recovery-aware E/t refit']
COLORS = ['#777777', '#d7952c', '#3683b5', '#ba4966']
BOOTSTRAP_SEED = 202610296701


def csv_file(name, rows):
    with (OUTPUT/name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def graph_flags(true_t, predicted_t, matches):
    """Original 10 ns proximity-graph definitions, using time-only matches.

    An unmatched truth near an assigned prediction flags a merge. An extra
    prediction near an assigned truth flags a split. Count deficits alone do
    not identify temporal merging. All present events are inside the gate.
    """
    close = np.abs(np.asarray(true_t)[:, None]-np.asarray(predicted_t)[None, :]) <= 10.
    used_t, used_p = {i for i, _ in matches}, {j for _, j in matches}
    merge = {j for i in range(len(true_t)) if i not in used_t for j in used_p if close[i, j]}
    split = {i for i in used_t if any(close[i, j] for j in range(len(predicted_t)) if j not in used_p)}
    return len(split), len(merge)


def validate_means(response):
    cases = load(OUTPUT/'inputs/development.cases.json')
    records = {r['root']: r for r in load(OUTPUT/'inputs/development.energy.json')}
    with np.load(OUTPUT/'inputs/development.truth.npz') as f:
        t, q, curves = f['times_ns'], f['charge_pe'], f['components_mV']
    singles = [i for i, c in enumerate(cases) if c['fold'] == 1 and c['kind'] == 'single']
    assert len(singles) == 280
    rows = []
    max_single_difference = 0.
    physical = dict(min_voltage_mV=0., max_negative_voltage_area_fraction=0.,
                    max_positive_charge_change_relative=0., max_full_area_discrepancy=0.,
                    max_negative_innovation_fraction=0., max_pre_neighbor_difference_mV=0.,
                    max_isolated_q_truncation_relative=0.)
    checked = singles + [i for i, c in enumerate(cases) if c['fold'] == 1 and
                        c['kind'] == 'pair' and c['gap_ns'] in GAPS]
    assert len(checked) == 1180
    for i in checked:
        c = cases[i]; k = len(c['mother_roots'])
        energies = [records[p]['edep_MeV'] for p in c['mother_roots']]
        z, h, q0 = response.drivers(t[i, :k], energies)
        a, qa = response.from_drivers(z, h, 0.)
        b, qb = response.from_drivers(z, h)
        physical['min_voltage_mV'] = min(physical['min_voltage_mV'], float(b.min()))
        physical['max_negative_voltage_area_fraction'] = max(physical['max_negative_voltage_area_fraction'],
            float(np.max(np.maximum(-b, 0).sum(-1)/np.maximum(np.abs(b).sum(-1), 1e-12))))
        physical['max_positive_charge_change_relative'] = max(physical['max_positive_charge_change_relative'], float(np.max(qb/qa-1)))
        physical['max_full_area_discrepancy'] = max(physical['max_full_area_discrepancy'],
            float(np.max(np.abs(b.sum(-1)*response.dt/response.area_per_pe/qb-1))))
        physical['max_negative_innovation_fraction'] = max(physical['max_negative_innovation_fraction'],
            float(np.max(np.maximum(-z, 0).sum(-1)/np.maximum(np.abs(z).sum(-1), 1e-12))))
        physical['max_isolated_q_truncation_relative'] = max(physical['max_isolated_q_truncation_relative'], float(np.max(np.abs(qa/q0-1))))
        if k == 1:
            max_single_difference = max(max_single_difference, float(np.max(np.abs(a-b))))
            rows.append(dict(index=i, condition=c['condition'], true_charge_pe=float(q[i, 0]),
                predicted_charge_pe=float(qb[0]), relative_charge_error=float(qb[0]/q[i, 0]-1),
                relative_voltage_rmse=float(np.linalg.norm(b[0, :1024]-curves[i, 0])/np.linalg.norm(curves[i, 0]))))
        else:
            for j in range(2):
                early = response.full_grid < t[i, 1-j]
                if early.any():
                    physical['max_pre_neighbor_difference_mV'] = max(physical['max_pre_neighbor_difference_mV'], float(np.max(np.abs(a[j, early]-b[j, early]))))
        assert np.isfinite(b).all() and (qb > 0).all()
    assert max_single_difference == 0.
    assert physical['max_pre_neighbor_difference_mV'] < 1e-10
    assert physical['max_negative_voltage_area_fraction'] < 1e-10
    assert physical['max_positive_charge_change_relative'] < 1e-10
    assert physical['max_full_area_discrepancy'] < 1e-5
    csv_file('heldout-single-means.csv', rows)
    summary = []
    for condition in ['weak', 'low', 'mid', 'high']:
        rr = [r for r in rows if r['condition'] == condition]
        summary.append(dict(condition=condition, n=len(rr), relative_charge_rmse=float(np.sqrt(np.mean([r['relative_charge_error']**2 for r in rr]))),
            median_relative_charge_error=float(np.median([r['relative_charge_error'] for r in rr]))))
    result = dict(passed=True, single_frames=280, pair_frames=900, exact_single_invariance=True,
        single_summary=summary, numerical_response_checks=physical,
        interpretation='Innovation proxy is not a microscopic avalanche count; no retuning after heldout validation')
    save(RECORD/'heldout-response-validation.json', result)
    return result


def scoring_tables():
    cases = load(OUTPUT/'inputs/evaluation.cases.json')
    scores = lines(OUTPUT/'scores.jsonl')
    predictions = {(r['method'], r['index']): r for r in lines(OUTPUT/'predictions.jsonl')}
    old_path = PROJECT/'jupyter/artifacts/sipm_energy_guardband_20261029/historical-energy.scores.jsonl'
    old = {r['index']: r for r in lines(old_path) if r['method'] == 'joint_conserved_nll'}
    with np.load(OUTPUT/'inputs/evaluation.labels.npz') as f:
        q, t, e = f['charge_pe'], f['times_ns'], f['energy_MeV']
    components = []
    for row in scores:
        i, method = row['index'], row['method']; p = predictions[method, i]
        mask = (q[i] >= 500.) & (e[i] > 0); tt, ee, qq = t[i, mask], e[i, mask], q[i, mask]
        pairs = [(m['target_index'], m['prediction_index']) for m in row['matches']]
        split, merge = graph_flags(tt, p['times_ns'], pairs)
        row['count_over_number_flag'], row['count_under_number_flag'] = row.pop('split'), row.pop('merge')
        row['temporal_split'], row['temporal_merge'] = split, merge
        tol = min(10., row['gap_ns']/4) if row['gap_ns'] is not None else 10.
        strict_matches = time_matches(tt, p['times_ns'], tol)
        row['failure_stage'] = ('wrong_count' if len(p['times_ns']) != len(tt) else
            'strict_time' if len(strict_matches) != len(tt) else
            'energy' if not row['frame_success_10'] else 'success')
        assert row['output_count'] == row['target']-row['missed']+row['extra']
        assert len(p['times_ns']) == len(predictions['static_nll', i]['times_ns'])
        if method == 'static_nll':
            for key in ['target', 'missed', 'extra', 'strict_10', 'frame_success_10', 'output_count']:
                assert row[key] == old[i][key], (i, key)
            assert (split, merge) == (old[i]['split'], old[i]['merge']), i
        if row['kind'] == 'pair':
            for m in row['matches']:
                a, b = m['target_index'], m['prediction_index']
                components.append(dict(index=i, method=method, condition=row['condition'], gap_ns=row['gap_ns'],
                    component=a, species=row['condition'].split('_')[a],
                    relative_energy_error=m['relative_energy_error'],
                    relative_charge_error=float(p['charge_pe'][b]/qq[a]-1), time_error_ns=m['time_error_ns']))
    write_lines(OUTPUT/'audited-scores.jsonl', scores)
    csv_file('matched-components.csv', components)
    grouped = defaultdict(list)
    for r in scores: grouped[r['method'], r['kind'], r['condition'], r['gap_ns']].append(r)
    summaries = []
    fields = ['target', 'missed', 'extra', 'output_count', 'strict_5', 'strict_10', 'strict_20',
              'frame_success_10', 'temporal_split', 'temporal_merge', 'count_over_number_flag',
              'count_under_number_flag', 'energy_out_of_domain']
    for (method, kind, condition, gap), rr in grouped.items():
        matches = [m for r in rr for m in r['matches'] if m['relative_energy_error'] is not None]
        errs = error_summary([m['true_energy_MeV'] for m in matches], [m['estimated_energy_MeV'] for m in matches]) if matches else {}
        summaries.append(dict(method=method, kind=kind, condition=condition, gap_ns=gap, frames=len(rr),
            **{key: sum(r[key] for r in rr) for key in fields},
            median_relative_energy=errs.get('median_relative'), sigma68_relative_energy=errs.get('sigma68_relative'),
            finite_time_matches=len(matches), boundary=sum(predictions[method, r['index']].get('boundary', False) for r in rr),
            not_converged=sum(not predictions[method, r['index']].get('converged', True) for r in rr)))
    csv_file('energy-accuracy.csv', summaries)
    failures, totals = [], []
    for method in VARIANTS:
        for group in ['weak', 'strong']:
            rr = [r for r in scores if r['method'] == method and r['kind'] == 'pair' and ('weak' in r['condition']) == (group == 'weak')]
            failures.append(dict(method=method, group=group, frames=len(rr), **{key: sum(r['failure_stage'] == key for r in rr)
                for key in ['wrong_count', 'strict_time', 'energy', 'success']}))
        for kind in ['pair', 'single', 'empty']:
            rr = [r for r in scores if r['method'] == method and r['kind'] == kind]
            totals.append(dict(method=method, kind=kind, frames=len(rr), **{key: sum(r[key] for r in rr) for key in fields}))
    csv_file('failure-stages.csv', failures); csv_file('event-counts.csv', totals)
    weights = np.zeros((2000, 40), int); rng = np.random.default_rng(BOOTSTRAP_SEED)
    for start in [0, 20]:
        draws = rng.integers(start, start+20, size=(2000, 20))
        for k in range(2000): weights[k] += np.bincount(draws[k], minlength=40)
    bootstrap, means = {}, {}
    for method in VARIANTS:
        for group in ['weak', 'strong']:
            vv, bb = [], []
            for (m, kind, condition, gap), rr in grouped.items():
                if m != method or kind != 'pair' or ('weak' in condition) != (group == 'weak'): continue
                ids = np.array([r['block'] for r in rr]); values = np.array([r['frame_success_10'] for r in rr])
                num = np.bincount(ids, weights=values, minlength=40); den = np.bincount(ids, minlength=40)
                vv.append(values.mean()); bb.append((weights@num)/(weights@den))
            means[method, group] = float(np.mean(vv)); bootstrap[method, group] = np.mean(bb, axis=0)
    contrasts = []
    for method, baseline in [(m, 'static_nll') for m in VARIANTS[1:]] + [('recovery_energy_refit', 'isolated_energy_refit')]:
        for group in ['weak', 'strong']:
            diff = bootstrap[method, group]-bootstrap[baseline, group]
            ci = np.quantile(diff, [.025, .975])
            contrasts.append(dict(method=method, baseline=baseline, group=group,
                difference=means[method, group]-means[baseline, group], low95=float(ci[0]), high95=float(ci[1])))
    csv_file('paired-comparisons.csv', contrasts)
    return scores, predictions, summaries, failures, totals, contrasts


def figures(response, summaries, predictions):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'figure.dpi': 110, 'savefig.dpi': 150, 'svg.fonttype': 'none'})
    saved = []
    def finish(fig, name):
        for suffix in ['.png', '.svg']: fig.savefig(OUTPUT/(name+suffix), bbox_inches='tight')
        plt.close(fig); saved.append(name)
    cases = load(OUTPUT/'inputs/evaluation.cases.json')
    with np.load(OUTPUT/'inputs/evaluation.labels.npz') as f: e, t, q = f['energy_MeV'], f['times_ns'], f['charge_pe']
    truth_path = PROJECT/'jupyter/artifacts/alpha_pair_joint_confirmation_20261026/labels.truth.npz'
    with np.load(truth_path) as f: components = f['components_mV']
    with np.load(OUTPUT/'inputs/observations.npz') as f: voltage = f['voltage_mV']
    # Fixed first case in the protocol pilot strata; no selection by outcome.
    example = next(i for i, c in enumerate(cases) if c['condition'] == 'weak_high' and c['gap_ns'] == 40.)
    z, h, _ = response.drivers(t[example], e[example])
    iso, _ = response.from_drivers(z, h, 0.); recovered, _ = response.from_drivers(z, h)
    fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True, layout='constrained')
    x = response.sample_grid
    ax[0].plot(x, h[0, :1024], color=COLORS[2]); ax[0].axvline(t[example, 1], color='black', linestyle=':')
    ax[0].set(ylabel='Neighbor activity proxy', title='(a) Causal history acting on the weak component')
    for yy, label, color in [(components[example, 0], 'Simulated shared weak component', '#bbb'),
            (iso[0, :1024], 'Isolated mean at true E/t', COLORS[2]),
            (recovered[0, :1024], 'Recovery mean at true E/t', COLORS[3])]:
        ax[1].plot(x, yy, label=label, color=color)
    ax[1].set(xlabel='Time from window start / ns', ylabel='Weak component / mV', xlim=(0, 2048),
              title='(b) Known-energy diagnostic; this is not a blind fit')
    ax[1].legend(fontsize=9)
    fig.suptitle('Fig. 0007-21-1 | Causal energy response, weak -> high, 40 ns\nFull observed window; first protocol-pilot example')
    finish(fig, '0007-21-1_causal_response')

    means = lines(OUTPUT/'mean-charge-diagnostics.jsonl')
    fig, axes = plt.subplots(1, 2, figsize=(11, 6.8), layout='constrained')
    for panel, (ax, method) in enumerate(zip(axes, ['isolated', 'recovery'])):
        values = np.array([[100*np.sqrt(np.mean([r['relative_charge_error']**2 for r in means
            if r['fold'] == 1 and r['method'] == method and r['condition'] == c and r['gap_ns'] == g]))
            for g in GAPS] for c in CONDITIONS])
        im = ax.imshow(values, vmin=0, vmax=15, cmap='YlOrRd', aspect='auto')
        for i in range(10):
            for j in range(3): ax.text(j, i, f'{values[i,j]:.1f}', ha='center', va='center', color='white' if values[i,j]>9 else 'black')
        ax.set_xticks(range(3), ['8', '40', '200']); ax.set_yticks(range(10), [c.replace('_', ' -> ') for c in CONDITIONS])
        ax.set(xlabel='Separation / ns', title=f'({chr(97+panel)}) {method.title()} response')
    fig.colorbar(im, ax=axes, label='Component charge relative RMSE / %', shrink=.85)
    fig.suptitle('Fig. 0007-22-1 | Held-out development mean-response check\n900 pairs; both components included; known true E/t')
    finish(fig, '0007-22-1_heldout_charge')

    fig, axes = plt.subplots(2, 2, figsize=(12.8, 12), layout='constrained')
    for panel, (ax, method, label) in enumerate(zip(axes.flat, VARIANTS, LABELS)):
        lookup = {(r['condition'], r['gap_ns']): r for r in summaries if r['method'] == method and r['kind'] == 'pair'}
        counts = np.array([[lookup[c, g]['frame_success_10'] for g in GAPS] for c in CONDITIONS])
        im = ax.imshow(counts/40*100, vmin=0, vmax=100, cmap='YlGnBu', aspect='auto')
        for i in range(10):
            for j in range(3): ax.text(j, i, f'{counts[i,j]}/40', ha='center', va='center', color='white' if counts[i,j]>23 else 'black')
        ax.set_xticks(range(3), ['8', '40', '200']); ax.set_yticks(range(10), [c.replace('_', ' -> ') for c in CONDITIONS])
        ax.set(xlabel='Separation / ns', title=f'({chr(97+panel)}) {label}')
    fig.colorbar(im, ax=axes, label='Strict full-pair success / %', shrink=.75)
    fig.suptitle('Fig. 0007-23-1 | Retrospective energy reconstruction; original blind K retained\nBoth E within 10%, time within min(10 ns, gap/4), no extra output')
    finish(fig, '0007-23-1_pair_accuracy')

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.7), layout='constrained')
    groups = ['weak', 'low', 'mid', 'high']; xx = np.arange(4)
    for j, (method, label, color) in enumerate(zip(VARIANTS, LABELS, COLORS)):
        rr = [next(r for r in summaries if r['method'] == method and r['kind'] == 'single' and r['condition'] == c) for c in groups]
        pos = xx+(j-1.5)*.17
        axes[0].errorbar(pos, [r['median_relative_energy']*100 for r in rr],
            yerr=[r['sigma68_relative_energy']*100 for r in rr], fmt='o', color=color, capsize=3, label=label)
        axes[1].bar(pos, [r['strict_10']/80*100 for r in rr], width=.16, color=color)
    axes[0].axhline(0, color='#aaa', linestyle=':'); axes[0].set(ylabel='Relative deposited-energy error / %', title='(a) Median +/- sigma68 (not a confidence interval)')
    axes[1].set(ylabel='Strict event energy recovery / %', ylim=(0, 105), title='(b) All 80 true events per group in denominator')
    for ax in axes: ax.set_xticks(xx, groups)
    axes[0].legend(fontsize=8, loc='lower right')
    fig.suptitle('Fig. 0007-24-1 | Historical single-event controls\nSame count selection; 80 empty controls separately retain zero output')
    finish(fig, '0007-24-1_controls')

    examples = [next(i for i, c in enumerate(cases) if c['condition'] == 'high_high' and c['gap_ns'] == 8.), example]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True, layout='constrained')
    for col, i in enumerate(examples):
        axes[0, col].plot(x, voltage[i]-response.baseline-response.background, color='#aaa', lw=.8, label='Observed voltage')
        chosen_component = 0
        axes[1, col].plot(x, components[i, chosen_component], color='#aaa', lw=.9, label='True first component (evaluation only)')
        for method, color in [('static_nll', COLORS[0]), ('isolated_energy_refit', COLORS[2]), ('recovery_energy_refit', COLORS[3])]:
            row = predictions[method, i]
            yy = (response.legacy_components(row['times_ns'], row['charge_pe']) if method == 'static_nll' else
                  response.components(row['times_ns'], row['energy_MeV'], 0. if method == 'isolated_energy_refit' else None)[0])
            label = LABELS[VARIANTS.index(method)]
            axes[0, col].plot(x, yy.sum(0), color=color, label=label, lw=1.)
            pair = dict(time_matches(t[i], row['times_ns']))
            if chosen_component in pair: axes[1, col].plot(x, yy[pair[chosen_component]], color=color, label=label, lw=1.)
        axes[0, col].set(title=f'({chr(97+col)}) {cases[i]["condition"].replace("_", " -> ")}, {cases[i]["gap_ns"]:g} ns; frame {i}', ylabel='Total signal / mV')
        axes[1, col].set(title=f'({chr(99+col)}) First component, ordinary 10 ns time match', xlabel='Time from window start / ns', ylabel='Component / mV', xlim=(0, 2048))
        axes[0, col].legend(fontsize=8); axes[1, col].legend(fontsize=8)
    fig.suptitle('Fig. 0007-25-1 | Full-window fits and component decomposition\nFirst predeclared pilot-stratum examples; true components never enter inference')
    finish(fig, '0007-25-1_waveforms')
    save(RECORD/'plot-inputs.json', dict(truth_path=str(truth_path.relative_to(PROJECT)), sha256=digest(truth_path),
        waveform_examples=examples, causal_example=example, selected_by='First case in protocol pilot strata, not score'))
    return saved


def main():
    verify(); verify('method-freeze.json'); verify('report-source-freeze.json')
    response = RecoveryEnergyResponse.load(EXPERIMENT)
    heldout = validate_means(response)
    scores, predictions, summary, failures, totals, contrasts = scoring_tables()
    images = figures(response, summary, predictions)
    report = dict(load(OUTPUT/'report.json'), audited_summary=summary, failures=failures, event_totals=totals,
        paired_comparisons=contrasts, heldout_response=heldout, figures=images,
        bootstrap=dict(seed=BOOTSTRAP_SEED, repetitions=2000, blocks=40, strata=[20, 20],
            scope='Post-result descriptive paired intervals, calibration and model selection excluded'),
        score_note='Raw scores split/merge are count-only flags. Audited scores use original temporal proximity graph; static baseline exactly reproduced.',
        default_changed=False, fresh_confirmation=False, new_ROOTs=0)
    save(OUTPUT/'audited-report.json', report)
    save(RECORD/'report-validation.json', dict(passed=True, frames=1600, methods=4, audited_score_rows=len(scores),
        original_baseline_reproduced=True, original_blind_count_retained=True, figures=len(images), heldout_singles=280))
    print('RECOVERY_REPORT', contrasts, totals, heldout, flush=True)


if __name__ == '__main__': main()
