"""Paired diagnostic tables and numbered figures; no model retuning."""
from __future__ import annotations

from collections import Counter, defaultdict
import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .alpha_pair_energy_oracle_study import (
    PROJECT, EXPERIMENT, OUTPUT, RECORD, METHODS, GAPS, load, lines, save, digest, verify,
)
from .alpha_pair_recovery_energy import RecoveryEnergyResponse
from .sipm_energy_calibration import time_matches

CONDITIONS = ['weak_weak', 'low_weak', 'weak_low', 'high_weak', 'weak_high',
              'low_low', 'mid_mid', 'high_high', 'low_high', 'high_low']
TITLES = ['Blind count + free E/t', 'Known count, blind time starts',
          'Known count + fixed true times', 'Known count + true time starts, then free E/t']
COLORS = ['#555555', '#dc9b2d', '#397cb3', '#b34f76']


def csv_file(name, rows):
    with (OUTPUT/name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def block_design(cases):
    """Connected optical-parent components, stratified by design composition.

    Shared ROOTs can cross original dataset names. Connecting on ROOT identity
    prevents treating repeated digitizations of one mother as independent.
    """
    ids = [i for i, c in enumerate(cases) if c['fold'] == 1 and c['mother_roots']]
    parent = {root: root for i in ids for root in cases[i]['mother_roots']}
    def find(root):
        while parent[root] != root:
            parent[root] = parent[parent[root]]; root = parent[root]
        return root
    for i in ids:
        roots = cases[i]['mother_roots']
        for other in roots[1:]: parent[find(other)] = find(roots[0])
    root_groups = defaultdict(list)
    for root in parent: root_groups[find(root)].append(root)
    groups = sorted([sorted(group) for group in root_groups.values()])
    lookup = {root: j for j, group in enumerate(groups) for root in group}
    frame_blocks = {i: lookup[cases[i]['mother_roots'][0]] for i in ids}
    composition = defaultdict(Counter)
    for i, j in frame_blocks.items():
        c = cases[i]; composition[j][(c['dataset'], c['kind'], c['condition'], c['gap_ns'])] += 1
        assert len({lookup[root] for root in c['mother_roots']}) == 1
    strata = defaultdict(list)
    for j in range(len(groups)):
        signature = json.dumps(sorted([(list(key), count) for key, count in composition[j].items()], key=str))
        strata[signature].append(j)
    result = dict(blocks=groups, frame_blocks={str(i): j for i, j in frame_blocks.items()},
                  strata=list(strata.values()), count=len(groups), rule='ROOT-connected components; same design composition within each bootstrap stratum')
    save(RECORD/'mother-block-design.json', result)
    return result


def tables(report, scores, predictions):
    cases = load(OUTPUT/'inputs/cases.json')
    lookup = {(r['method'], r['index']): r for r in scores}
    rows = []
    for r in report['summary']:
        row = {k: v for k, v in r.items() if k != 'matched_energy'}
        error = r['matched_energy'] or {}
        row.update(matched_n=error.get('n', 0), median_energy_relative=error.get('median_relative'),
                   sigma68_energy_relative=error.get('sigma68_relative'), rmse_MeV=error.get('rmse_MeV'))
        rows.append(row)
    csv_file('energy-accuracy.csv', rows); csv_file('aggregate.csv', report['aggregate'])
    transitions = []
    for group in ['weak', 'strong']:
        ids = [i for i, c in enumerate(cases) if c['fold'] == 1 and c['kind'] == 'pair' and ('weak' in c['condition']) == (group == 'weak')]
        for after, before in [('true_count_blind_start', 'blind_energy_gate'),
                              ('fixed_true_time', 'true_count_blind_start'),
                              ('free_from_true_time', 'fixed_true_time'),
                              ('free_from_true_time', 'true_count_blind_start')]:
            a = np.array([lookup[after, i]['frame_success_10'] for i in ids], bool)
            b = np.array([lookup[before, i]['frame_success_10'] for i in ids], bool)
            transitions.append(dict(group=group, after=after, before=before, frames=len(ids),
                both_success=int(np.sum(a & b)), gained=int(np.sum(a & ~b)), lost=int(np.sum(~a & b)),
                both_failed=int(np.sum(~a & ~b)), difference=float(a.mean()-b.mean())))
    csv_file('paired-transitions.csv', transitions)
    components = []
    for s in scores:
        if s['fold'] != 1 or s['kind'] != 'pair': continue
        for m in s['matches']:
            components.append(dict(index=s['index'], method=s['method'], condition=s['condition'], gap_ns=s['gap_ns'],
                component=m['target_index'], species=s['condition'].split('_')[m['target_index']],
                relative_energy_error=m['relative_energy_error'], time_error_ns=m['time_error_ns']))
    csv_file('matched-components.csv', components)
    weak = []
    for method in METHODS:
        for condition in CONDITIONS[:5]:
            for gap in GAPS:
                rr = [r for r in components if r['method'] == method and r['condition'] == condition and r['gap_ns'] == gap and r['species'] == 'weak']
                values = np.array([r['relative_energy_error'] for r in rr], float)
                ids = [i for i, c in enumerate(cases) if c['fold'] == 1 and c['kind'] == 'pair' and c['condition'] == condition and c['gap_ns'] == gap]
                target = len(ids)*(2 if condition == 'weak_weak' else 1)
                weak.append(dict(method=method, condition=condition, gap_ns=gap, true_weak_components=target,
                    matched=len(rr), unmatched=target-len(rr), within_energy10=int(np.sum(np.abs(values) <= .1)),
                    mean_relative=float(np.mean(values)) if len(values) else None,
                    median_relative=float(np.median(values)) if len(values) else None,
                    sigma68_relative=float(np.diff(np.quantile(values, [.16, .84]))[0]/2) if len(values) else None))
    csv_file('weak-energy.csv', weak)
    search = lines(OUTPUT/'search-diagnostics.jsonl'); search_summary = []
    for condition in CONDITIONS:
        for gap in GAPS:
            rr = [r for r in search if r['fold'] == 1 and r['condition'] == condition and r['gap_ns'] == gap]
            values = np.array([r['blind_search_minus_true_start_cost'] for r in rr])
            search_summary.append(dict(condition=condition, gap_ns=gap, frames=len(rr),
                lower_cost_from_true_start=int(np.sum(values > 1.)), median_blind_minus_true_start_cost=float(np.median(values)),
                fixed_to_free_cost_reduction_median=float(np.median([r['fixed_minus_released_cost'] for r in rr])),
                true_start_success_count=sum(lookup['free_from_true_time', r['index']]['frame_success_10'] for r in rr),
                known_count_success_count=sum(lookup['true_count_blind_start', r['index']]['frame_success_10'] for r in rr)))
    csv_file('search-gaps.csv', search_summary)
    blocks = block_design(cases); nblocks = blocks['count']
    rng = np.random.default_rng(202610297301); weights = np.zeros((2000, nblocks), int)
    for group in blocks['strata']:
        draws = rng.choice(group, (2000, len(group)), replace=True)
        for i, draw in enumerate(draws): weights[i] += np.bincount(draw, minlength=nblocks)
    boot, points = {}, {}
    frame_block = {int(i): j for i, j in blocks['frame_blocks'].items()}
    for method in METHODS:
        for weak_group, label in [(True, 'weak'), (False, 'strong')]:
            values, replicates = [], []
            for condition in CONDITIONS:
                if ('weak' in condition) != weak_group: continue
                for gap in GAPS:
                    rr = [s for s in scores if s['fold'] == 1 and s['method'] == method and s['kind'] == 'pair' and s['condition'] == condition and s['gap_ns'] == gap]
                    ids = np.array([frame_block[r['index']] for r in rr]); success = np.array([r['frame_success_10'] for r in rr], float)
                    den = np.bincount(ids, minlength=nblocks); num = np.bincount(ids, weights=success, minlength=nblocks)
                    assert np.all(weights@den > 0)
                    replicates.append((weights@num)/(weights@den)); values.append(success.mean())
            boot[method, label] = np.mean(replicates, axis=0); points[method, label] = float(np.mean(values))
    intervals = []
    for transition in transitions:
        a, b, group = transition['after'], transition['before'], transition['group']
        low, high = np.quantile(boot[a, group]-boot[b, group], [.025, .975])
        intervals.append(dict(transition, low95=float(low), high95=float(high)))
    csv_file('paired-intervals.csv', intervals)
    return rows, weak, search_summary, intervals, blocks


def figures(report, metrics, weak, intervals, predictions):
    plt.rcParams.update({'font.size': 10, 'figure.dpi': 110, 'savefig.dpi': 150,
                        'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
    names = []
    def finish(fig, name):
        for suffix in ['.png', '.svg']: fig.savefig(OUTPUT/(name+suffix), bbox_inches='tight')
        plt.close(fig); names.append(name)
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 12), layout='constrained')
    for panel, (ax, method, title) in enumerate(zip(axes.flat, METHODS, TITLES)):
        lookup = {(r['condition'], r['gap_ns']): r for r in metrics if r['fold'] == 1 and r['method'] == method and r['kind'] == 'pair'}
        values = np.array([[lookup[c, g]['frame_success_10']/lookup[c, g]['frames']*100 for g in GAPS] for c in CONDITIONS])
        im = ax.imshow(values, vmin=0, vmax=100, cmap='YlGnBu', aspect='auto')
        for i, c in enumerate(CONDITIONS):
            for j, g in enumerate(GAPS):
                row = lookup[c, g]; ax.text(j, i, f"{row['frame_success_10']}/{row['frames']}", ha='center', va='center', color='white' if values[i, j] > 60 else 'black')
        ax.set_xticks(range(3), ['8', '40', '200']); ax.set_yticks(range(10), [c.replace('_', ' -> ') for c in CONDITIONS])
        ax.set(xlabel='Separation / ns', title=f'({chr(97+panel)}) {title}')
    fig.colorbar(im, ax=axes, label='Strict full-pair energy success / %', shrink=.8)
    fig.suptitle('Fig. 0007-29-1 | Same 900 development-validation pairs and common weighting\nE within 10%, time within min(10 ns, gap/4); panels (b-d) use truth assistance')
    finish(fig, '0007-29-1_oracle_comparison')

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.8), layout='constrained')
    short = ['Blind', 'Known K', 'Known K,t', 'True-t start, free E/t']
    for group_index, (group, label) in enumerate([('weak', 'With weak alpha'), ('strong', 'Without weak alpha')]):
        values = [next(r['strict_pair_success'] for r in report['aggregate'] if r['method'] == m and r['fold'] == 1 and r['group'] == group) for m in METHODS]
        axes[0].plot(np.arange(4), np.asarray(values)*100, 'o-', label=label)
        for i, v in enumerate(values): axes[0].text(i+.02, 100*v+(-2 if group_index else 2), f'{v*100:.1f}', fontsize=9)
    axes[0].set_xticks(range(4), short, rotation=20, ha='right'); axes[0].set(ylabel='Strict full-pair success / %', ylim=(0, 105), title='(a) Truth assistance changes the task')
    axes[0].legend(fontsize=9)
    pairs = [('true_count_blind_start', 'blind_energy_gate'), ('fixed_true_time', 'true_count_blind_start'),
             ('free_from_true_time', 'fixed_true_time'), ('free_from_true_time', 'true_count_blind_start')]
    for i, pair in enumerate(pairs):
        for offset, group, color in [(-.12, 'weak', '#3b78a3'), (.12, 'strong', '#b4772d')]:
            r = next(r for r in intervals if (r['after'], r['before']) == pair and r['group'] == group)
            x = r['difference']*100
            axes[1].errorbar(x, i+offset, xerr=[[x-r['low95']*100], [r['high95']*100-x]], fmt='o', color=color, capsize=3)
    axes[1].set_yticks(range(4), ['Known K - blind', 'Fixed true t - known K', 'Release t - fixed true t', 'True-t start - blind starts'])
    axes[1].axvline(0, color='#aaa', linestyle=':'); axes[1].invert_yaxis()
    axes[1].set(xlabel='Paired change / percentage points', title='(b) Mother-component bootstrap 95% intervals')
    fig.suptitle('Fig. 0007-30-1 | Count, timing and search contrasts\nBlue: with weak alpha; brown: without weak alpha; conditional development intervals')
    finish(fig, '0007-30-1_paired_diagnostics')

    fig, axes = plt.subplots(2, 1, figsize=(13, 8.8), sharex=True, layout='constrained')
    keys = [(c, g) for c in CONDITIONS[:5] for g in GAPS]
    for j, method in enumerate(METHODS):
        rr = [next(r for r in weak if r['method'] == method and (r['condition'], r['gap_ns']) == key) for key in keys]
        x = np.arange(len(keys))+(j-1.5)*.16
        y = np.array([np.nan if r['median_relative'] is None else r['median_relative']*100 for r in rr])
        spread = np.array([np.nan if r['sigma68_relative'] is None else r['sigma68_relative']*100 for r in rr])
        axes[0].errorbar(x, y, yerr=spread, fmt='o', color=COLORS[j], capsize=2, ms=3, label=short[j])
        axes[1].bar(x, [r['within_energy10']/r['true_weak_components']*100 for r in rr], width=.15, color=COLORS[j])
    axes[0].axhline(0, color='#aaa', linestyle=':'); axes[0].set(ylabel='Relative weak-energy error / %', title='(a) Time-matched median +/- sigma68; not a confidence interval')
    axes[0].legend(fontsize=9, ncol=2)
    axes[1].set_xticks(range(len(keys)), [f'{c.replace("_", " -> ")}\n{g:g} ns' for c, g in keys], rotation=40, ha='right')
    axes[1].set(ylabel='Weak E within 10% / all true weak events / %', ylim=(0, 105), title='(b) Ordinary 10 ns time match; this is not the strict pair endpoint')
    fig.suptitle('Fig. 0007-31-1 | Weak component energy and retained-event denominators\nNo signal removal or normalization to recovered counts')
    finish(fig, '0007-31-1_weak_energy')

    cases = load(OUTPUT/'inputs/cases.json'); response = RecoveryEnergyResponse.load(EXPERIMENT)
    with np.load(OUTPUT/'inputs/observations.npz') as f: voltage = f['voltage_mV']
    with np.load(OUTPUT/'inputs/evaluation.labels.npz') as f: components, times = f['components_mV'], f['times_ns']
    examples = [next(i for i, c in enumerate(cases) if c['fold'] == 1 and c['condition'] == condition and c['gap_ns'] == 40.) for condition in ['high_weak', 'weak_high']]
    fig, axes = plt.subplots(3, 2, figsize=(13.2, 11), sharex=True, layout='constrained'); x = response.sample_grid
    for col, i in enumerate(examples):
        weak_position = cases[i]['condition'].split('_').index('weak'); strong_position = 1-weak_position
        axes[0, col].plot(x, voltage[i]-response.baseline-response.background, color='#bbb', lw=.7, label='Observed')
        for row, position in [(1, weak_position), (2, strong_position)]: axes[row, col].plot(x, components[i, position], color='#aaa', lw=.9, label='True component')
        for method, color, label in zip(METHODS, COLORS, short):
            p = predictions[method, i]; curves, _ = response.components(p['times_ns'], p['energy_MeV'])
            axes[0, col].plot(x, curves.sum(0), color=color, lw=1., label=label)
            matched = dict(time_matches(times[i], p['times_ns']))
            for row, position in [(1, weak_position), (2, strong_position)]:
                if position in matched: axes[row, col].plot(x, curves[matched[position]], color=color, lw=1., label=label)
        for row, label in enumerate(['Total', 'Weak component', 'Strong component']):
            axes[row, col].set(ylabel=f'{label} / mV', xlim=(0, 2048), title=f'({chr(97+row*2+col)}) {label}')
            axes[row, col].legend(fontsize=7, loc='upper right')
        axes[0, col].set_title(f'({chr(97+col)}) {cases[i]["condition"].replace("_", " -> ")}, 40 ns; frame {i}')
        axes[2, col].set_xlabel('Time from observed-window start / ns')
    fig.suptitle('Fig. 0007-32-1 | Full-window total and separated-component diagnostics\nFirst held-out high/weak and weak/high 40 ns cases; true components used only for plotting')
    finish(fig, '0007-32-1_full_waveforms')
    save(RECORD/'plot-provenance.json', dict(indices=examples, rule='First fold1 case in each declared high_weak/weak_high40ns stratum, not selected by outcome'))
    return names


def main():
    verify(); verify('execution-freeze.json'); verify('selection-freeze.json'); verify('report-freeze.json')
    report = load(OUTPUT/'report.json'); scores = lines(OUTPUT/'scores.jsonl')
    predictions = {(p['method'], p['index']): p for p in lines(OUTPUT/'predictions.jsonl')}
    assert len(scores) == len(predictions) == 6160
    metrics, weak, search, intervals, blocks = tables(report, scores, predictions)
    names = figures(report, metrics, weak, intervals, predictions)
    result = dict(report, weak_components=weak, search_summary=search, paired_intervals=intervals, figures=names,
        bootstrap=dict(seed=202610297301, repetitions=2000, blocks=blocks['count'], stratum_sizes=list(map(len, blocks['strata'])),
            note='Optical ROOT-connected components, repeated digitizations linked; fixed model/gate conditional intervals, no selection or parameter uncertainty'),
        interpretation='Only the blind_energy_gate row is blind. Other rows use declared true count/time information. Current development weights differ from historical replay.')
    save(OUTPUT/'audited-report.json', result)
    save(RECORD/'report-validation.json', dict(passed=True, score_rows=6160, frames=1540, validation_pairs=900,
        validation_mothers=280, figures=len(names), parent_blocks=blocks['count'], parameter_retuning=False))
    print('ORACLE_REPORT_COMPLETE', report['aggregate'], intervals, 'blocks', blocks['count'], flush=True)


if __name__ == '__main__': main()
