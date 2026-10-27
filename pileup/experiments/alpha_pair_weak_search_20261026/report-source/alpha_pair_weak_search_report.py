"""Checks, paired summaries and figures for the frozen weak-search study."""
from __future__ import annotations

import argparse
import csv
import json

import numpy as np

from .alpha_pair_weak_search_study import PROJECT, EXPERIMENT, OUTPUT, RECORD, OLDER, LATEST, VARIANTS, rows, check_freeze
from .alpha_pair_resolution_data import load, save, digest

LABELS = dict(baseline='Original joint NLL', legacy_extra='Extra original-share starts', weak_extra='Extra weak-share starts')
COLORS = dict(baseline='#555555', legacy_extra='#0077bb', weak_extra='#009988')
CONDITIONS = ['low_low', 'mid_mid', 'high_high', 'low_high', 'high_low',
              'low_weak', 'weak_low', 'high_weak', 'weak_high', 'weak_weak']


def validate():
    from .alpha_pair_joint_likelihood import ConservingResponse, JointNoise, calibrate_gate, select
    from .alpha_pair_profile_diagnostics import cost_at
    from .alpha_pair_conditional_study import evaluate_fit
    from .alpha_pair_null_report import summarize
    check_freeze(); provenance = load(RECORD / 'input-provenance.json')
    for path, expected in provenance['files'].items(): assert digest(PROJECT / path) == expected, path
    for filename, expected in provenance['output_hashes'].items(): assert digest(OUTPUT / filename) == expected, filename
    groups = [set(v) for v in provenance['mother_groups'].values()]
    assert all(not groups[i] & groups[j] for i in range(len(groups)) for j in range(i))
    assert load(RECORD / 'numerical-tests.json')['passed']
    for p, sha in load(RECORD / 'numerical-tests.json')['files'].items(): assert digest(PROJECT / p) == sha
    cases = load(OUTPUT / 'labels.cases.json'); protocol = load(EXPERIMENT / 'protocol.json')
    with np.load(OUTPUT / 'labels.truth.npz') as f: q = f['charge_pe']; t = f['times_ns']
    with np.load(OUTPUT / 'observations.npz') as f: voltage = f['voltage_mV']
    response = ConservingResponse.load(EXPERIMENT / 'joint-response.npz')
    noise = JointNoise.load(EXPERIMENT / 'joint_conserved_nll-noise.npz')
    predictions = {name: {r['index']: r for r in rows(OUTPUT / f'predictions.{name}.jsonl')} for name in LABELS}
    assert all(len(p) == 3480 for p in predictions.values())
    gates = load(EXPERIMENT / 'gate-calibrations.json'); scores = rows(OUTPUT / 'scores.jsonl')
    saved = {(r['method'], r['index']): r for r in scores}; assert len(saved) == len(scores) == 9600
    report = load(OUTPUT / 'report.json'); summary = []; convergence = {}; runtimes = {}; cost_replay_error = 0.
    for name, prediction in predictions.items():
        assert gates[name] == calibrate_gate([prediction[i]['candidates'] for i in range(280)])
        if name != 'baseline':
            run = load(RECORD / (name + '-inference.json'))
            assert run['worker_labels_blocked'] and run['truth_inputs'] == []
            assert run['prediction_sha256'] == digest(OUTPUT / f'predictions.{name}.jsonl')
            assert run['waveform_sha256'] == digest(OUTPUT / 'observations.npz')
            samples = np.array([r['seconds'] for r in prediction.values()])
            runtimes[name] = dict(incremental_seconds_quantiles=np.quantile(samples, [.16, .5, .84, 1.]).tolist(),
                                 sum_worker_elapsed_seconds=float(samples.sum()),
                                 extra_evaluations=sum(r['search']['evaluations'] for r in prediction.values()))
        for i, r in prediction.items():
            candidates = r['candidates']; old = predictions['baseline'][i]['candidates']
            assert candidates[:2] == old[:2]
            if len(candidates) == 3:
                assert candidates[2]['cost'] <= old[2]['cost'] + 1e-10
                replay = cost_at(response, noise, voltage[i], candidates[2]['times_ns'], candidates[2]['charge_pe'])['cost']
                error = abs(replay - candidates[2]['cost']) / max(1., abs(replay))
                cost_replay_error = max(cost_replay_error, error); assert error < 1e-7
        for dataset in ('older', 'latest'):
            grouped = {}; failures = 0
            for i, c in enumerate(cases):
                if c['dataset_role'] != dataset: continue
                fit = select(prediction[i]['candidates'], gates[name]); score = evaluate_fit(fit, t[i], q[i], c)
                stored = saved[name, i]
                assert fit == stored['fit'] and all(stored[k] == v for k, v in score.items())
                failures += not fit.get('converged', True)
                grouped.setdefault((c['kind'], c['condition'], c['gap_ns']), []).append(stored)
            summary += [dict(summarize(rr), dataset_role=dataset) for rr in grouped.values()]
            convergence[name + '/' + dataset] = failures
    assert summary == report['summary'] and len(summary) == 210
    # Baseline replay must reproduce both prior public reports exactly in key metrics.
    prior = {'older': load(OLDER / 'development-report.json')['results']['joint_conserved_nll'],
             'latest': load(LATEST / 'report.json')['results']['joint_conserved_nll']}
    for dataset, old in prior.items():
        for key in ('weak_mean', 'strong_mean', 'controls_pass'):
            assert report['results']['baseline'][dataset][key] == old[key]
    with (OUTPUT / 'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0])); writer.writeheader(); writer.writerows(summary)
    save(OUTPUT / 'runtime-and-convergence.json', dict(runtimes=runtimes, selected_not_converged=convergence,
        interpretation='All selected non-converged results remain scored; elapsed time is extra to the cached original search, includes consistency checks, and mixes different concurrent execution phases; not a controlled speed ranking'))

    selection_rows = []
    for name, prediction in predictions.items():
        for i, c in enumerate(cases):
            if c['dataset_role'] == 'gate' or c['kind'] != 'pair': continue
            old = predictions['baseline'][i]['candidates'][2]; new = prediction[i]['candidates'][2]
            known = evaluate_fit(new, t[i], q[i], c); old_known = evaluate_fit(old, t[i], q[i], c)
            below = min(new['charge_pe']) < 500.; accepted = saved[name, i]['fit']['birth_accepted']
            selection_rows.append(dict(method=name, index=i, dataset_role=c['dataset_role'], condition=c['condition'],
                gap_ns=c['gap_ns'], eligible=known['target'] == 2, known_K_success=known['success'],
                blind_success=saved[name, i]['success'], candidate_min_Q_below500=below,
                birth_accepted=accepted, successful_K2_rejected=known['success'] and not accepted,
                candidate_success_lost_vs_baseline=old_known['success'] and not known['success'],
                old_eligible_new_belowcut=min(old['charge_pe']) >= 500. and below,
                baseline_success_lost_with_new_belowcut=saved['baseline', i]['success'] and not saved[name, i]['success'] and below))
    selection_summary = []
    for dataset in ('older', 'latest'):
        for name in LABELS:
            for condition in CONDITIONS:
                for gap in (8., 40., 200.):
                    rr = [r for r in selection_rows if r['dataset_role'] == dataset and r['method'] == name
                          and r['condition'] == condition and r['gap_ns'] == gap]
                    assert len(rr) == 40
                    keys = [key for key in rr[0] if key not in ('method', 'index', 'dataset_role', 'condition', 'gap_ns')]
                    selection_summary.append(dict(dataset_role=dataset, method=name, condition=condition, gap_ns=gap,
                        frames=len(rr), **{key: sum(r[key] for r in rr) for key in keys}))
    save(OUTPUT / 'candidate-selection-audit.json', dict(summary=selection_summary, rows=selection_rows,
        role='Posthoc description of stored candidates and frozen gate decisions; no new fits or selection changes',
        warning='Truth Q>=500 defines the evaluation target; fitted Q>=500 is a separate existing selection guard. Known-K results use true count and are diagnostic.'))

    diagnostic = rows(OUTPUT / 'diagnostics.truth.jsonl')
    assert len({(r['index'], r['kind']) for r in diagnostic}) == len(diagnostic) == 247
    diagnostic_by_id = {r['index']: r for r in diagnostic if r['kind'] == 'parameters'}
    assert set(diagnostic_by_id) == set(protocol['diagnostics']['indices'])
    diagnostic_scores = []; comparisons = []; rank_rows = []
    for i, d in diagnostic_by_id.items():
        c = cases[i]; truth_cost = cost_at(response, noise, voltage[i], t[i], q[i])['cost']
        assert abs(truth_cost - d['truth_point']['cost']) < 1e-6
        assert d['truth_start']['cost'] <= truth_cost + 1e-6
        assert d['fixed_times']['times_ns'] == t[i].tolist()
        fits = {name + '_known_K': prediction[i]['candidates'][2] for name, prediction in predictions.items()}
        fits.update(truth_start=d['truth_start'], fixed_times=d['fixed_times'])
        for method, fit in fits.items():
            score = evaluate_fit(fit, t[i], q[i], c)
            diagnostic_scores.append(dict(score, method=method, index=i, condition=c['condition'], gap_ns=c['gap_ns'],
                                          fit=fit, truth_used=True))
        comparisons.append(dict(index=i, condition=c['condition'], gap_ns=c['gap_ns'],
            original_minus_truth_start=d['baseline_double']['cost'] - d['truth_start']['cost'],
            weak_search_minus_truth_start=predictions['weak_extra'][i]['candidates'][2]['cost'] - d['truth_start']['cost'],
            exact_truth_minus_truth_start=truth_cost - d['truth_start']['cost']))
        rank_bank = [d['baseline_double'], d['truth_start'], d['fixed_times']]
        rank_rows.append(dict(index=i, condition=c['condition'], gap_ns=c['gap_ns'],
            original=evaluate_fit(rank_bank[0], t[i], q[i], c)['success'],
            with_true_start_ranked_by_nll=evaluate_fit(min(rank_bank[:2], key=lambda f: f['cost']), t[i], q[i], c)['success'],
            with_true_start_and_fixed_times_ranked_by_nll=evaluate_fit(min(rank_bank, key=lambda f: f['cost']), t[i], q[i], c)['success']))
    diagnostic_summary = []
    for condition in ('high_weak', 'weak_high'):
        for gap in (8., 40., 200.):
            for method in fits:
                rr = [r for r in diagnostic_scores if r['condition'] == condition and r['gap_ns'] == gap and r['method'] == method]
                diagnostic_summary.append(dict(method=method, condition=condition, gap_ns=gap, frames=len(rr),
                    success=sum(r['success'] for r in rr), matched_two=sum(r['matched'] == 2 for r in rr),
                    not_converged=sum(not r['fit'].get('converged', True) for r in rr)))
    profiles = [r for r in diagnostic if r['kind'] == 'profile']; dense = next(r for r in diagnostic if r['kind'] == 'dense')
    assert {r['index'] for r in profiles} == set(protocol['diagnostics']['profile_indices'])
    coarse = next(r for r in profiles if r['index'] == dense['index'])
    ti = [int(np.argmin(abs(np.asarray(dense['weak_times_ns']) - v))) for v in coarse['weak_times_ns']]
    qi = [int(np.argmin(abs(np.asarray(dense['weak_charge_pe']) - v))) for v in coarse['weak_charge_pe']]
    common_error = float(np.max(abs(np.array(coarse['costs']) - np.array(dense['costs'])[np.ix_(qi, ti)])))
    assert common_error < 1e-6 and dense['grid_min_cost'] <= coarse['grid_min_cost'] + 1e-6
    profile_summary = []
    for p in profiles + [dense]:
        costs = np.array(p['costs']); j, i = np.unravel_index(np.argmin(costs), costs.shape)
        index = p['index']; c = cases[index]
        score = evaluate_fit(p['polished'], t[index], q[index], c)
        profile_summary.append(dict(index=index, gap_ns=c['gap_ns'], replicate=c['replicate'], kind=p['kind'],
            grid_min_cost=p['grid_min_cost'], boundary_minimum=bool(j in (0, costs.shape[0] - 1) or i in (0, costs.shape[1] - 1)),
            nuisance_nonconverged=int((~np.array(p['converged'])).sum()),
            polished_success=score['success'], polished_not_converged=not p['polished']['converged'],
            baseline_minus_polished=predictions['baseline'][index]['candidates'][2]['cost'] - p['polished']['cost'],
            weak_search_minus_polished=predictions['weak_extra'][index]['candidates'][2]['cost'] - p['polished']['cost']))
    save(OUTPUT / 'diagnostic-summary.json', dict(summary=diagnostic_summary, cost_comparisons=comparisons,
        profiles=profile_summary, dense_grid_common_max_cost_error=common_error,
        dense_grid_gain=coarse['grid_min_cost'] - dense['grid_min_cost'], rows=diagnostic_scores,
        warning='Known count and true-centred profiles are diagnostics; neither local optimization nor noisy truth costs prove a fundamental limit'))
    rank_summary = []
    for condition in ('high_weak', 'weak_high'):
        for gap in (8., 40., 200.):
            rr = [r for r in rank_rows if r['condition'] == condition and r['gap_ns'] == gap]
            rank_summary.append(dict(condition=condition, gap_ns=gap, frames=len(rr),
                **{key: sum(r[key] for r in rr) for key in ('original', 'with_true_start_ranked_by_nll', 'with_true_start_and_fixed_times_ranked_by_nll')}))
    audit_path = OUTPUT / 'exploratory-rank-audit.json'
    if audit_path.exists(): assert load(audit_path)['summary'] == rank_summary
    save(audit_path, dict(summary=rank_summary, rows=rank_rows, truth_used=True,
        role='Posthoc diagnostic of stored candidates only; no new fitting and no change to the declared blind experiment',
        warning='Not a blind performance bound; lower-cost noisy-data solutions need not minimize parameter error'))

    block = np.asarray(load(OUTPUT / 'mother-blocks.truth.json')); bootstrap = protocol['bootstrap']
    rng = np.random.default_rng(bootstrap['seed']); B = bootstrap['replicates']; paired = []; aggregate = []; point_intervals = []
    for dataset in ('older', 'latest'):
        weights = np.zeros((B, 40), int)
        for start in (0, 20):
            draws = rng.integers(start, start + 20, size=(B, 20))
            for j in range(B): weights[j] += np.bincount(draws[j], minlength=40)
        estimates = {}; point = {}
        for name in LABELS:
            grouped = {True: [], False: []}
            for condition in CONDITIONS:
                for gap in (8., 40., 200.):
                    ids = [i for i, c in enumerate(cases) if c['dataset_role'] == dataset and c['kind'] == 'pair'
                           and c['condition'] == condition and c['gap_ns'] == gap]
                    den = np.bincount(block[ids], weights=[saved[name, i]['target'] == 2 for i in ids], minlength=40)
                    num = np.bincount(block[ids], weights=[saved[name, i]['target'] == 2 and saved[name, i]['success'] for i in ids], minlength=40)
                    assert np.all(weights @ den > 0)
                    estimate = (weights @ num) / (weights @ den)
                    grouped['weak' in condition].append(estimate); point[name, condition, gap] = estimate
                    if name != 'baseline':
                        changes = [int(saved[name, i]['success']) - int(saved['baseline', i]['success']) for i in ids if saved[name, i]['target'] == 2]
                        paired.append(dict(dataset=dataset, method=name, condition=condition, gap_ns=gap, frames=len(changes),
                            wins=sum(v > 0 for v in changes), losses=sum(v < 0 for v in changes), ties=sum(v == 0 for v in changes)))
            estimates[name] = {kind: np.mean(grouped[weak], axis=0) for weak, kind in [(True, 'weak'), (False, 'strong')]}
        for method, base in [('legacy_extra', 'baseline'), ('weak_extra', 'baseline'), ('weak_extra', 'legacy_extra')]:
            for group in ('weak', 'strong'):
                delta = estimates[method][group] - estimates[base][group]
                aggregate.append(dict(dataset=dataset, method=method, baseline=base, group=group,
                    difference=report['results'][method][dataset][group + '_mean'] - report['results'][base][dataset][group + '_mean'],
                    paired_block95=np.quantile(delta, [.025, .975]).tolist()))
            for condition in ('high_weak', 'weak_high'):
                for gap in (8., 40., 200.):
                    delta = point[method, condition, gap] - point[base, condition, gap]
                    point_intervals.append(dict(dataset=dataset, method=method, baseline=base, condition=condition, gap_ns=gap,
                        paired_block95=np.quantile(delta, [.025, .975]).tolist()))
    save(OUTPUT / 'paired-development.json', dict(aggregate=aggregate, point_intervals=point_intervals, pairs=paired,
        interpretation='Exposed development; mother-block intervals do not undo method-selection bias', **bootstrap))
    save(RECORD / 'science-validation.json', dict(passed=True, scores_replayed=len(scores), summary_rows=len(summary),
        diagnostic_scores=len(diagnostic_scores), diagnostic_records=len(diagnostic), native_mother_groups_disjoint=True,
        original_inputs_unchanged=True, freeze_unchanged=True, gate_recalibrated_for_each_search=True,
        maximum_candidate_relative_cost_replay_error=cost_replay_error,
        original_zero_and_one_candidates_unchanged=True, selected_not_converged=convergence,
        dense_grid_common_max_cost_error=common_error, chosen=report['chosen'], new_roots=0,
        figures_pending=True, notebook_pending=True))
    print('WEAK_SEARCH_VALIDATED', len(scores), len(summary), len(diagnostic), flush=True)


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from matplotlib.lines import Line2D
    from .alpha_pair_joint_likelihood import ConservingResponse, JointNoise
    assert load(RECORD / 'science-validation.json')['passed']
    folder = OUTPUT / 'figures'; folder.mkdir(exist_ok=True)
    report = load(OUTPUT / 'report.json'); summary = report['summary']
    cases = load(OUTPUT / 'labels.cases.json'); scored = rows(OUTPUT / 'scores.jsonl')
    scores = {(r['method'], r['index']): r for r in scored}
    predictions = {name: {r['index']: r['candidates'] for r in rows(OUTPUT / f'predictions.{name}.jsonl')} for name in LABELS}
    diagnostic = rows(OUTPUT / 'diagnostics.truth.jsonl')
    profiles = [r for r in diagnostic if r['kind'] == 'profile']
    profiles.sort(key=lambda r: (cases[r['index']]['replicate'], cases[r['index']]['gap_ns']))
    norm = SymLogNorm(linthresh=1., vmin=0., vmax=300.)

    def profile_panel(ax, profile, cost_reference=None):
        weak = profile['weak_truth_index']; strong = 1 - weak
        offset = profile['true_t'][strong]; charge = profile['true_q'][weak]
        x = np.array(profile['weak_times_ns']) - offset; y = np.array(profile['weak_charge_pe']) / charge
        cost = np.array(profile['costs']); reference = float(cost.min()) if cost_reference is None else cost_reference
        mesh = ax.pcolormesh(x, y, np.maximum(cost - reference, 0.), shading='nearest', norm=norm, cmap='viridis')
        ax.plot(profile['true_t'][weak] - offset, 1., '*', color='#ee7733', markeredgecolor='black', markersize=13)
        outside = []
        for name, marker, color in [('baseline', 'o', 'white'), ('weak_extra', 'x', '#ee3377')]:
            fit = predictions[name][profile['index']][2]; smaller = int(np.argmin(fit['charge_pe']))
            xx = fit['times_ns'][smaller] - offset; yy = fit['charge_pe'][smaller] / charge
            if x[0] <= xx <= x[-1] and y[0] <= yy <= y[-1]:
                ax.plot(xx, yy, marker, color=color, markeredgecolor=color, markersize=7, fillstyle='none', mew=1.7)
            else: outside.append('original' if name == 'baseline' else 'weak-start')
        if outside: ax.text(.02, .02, 'Outside: ' + ', '.join(outside), transform=ax.transAxes, fontsize=7, color='white',
                            bbox=dict(facecolor='black', alpha=.5, edgecolor='none'))
        ax.set(yscale='log', ylim=(y[0], y[-1]), xlim=(x[0], x[-1]),
               xlabel='Weak time - true strong time (ns)', ylabel='Weak Q / true weak Q')
        ax.set_yticks([.25, .5, 1., 2., 3.], labels=['.25', '.5', '1', '2', '3'])
        return mesh

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for ax, p, panel in zip(axes.flat, profiles, 'abcdef'):
        mesh = profile_panel(ax, p); c = cases[p['index']]
        ax.set_title(f"({panel}) Gap {c['gap_ns']:g} ns, pair {c['replicate']}")
    color_axis = fig.add_axes([.89, .21, .016, .59])
    fig.colorbar(mesh, cax=color_axis, extend='max', label='NLL - local grid minimum (not a confidence level)')
    handles = [Line2D([], [], marker='*', color='#ee7733', linestyle='none', label='True weak event (diagnostic)'),
               Line2D([], [], marker='o', color='#555555', fillstyle='none', linestyle='none', label='Original smaller-Q K=2 output'),
               Line2D([], [], marker='x', color='#ee3377', linestyle='none', label='Weak-start smaller-Q K=2 output')]
    fig.legend(handles=handles, loc='lower center', ncol=3, fontsize=9)
    fig.suptitle('Figure 0004-66-1 | Six prespecified, truth-centred local profiles; strong Q/t reoptimized at every point')
    fig.subplots_adjust(left=.07, right=.83, bottom=.11, top=.9, hspace=.34, wspace=.3)
    fig.savefig(folder / '0004-66-1.png', dpi=150); plt.close(fig)

    dense = next(r for r in diagnostic if r['kind'] == 'dense')
    coarse = next(r for r in profiles if r['index'] == dense['index'])
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    for ax, p, label in zip(axes, (coarse, dense), ('(a) 13 x 9 grid', '(b) 25 x 17 grid')):
        mesh = profile_panel(ax, p, dense['grid_min_cost']); ax.set_title(label)
    color_axis = fig.add_axes([.87, .18, .02, .64])
    fig.colorbar(mesh, cax=color_axis, extend='max', label='NLL - refined grid minimum')
    fig.suptitle('Figure 0004-66-2 | Prespecified 40 ns pair 0: grid-density check; local diagnostic only')
    fig.subplots_adjust(left=.08, right=.81, bottom=.15, top=.84, wspace=.34)
    fig.savefig(folder / '0004-66-2.png', dpi=150); plt.close(fig)

    ds = load(OUTPUT / 'diagnostic-summary.json')['summary']
    methods = ['baseline_known_K', 'legacy_extra_known_K', 'weak_extra_known_K', 'truth_start', 'fixed_times']
    labels = ['Known K, original search', 'Known K, extra original starts', 'Known K, extra weak starts',
              'True Q/t initialization', 'Fixed true times, fit Q']
    colors = ['#555555', '#0077bb', '#009988', '#aa4499', '#ee7733']
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6))
    for ax, condition, panel in zip(axes, ('high_weak', 'weak_high'), 'ab'):
        for j, (method, label, color) in enumerate(zip(methods, labels, colors)):
            rr = sorted([r for r in ds if r['method'] == method and r['condition'] == condition], key=lambda r: r['gap_ns'])
            ax.bar(np.arange(3) + (j - 2) * .15, [100 * r['success'] / r['frames'] for r in rr], width=.145, label=label, color=color)
        ax.set(title=f'({panel}) {condition}', xticks=np.arange(3), xticklabels=['8', '40', '200'],
               xlabel='Separation (ns)', ylabel='Strict pair recovery (%)', ylim=(0, 103)); ax.grid(axis='y', alpha=.2)
    handles, labels = axes[0].get_legend_handles_labels(); fig.legend(handles, labels, loc='lower center', ncol=3, fontsize=9)
    fig.suptitle('Figure 0004-67-1 | Truth-assisted diagnosis on all latest high/weak pairs; not blind performance')
    fig.tight_layout(rect=(0, .12, 1, .94)); fig.savefig(folder / '0004-67-1.png', dpi=150); plt.close(fig)

    fig, axes = plt.subplots(4, 5, figsize=(18, 14), sharex=True, sharey=True)
    for dataset, start in [('older', 0), ('latest', 2)]:
        for j, condition in enumerate(CONDITIONS):
            ax = axes[start + j // 5, j % 5]
            for name in LABELS:
                rr = sorted([r for r in summary if r['method'] == name and r['dataset_role'] == dataset and r['condition'] == condition], key=lambda r: r['gap_ns'])
                ax.plot([r['gap_ns'] for r in rr], [100 * r['efficiency'] for r in rr], 'o-', color=COLORS[name], label=LABELS[name], lw=1.5)
            ax.set(title=f'{dataset}: {condition}', xscale='log', xticks=[8, 40, 200], xticklabels=['8', '40', '200'], ylim=(-3, 103))
            ax.grid(alpha=.2)
            if j % 5 == 0: ax.set_ylabel('Strict pair recovery (%)')
            if start + j // 5 == 3: ax.set_xlabel('Separation (ns)')
    handles, labels = axes[0, 0].get_legend_handles_labels(); fig.legend(handles, labels, loc='lower center', ncol=3)
    fig.suptitle('Figure 0004-68-1 | All ten combinations, two exposed development sets; separate recalibrated gates')
    fig.tight_layout(rect=(0, .04, 1, .96)); fig.savefig(folder / '0004-68-1.png', dpi=150); plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(17, 9))
    controls = ['all', 'weak', 'low', 'mid', 'high']
    for row, dataset in enumerate(('older', 'latest')):
        for j, name in enumerate(LABELS):
            rr = [next(r for r in summary if r['method'] == name and r['dataset_role'] == dataset and r['kind'] != 'pair' and r['condition'] == c) for c in controls]
            x = np.arange(5) + (j - 1) * .2
            axes[row, 0].bar(x, [100 * r['extra_fraction'] for r in rr], width=.19, color=COLORS[name], label=LABELS[name])
            axes[row, 1].plot(x, [100 * r['efficiency'] for r in rr], 'o-', color=COLORS[name])
            axes[row, 2].plot(x[1:], [100 * r['q_error_p16_p50_p84'][1] for r in rr[1:]], 'o-', color=COLORS[name])
        for col, title in enumerate(('Extra-output frames (%)', 'Strict success (%)', 'Matched median Q error (%)')):
            axes[row, col].set(title=f'{dataset}: {title}', xticks=np.arange(5), xticklabels=['empty', 'weak', 'low', 'mid', 'high']); axes[row, col].grid(alpha=.2)
        axes[row, 0].axhline(5., color='k', ls='--', lw=.8)
        axes[row, 1].axhline(95., color='k', ls='--', lw=.8)
        for v in (-5., 5.): axes[row, 2].axhline(v, color='k', ls='--', lw=.8)
    handles, labels = axes[0, 0].get_legend_handles_labels(); fig.legend(handles, labels, loc='lower center', ncol=3)
    fig.suptitle('Figure 0004-68-2 | Empty/native-single controls: 80 frames per group in each exposed set')
    fig.tight_layout(rect=(0, .05, 1, .94)); fig.savefig(folder / '0004-68-2.png', dpi=150); plt.close(fig)

    response = ConservingResponse.load(EXPERIMENT / 'joint-response.npz')
    noise = JointNoise.load(EXPERIMENT / 'joint_conserved_nll-noise.npz')
    with np.load(OUTPUT / 'observations.npz') as f: voltage = f['voltage_mV']
    with np.load(LATEST / 'labels.truth.npz') as f: components = f['components_mV']
    examples = []; fig, axes = plt.subplots(3, 2, figsize=(16, 10)); sample = response.sample_grid
    for column, gap in enumerate((40., 200.)):
        index = next(i for i, c in enumerate(cases) if c['dataset_role'] == 'latest' and c['condition'] == 'high_weak' and c['gap_ns'] == gap and c['replicate'] == 0)
        case = cases[index]; all_methods = {name: scores[name, index] for name in LABELS}
        examples.append(dict(index=index, case=case, methods=all_methods))
        axes[0, column].plot(sample, voltage[index], color='#aaaaaa', lw=.8, label='Observed voltage')
        weak = int(np.argmin(all_methods['baseline']['true_q']))
        axes[1, column].plot(sample, components[case['original_index'], weak], color='black', lw=1., label='True weak component (diagnostic)')
        for name in LABELS:
            fit = all_methods[name]['fit']; count = len(fit['times_ns'])
            curves = response.components(fit['times_ns'], fit['charge_pe']) if count else np.empty((0, len(sample)))
            predicted = curves.sum(0) + response.baseline + noise.background_mV
            weak_curve = curves[int(np.argmin(fit['charge_pe']))] if count == 2 else np.zeros_like(sample)
            axes[0, column].plot(sample, predicted, color=COLORS[name], lw=1.1, label=LABELS[name])
            axes[1, column].plot(sample, weak_curve, color=COLORS[name], lw=1.1, label=LABELS[name])
            axes[2, column].plot(sample, voltage[index] - predicted, color=COLORS[name], lw=.7, label=LABELS[name])
        for row, title in enumerate(('Total voltage', 'Weak / smaller-Q output component', 'Observed - fitted voltage')):
            ax = axes[row, column]; ax.set(title=f'({"abcdef"[2*row+column]}) {gap:g} ns: {title}', xlim=(0, 2048), ylabel='mV'); ax.grid(alpha=.2)
            for onset in all_methods['baseline']['true_t']: ax.axvline(onset, ls=':', color='k', alpha=.3, lw=.8)
            if row == 2: ax.set_xlabel('Record time (ns)')
        axes[1, column].legend(fontsize=8)
    handles, labels = axes[0, 0].get_legend_handles_labels(); fig.legend(handles, labels, loc='lower center', ncol=4, fontsize=9)
    fig.suptitle('Figure 0004-69-1 | Prespecified latest high_weak pair 0, full 0-2048 ns; every fitted curve comes from blind selection')
    fig.tight_layout(rect=(0, .05, 1, .95)); fig.savefig(folder / '0004-69-1.png', dpi=150); plt.close(fig)
    save(OUTPUT / 'examples.json', examples)
    save(RECORD / 'figure-delivery.json', dict(passed=True, figures={str(p.relative_to(PROJECT)): digest(p) for p in sorted(folder.glob('*.png'))}))
    print('WEAK_SEARCH_FIGURES', len(list(folder.glob('*.png'))), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('validate', 'figures'))
    args = parser.parse_args(); globals()[args.command]()
