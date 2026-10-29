"""Descriptive response-history diagnostics; never a blind-reconstruction score.

Historical error terms are added on the same finite time-matched subset. The
offline control compares shared and separate cell histories, preserving optical
parents. Recovery-disabled results are reported in charge only.
"""
from __future__ import annotations

import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .sipm_energy_calibration import invert_charge
from .sipm_energy_history_study import (
    EXPERIMENT, OUTPUT, RECORD, BLOCKS, GAPS, METHODS, check_freeze,
    digest, lines, load, save,
)

CONDITIONS = ['weak_weak', 'low_weak', 'weak_low', 'high_weak', 'weak_high',
              'low_low', 'mid_mid', 'high_high', 'low_high', 'high_low']
LABELS = {'model_relative_residual': 'Model GLS', 'joint_conserved_nll': 'Joint NLL'}
ORDERS = ['high_weak', 'weak_high']


def csv_file(name, rows):
    with (OUTPUT / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def finite_number(value):
    return float(value) if np.isfinite(value) else None


def mean_or_none(values):
    return float(np.mean(values)) if len(values) else None


def interval(values):
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    return [float(x) for x in np.quantile(values, [.025, .975])] if len(values) else [None, None]


def finite_mean(values, axes):
    """Conditional mean, retaining an explicit zero-denominator NaN."""
    values = np.asarray(values, float)
    count = np.isfinite(values).sum(axis=axes)
    return np.divide(np.nansum(values, axis=axes), count,
                     out=np.full(np.shape(count), np.nan), where=count > 0)


def historical_tables():
    components = lines(OUTPUT / 'historical-components.jsonl')
    fitted = lines(OUTPUT / 'historical-fitted-components.jsonl')
    cases = load(OUTPUT / 'inputs/cases.json')
    scores = lines(OUTPUT / 'inputs/scores.jsonl')
    response, decomposition, counts = [], [], []
    for condition in CONDITIONS:
        for gap in GAPS:
            ids = {i for i, c in enumerate(cases)
                   if c['kind'] == 'pair' and c['condition'] == condition and c['gap_ns'] == gap}
            assert len(ids) == 40
            for position in ['first', 'second']:
                rr = [r for r in components if r['index'] in ids and r['position'] == position]
                assert len(rr) == 40
                response.append(dict(condition=condition, gap_ns=gap, position=position, n=len(rr),
                    species=rr[0]['species'],
                    median_charge_change=float(np.median([r['charge_change_relative'] for r in rr])),
                    mean_charge_change=float(np.mean([r['charge_change_relative'] for r in rr])),
                    mean_native_energy_error=float(np.mean([r['native_energy_relative'] for r in rr])),
                    mean_shared_energy_error=float(np.mean([r['shared_energy_relative'] for r in rr]))))
                for method in METHODS:
                    ff = [r for r in fitted if r['index'] in ids and r['position'] == position and r['method'] == method]
                    matched = [r for r in ff if r['matched']]
                    valid = [r for r in matched if r['fit_energy_relative'] is not None]
                    terms = [mean_or_none([r[k] for r in valid]) for k in
                             ['native_energy_relative', 'response_energy_change', 'fit_energy_change', 'fit_energy_relative']]
                    if valid:
                        assert abs(sum(terms[:3]) - terms[3]) < 1e-12
                    decomposition.append(dict(condition=condition, gap_ns=gap, position=position,
                        species=rr[0]['species'], method=method, n=len(ff), time_matched=len(matched),
                        unmatched=len(ff)-len(matched), matched_outside_energy_domain=len(matched)-len(valid),
                        finite_matched=len(valid), native_mean=terms[0], response_mean=terms[1],
                        fit_mean=terms[2], total_mean=terms[3],
                        energy10_among_all=sum(abs(r['fit_energy_relative']) <= .1 for r in valid),
                        cancellation_count=sum(r['cancellation'] for r in valid),
                        sigma68_fit_energy=None if not valid else float(np.diff(np.quantile(
                            [r['fit_energy_relative'] for r in valid], [.16, .84]))[0]/2),
                        median_fit_charge_error=None if not valid else float(np.median([
                            r['fit_charge_pe']/r['shared_charge_pe']-1 for r in valid]))))
            for method in METHODS:
                ss = [s for s in scores if s['method'] == method and s['index'] in ids]
                assert len(ss) == 40
                counts.append(dict(condition=condition, gap_ns=gap, method=method, frames=40,
                    true_components=80, outputs=sum(s['output_count'] for s in ss),
                    missed=sum(s['missed'] for s in ss), extra=sum(s['extra'] for s in ss),
                    split=sum(s['split'] for s in ss), merge=sum(s['merge'] for s in ss),
                    strict_energy10=sum(s['strict_10'] for s in ss),
                    successful_pairs=sum(s['frame_success_10'] for s in ss)))
    csv_file('historical-response-summary.csv', response)
    csv_file('historical-decomposition-summary.csv', decomposition)
    csv_file('historical-frame-counts.csv', counts)
    return response, decomposition, counts


def counterfactual_tables():
    plan = load(EXPERIMENT / 'simulation-plan.json')
    results = {}
    for task in plan:
        result = load(OUTPUT / 'simulations' / (task['id'] + '.json'))
        assert result['task'] == task and result['charge_closure_error_pe'] < 1e-8
        results[task['id']] = result
    assert len(results) == 912
    cases = load(OUTPUT / 'inputs/cases.json')
    with np.load(OUTPUT / 'inputs/labels.npz') as file:
        energies = file['energy_MeV']
    parent_energy = {c['parents'][0]: float(energies[i, 0])
                     for i, c in enumerate(cases) if c['kind'] == 'single'}
    model = load(EXPERIMENT / 'energy-model.json')
    protocol = load(EXPERIMENT / 'protocol.json')
    rng = np.random.default_rng(protocol['bootstrap_seed'])
    # All four columns use identical parent and replicate draws. Each sampled
    # parent contributes a ratio of four-replicate mean charges, with equal
    # optical-parent weight. These intervals describe the eight selected pairs.
    boot_parent = rng.integers(0, 8, (2000, 8))
    boot_mc = rng.integers(0, 4, (2000, 8, 4))
    raw, charge_summary, contrast_summary, energy_summary = [], [], [], []
    for condition in ORDERS:
        for gap in GAPS:
            for species in ['weak', 'high']:
                position = condition.split('_').index(species)
                charge = np.zeros((8, 4, 4))  # joint on, single on, joint off, single off
                energy = np.zeros((8, 4, 2))  # joint on, single on
                for slot, block in enumerate(BLOCKS):
                    parent = block if species == 'weak' else 240 + block
                    edep = parent_energy[parent]
                    for mc in range(4):
                        for mode, recovery in enumerate([True, False]):
                            prefix = f'p{slot}_m{mc}_r{int(recovery)}'
                            pair = results[f'{prefix}_{condition}_{int(gap)}']['components'][position]
                            offset = gap if position else 0
                            single = results[f'{prefix}_single_{species}_{int(offset)}']['components'][0]
                            assert pair['parent'] == single['parent'] == parent
                            qj, qs = pair['charge_pe'], single['charge_pe']
                            charge[slot, mc, mode*2:mode*2+2] = [qj, qs]
                            ej, es = (np.asarray(invert_charge(model, [qj, qs])) / edep - 1
                                      if recovery else np.array([np.nan, np.nan]))
                            if recovery:
                                energy[slot, mc] = [ej, es]
                            raw.append(dict(condition=condition, gap_ns=gap, species=species,
                                position='first' if position == 0 else 'second', block=block, mc=mc,
                                parent=parent, recovery=recovery, edep_MeV=edep,
                                joint_charge_pe=qj, isolated_charge_pe=qs, relative_charge_change=qj/qs-1,
                                joint_direct=pair['direct_avalanches'], isolated_direct=single['direct_avalanches'],
                                joint_mean_recovery=pair['mean_recovery'], isolated_mean_recovery=single['mean_recovery'],
                                joint_energy_relative=finite_number(ej), isolated_energy_relative=finite_number(es)))
                mean_q = charge.mean(axis=1)
                parent_ratio = mean_q[:, [0, 2]] / mean_q[:, [1, 3]] - 1
                sample_q = charge[boot_parent[:, :, None], boot_mc].mean(axis=2)
                boot_ratio = (sample_q[:, :, [0, 2]] / sample_q[:, :, [1, 3]] - 1).mean(axis=1)
                for mode, recovery in enumerate([True, False]):
                    low, high = interval(boot_ratio[:, mode])
                    charge_summary.append(dict(condition=condition, gap_ns=gap, species=species,
                        recovery=recovery, optical_pairs=8, replicates_per_pair=4,
                        mean_relative_charge_change=float(parent_ratio[:, mode].mean()),
                        bootstrap_low95=low, bootstrap_high95=high))
                low, high = interval(boot_ratio[:, 0] - boot_ratio[:, 1])
                contrast_summary.append(dict(condition=condition, gap_ns=gap, species=species,
                    mean_on_minus_off=float(np.mean(parent_ratio[:, 0] - parent_ratio[:, 1])),
                    bootstrap_low95=low, bootstrap_high95=high, optical_pairs=8, replicates_per_pair=4))
                eb = energy[boot_parent[:, :, None], boot_mc]
                for k, kind in enumerate(['joint', 'isolated']):
                    values = energy[:, :, k]
                    low, high = interval(finite_mean(eb[:, :, :, k], (1, 2)))
                    energy_summary.append(dict(condition=condition, gap_ns=gap, species=species, kind=kind,
                        n=32, valid=int(np.isfinite(values).sum()),
                        within_energy10=int(np.sum(np.isfinite(values) & (np.abs(values) <= .1))),
                        mean_energy_relative=finite_number(finite_mean(values, (0, 1))),
                        median_energy_relative=finite_number(np.nanmedian(values)),
                        sigma68_relative=finite_number(np.diff(np.nanquantile(values, [.16, .84]))[0]/2),
                        bootstrap_low95=low, bootstrap_high95=high))
    window = []
    for slot, block in enumerate(BLOCKS):
        for species in ['weak', 'high']:
            long = results[f'p{slot}_long_{species}']['components'][0]
            short = results[f'p{slot}_m0_r1_single_{species}_0']['components'][0]
            window.append(dict(block=block, species=species, parent=long['parent'],
                charge_20us_pe=long['charge_pe'], charge_4096ns_pe=short['charge_pe'],
                long_run_fraction_after_2048ns=1-long['q_before_2048_pe']/long['charge_pe'],
                long_run_fraction_after_4096ns=1-long['q_before_4096_pe']/long['charge_pe'],
                prefix_difference_relative=long['q_before_4096_pe']/short['charge_pe']-1,
                photons_dropped_4096ns=short['photons_dropped']))
    # Extending the window changes the Poisson dark-event schedule, even at the
    # same seed. Prefix identity is therefore NOT assumed or asserted.
    csv_file('counterfactual-components.csv', raw)
    csv_file('counterfactual-summary.csv', charge_summary)
    csv_file('recovery-contrast.csv', contrast_summary)
    csv_file('known-charge-energy-summary.csv', energy_summary)
    csv_file('window-control.csv', window)
    return charge_summary, contrast_summary, energy_summary, window


def figures(response, decomposition, charge, energy):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'figure.dpi': 120, 'savefig.dpi': 160, 'svg.fonttype': 'none'})
    saved = []

    def finish(fig, name):
        for suffix in ['.png', '.svg']:
            fig.savefig(OUTPUT / (name+suffix), bbox_inches='tight')
        plt.close(fig)
        saved.append(name)

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 6.5), layout='constrained')
    lookup = {(r['condition'], r['gap_ns'], r['position']): r for r in response}
    limit = max(5, np.ceil(max(abs(r['median_charge_change'])*100 for r in response)/5)*5)
    for panel, (ax, position) in enumerate(zip(axes, ['first', 'second'])):
        values = np.array([[lookup[c, g, position]['median_charge_change']*100 for g in GAPS]
                           for c in CONDITIONS])
        im = ax.imshow(values, vmin=-limit, vmax=limit, cmap='RdBu_r', aspect='auto')
        for i in range(10):
            for j in range(3):
                ax.text(j, i, f'{values[i, j]:+.1f}%', ha='center', va='center',
                        color='white' if abs(values[i, j]) > .65*limit else 'black')
        ax.set_xticks(range(3), [f'{g:g}' for g in GAPS])
        ax.set_yticks(range(10), [c.replace('_', ' -> ') for c in CONDITIONS])
        ax.set(xlabel='Pulse separation / ns', title=f'({chr(97+panel)}) {position.title()} component; n=40 / cell')
    fig.colorbar(im, ax=axes, label='Median shared Q / native isolated Q - 1 / %', shrink=.85)
    fig.suptitle('Fig. 0007-15-1 | Historical charge shift, matched optical parents\nDifferent digitizers, windows and RNG: not a pure recovery effect')
    finish(fig, '0007-15-1_historical_charge')

    fig, axes = plt.subplots(2, 1, figsize=(11.6, 9), layout='constrained')
    colors = ['#929292', '#d19334', '#477eaa', '#222222']
    for panel, (ax, method) in enumerate(zip(axes, METHODS)):
        rr = [next(r for r in decomposition if r['method'] == method and r['condition'] == c
                   and r['gap_ns'] == g and r['species'] == 'weak') for c in ORDERS for g in GAPS]
        for k, (field, label) in enumerate([('native_mean', 'Isolated calibration'), ('response_mean', 'Response change'),
                                           ('fit_mean', 'Fit change'), ('total_mean', 'Total error')]):
            vals = [np.nan if r[field] is None else r[field]*100 for r in rr]
            ax.bar(np.arange(6)+(k-1.5)*.18, vals, width=.18, color=colors[k], label=label)
        ax.axhline(0, lw=.8, color='gray')
        ax.set_xticks(range(6), [f"{r['condition'].replace('_', ' -> ')}\n{r['gap_ns']:g} ns; n={r['finite_matched']}/40" for r in rr])
        ax.set(ylabel='Mean relative energy error / %', title=f'({chr(97+panel)}) {LABELS[method]}: weak component')
        ax.grid(axis='y', alpha=.13)
        ax.legend(ncol=4, fontsize=8)
    fig.suptitle('Fig. 0007-16-1 | Exact additive energy-error decomposition\nSame finite time-matched subset for every term; unmatched events excluded from means')
    finish(fig, '0007-16-1_error_decomposition')

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.8), layout='constrained', sharey=True)
    for panel, (ax, species) in enumerate(zip(axes, ['weak', 'high'])):
        for order, color in zip(ORDERS, ['#ba5540', '#3b789c']):
            for recovery, style in [(True, '-'), (False, '--')]:
                rr = [next(r for r in charge if r['condition'] == order and r['gap_ns'] == g
                           and r['species'] == species and r['recovery'] == recovery) for g in GAPS]
                values = np.array([r['mean_relative_charge_change'] for r in rr])*100
                low = np.array([r['bootstrap_low95'] for r in rr])*100
                high = np.array([r['bootstrap_high95'] for r in rr])*100
                # Plot intervals directly; percentile intervals need not contain a point estimate.
                ax.vlines(GAPS, low, high, color=color, alpha=.65)
                ax.plot(GAPS, values, style, color=color, marker='o' if recovery else 's',
                        label=f"{order.replace('_', ' -> ')}; recovery {'on' if recovery else 'off'}")
        ax.axhline(0, color='gray', lw=.8)
        ax.set_xscale('log'); ax.set_xticks(GAPS, ['8', '40', '200'])
        ax.set(xlabel='Pulse separation / ns', ylabel='Mean joint Q / isolated Q - 1 / %',
               title=f'({chr(97+panel)}) {species.title()} component')
        ax.grid(alpha=.15)
    axes[0].legend(fontsize=8)
    fig.suptitle('Fig. 0007-17-1 | Shared versus separate cell histories\n8 optical pairs x 4 digitizer replicates; paired nested bootstrap 95% intervals')
    finish(fig, '0007-17-1_recovery_control')

    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8), layout='constrained', sharex=True)
    for i, order in enumerate(ORDERS):
        for j, species in enumerate(['weak', 'high']):
            ax = axes[i, j]
            for kind, color, marker in [('joint', '#ba5540', 'o'), ('isolated', '#3b789c', 's')]:
                rr = [next(r for r in energy if r['condition'] == order and r['gap_ns'] == g
                           and r['species'] == species and r['kind'] == kind) for g in GAPS]
                values = [np.nan if r['mean_energy_relative'] is None else r['mean_energy_relative']*100 for r in rr]
                low = [np.nan if r['bootstrap_low95'] is None else r['bootstrap_low95']*100 for r in rr]
                high = [np.nan if r['bootstrap_high95'] is None else r['bootstrap_high95']*100 for r in rr]
                ax.vlines(GAPS, low, high, color=color, alpha=.65)
                ax.plot(GAPS, values, marker=marker, color=color, label=kind.title()+' true Q')
                for x, y, r in zip(GAPS, values, rr):
                    if r['valid'] != r['n']:
                        ax.annotate(f"{r['valid']}/{r['n']} finite", (x, y), fontsize=8)
            ax.axhline(0, color='gray', lw=.8); ax.axhspan(-10, 10, color='gray', alpha=.08)
            ax.set_xscale('log'); ax.set_xticks(GAPS, ['8', '40', '200']); ax.grid(alpha=.12)
            ax.set(xlabel='Pulse separation / ns', ylabel='Mean relative E error / %',
                   title=f"({chr(97+i*2+j)}) {order.replace('_', ' -> ')}: {species}")
            ax.legend(fontsize=8)
    fig.suptitle('Fig. 0007-18-1 | Frozen isolated calibration applied to known charges\nRecovery on only; truth-assisted diagnostic, not blind event recovery')
    finish(fig, '0007-18-1_known_charge_energy')
    return saved


def main():
    check_freeze()
    assert load(RECORD/'historical-validation.json')['passed']
    assert load(RECORD/'simulation-validation.json')['passed']
    for name, sha in load(EXPERIMENT/'input-hashes.json').items():
        assert digest(OUTPUT/name) == sha, name
    response, decomposition, counts = historical_tables()
    charge, contrast, energy, window = counterfactual_tables()
    names = figures(response, decomposition, charge, energy)
    report = dict(kind='Retrospective mechanistic diagnostic; no new fits, transport or calibration',
        historical_response=response, historical_decomposition=decomposition, historical_counts=counts,
        counterfactual_charge=charge, recovery_contrast=contrast, known_charge_energy=energy, window_controls=window,
        figures=names, figures_sha256={n: digest(OUTPUT/(n+'.png')) for n in names},
        window_note='The longer window changes dark-event scheduling; equal-seed prefixes need not be identical.',
        independent_sampling='Eight reused optical pairs, four new digitizer seeds each; not 912 independent optical events.',
        calibration_unchanged=True, fits_unchanged=True, default_recovery_unchanged=True)
    save(OUTPUT/'report.json', report)
    save(RECORD/'report-validation.json', dict(passed=True, historical_components=2400,
        fitted_component_rows=4800, historical_decomposition_groups=len(decomposition),
        historical_frame_groups=len(counts), counterfactual_component_rows=768,
        counterfactual_groups=len(charge), recovery_contrasts=len(contrast),
        known_charge_energy_groups=len(energy), window_controls=len(window),
        figures=len(names), report_sha256=digest(OUTPUT/'report.json'),
        scope='Truth-assisted diagnostic only; all unmatched and extra counts retained'))
    print('HISTORY_REPORT_COMPLETE', names, flush=True)


if __name__ == '__main__':
    main()
