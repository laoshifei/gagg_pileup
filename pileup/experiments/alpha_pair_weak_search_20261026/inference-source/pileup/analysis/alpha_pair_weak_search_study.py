"""Frozen development study: bounded profiles and weak-charge search starts."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from importlib.metadata import version
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

from .alpha_pair_resolution_data import load, save, digest

PROJECT = Path(os.environ.get('GAGG_SEARCH_PROJECT', Path(__file__).resolve().parents[2])).resolve()
NAME = 'alpha_pair_weak_search_20261026'
EXPERIMENT = PROJECT / 'pileup/experiments' / NAME
OUTPUT = PROJECT / 'jupyter/artifacts' / NAME
RECORD = PROJECT / 'jupyter/run_records' / NAME
PREVIOUS = PROJECT / 'pileup/experiments/alpha_pair_joint_confirmation_20261026'
OLDER = PROJECT / 'jupyter/artifacts/alpha_pair_joint_likelihood_20261026'
LATEST = PROJECT / 'jupyter/artifacts/alpha_pair_joint_confirmation_20261026'
VARIANTS = ('legacy_extra', 'weak_extra')


def rows(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def check_freeze():
    for relative, sha in load(EXPERIMENT / 'freeze.json')['files'].items():
        assert digest(EXPERIMENT / relative) == sha, relative


def prepare():
    from .alpha_pair_weak_search import SHARES, GAPS
    for filename in ('joint-response.npz', 'joint_conserved_nll-noise.npz'):
        assert digest(PREVIOUS / filename) == digest(PROJECT / 'pileup/experiments/alpha_pair_joint_likelihood_20261026' / filename)
    old_cases = load(OLDER / 'labels.cases.json'); latest_cases = load(LATEST / 'labels.cases.json')
    assert len(old_cases) == 1880 and len(latest_cases) == 1600
    assert all(c['role'] == 'calibration' and c['kind'] == 'single' for c in old_cases[:280])
    assert all(c['role'] == 'development' for c in old_cases[280:])
    cases = [dict(c, dataset_role='gate' if i < 280 else 'older', original_index=i)
             for i, c in enumerate(old_cases)]
    cases += [dict(c, dataset_role='latest', original_index=i) for i, c in enumerate(latest_cases)]
    arrays = {}
    for key, filename in [('voltage_mV', 'observations.npz'), ('charge_pe', 'labels.truth.npz'), ('times_ns', 'labels.truth.npz')]:
        with np.load(OLDER / filename) as f: a = f[key]
        with np.load(LATEST / filename) as f: b = f[key]
        arrays[key] = np.concatenate((a, b))
    old_predictions = {r['index']: r['candidates'] for r in rows(OLDER / 'predictions.joint_conserved_nll.jsonl')}
    latest_predictions = {r['index']: r['candidates'] for r in rows(LATEST / 'predictions.joint.jsonl')}
    predictions = [old_predictions[i] for i in range(1880)] + [latest_predictions[i] for i in range(1600)]
    provenance = load(PROJECT / 'jupyter/run_records/alpha_pair_joint_likelihood_20261026/data-provenance.json')
    training = set(provenance['fold_mothers'][0]); gates = set(provenance['fold_mothers'][1])
    older_mothers = set(provenance['validation_mothers'])
    latest_native = load(LATEST / 'native-records.json')['test_single']
    latest_mothers = {r['root'] for r in latest_native}
    sets = [training, gates, older_mothers, latest_mothers]
    assert [len(s) for s in sets] == [280, 280, 320, 320]
    assert all(not sets[i] & sets[j] for i in range(4) for j in range(i))
    old_native = load(PROJECT / 'jupyter/artifacts/alpha_pair_null_gate_20261023/test/native-records.json')['test_single']
    blocks = np.full(len(cases), -1, int)
    for i, c in enumerate(cases):
        if c['dataset_role'] == 'gate' or c['kind'] != 'pair': continue
        native = old_native if c['dataset_role'] == 'older' else latest_native
        pair_blocks = {native[j]['replicate'] // 2 for j in c['parents']}
        assert len(pair_blocks) == 1, (i, c)
        blocks[i] = pair_blocks.pop()
    diagnostic_ids = [i for i, c in enumerate(cases) if c['dataset_role'] == 'latest'
                      and c['kind'] == 'pair' and c['condition'] in ('high_weak', 'weak_high')]
    assert len(diagnostic_ids) == 240
    profile_ids = [i for gap in (8., 40., 200.) for i in [j for j in diagnostic_ids
                   if cases[j]['condition'] == 'high_weak' and cases[j]['gap_ns'] == gap][:2]]
    assert len(profile_ids) == 6
    dense_id = next(i for i in profile_ids if cases[i]['gap_ns'] == 40.)
    pilot_groups = {}
    for i, c in enumerate(cases):
        pilot_groups.setdefault((c['dataset_role'], c['kind'], c['condition'], c['gap_ns']), i)
    for folder in (EXPERIMENT, OUTPUT, RECORD): folder.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(PROJECT / 'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb', RECORD / 'notebook.before.ipynb')
    for filename in ('joint-response.npz', 'joint_conserved_nll-noise.npz', 'settings.json'):
        shutil.copyfile(PREVIOUS / filename, EXPERIMENT / filename)
    save(EXPERIMENT / 'baseline-gate.json', load(PREVIOUS / 'gate-calibrations.json')['joint_conserved_nll'])
    np.savez_compressed(OUTPUT / 'observations.npz', voltage_mV=arrays['voltage_mV'])
    np.savez_compressed(OUTPUT / 'labels.truth.npz', charge_pe=arrays['charge_pe'], times_ns=arrays['times_ns'])
    save(OUTPUT / 'labels.cases.json', cases)
    save(OUTPUT / 'mother-blocks.truth.json', blocks.tolist())
    with (OUTPUT / 'predictions.baseline.jsonl').open('w') as stream:
        for i, candidates in enumerate(predictions): stream.write(json.dumps(dict(index=i, candidates=candidates)) + '\n')
    protocol = dict(created=datetime.now().astimezone().isoformat(),
        scope='Pure alpha; fixed95/351ns; local unknown K<=2; 2ns x1024 full analog samples; sufficient dynamic range',
        data_roles='280 native gate mothers; two previously exposed 1600-frame sets, both development; no new independent test yet',
        unchanged='Conditional mean, area normalization, covariance, logdet, K0/K1 models and raw-presence guard',
        variants={k: list(v) for k, v in SHARES.items()}, gap_bank_ns=list(GAPS),
        search='Best share per gap/direction; best4 per direction; eight20-iteration screens and two70-iteration polishes; retain original double candidate',
        budget='Same extra nonlinear optimization limits for both variants; weak bank has more inexpensive initial scores; all timing is incremental to cached baseline',
        gate='Recalibrate each searched K2 bank on the same280 independent native singles; plus-one/tie upper-tail rank<=.025, fitted-Q edges1500/5000/12500pe',
        inference_frames=3480, pilot_indices=list(pilot_groups.values()),
        diagnostics=dict(indices=diagnostic_ids, profile_indices=profile_ids, dense_profile_index=dense_id,
            profiles='13 weak-time offsets x9 weak-Q factors, two strong nuisance starts; true-centred region and labels explicitly diagnostic',
            strong_time_bounds_ns=48., dense_check='One prespecified40ns frame; midpoint-refined25x17 grid',
            scalar_warning='Lower noisy-data cost than exact truth does not by itself prove bias; profiles are local and not Wilks confidence regions'),
        success='Both trueQ>=500pe; eachQ error<=20%; each t error<=min(10ns,gap/4); no extra output',
        selection=dict(control_extra_max=.05, single_success_min=.95, single_median_absolute_q_bias_max=.05,
            latest_high_weak_40_gain_min=.10, older_high_weak_40_must_not_decrease=True,
            high_weak_200_must_not_decrease_in_either_set=True,
            weak_high_8_max_drop=.05, weak_average_must_not_decrease=True, strong_average_max_drop=.02,
            tie_break='Latest high_weak40+200 successes, then latest weak mean, then older weak mean; fixed variant order'),
        bootstrap=dict(replicates=2000, seed=276990226,
            design='Per dataset, 40 mother blocks; first20/last20 stratified; common weights across conditions/gaps/methods; descriptive development intervals'),
        next='Only a candidate satisfying declared development constraints proceeds to frozen fresh-parent confirmation',
        literature='Reuse Fixsen2004 existing note, TES simulation and AI pending review; profiles and search bank are project diagnostics, not a result of that paper')
    save(EXPERIMENT / 'protocol.json', protocol)
    input_paths = [folder / filename for folder, names in [(OLDER, ['observations.npz', 'labels.truth.npz', 'labels.cases.json', 'predictions.joint_conserved_nll.jsonl']),
                    (LATEST, ['observations.npz', 'labels.truth.npz', 'labels.cases.json', 'predictions.joint.jsonl', 'native-records.json'])] for filename in names]
    save(RECORD / 'input-provenance.json', dict(files={str(p.relative_to(PROJECT)): digest(p) for p in input_paths},
        mother_groups={k: sorted(v) for k, v in zip(('training', 'gate', 'older', 'latest'), sets)}, mothers_disjoint=True,
        output_hashes={p.name: digest(p) for p in (OUTPUT / 'observations.npz', OUTPUT / 'labels.truth.npz', OUTPUT / 'labels.cases.json', OUTPUT / 'predictions.baseline.jsonl', OUTPUT / 'mother-blocks.truth.json')}))
    save(RECORD / 'environment.json', dict(python=sys.version, executable=sys.executable,
        dependencies={name: version(name) for name in ('numpy', 'scipy', 'matplotlib', 'nbformat', 'pandas')},
        max_workers=8, blas_threads=1,
        git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=PROJECT, text=True).strip(),
        uncommitted='New weak-search, profile-diagnostic and study/test modules; exact executed source hashes in freeze.json',
        no_physics_parameter_changes=True, no_new_roots=True))
    (EXPERIMENT / 'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for name in ('SEARCH','JOINT','CONFIRM','NULL','CONDITIONAL','RESOLUTION','PAIR','COMPONENT'):
    os.environ['GAGG_'+name+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_weak_search_study import main
    main()
''')
    for namespace in ('analysis', 'pileup/analysis'):
        destination = EXPERIMENT / 'inference-source' / namespace; destination.mkdir(parents=True)
        for p in (PROJECT / namespace).glob('*.py'): shutil.copyfile(p, destination / p.name)
    if (PROJECT / 'pileup/__init__.py').exists():
        shutil.copyfile(PROJECT / 'pileup/__init__.py', EXPERIMENT / 'inference-source/pileup/__init__.py')
    shutil.copytree(PREVIOUS / 'source', EXPERIMENT / 'source')
    files = [p for directory in ('inference-source', 'source') for p in (EXPERIMENT / directory).rglob('*') if p.is_file()]
    files += [EXPERIMENT / filename for filename in ('joint-response.npz', 'joint_conserved_nll-noise.npz', 'settings.json', 'baseline-gate.json', 'protocol.json', 'run.py')]
    save(EXPERIMENT / 'freeze.json', dict(created=datetime.now().astimezone().isoformat(), files={str(p.relative_to(EXPERIMENT)): digest(p) for p in files}))
    print('WEAK_SEARCH_PREPARED', len(cases), 'frames;', len(files), 'frozen files', flush=True)


def initialize(blind=True):
    global RESPONSE, NOISE
    from .alpha_pair_joint_likelihood import ConservingResponse, JointNoise
    RESPONSE = ConservingResponse.load(EXPERIMENT / 'joint-response.npz')
    NOISE = JointNoise.load(EXPERIMENT / 'joint_conserved_nll-noise.npz')
    if blind:
        from .alpha_pair_conditional_study import forbid_labels
        sys.addaudithook(forbid_labels)


def fit_one(task):
    from .alpha_pair_weak_search import extend_candidates
    from .alpha_pair_profile_diagnostics import cost_at
    i, voltage, baseline, variant = task; started = time.monotonic()
    result = extend_candidates(RESPONSE, NOISE, voltage, variant, baseline)
    for candidate in result['candidates'][1:]:
        replay = cost_at(RESPONSE, NOISE, voltage, candidate['times_ns'], candidate['charge_pe'])['cost']
        assert abs(replay - candidate['cost']) <= 1e-7 * max(1., abs(replay)), (i, replay, candidate['cost'])
    if len(baseline) == 3:
        assert result['candidates'][2]['cost'] <= baseline[2]['cost'] + 1e-10
    assert result['candidates'][:2] == baseline[:2]
    return dict(index=i, seconds=time.monotonic() - started, **result)


def infer(variant, workers=4, pilot=False):
    check_freeze(); assert 1 <= workers <= 8
    with np.load(OUTPUT / 'observations.npz') as f: voltage = f['voltage_mV']
    baseline = {r['index']: r['candidates'] for r in rows(OUTPUT / 'predictions.baseline.jsonl')}
    indices = load(EXPERIMENT / 'protocol.json')['pilot_indices'] if pilot else list(range(len(voltage)))
    path = OUTPUT / f'predictions.{variant}.jsonl'
    old = rows(path) if path.exists() else []; done = {r['index'] for r in old}; assert len(done) == len(old)
    started = time.monotonic()
    with path.open('a') as stream:
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn'), initializer=initialize) as pool:
            futures = [pool.submit(fit_one, (i, voltage[i], baseline[i], variant)) for i in indices if i not in done]
            for n, future in enumerate(as_completed(futures), 1):
                result = future.result(); stream.write(json.dumps(result, allow_nan=False) + '\n'); stream.flush()
                if n % 50 == 0: os.fsync(stream.fileno()); print('WEAK_SEARCH_INFER', variant, n, '/', len(futures), flush=True)
        os.fsync(stream.fileno())
    save(RECORD / (variant + ('.pilot' if pilot else '') + '-inference.json'), dict(passed=True,
        requested=len(indices), new_frames=len(futures), seconds=time.monotonic() - started, workers=workers,
        worker_labels_blocked=True, truth_inputs=[], waveform_sha256=digest(OUTPUT / 'observations.npz'),
        baseline_sha256=digest(OUTPUT / 'predictions.baseline.jsonl'), prediction_sha256=digest(path),
        freeze_sha256=digest(EXPERIMENT / 'freeze.json')))
    print('WEAK_SEARCH_COMPLETE', variant, len(indices), flush=True)


def diagnose_one(task):
    from .alpha_pair_profile_diagnostics import truth_diagnosis, profile_map
    i, voltage, times, charge, baseline, kind = task; started = time.monotonic()
    result = (truth_diagnosis(RESPONSE, NOISE, voltage, times, charge, baseline) if kind == 'parameters'
              else profile_map(RESPONSE, NOISE, voltage, times, charge, baseline, dense=kind == 'dense'))
    return dict(index=i, kind=kind, seconds=time.monotonic() - started, **result)


def diagnostics(workers=4, pilot=False):
    check_freeze(); protocol = load(EXPERIMENT / 'protocol.json')['diagnostics']
    with np.load(OUTPUT / 'observations.npz') as f: y = f['voltage_mV']
    with np.load(OUTPUT / 'labels.truth.npz') as f: q = f['charge_pe']; t = f['times_ns']
    baseline = {r['index']: r['candidates'] for r in rows(OUTPUT / 'predictions.baseline.jsonl')}
    tasks = [(i, 'parameters') for i in protocol['indices']]
    tasks += [(i, 'profile') for i in protocol['profile_indices']] + [(protocol['dense_profile_index'], 'dense')]
    if pilot: tasks = tasks[:1] + [(protocol['profile_indices'][0], 'profile')]
    path = OUTPUT / 'diagnostics.truth.jsonl'; previous = rows(path) if path.exists() else []
    done = {(r['index'], r['kind']) for r in previous}; assert len(done) == len(previous)
    started = time.monotonic()
    with path.open('a') as stream:
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn'), initializer=initialize, initargs=(False,)) as pool:
            futures = [pool.submit(diagnose_one, (i, y[i], t[i], q[i], baseline[i], kind)) for i, kind in tasks if (i, kind) not in done]
            for n, future in enumerate(as_completed(futures), 1):
                result = future.result(); stream.write(json.dumps(result, allow_nan=False) + '\n'); stream.flush()
                if n % 20 == 0 or result['kind'] != 'parameters':
                    print('WEAK_PROFILE_DIAGNOSTIC', n, '/', len(futures), result['kind'], result['index'], flush=True)
        os.fsync(stream.fileno())
    save(RECORD / ('diagnostic-pilot.json' if pilot else 'diagnostic-run.json'), dict(passed=True, records=len(tasks),
        new_records=len(futures), workers=workers, seconds=time.monotonic() - started, truth_used=True,
        outputs_sha256=digest(path), freeze_sha256=digest(EXPERIMENT / 'freeze.json')))
    print('WEAK_DIAGNOSTICS_COMPLETE', len(tasks), flush=True)


def select():
    from .alpha_pair_joint_likelihood import calibrate_gate, select as choose
    from .alpha_pair_conditional_study import evaluate_fit
    from .alpha_pair_null_report import summarize
    check_freeze(); cases = load(OUTPUT / 'labels.cases.json')
    with np.load(OUTPUT / 'labels.truth.npz') as f: q = f['charge_pe']; t = f['times_ns']
    reports = {}; gates = {}; all_summary = []; all_scores = []; predictions = {}
    for name in ('baseline', *VARIANTS):
        predictions[name] = {r['index']: r['candidates'] for r in rows(OUTPUT / f'predictions.{name}.jsonl')}
        assert len(predictions[name]) == len(cases) == 3480
        gate = calibrate_gate([predictions[name][i] for i in range(280)])
        if name == 'baseline': assert gate == load(EXPERIMENT / 'baseline-gate.json')
        gates[name] = gate; reports[name] = {}
        for dataset in ('older', 'latest'):
            groups = {}
            for i, case in enumerate(cases):
                if case['dataset_role'] != dataset: continue
                fit = choose(predictions[name][i], gate)
                row = dict(evaluate_fit(fit, t[i], q[i], case), method=name, dataset_role=dataset,
                           kind=case['kind'], condition=case['condition'], gap_ns=case['gap_ns'], index=i,
                           fit=fit, true_q=q[i].tolist(), true_t=t[i].tolist())
                all_scores.append(row); groups.setdefault((case['kind'], case['condition'], case['gap_ns']), []).append(row)
            summary = [dict(summarize(rr), dataset_role=dataset) for rr in groups.values()]; all_summary += summary
            controls = [r for r in summary if r['kind'] != 'pair']
            control_pass = all(r['extra_fraction'] <= .05 for r in controls) and all(
                r['efficiency'] >= .95 and abs(r['q_error_p16_p50_p84'][1]) <= .05 for r in controls if r['kind'] == 'single')
            reports[name][dataset] = dict(controls_pass=control_pass,
                weak_mean=float(np.mean([r['efficiency'] for r in summary if r['kind'] == 'pair' and 'weak' in r['condition']])),
                strong_mean=float(np.mean([r['efficiency'] for r in summary if r['kind'] == 'pair' and 'weak' not in r['condition']])),
                high_weak={str(int(r['gap_ns'])): r['efficiency'] for r in summary if r['condition'] == 'high_weak'},
                weak_high={str(int(r['gap_ns'])): r['efficiency'] for r in summary if r['condition'] == 'weak_high'},
                controls={r['condition']: dict(extra=r['extra_frames'], success=r['success'], frames=r['frames'],
                    median_q_bias=None if r['q_error_p16_p50_p84'] is None else r['q_error_p16_p50_p84'][1]) for r in controls})
            print('WEAK_SEARCH_RESULT', name, dataset, json.dumps(reports[name][dataset]), flush=True)
    eligible = []; checks = {}
    for name in VARIANTS:
        latest = reports[name]['latest']; base_latest = reports['baseline']['latest']
        checks[name] = dict(latest_high_weak_40_gain=latest['high_weak']['40'] >= base_latest['high_weak']['40'] + .10 - 1e-12)
        for dataset in ('older', 'latest'):
            result = reports[name][dataset]; baseline = reports['baseline'][dataset]
            checks[name].update({dataset + '_controls': result['controls_pass'],
                dataset + '_high_weak_40': result['high_weak']['40'] >= baseline['high_weak']['40'] - 1e-12,
                dataset + '_high_weak_200': result['high_weak']['200'] >= baseline['high_weak']['200'] - 1e-12,
                dataset + '_weak_high_8': result['weak_high']['8'] >= baseline['weak_high']['8'] - .05 - 1e-12,
                dataset + '_weak_mean': result['weak_mean'] >= baseline['weak_mean'] - 1e-12,
                dataset + '_strong_mean': result['strong_mean'] >= baseline['strong_mean'] - .02 - 1e-12})
        if all(checks[name].values()): eligible.append(name)
    chosen = max(eligible, key=lambda name: (reports[name]['latest']['high_weak']['40'] + reports[name]['latest']['high_weak']['200'],
                 reports[name]['latest']['weak_mean'], reports[name]['older']['weak_mean'])) if eligible else None
    save(EXPERIMENT / 'gate-calibrations.json', gates)
    result = dict(chosen=chosen, eligibility=checks, results=reports, summary=all_summary,
                  data_role='Both datasets are exposed development; not new independent performance', independent_test_started=False)
    save(EXPERIMENT / 'development-selection.json', result)
    save(OUTPUT / 'report.json', result)
    with (OUTPUT / 'scores.jsonl').open('w') as stream:
        for row in all_scores: stream.write(json.dumps(row, allow_nan=False) + '\n')
    print('WEAK_SEARCH_CHOSEN', chosen, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare', 'infer', 'diagnostics', 'select'))
    parser.add_argument('--variant', choices=VARIANTS, default='weak_extra')
    parser.add_argument('--workers', type=int, default=4); parser.add_argument('--pilot', action='store_true')
    args = parser.parse_args()
    if args.command == 'infer': infer(args.variant, args.workers, args.pilot)
    elif args.command == 'diagnostics': diagnostics(args.workers, args.pilot)
    else: globals()[args.command]()


if __name__ == '__main__': main()
