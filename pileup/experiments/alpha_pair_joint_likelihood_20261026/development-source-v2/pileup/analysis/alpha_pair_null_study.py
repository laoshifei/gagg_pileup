"""Physical weak-alpha controls and calibration of residual pulse births."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

from .alpha_pair_resolution_data import load, save, digest

PROJECT = Path(os.environ.get('GAGG_NULL_PROJECT', Path(__file__).resolve().parents[2])).resolve()
NAME = 'alpha_pair_null_gate_20261023'
EXPERIMENT = PROJECT / 'pileup/experiments' / NAME
OUTPUT = PROJECT / 'jupyter/artifacts' / NAME
RECORD = PROJECT / 'jupyter/run_records' / NAME
PARENT = PROJECT / 'pileup/experiments/alpha_pair_residual_20261023'
GAPS = [4., 8., 12., 20., 40., 80., 200.]
CONDITIONS = ['low_low', 'mid_mid', 'high_high', 'low_high', 'high_low',
              'low_weak', 'weak_low', 'high_weak', 'weak_high', 'weak_weak']


def prepare():
    for folder in (EXPERIMENT, OUTPUT, RECORD):
        folder.mkdir(parents=True, exist_ok=False)
    for name in ('codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb', 'Run_Catalog.ipynb'):
        shutil.copyfile(PROJECT/'jupyter'/name, RECORD/(name+'.before'))
    shutil.copyfile(PROJECT/'docs/literature-use.md', RECORD/'literature-use.before.md')
    for name in ('source', 'bin'):
        shutil.copytree(PARENT/name, EXPERIMENT/name)
    for name in ('settings.json', 'response.npz', 'calibration.json'):
        shutil.copyfile(PARENT/name, EXPERIMENT/name)
    physics = {str(p.relative_to(EXPERIMENT)): digest(p)
               for part in ('source', 'bin') for p in (EXPERIMENT/part).rglob('*') if p.is_file()}
    save(EXPERIMENT/'physics-freeze.json', dict(created=datetime.now().astimezone().isoformat(),
        files=physics, parent=str(PARENT.relative_to(PROJECT)), physics_unchanged=True))
    save(RECORD/'environment.json', dict(python=sys.version, executable=sys.executable,
        numpy=np.__version__, cpu_affinity=sorted(os.sched_getaffinity(0)),
        head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=PROJECT, text=True).strip()))
    save(EXPERIMENT/'development-plan.json', dict(
        question='Reduce false births without deleting true 500-1000pe second events',
        inherited_response_and_search=True, neural_networks_retrained=False,
        old_tests_exposed_for_development=True,
        pilot_source_energies_MeV=[1.5, 1.75, 2.], pilot_per_energy=4,
        pilot_seed_start=264000000,
        gate_scores=['relative_residual', 'relative_signal'], gate_alpha=.025,
        observable_single_fit_charge_edges_pe=[1500., 5000., 12500.],
        calibration='All single controls from two exposed parent campaigns plus new weak singles',
        source_energy_choice='Choose pilot energy with median native Q closest to750pe; no photon scaling',
        validation='Exposed old pairs plus newly simulated weak-pair development; independent test after freezing',
        controls=dict(extra_fraction_max=.05, single_success_min=.95, abs_median_charge_bias_max=.05),
        scope='Pure alpha, fixed95/351ns, local unknown K<=2, analog mV, sufficient dynamic range'))
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os, sys
folder = Path(__file__).resolve().parent
os.environ['GAGG_NULL_PROJECT'] = str(folder.parents[2])
source = folder / 'inference-source'
sys.path.insert(0, str(source if source.exists() else folder.parents[2]))
if __name__ == '__main__':
    from pileup.analysis.alpha_pair_null_study import main
    main()
''')
    print('PREPARED', EXPERIMENT, flush=True)


def check_physics():
    for name, sha in load(EXPERIMENT/'physics-freeze.json')['files'].items():
        assert digest(EXPERIMENT/name) == sha, name


def create_stage(stage, energies, repetitions, seed):
    folder = EXPERIMENT/stage
    assert not folder.exists(), 'Stage already exists; resume it instead'
    for part in ('macros', 'logs', 'bin'):
        (folder/part).mkdir(parents=True, exist_ok=True)
    shutil.copyfile(EXPERIMENT/'bin/gagg-brems', folder/'bin/gagg-brems')
    plans = [dict(label=f'{stage}_alpha_{energy:g}_{rep:03d}', particle='alpha', energy_MeV=energy,
        x_mm=0., y_mm=0., replicate=rep, role=stage, gamma_model='brems')
        for energy in energies for rep in range(repetitions)]
    for i, row in enumerate(plans):
        row['seed'] = seed+i
    new = {r['seed'] for r in plans}
    audit = []
    paths = list((PROJECT/'pileup/experiments').glob('*/manifest.json'))
    paths += list(EXPERIMENT.glob('*/manifest.json'))
    for path in sorted(paths):
        # Photon/model manifests can lack runs; they are harmless here.
        old = load(path)
        if not isinstance(old, dict):
            continue
        used = {r.get('seed') for key in ('planned_runs', 'runs') for r in old.get(key, [])}
        assert not new & used, str(path)
        if used:
            audit.append(dict(path=str(path.relative_to(PROJECT)), sha256=digest(path)))
    save(RECORD/stage/'seed-independence.json', dict(passed=True, seeds=sorted(new), checked=audit))
    previous = load(PARENT/'manifest.json')
    save(folder/'manifest.json', dict(created=datetime.now().astimezone().isoformat(),
        planned_runs=plans, runs=[], preflight_validated=False, max_workers=4, window_us=20.,
        environment=previous['environment'], binary_sha256=digest(folder/'bin/gagg-brems')))
    save(folder/'protocol.json', dict(stage=stage, energies_MeV=energies,
        repetitions=repetitions, source_seed_start=seed, snapshot='../source',
        inference_truth_inputs=[], physics_sha256=digest(EXPERIMENT/'physics-freeze.json')))
    print('STAGE_PREPARED', stage, len(plans), flush=True)


def adapt_data(stage):
    from . import alpha_pair_resolution_data as data
    data.PROJECT = PROJECT
    data.EXPERIMENT = EXPERIMENT/stage
    data.OUTPUT = OUTPUT/stage
    data.RECORD = RECORD/stage
    data.check_freeze = check_physics
    return data


def native(stage, workers=4):
    if stage == 'test':
        check_freeze()
        assert load(RECORD/'test/digitizer-independence.json')['passed']
    from . import alpha_gamma_regularized_campaign as campaign
    original = campaign.macro
    # All stages point to this new experiment's preserved physical snapshot.
    def actual_macro(folder, row, output):
        body = original(EXPERIMENT, row, output)
        save(RECORD/stage/'attempts'/(row['label']+'.json'), dict(plan=row,
            root=str(output.relative_to(PROJECT)), created=datetime.now().astimezone().isoformat(),
            boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
            macro_contents=body))
        return body
    campaign.macro = actual_macro
    try:
        adapt_data(stage).native(workers)
    finally:
        campaign.macro = original


def pilot_report():
    records = load(OUTPUT/'pilot/native-records.json')['test_single']
    extra = OUTPUT/'pilot-refinement/native-records.json'
    if extra.exists():
        records += load(extra)['test_single']
    rows = []
    for r in records:
        with np.load(PROJECT/r['truth']) as f:
            q = float(np.sum(f['charge_pe']))
        rows.append(dict(energy_MeV=r['energy_MeV'], charge_pe=q, seed=r['seed'], root=r['root']))
    summaries = [dict(energy_MeV=e, median_charge_pe=float(np.median([r['charge_pe'] for r in rows if r['energy_MeV']==e])))
                 for e in sorted({r['energy_MeV'] for r in rows})]
    chosen = min(summaries, key=lambda r: abs(r['median_charge_pe']-750.))
    destination = OUTPUT/'pilot-report.json'
    if destination.exists() and not (OUTPUT/'pilot-report-initial.json').exists():
        shutil.copyfile(destination, OUTPUT/'pilot-report-initial.json')
    save(destination, dict(rows=rows, summary=summaries, suggested=chosen,
        note='Source energy, not incident/deposited/visible energy; charge is actual avalanche-equivalent pe'))
    print('PILOT', summaries, 'SUGGESTED', chosen, flush=True)


def cache(stage):
    pointer = RECORD/stage/'local-cache.json'
    if pointer.exists():
        path = Path(load(pointer)['directory'])
        if path.exists():
            return path
    path = Path(tempfile.mkdtemp(prefix='gagg-null-'+stage+'-'))
    save(pointer, dict(directory=str(path)))
    return path


def create_development():
    choice = load(OUTPUT/'pilot-report.json')['suggested']
    assert 500. <= choice['median_charge_pe'] <= 1000., 'Pilot did not bracket the intended weak charge'
    create_stage('development', [choice['energy_MeV']], 80, 264010000)
    save(EXPERIMENT/'weak-energy.json', dict(**choice, pilot_sha256=digest(OUTPUT/'pilot-report.json'),
        chosen_at=datetime.now().astimezone().isoformat(), photon_scaling=False))


def initialize_physical(local, count):
    from . import alpha_pair_resolution_data as data
    from . import sipm_response
    from .alpha_pair_digitizer_optimized import simulate_avalanches
    local = Path(local)
    data._DONORS = []
    for i in range(count):
        with np.load(local/'photons'/f'{i:04}.npz') as f:
            data._DONORS.append({k: f[k] for k in f.files})
    data._CONFIG = load(local/'sensor-config.json')
    sipm_response.simulate_avalanches = simulate_avalanches


def physical_one(task):
    from .alpha_pair_resolution_data import joint_case
    i, case, output = task
    path = Path(output)/'physical-cases'/f'{i:05}.npz'
    if path.exists():
        with np.load(path) as f:
            return i, (f['y'], f['q'], f['t'], f['components'], json.loads(str(f['case'])))
    row = joint_case(case)
    temporary = path.with_suffix('.partial')
    with temporary.open('wb') as stream:
        np.savez_compressed(stream, y=row[0], q=row[1], t=row[2], components=row[3], case=json.dumps(row[4]))
    temporary.replace(path)
    return i, row


def physical(stage, workers=4):
    if stage == 'test': check_freeze()
    data = adapt_data(stage)
    output = OUTPUT/stage
    (output/'photons').mkdir(exist_ok=True)
    records = data.donor_cache()
    local = cache(stage)
    (local/'photons').mkdir(exist_ok=True)
    (output/'physical-cases').mkdir(exist_ok=True)
    photon_sources = [(output/'photons'/f'{i:04}.npz') for i in range(len(records))]
    if stage == 'development':
        parent_output = PROJECT/'jupyter/artifacts'/PARENT.name
        old = load(parent_output/'native-records.json')['test_single']
        records += old
        photon_sources += [(parent_output/'photons'/f'{i:04}.npz') for i in range(len(old))]
    for i, path in enumerate(photon_sources):
        target = local/'photons'/f'{i:04}.npz'
        if not target.exists(): shutil.copyfile(path, target)
        assert digest(target) == digest(path)
    save(output/'physical-parents.json', dict(records=records,
        photon_files=[dict(path=str(p.relative_to(PROJECT)), sha256=digest(p)) for p in photon_sources]))
    save(local/'sensor-config.json', load(PROJECT/records[0]['metadata'])['config'])
    weak_energy = load(EXPERIMENT/'weak-energy.json')['energy_MeV']
    energies = dict(weak=weak_energy, low=2.5, mid=5.486, high=8.)
    groups = {name: [i for i, r in enumerate(records) if r['energy_MeV'] == energy]
              for name, energy in energies.items()}
    assert all(len(g) == 80 for g in groups.values()), {k:len(v) for k,v in groups.items()}
    conditions = CONDITIONS if stage == 'test' else CONDITIONS[5:]
    base = 265500000 if stage == 'test' else 264500000
    cases = []
    for ci, condition in enumerate(conditions):
        first, second = condition.split('_')
        pairs = (list(zip(groups[first][::2], groups[second][1::2])) if first == second
                 else list(zip(groups[first][:40], groups[second][:40])))
        for gap in GAPS:
            for rep, parents in enumerate(pairs):
                seed = base+ci*1000+rep
                cases.append(dict(kind='pair', condition=condition, gap_ns=gap, replicate=rep,
                    parents=list(parents), seed=seed,
                    first_ns=float(np.random.default_rng(seed).uniform(96., 256.)),
                    cluster='parents-'+'-'.join(map(str, sorted(parents)))))
    if stage == 'test':
        for rep in range(80):
            cases.append(dict(kind='empty', condition='all', gap_ns=None, replicate=rep, parents=[],
                seed=base+50000+rep, first_ns=128., cluster=f'empty-{rep}'))
    # Source labels belong exclusively to generation and later scoring.
    save(output/'generation-plan.json', dict(cases=cases, digitizer_base=base,
        common_random_numbers_across_gaps=True, photon_scaling=False, direct_voltage_addition=False))
    rows = [None]*len(cases)
    initialize_physical(local, len(records))
    started = time.monotonic()
    i, row = physical_one((0, cases[0], output)); rows[i] = row
    # Check the optimized digitizer against the original scalar algorithm on
    # the first predetermined case, before processing the full sample.
    import importlib
    from . import sipm_response
    from .alpha_pair_resolution_data import joint_case
    from .alpha_pair_digitizer_optimized import simulate_avalanches
    importlib.reload(sipm_response)
    try:
        reference = joint_case(cases[0])
    finally:
        sipm_response.simulate_avalanches = simulate_avalanches
    for actual, expected in zip(row[:4], reference[:4]):
        np.testing.assert_array_equal(actual, expected)
    save(RECORD/stage/'physical-preflight.json', dict(seconds=time.monotonic()-started,
        closure_error_mV=row[4]['component_sum_max_error_mV'], scalar_digitizer_exact=True, passed=True))
    with ProcessPoolExecutor(max_workers=workers, initializer=initialize_physical,
            initargs=(local, len(records)), mp_context=mp.get_context('spawn')) as pool:
        futures = [pool.submit(physical_one, (i, case, output)) for i, case in enumerate(cases) if i]
        for n, future in enumerate(as_completed(futures), 2):
            i, row = future.result(); rows[i] = row
            if n % 100 == 0: print('PHYSICAL', stage, n, '/', len(cases), flush=True)
    native_count = 320 if stage == 'test' else 80
    for i, r in enumerate(records[:native_count]):
        with np.load(PROJECT/r['waveform']) as f:
            timegrid = f['time_ns']; y = f['analog_mV']
        with np.load(PROJECT/r['truth']) as f:
            q = float(f['charge_pe'][0]); first = float(f['first_active_ns'][0])
        pre = float(np.random.default_rng(base+80000+i).uniform(96., 256.))
        start = int(np.floor((5000.-pre-timegrid[0])/2.))
        origin = float(timegrid[start]-.25)
        name = next(name for name, e in energies.items() if r['energy_MeV'] == e)
        case = dict(kind='single', condition=name, gap_ns=None, replicate=r['replicate'], parents=[i],
            seed=r['seed'], cluster=f'parent-{i}', crop_first_index=start, origin_ns=origin,
            crop='Nominal source trigger plus independent jitter; not true onset')
        rows.append((y[start:start+1024], np.array([q, 0.]), np.array([first-origin, 0.]),
                     np.zeros((2, 1024)), case))
    np.savez_compressed(output/'observations.npz', voltage_mV=np.asarray([r[0] for r in rows], np.float32))
    np.savez_compressed(output/'labels.truth.npz', charge_pe=np.asarray([r[1] for r in rows]),
        times_ns=np.asarray([r[2] for r in rows]), components_mV=np.asarray([r[3] for r in rows], np.float32))
    save(output/'labels.cases.json', [r[4] for r in rows])
    save(RECORD/stage/'physical-validation.json', dict(passed=True, frames=len(rows),
        native_parents=native_count, total_parents=len(records), workers=workers,
        max_closure_error_mV=max(r[4].get('component_sum_max_error_mV', 0.) for r in rows),
        observation_sha256=digest(output/'observations.npz'), photon_scaling=False,
        description='Joint offline microcell digitization; singles are native voltages'))
    print('PHYSICAL_COMPLETE', stage, len(rows), flush=True)


_RESPONSE = None
_PILOT = 'observed'


def initialize_inference(response_path, pilot='observed'):
    from .alpha_pair_conditional import ConditionalResponse
    from .alpha_pair_conditional_study import forbid_labels
    global _RESPONSE, _PILOT
    sys.addaudithook(forbid_labels)
    _RESPONSE = ConditionalResponse.load(response_path)
    _PILOT = pilot


def infer_one(task):
    from .alpha_pair_null_gate import waveform_candidates
    i, voltage = task
    started = time.monotonic()
    candidates = waveform_candidates(_RESPONSE, voltage, pilot=_PILOT)
    return dict(index=i, candidates=candidates, seconds=time.monotonic()-started)


def infer(stage, workers=4, limit=None, variant='observed'):
    if stage == 'test': check_freeze()
    output = OUTPUT/stage; local = cache(stage)
    source = output/'observations.npz'
    target = local/'observations.npz'
    if not target.exists(): shutil.copyfile(source, target)
    assert digest(target) == digest(source)
    with np.load(target) as f:
        assert f.files == ['voltage_mV']; y = f['voltage_mV']
    filename = 'predictions.jsonl' if variant == 'observed' else 'predictions.model.jsonl'
    dest = local/filename; durable = output/filename
    if durable.exists() and not dest.exists(): shutil.copyfile(durable, dest)
    done = {}
    if dest.exists():
        for line in dest.read_text().splitlines():
            row = json.loads(line); assert row['index'] not in done; done[row['index']] = row
    count = len(y) if limit is None else min(limit, len(y))
    def checkpoint():
        temporary = durable.with_suffix('.partial')
        shutil.copyfile(dest, temporary); temporary.replace(durable)
    with ProcessPoolExecutor(max_workers=workers, initializer=initialize_inference,
            initargs=(EXPERIMENT/'response.npz', variant), mp_context=mp.get_context('spawn')) as pool:
        futures = [pool.submit(infer_one, (i, y[i])) for i in range(count) if i not in done]
        with dest.open('a') as stream:
            for future in as_completed(futures):
                row = future.result(); done[row['index']] = row
                stream.write(json.dumps(row, allow_nan=False)+'\n'); stream.flush()
                if len(done) % 50 == 0 or len(done) == count:
                    checkpoint(); print('INFER', stage, len(done), '/', count, flush=True)
    checkpoint()
    record_name = 'inference.json' if variant == 'observed' else 'model-inference.json'
    save(RECORD/stage/record_name, dict(passed=True, frames=len(done), requested=count, pilot=variant,
        waveform_sha256=digest(source), prediction_sha256=digest(durable),
        mean_seconds=float(np.mean([r['seconds'] for r in done.values()])),
        input_files=[str(source.relative_to(PROJECT)), str((EXPERIMENT/'response.npz').relative_to(PROJECT))],
        worker_label_access_blocked=True, truth_inputs=[], neural_inputs=[]))


def historical_rows():
    parent_output = PROJECT/'jupyter/artifacts'/PARENT.name
    rows = []
    for dataset in ('validation', 'test'):
        cases = load(parent_output/(dataset+'.cases.json'))
        with np.load(parent_output/(dataset+'.truth.npz')) as f:
            charge = f['charge_pe']; times = f['times_ns']
        predictions = sorted(map(json.loads, (parent_output/(dataset+'.predictions.jsonl')).read_text().splitlines()),
                             key=lambda r: r['index'])
        assert len(predictions) == len(cases)
        for i, (prediction, case) in enumerate(zip(predictions, cases)):
            assert prediction['index'] == i
            rows.append(dict(case=case, charge=charge[i], times=times[i], dataset='exposed_'+dataset,
                             index=i, candidates=prediction['policies']['rho90']['waveform_only']))
    return rows


def new_rows(stage, variant='observed'):
    output = OUTPUT/stage
    cases = load(output/'labels.cases.json')
    with np.load(output/'labels.truth.npz') as f:
        charge = f['charge_pe']; times = f['times_ns']
    filename = 'predictions.jsonl' if variant == 'observed' else 'predictions.model.jsonl'
    predictions = sorted(map(json.loads, (output/filename).read_text().splitlines()), key=lambda r:r['index'])
    assert len(predictions) == len(cases)
    assert [r['index'] for r in predictions] == list(range(len(cases)))
    return [dict(case=cases[i], charge=charge[i], times=times[i], dataset=stage,
                 index=i, candidates=predictions[i]['candidates']) for i in range(len(cases))]


def selection_summary(rows, fits):
    from .alpha_pair_conditional_study import evaluate_fit
    groups = {}; controls = {}; scores = []; excluded = {}
    for row, fit in zip(rows, fits):
        case = row['case']; score = evaluate_fit(fit, row['times'], row['charge'], case)
        scores.append(score)
        if case['kind'] == 'pair':
            key = (case['condition'], case['gap_ns'])
            groups.setdefault(key, [])
            if score['target'] == 2:
                groups[key].append(score['success'])
            else:
                excluded[key] = excluded.get(key, 0)+1
        else:
            controls.setdefault(case['condition'], []).append(score)
    control_results = {}
    for name, ss in controls.items():
        errors = [e for s in ss for e in s['q_errors']]
        control_results[name] = dict(n=len(ss), success=float(np.mean([s['success'] for s in ss])),
            extra_fraction=float(np.mean([s['extra'] > 0 for s in ss])),
            q_median_relative_error=float(np.median(errors)) if errors else None)
    strong = [np.mean(v) if v else 0. for (condition, gap), v in groups.items() if 'weak' not in condition]
    weak = [np.mean(v) if v else 0. for (condition, gap), v in groups.items() if 'weak' in condition]
    passed = all(r['extra_fraction'] <= .05 for r in control_results.values()) and all(
        r['success'] >= .95 and r['q_median_relative_error'] is not None
        and abs(r['q_median_relative_error']) <= .05 for k, r in control_results.items() if k != 'all')
    return dict(controls=control_results, eligible=passed,
        strong_pair_success=float(np.mean(strong)), weak_pair_success=float(np.mean(weak)),
        utility=float((np.mean(strong)+np.mean(weak))/2.),
        points=[dict(condition=k[0], gap_ns=k[1], n=len(v),
                     below_two_target_threshold=excluded.get(k, 0),
                     success=float(np.mean(v)) if v else None) for k,v in groups.items()])


def select_gates(variant='observed'):
    from .alpha_pair_null_gate import calibrate, select
    from .alpha_pair_residual import choose
    rows = historical_rows()+new_rows('development') if variant == 'observed' else new_rows('model-development', 'model')
    singles = [r for r in rows if r['case']['kind'] == 'single']
    plan = load(EXPERIMENT/'development-plan.json')
    calibration = {}; results = {}
    for mode in plan['gate_scores']:
        cal = calibrate([r['candidates'] for r in singles], edges=plan['observable_single_fit_charge_edges_pe'],
                        mode=mode, alpha=plan['gate_alpha'], require_birth_penalty=variant == 'observed')
        calibration[mode] = cal
        results[mode] = selection_summary(rows, [select(r['candidates'], cal) for r in rows])
    eligible = [name for name, result in results.items() if result['eligible']]
    if eligible:
        chosen = max(eligible, key=lambda name: (results[name]['utility'], name == 'relative_residual'))
    else:
        chosen = None
    baseline = selection_summary(rows, [choose(r['candidates'], .001) for r in rows])
    prefix = '' if variant == 'observed' else 'model-'
    save(EXPERIMENT/(prefix+'null-calibration.json'), calibration)
    save(EXPERIMENT/(prefix+'null-selection.json'), dict(chosen=chosen, results=results, baseline=baseline,
        calibration_rows=len(singles), calibration_sources=[dict(dataset=r['dataset'], index=r['index'],
            parents=r['case']['parents'], condition=r['case']['condition']) for r in singles],
        criterion='Pass every control and bias condition; maximize equal-weight strong/weak pair means; tie prefers residual ratio',
        birth_selection='Old penalty then rank' if variant=='observed' else 'Rank replaces K1-to-K2 penalty; zero/single detection unchanged',
        controls_in_selection_are_calibration_not_independent=True,
        gate_alpha=.025, inference_truth_inputs=[], new_test_not_started=True))
    print('GATE_SELECTION', chosen, flush=True)
    for name, result in results.items():
        print(name, {k:v for k,v in result.items() if k != 'points'}, flush=True)
    print('BASELINE', {k:v for k,v in baseline.items() if k != 'points'}, flush=True)


def prepare_model_development():
    stage = 'model-development'; output = OUTPUT/stage
    output.mkdir(exist_ok=False)
    parent_output = PROJECT/'jupyter/artifacts'/PARENT.name
    observations = []; cases = []; charge = []; times = []; predictions = []; provenance = []
    inputs = [('exposed_'+part, parent_output, part+'.waveforms.npz', part+'.truth.npz',
               part+'.cases.json', part+'.predictions.jsonl') for part in ('validation', 'test')]
    inputs.append(('development', OUTPUT/'development', 'observations.npz', 'labels.truth.npz',
                   'labels.cases.json', 'predictions.jsonl'))
    for dataset, folder, wave, truth, casefile, predictionfile in inputs:
        source_cases = load(folder/casefile)
        with np.load(folder/wave) as f: y = f['voltage_mV']
        with np.load(folder/truth) as f: q = f['charge_pe']; t = f['times_ns']
        pred = {r['index']:r for r in map(json.loads, (folder/predictionfile).read_text().splitlines())}
        assert len(pred) == len(source_cases)
        for i, case in enumerate(source_cases):
            if dataset != 'development' and case['kind'] == 'pair' and case['replicate'] >= 10:
                continue
            j = len(cases)
            observations.append(y[i]); charge.append(q[i]); times.append(t[i])
            cases.append(dict(case, origin_dataset=dataset, origin_index=i))
            candidates = pred[i]['candidates'] if dataset == 'development' else pred[i]['policies']['rho90']['waveform_only']
            predictions.append(dict(index=j, candidates=candidates))
        provenance.append(dict(dataset=dataset, files={str((folder/p).relative_to(PROJECT)):digest(folder/p)
            for p in (wave, truth, casefile, predictionfile)}))
    np.savez_compressed(output/'observations.npz', voltage_mV=np.asarray(observations, np.float32))
    np.savez_compressed(output/'labels.truth.npz', charge_pe=np.asarray(charge), times_ns=np.asarray(times))
    save(output/'labels.cases.json', cases)
    (output/'predictions.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in predictions))
    save(EXPERIMENT/'model-pilot-development.json', dict(created=datetime.now().astimezone().isoformat(),
        reason='Observed-pilot GLS gives about-10% weak-single charge bias; model pilot fixes the exploratory bias',
        pilot='Conditional-mean OLS single fit from observation; one fixed precision for every K trial',
        changed_detector_parameters=False, inherited_mean=True, rho=.9,
        frames=len(cases), old_pairs_per_point_per_campaign=10, sources=provenance,
        gate_candidates=['relative_residual', 'relative_signal'], alpha=.025,
        selection='All control constraints; equal weight for strong-pair and weak-pair mean success',
        no_new_test_started=True))
    print('MODEL_DEVELOPMENT_PREPARED', len(cases), flush=True)


def check_freeze():
    for name, sha in load(EXPERIMENT/'freeze.json')['files'].items():
        assert digest(EXPERIMENT/name) == sha, name


def freeze():
    assert not (EXPERIMENT/'freeze.json').exists()
    choice = load(EXPERIMENT/'model-null-selection.json')
    assert choice['chosen'] is not None, 'A candidate must pass development requirements before independent production'
    assert load(RECORD/'waveform-search-replay.json')['passed']
    assert load(RECORD/'development/physical-preflight.json')['scalar_digitizer_exact']
    check_physics()
    weak = load(EXPERIMENT/'weak-energy.json')['energy_MeV']
    save(EXPERIMENT/'test-plan.json', dict(created=datetime.now().astimezone().isoformat(),
        source_energies_MeV=[weak,2.5,5.486,8.], native_roots_per_energy=80, native_roots=320,
        source_seed_start=265000000, digitizer_seed_start=265500000,
        conditions=CONDITIONS, gaps_ns=GAPS, pairs_per_point=40, empty_frames=80,
        expected_frames=3200, fixed_taus_ns=[95.,351.], local_maximum_K=2,
        observation='Analog mV, 1024 samples at2ns; sufficient dynamic range',
        selected_method='model_'+choice['chosen'], calibration_alpha=.025,
        birth_selection='Empirical post-search rank replaces fixed birth penalty; unchanged zero/single detection',
        metrics='Both true Q>=500pe, each Q within20%, each t within min(10ns,gap/4), no extra output',
        below_threshold='Retain all generated frames; separately report how many do not have two Q>=500pe events',
        controls=dict(extra_fraction_max=.05,single_success_min=.95,abs_median_charge_bias_max=.05),
        time_definition='First photon inside active SiPM plane; not source emission or pulse maximum',
        examples_predeclared=[dict(condition='high_weak',replicate=0,gap_ns=gap) for gap in (40.,200.)],
        inference_truth_inputs=[], models_retrained=False, new_test_generation_started=False,
        independence='Fresh optical parents; reused across gaps/conditions; no pooled independent exposure claim',
        scope='Pure central-incidence alpha and fixed lifetimes; no alpha/gamma or random high-rate claim'))
    for namespace in ('analysis','pileup/analysis'):
        dest = EXPERIMENT/'inference-source'/namespace; dest.mkdir(parents=True)
        for path in (PROJECT/namespace).glob('*.py'): shutil.copyfile(path,dest/path.name)
    (EXPERIMENT/'inference-source/pileup/__init__.py').write_text('')
    files = dict(load(EXPERIMENT/'physics-freeze.json')['files'])
    files.update({str(p.relative_to(EXPERIMENT)):digest(p) for p in (EXPERIMENT/'inference-source').rglob('*.py')})
    for name in ('physics-freeze.json','settings.json','response.npz','calibration.json','development-plan.json',
                 'weak-energy.json','model-pilot-development.json','model-null-calibration.json','model-null-selection.json',
                 'test-plan.json','run.py','report.py','production.py'):
        files[name] = digest(EXPERIMENT/name)
    save(EXPERIMENT/'freeze.json', dict(frozen_at=datetime.now().astimezone().isoformat(), files=files,
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        development_observations_sha256=digest(OUTPUT/'model-development/observations.npz'),
        development_predictions_sha256=digest(OUTPUT/'model-development/predictions.model.jsonl'),
        new_test_generation_started=False, inference_truth_inputs=[]))
    (RECORD/'tracked-status-at-freeze.txt').write_text(subprocess.check_output(
        ['git','--no-optional-locks','-c','core.fileMode=false','status','--short','--untracked-files=no'],cwd=PROJECT,text=True))
    save(RECORD/'source-snapshot-note.json',dict(
        note='All reusable analysis Python and physical source are copied, including uncommitted files; do not infer source solely from Git HEAD',
        frozen_file_count=len(files)))
    print('FROZEN',len(files),flush=True)


def create_test():
    check_freeze(); p=load(EXPERIMENT/'test-plan.json')
    create_stage('test',p['source_energies_MeV'],p['native_roots_per_energy'],p['source_seed_start'])
    # Check every new avalanche/readout stream against earlier G4 streams.
    source = {p['source_seed_start']+i for i in range(p['native_roots'])}
    joint = {p['digitizer_seed_start']+ci*1000+rep for ci in range(len(CONDITIONS)) for rep in range(40)}
    joint.update(p['digitizer_seed_start']+50000+rep for rep in range(80))
    new_streams = {s+offset for s in source|joint for offset in (0,100000,200000)}
    compared = []
    paths = list((PROJECT/'pileup/experiments').glob('*/manifest.json'))+list(EXPERIMENT.glob('*/manifest.json'))
    for path in sorted(paths):
        if path == EXPERIMENT/'test/manifest.json': continue
        old=load(path); bases={r['seed'] for key in ('runs','planned_runs') for r in old.get(key,[]) if 'seed' in r}
        assert not new_streams & {s+offset for s in bases for offset in (0,100000,200000)}, str(path)
        compared.append(str(path.relative_to(PROJECT)))
    case_paths = list((PROJECT/'jupyter/artifacts').glob('*/*.cases.json'))
    case_paths += [OUTPUT/'development/labels.cases.json']
    for path in case_paths:
        old=load(path)
        if isinstance(old,dict): old=old.get('cases',[])
        bases={r['seed'] for r in old if 'seed' in r}
        assert not new_streams & {s+offset for s in bases for offset in (0,100000,200000)},str(path)
    save(RECORD/'test/digitizer-independence.json',dict(passed=True,new_streams=len(new_streams),
        compared_source_manifests=compared,compared_cases=[str(p.relative_to(PROJECT)) for p in case_paths],
        cross_gap_common_random_numbers=True))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('prepare', 'create-pilot', 'create-pilot-refinement', 'native', 'pilot-report',
                                    'create-development', 'physical', 'infer', 'select', 'prepare-model-development',
                                    'freeze','create-test'))
    p.add_argument('--stage', default='pilot')
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--limit', type=int)
    p.add_argument('--variant', choices=('observed','model'), default='observed')
    a = p.parse_args()
    if a.action == 'prepare': prepare()
    elif a.action == 'create-pilot': create_stage('pilot', [1.5, 1.75, 2.], 4, 264000000)
    elif a.action == 'create-pilot-refinement':
        save(EXPERIMENT/'pilot-refinement-plan.json', dict(source_energy_MeV=1.625,
            reason='Initial1.5MeV median422.8pe and1.75MeV median1034.1pe bracket750pe; add midpoint',
            initial_report_sha256=digest(OUTPUT/'pilot-report.json')))
        create_stage('pilot-refinement', [1.625], 4, 264000100)
    elif a.action == 'native': native(a.stage, a.workers)
    elif a.action == 'pilot-report': pilot_report()
    elif a.action == 'create-development': create_development()
    elif a.action == 'physical': physical(a.stage, a.workers)
    elif a.action == 'infer': infer(a.stage, a.workers, a.limit, a.variant)
    elif a.action == 'select': select_gates(a.variant)
    elif a.action == 'prepare-model-development': prepare_model_development()
    elif a.action == 'freeze': freeze()
    elif a.action == 'create-test': create_test()


if __name__ == '__main__':
    main()
