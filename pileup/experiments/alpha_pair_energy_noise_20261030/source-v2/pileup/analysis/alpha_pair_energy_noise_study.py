"""Independent-repeat conditional mean/noise diagnosis with frozen E response."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
from datetime import datetime
from importlib.metadata import version
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

from .alpha_pair_digitizer_optimized import simulate_avalanches
from .sipm_response import PhotonHits, SensorConfig, AvalancheResult, map_cells
from .sipm_readout import ReadoutConfig, synthesize_readout
from .alpha_pair_recovery_energy import RecoveryEnergyResponse
from .alpha_pair_energy_oracle import VoltageSeeds
from .alpha_pair_energy_noise import (diagonal_factor, fit_known_times,
    replace_conditional_mean, gaussian_waveforms, local_energy_projection)

PROJECT = Path(os.environ.get('GAGG_ENERGY_NOISE_PROJECT', Path(__file__).resolve().parents[2])).resolve()
NAME = 'alpha_pair_energy_noise_20261030'
EXPERIMENT = PROJECT/'pileup/experiments'/NAME
OUTPUT = PROJECT/'jupyter/artifacts'/NAME
RECORD = PROJECT/'jupyter/run_records'/NAME
OLD = PROJECT/'pileup/experiments/alpha_pair_energy_oracle_20261029'
HISTORY = PROJECT/'jupyter/artifacts/sipm_energy_history_20261029'
BLOCKS = [0, 10, 20, 30]
NREF, NEVAL = 16, 16


def load(p): return json.loads(Path(p).read_text())
def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p, obj):
    p = Path(p); p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(p.name+'.partial')
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    tmp.replace(p)


def prepare():
    # Permit resuming preparation only before an immutable freeze or any run.
    assert not (EXPERIMENT/'freeze.json').exists()
    assert not (OUTPUT/'cases').exists() and not (OUTPUT/'pilot').exists()
    for p in [EXPERIMENT, OUTPUT, RECORD]: p.mkdir(parents=True, exist_ok=True)
    inputs = OUTPUT/'inputs'; (inputs/'photons').mkdir(parents=True, exist_ok=True)
    for name in ['template.npz', 'energy-model.json', 'sensor-config.json', 'response-parameters.json']:
        shutil.copyfile(OLD/name, EXPERIMENT/name)
    assert load(EXPERIMENT/'sensor-config.json') == load(PROJECT/'pileup/experiments/sipm_energy_history_20261029/sensor-config.json')
    parents = load(HISTORY/'inputs/parents.json')['records']
    cases = load(HISTORY/'inputs/cases.json')
    with np.load(HISTORY/'inputs/labels.npz') as f: energies = f['energy_MeV'].copy()
    truth = {c['parents'][0]: float(energies[i, 0]) for i, c in enumerate(cases) if c['kind'] == 'single'}
    evidence = {r['parent']: r for r in load(PROJECT/'pileup/experiments/sipm_energy_history_20261029/photon-inputs.json')}
    index = []
    for parent in BLOCKS+[240+b for b in BLOCKS]:
        src = HISTORY/'inputs/photons'/f'{parent:04}.npz'
        assert digest(src) == evidence[parent]['sha256']
        shutil.copyfile(src, inputs/'photons'/src.name)
        shutil.copyfile(PROJECT/parents[parent]['metadata'], inputs/f'parent-{parent:04}-metadata.json')
        index.append(dict(parent=parent, edep_MeV=truth[parent], **parents[parent], photon_sha256=digest(src),
                          source_photon_file=str(src.relative_to(PROJECT))))
    save(inputs/'parents.json', index)
    tasks = []
    for block in BLOCKS:
        for condition, gap in [('weak', None), ('high', None), ('high_weak', 8.),
                               ('high_weak', 40.), ('high_weak', 200.), ('weak_high', 40.)]:
            ps = [block if label == 'weak' else block+240 for label in condition.split('_')]
            times = [128.] if len(ps) == 1 else [128., 128.+gap]
            tasks.append(dict(id=len(tasks), block=block, condition=condition, gap_ns=gap,
                parents=ps, roots=[parents[p]['root'] for p in ps], times_ns=times,
                energy_MeV=[truth[p] for p in ps], labels=condition.split('_'),
                seed_base=703000000+len(tasks)*1000, window_ns=4096., n_reference=NREF, n_evaluation=NEVAL))
    save(EXPERIMENT/'plan.json', tasks)
    try:
        status = subprocess.run(['git', '-c', 'core.fileMode=false', 'diff', '--name-only', '--', 'README.md', 'docs/literature-use.md'],
                                capture_output=True, text=True, timeout=20)
        git_paths = status.stdout.splitlines()
    except subprocess.TimeoutExpired:
        git_paths = ['README.md', 'docs/literature-use.md']
        save(RECORD/'git-status-limitation.json', dict(timeout=True, paths_from_prior_completed_scoped_check=git_paths,
             note='Shared-filesystem git scan timed out; actual analysis source is captured and hashed below.'))
    save(EXPERIMENT/'protocol.json', dict(created=datetime.now().astimezone().isoformat(),
        git_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=PROJECT, text=True).strip(),
        preexisting_document_diff_paths=git_paths, implementation_source_files='See source snapshot and SHA256 freeze',
        scope='Known K/t pure-alpha conditional diagnostic, exposed existing optical parents; no blind improvement or fresh optical confirmation',
        blocks=BLOCKS, optical_parents=8, conditions_per_block=6, reference_repeats=NREF, evaluation_repeats=NEVAL,
        total_redigitizations=len(tasks)*(NREF+NEVAL), gaussian_evaluation_draws=len(tasks)*NEVAL,
        sample_spacing_ns=2., samples=1024, first_active_photon_ns=128., full_digitization_ns=4096.,
        fixed_scintillation_taus_ns=[95., 351.], sensor='Frozen S13360-6050PE proxy:14400 cells,29ns recovery',
        noise='Shared PDE/XT/dark/recovery, then0.1mV independent electronics; photons fixed within each optical mother',
        rate='Controlled singles/pairs, no Poisson source rate; cannot infer MHz performance from these controls',
        mean_model='Frozen recovery-energy strength2.515148098, existing isolated energy calibration, no relearning',
        domains=['physical_repeats', 'independent_reference_mean_replaced', 'matched_gaussian'],
        mean_replacement='y_evaluation - mean(y_reference) + model(true E,t); oracle only; finite reference error remains and is shared',
        covariance='One fixed factor from voltage-only one-event fit to independent reference mean; reused across all domains/weights',
        weights='GLS versus identical diagonal marginal variances; removes correlations only, not an exhaustive test of weights',
        gaussian='Exact frozen covariance in fitted FIR domain; voltage prefix noise fixed to zero; not native photon statistics',
        fitting='Known K and fixed true t only; no E/Q truth in fitting API; original linear +5 charge-share starts,120eval each',
        endpoints='All injected components E±5/10/20%, paired success, signed bias,sigma68,RMSE, boundaries,convergence; actual Q>=500 counts also reported',
        count_metrics='K/t fixed: missed/extra/split/merge/PSD not estimated; never present fixed counts as blind detection performance',
        selection='Four evenly spaced old mother blocks, fixed before outcomes, no selection by previous success/failure',
        limitations=['Conditional redigitization omits new transport fluctuations','Only4 optical pairs, repeated draws not independent optical events',
                     'Reference-mean differences include fixed optical realization as well as response mismatch',
                     'Local linear projection is no bound or confidence calibration','No alpha/gamma or free-tau/rate/real-detector claim'],
        max_workers=2, blas_threads=1, pilot_cost_limit_seconds=1800., pilot_cost_safety_factor=2.,
        no_new_G4=True, no_new_ROOT=True, no_new_macros=True, no_ADC_rejection=True, default_changed=False))
    source = EXPERIMENT/'source/pileup/analysis'; source.mkdir(parents=True)
    (source.parent/'__init__.py').write_text(''); (source/'__init__.py').write_text('')
    names = ['alpha_pair_energy_noise.py', 'alpha_pair_energy_noise_study.py', 'test_alpha_pair_energy_noise.py',
             'test_alpha_pair_recovery_energy.py', 'alpha_pair_recovery_energy.py', 'alpha_pair_energy_oracle.py',
             'alpha_pair_model_sampler.py', 'alpha_pair_digitizer_optimized.py', 'sipm_response.py', 'sipm_readout.py']
    for name in names: shutil.copyfile(PROJECT/'pileup/analysis'/name, source/name)
    (EXPERIMENT/'run.py').write_text("from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\nos.environ['GAGG_ENERGY_NOISE_PROJECT']=str(f.parents[2])\nsys.path.insert(0,str(f/'source'))\nif __name__=='__main__':\n    from pileup.analysis.alpha_pair_energy_noise_study import main\n    main()\n")
    save(EXPERIMENT/'input-hashes.json', {str(p.relative_to(OUTPUT)): digest(p) for p in inputs.rglob('*') if p.is_file()})
    save(EXPERIMENT/'freeze.json', dict(files={str(p.relative_to(EXPERIMENT)): digest(p)
         for p in EXPERIMENT.rglob('*') if p.is_file()}, before_simulation_and_outcomes=True))
    for name in ['codex_0007_SiPM_Energy_Calibration.ipynb', 'Research_Plan_20260921.ipynb']:
        shutil.copyfile(PROJECT/'jupyter'/name, RECORD/(name+'.before'))
    save(RECORD/'environment.json', dict(python=sys.version, executable=sys.executable,
        packages={n: version(n) for n in ['numpy','scipy','matplotlib','nbformat','uproot']}, workers=2, blas_threads=1))
    print('PREPARED', len(tasks), 'conditions', len(tasks)*(NREF+NEVAL), 'redigitizations', flush=True)


def verify():
    for path, expected in load(EXPERIMENT/'freeze.json')['files'].items():
        assert digest(EXPERIMENT/path) == expected, path
    for path, expected in load(EXPERIMENT/'input-hashes.json').items():
        assert digest(OUTPUT/path) == expected, path


def initialize():
    global DONORS, SENSOR, READOUT, RESPONSE
    DONORS = {}
    for p in (OUTPUT/'inputs/photons').glob('*.npz'):
        with np.load(p) as f: DONORS[int(p.stem)] = {k: f[k].copy() for k in f.files}
    config = load(EXPERIMENT/'sensor-config.json')
    SENSOR, READOUT = SensorConfig(**config['sensor']), ReadoutConfig(**config['readout'])
    RESPONSE = RecoveryEnergyResponse.load(EXPERIMENT)


def digitize_case(task, seed):
    columns = {k: [] for k in ['time_ns','x_mm','y_mm','energy_eV']}; owners = []; dropped = []
    for j, (parent, first) in enumerate(zip(task['parents'], task['times_ns'])):
        donor = DONORS[parent]; t = donor['time_ns']+first
        keep = (t >= 0)&(t < task['window_ns']); dropped.append(int((~keep).sum()))
        for key in columns: columns[key].append((t if key == 'time_ns' else donor[key])[keep])
        owners.extend([j]*int(keep.sum()))
    arrays = {k: np.concatenate(v) for k, v in columns.items()}; owner = np.asarray(owners, int)
    order = np.argsort(arrays['time_ns'], kind='stable'); arrays = {k: v[order] for k, v in arrays.items()}; owner = owner[order]
    hits = PhotonHits(**arrays, window_ns=task['window_ns'])
    cells = map_cells(hits, SENSOR)
    for j, expected in enumerate(task['times_ns']):
        assert np.isclose(hits.time_ns[(owner == j)&(cells >= 0)].min(), expected, atol=1e-9)
    av = simulate_avalanches(hits, SENSOR, seed)
    out = synthesize_readout(av, READOUT, seed+1000000)
    np.testing.assert_allclose(out.arrays['time_ns'][:1024], RESPONSE.sample_grid, atol=1e-12)
    rec = av.records; assigned = np.full(len(rec), -1, int); direct = rec['kind'] == 0
    assigned[direct] = owner[rec['photon_index'][direct]]
    for j in np.flatnonzero(rec['parent_id'] >= 0):
        assert rec['parent_id'][j] < j
        assigned[j] = assigned[rec['parent_id'][j]]
    q, components = [], []; signal = np.zeros_like(out.arrays['signal_mV'])
    for j in range(-1, len(task['parents'])):
        block = rec[assigned == j]
        contribution = synthesize_readout(AvalancheResult(block, task['window_ns'], av.full_charge_C, np.array([]), {}),
                          replace(READOUT, noise_rms_mV=0., baseline_mV=0.), 0).arrays['signal_mV']
        signal += contribution
        if j >= 0:
            q.append(float(block['charge_C'].sum()/av.full_charge_C)); components.append(contribution[:1024])
    closure = float(np.max(abs(signal-out.arrays['signal_mV'])))
    assert closure < 1e-8
    return out.arrays['analog_mV'][:1024], np.asarray(q), np.asarray(components), dict(closure_mV=closure, dropped_photons=dropped,
        avalanche_seed=seed, readout_seed=seed+1000000)


def one(task, pilot=False):
    started = time.monotonic(); nr = 2 if pilot else task['n_reference']; ne = 2 if pilot else task['n_evaluation']
    seed = task['seed_base']+(9000000 if pilot else 0)
    native = [digitize_case(task, seed+i) for i in range(nr+ne)]
    physical = np.array([r[0] for r in native]); charges = np.array([r[1] for r in native]); components = np.array([r[2] for r in native])
    generated_seconds = time.monotonic()-started
    ref, evaluation = physical[:nr], physical[nr:]
    y = ref.mean(0)-RESPONSE.baseline-RESPONSE.background
    crossing = np.flatnonzero(y > max(.5, .03*y.max()))
    onset = max(0., float(RESPONSE.sample_grid[crossing[0]])-3.)
    initial = VoltageSeeds(RESPONSE).refine(y, [onset])
    factor = RESPONSE.fixed_metric(initial)
    mean = RESPONSE.components(task['times_ns'], task['energy_MeV'])[0].sum(0)+RESPONSE.baseline+RESPONSE.background
    domain = dict(physical=evaluation, mean_replaced=replace_conditional_mean(evaluation, ref, mean),
                  gaussian=gaussian_waveforms(RESPONSE, mean, factor, seed+2000000, ne))
    factors = dict(gls=factor, diagonal=diagonal_factor(factor))
    rows, mean_fits = [], []
    for weight, metric in factors.items():
        for label, voltage in [('model_noiseless', mean), ('physical_reference_mean', ref.mean(0))]:
            fit = fit_known_times(RESPONSE, voltage, task['times_ns'], metric, initial)
            mean_fits.append(dict(domain=label, weight=weight, fit=fit))
            if label == 'model_noiseless': np.testing.assert_allclose(fit['energy_MeV'], task['energy_MeV'], rtol=1e-5, atol=1e-7)
        for name, observations in domain.items():
            for repeat, voltage in enumerate(observations):
                fit = fit_known_times(RESPONSE, voltage, task['times_ns'], metric, initial)
                rows.append(dict(domain=name, weight=weight, repeat=repeat, fit=fit))
    projection = local_energy_projection(RESPONSE, task['times_ns'], task['energy_MeV'], factor, ref, evaluation)
    folder = OUTPUT/('pilot' if pilot else 'cases'); folder.mkdir(exist_ok=True)
    np.savez_compressed(folder/f"{task['id']:02}.npz", physical=physical, components=components, charges=charges,
        gaussian=domain['gaussian'], model_mean=mean, factor=factor, diagonal_factor=factors['diagonal'])
    result = dict(task=task, pilot=pilot, n_reference=nr, n_evaluation=ne, initial_from_reference_voltage=initial,
                  fits=rows, mean_fits=mean_fits, projection=projection, simulation_records=[r[3] for r in native],
                  generation_seconds=generated_seconds, fitting_seconds=time.monotonic()-started-generated_seconds,
                  seconds=time.monotonic()-started, output_sha256=digest(folder/f"{task['id']:02}.npz"))
    save(folder/f"{task['id']:02}.json", result)
    return dict(id=task['id'], seconds=result['seconds'], generation_seconds=generated_seconds,
                fitting_seconds=result['fitting_seconds'], n_reference=nr, n_evaluation=ne)


def pilot():
    verify(); initialize()
    tasks = load(EXPERIMENT/'plan.json'); ids = [0, 3]
    timings = [one(tasks[i], pilot=True) for i in ids]
    # A pair is conservatively used for every case; the full run includes cheaper singles.
    worst = max(timings, key=lambda r: r['seconds'])
    projected = 2.*len(tasks)/2*(worst['generation_seconds']*(NREF+NEVAL)/4
                                           +worst['fitting_seconds']*NEVAL/2)
    record = dict(passed=projected <= 1800., projected_seconds=projected, limit_seconds=1800.,
                  safety_factor=2., workers=2, results=timings, criterion='cost and numerical integrity only; no accuracy-based selection')
    save(RECORD/'pilot-validation.json', record); print('PILOT', json.dumps(record), flush=True)
    assert record['passed'], 'Keep pilot and revise a separately versioned execution plan based on cost.'


def prepare_execution():
    """Cost-only reduction; retain the original frozen source and failed gate."""
    verify(); assert not (EXPERIMENT/'execution-plan.json').exists()
    timings = load(RECORD/'pilot-validation.json')['results']
    worst = max(timings, key=lambda r: r['seconds'])
    tasks = [dict(t, n_evaluation=8) for t in load(EXPERIMENT/'plan.json')]
    projected = 2.*len(tasks)/2*(worst['generation_seconds']*24/4+worst['fitting_seconds']*8/2)
    assert projected <= 1800.
    save(EXPERIMENT/'execution-plan.json', tasks)
    save(RECORD/'execution-plan-validation.json', dict(passed=True, projected_seconds=projected,limit_seconds=1800.,
        original_reference=16,original_evaluation=16,actual_reference=16,actual_evaluation=8,
        simulations=576,evaluation_frames=192,case_selection_unchanged=True,
        reason='Only conservative pilot wall time; no outcome-based mother/condition selection',
        original_failed_gate_preserved=True))
    src = EXPERIMENT/'source-v2'
    shutil.copytree(EXPERIMENT/'source', src)
    shutil.copyfile(PROJECT/'pileup/analysis/alpha_pair_energy_noise_study.py',src/'pileup/analysis/alpha_pair_energy_noise_study.py')
    (EXPERIMENT/'run_v2.py').write_text("from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\nos.environ['GAGG_ENERGY_NOISE_PROJECT']=str(f.parents[2])\nsys.path.insert(0,str(f/'source-v2'))\nif __name__=='__main__':\n    from pileup.analysis.alpha_pair_energy_noise_study import main\n    main()\n")
    files=[p for p in src.rglob('*') if p.is_file()]+[EXPERIMENT/'execution-plan.json',EXPERIMENT/'run_v2.py']
    save(EXPERIMENT/'execution-freeze.json',dict(files={str(p.relative_to(EXPERIMENT)):digest(p) for p in files},
        algorithm_unchanged=True,n_reference=16,n_evaluation=8,full_evaluation_not_started=True))
    print('EXECUTION_PLAN',projected,'seconds;576 digitizations,192 evaluation frames',flush=True)


def run(workers):
    assert 1 <= workers <= 2; verify(); assert load(RECORD/'execution-plan-validation.json')['passed']
    for path,sha in load(EXPERIMENT/'execution-freeze.json')['files'].items():assert digest(EXPERIMENT/path)==sha,path
    tasks = load(EXPERIMENT/'execution-plan.json'); pending = []
    for task in tasks:
        p = OUTPUT/'cases'/f"{task['id']:02}.json"
        if p.exists():
            old = load(p); assert old['task'] == task and digest(p.with_suffix('.npz')) == old['output_sha256']
        else: pending.append(task)
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn'), initializer=initialize) as pool:
        futures = [pool.submit(one, task) for task in pending]
        for count, future in enumerate(as_completed(futures), 1):
            result = future.result(); print('CASE_COMPLETE', count, '/', len(pending), json.dumps(result), flush=True)
    paths = sorted((OUTPUT/'cases').glob('*'))
    assert len(paths) == 2*len(tasks)
    save(RECORD/'result-freeze.json', dict(files={str(p.relative_to(OUTPUT)): digest(p) for p in paths}))
    save(RECORD/'run-validation.json', dict(passed=True, cases=len(tasks), new_cases=len(pending), seconds=time.monotonic()-started,
        workers=workers, simulations=sum(t['n_reference']+t['n_evaluation'] for t in tasks),
        evaluation_draws=sum(t['n_evaluation'] for t in tasks),
        fit_rows=sum(t['n_evaluation'] for t in tasks)*3*2, known_count_and_times=True, new_ROOTs=0, default_method_changed=False))
    print('RUN_COMPLETE', flush=True)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('command', choices=['prepare','pilot','plan','run']); parser.add_argument('--workers',type=int,default=2)
    args = parser.parse_args()
    if args.command == 'prepare': prepare()
    elif args.command == 'pilot': pilot()
    elif args.command == 'plan': prepare_execution()
    else: run(args.workers)


if __name__ == '__main__': main()
