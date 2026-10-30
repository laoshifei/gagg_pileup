"""Frozen development study separating energy, timing, search and count errors."""
from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
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

from .alpha_pair_recovery_energy import RecoveryEnergyResponse
from .alpha_pair_energy_oracle import blind_candidates, oracle_diagnostics, calibrate_gate, select
from .sipm_energy_calibration import score_events, error_summary

PROJECT = Path(os.environ.get('GAGG_ENERGY_ORACLE_PROJECT', Path(__file__).resolve().parents[2])).resolve()
NAME = 'alpha_pair_energy_oracle_20261029'
EXPERIMENT = PROJECT/'pileup/experiments'/NAME
OUTPUT = PROJECT/'jupyter/artifacts'/NAME
RECORD = PROJECT/'jupyter/run_records'/NAME
PREVIOUS = PROJECT/'pileup/experiments/alpha_pair_recovery_energy_20261029'
METHODS = ['blind_energy_gate', 'true_count_blind_start', 'fixed_true_time', 'free_from_true_time']
GAPS = [8., 40., 200.]


def load(p): return json.loads(Path(p).read_text())
def lines(p): return [json.loads(s) for s in Path(p).read_text().splitlines()]
def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p, value):
    p = Path(p); p.parent.mkdir(parents=True, exist_ok=True)
    temp = p.with_name(p.name+'.partial')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n'); temp.replace(p)
def write_lines(p, rows):
    with Path(p).open('w') as stream:
        for row in rows: stream.write(json.dumps(row, ensure_ascii=False, allow_nan=False)+'\n')
def freeze(name, paths, **metadata):
    save(EXPERIMENT/name, dict(created=datetime.now().astimezone().isoformat(),
        files={str(p.relative_to(EXPERIMENT)): digest(p) for p in paths}, **metadata))
def verify(name='freeze.json'):
    for path, sha in load(EXPERIMENT/name)['files'].items(): assert digest(EXPERIMENT/path) == sha, path
def verify_predictions(stage):
    for path, sha in load(RECORD/(stage+'-freeze.json'))['files'].items(): assert digest(OUTPUT/path) == sha, path


def prepare():
    assert not any(p.exists() for p in [EXPERIMENT, OUTPUT, RECORD])
    for name in ['freeze.json', 'method-freeze.json']:
        for path, sha in load(PREVIOUS/name)['files'].items(): assert digest(PREVIOUS/path) == sha, path
    for p in [EXPERIMENT, OUTPUT, RECORD]: p.mkdir(parents=True)
    pool = PROJECT/'jupyter/artifacts/alpha_pair_joint_likelihood_20261026'
    energy_path = PROJECT/'jupyter/artifacts/alpha_pair_recovery_energy_20261029/inputs/development.energy.json'
    cases = load(pool/'calibration.cases.json')
    indices = [i for i, c in enumerate(cases) if c['kind'] != 'pair' or c['gap_ns'] in GAPS]
    assert len(indices) == 2520
    selected = [dict(cases[i], source_index_in_pool=i) for i in indices]
    energy = {r['root']: r for r in load(energy_path)}
    with np.load(pool/'calibration.observations.npz') as f: y = f['voltage_mV'][indices]
    with np.load(pool/'calibration-ancestry.truth.npz') as f:
        q = f['charge_pe'][indices]; t = f['times_ns'][indices]; comp = f['components_mV'][indices]
    e = np.zeros_like(q); oracle_times = []
    for j, c in enumerate(selected):
        for k, root in enumerate(c['mother_roots']): e[j, k] = energy[root]['edep_MeV']
        oracle_times.append(t[j, q[j] > 0].tolist())
        assert len(oracle_times[-1]) == len(c['mother_roots'])
    inputs = OUTPUT/'inputs'; inputs.mkdir()
    np.savez_compressed(inputs/'observations.npz', voltage_mV=y)
    np.savez_compressed(inputs/'evaluation.labels.npz', charge_pe=q, times_ns=t, energy_MeV=e, components_mV=comp)
    save(inputs/'cases.json', selected); save(inputs/'oracle-times.truth.json', oracle_times)
    roots = [{p for c in selected if c['fold'] == f for p in c['mother_roots']} for f in [0, 1]]
    assert list(map(len, roots)) == [280, 280] and not roots[0]&roots[1]
    for name in ['template.npz', 'energy-model.json', 'sensor-config.json', 'response-parameters.json']:
        shutil.copyfile(PREVIOUS/name, EXPERIMENT/name)
    pilot = []
    for kind, condition, gap in [('empty', 'all', None), ('single', 'weak', None), ('single', 'high', None),
                                  ('pair', 'high_weak', 8.), ('pair', 'weak_high', 40.), ('pair', 'high_high', 8.)]:
        pilot.append(next(i for i, c in enumerate(selected) if c['fold'] == 0 and
            (c['kind'], c['condition'], c['gap_ns']) == (kind, condition, gap)))
    source_paths = [pool/'calibration.observations.npz', pool/'calibration-ancestry.truth.npz', pool/'calibration.cases.json', energy_path]
    protocol = dict(created=datetime.now().astimezone().isoformat(), git_baseline=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        purpose='Separate E-response, E/t coupling, finite search and count selection failures',
        scope='Existing development optical mothers only; no new transport or independent confirmation',
        data='2520 frames: each fold 900 pairs,280 singles,80 empties. Pair gaps8/40/200ns; weak groups20 pairs/cell/fold, other groups40.',
        train_mothers=280, validation_mothers=280, fixed_taus_ns=[95., 351.], sample_spacing_ns=2., observed_samples=1024,
        full_analog_range=True, new_ROOTs=0, particles='Pure central alpha; local K<=2; source E is not deposited E',
        model='Exactly inherited frozen recovery mean, energy calibration,29ns,14400cells,strength2.515148; no re-learning',
        time_convention='True first-active times in existing aligned records; inherited template support starts -0.25ns relative to nominal t, retained explicitly',
        common_weight='One frame-fixed unconditioned banded weight from voltage-only variable-projection single-pulse Q/t. Same for every K and oracle; not old NLL weighting.',
        objective='Whitened residual sum of squares; triangular solve equals r^T C^-1 r. No new likelihood/CI calibration.',
        seed_search='Adapt existing PairFitter variable-projection voltage starts; top4 recovery gap/share starts plus raw two-pulse start;24eval screening,best2 at120eval',
        bounds='Every free time spans observed window; logE constrained to frozen calibration domain; true E is never an initializer or bound',
        fixed_true_time='Known physical count and fixed true times; linear voltage E start plus5 share starts,120eval each',
        true_count_blind_start='True count only selects previously frozen voltage-only candidate; no true t/Q/E in candidate search',
        free_from_true_time='Extra diagnostic: release true times after fixed-time E fit; true time used as initializer, full-window bounds',
        blind_energy_gate='Select0/1/2 from same frozen candidates. Peak>1mV trigger; predictedQ>=500pe; positive score improvement; single-event plus-one rank<=.025 in Q bins1500/5000/12500.',
        gate_training='280 fold0 singles only; freeze before fold1 scoring. Response used same fold0 mothers, so empirical development gate, not a coverage guarantee.',
        endpoints='All truthQ>=500pe; E±5/10/20%, strict t<=min(10ns,gap/4), no extras for frame success; bias,sigma68,misses,extra,temporal split/merge',
        comparison='Only within this common-data/common-weight experiment. Previous38% historical score uses different mothers,weights,count rules and is not a paired baseline.',
        diagnostic_not_bound='Truth-assisted scores are neither blind efficiencies nor rigorous information bounds',
        pilot_indices=pilot, max_workers=2, cost_gate='2*mean pilot per-frame wall time*2520/2 <=3600seconds for combined stages',
        followup_gate='No fresh native confirmation unless fold1 weak results and single/empty controls support a frozen blind method',
        declared_before_outputs=True, statistical_seed=202610297301)
    save(EXPERIMENT/'protocol.json', protocol)
    save(EXPERIMENT/'input-hashes.json', {str(p.relative_to(OUTPUT)): digest(p) for p in inputs.iterdir()})
    save(RECORD/'provenance.json', dict(sources={str(p.relative_to(PROJECT)): digest(p) for p in source_paths},
        source_indices=indices, fold0_roots=sorted(roots[0]), fold1_roots=sorted(roots[1]), mother_disjoint=True,
        inherited_model={name: digest(PREVIOUS/name) for name in ['template.npz', 'energy-model.json', 'sensor-config.json', 'response-parameters.json']},
        previous_uncommitted_files_retained=True))
    save(RECORD/'environment.json', dict(python=sys.version, executable=sys.executable,
        packages={name: version(name) for name in ['numpy', 'scipy', 'matplotlib', 'nbformat']}, workers=2, blas_threads=1))
    source = EXPERIMENT/'source/pileup/analysis'; source.mkdir(parents=True)
    (source.parent/'__init__.py').write_text(''); (source/'__init__.py').write_text('')
    for name in ['alpha_pair_recovery_energy.py', 'sipm_energy_calibration.py']:
        shutil.copyfile(PREVIOUS/'source/pileup/analysis'/name, source/name)
    for name in ['alpha_pair_energy_oracle.py', 'alpha_pair_energy_oracle_study.py',
                 'test_alpha_pair_energy_oracle.py', 'test_alpha_pair_recovery_energy.py']:
        shutil.copyfile(PROJECT/'pileup/analysis'/name, source/name)
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
os.environ['GAGG_ENERGY_ORACLE_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_energy_oracle_study import main
    main()
''')
    freeze('freeze.json', [p for p in EXPERIMENT.rglob('*') if p.is_file()], inference_started=False)
    for name in ['codex_0007_SiPM_Energy_Calibration.ipynb', 'Research_Plan_20260921.ipynb']:
        shutil.copyfile(PROJECT/'jupyter'/name, RECORD/(name+'.before'))
    print('ORACLE_PREPARED', len(indices), flush=True)


def initialize():
    global RESPONSE
    RESPONSE = RecoveryEnergyResponse.load(EXPERIMENT)
    def guard(event, args):
        if event == 'open' and args and any(w in str(args[0]).lower() for w in ['.truth', '.labels', '.root', '/truth/', '/labels/']):
            raise PermissionError('Worker may not read truth files')
    sys.addaudithook(guard)
    for path in ['forbidden.truth.json', 'forbidden.labels.npz', 'forbidden.root']:
        try: open(path)
        except PermissionError: pass
        else: raise AssertionError('Truth-read guard failed')


def one(task):
    stage, index, voltage, *extra = task; start = time.monotonic()
    result = (blind_candidates(RESPONSE, voltage) if stage == 'blind' else
              oracle_diagnostics(RESPONSE, voltage, extra[0], extra[1]))
    return dict(index=index, seconds=time.monotonic()-start, **result)


def run_stage(stage, workers=2, indices=None):
    verify(); assert 1 <= workers <= 2
    with np.load(OUTPUT/'inputs/observations.npz') as f: y = f['voltage_mV']
    if indices is None: indices = list(range(len(y)))
    folder = OUTPUT/stage; folder.mkdir(exist_ok=True)
    pending = [i for i in indices if not (folder/f'{i:04}.json').exists()]
    truth_times = load(OUTPUT/'inputs/oracle-times.truth.json') if stage == 'oracle' else None
    tasks = [(stage, i, y[i], *([load(OUTPUT/'blind'/f'{i:04}.json'), truth_times[i]] if stage == 'oracle' else [])) for i in pending]
    start = time.monotonic()
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn'), initializer=initialize) as pool:
        futures = [pool.submit(one, task) for task in tasks]
        for count, f in enumerate(as_completed(futures), 1):
            r = f.result(); save(folder/f"{r['index']:04}.json", r)
            if count % 60 == 0: print('ORACLE_PROGRESS', stage, count, '/', len(pending), flush=True)
    return dict(stage=stage, new_frames=len(pending), seconds=time.monotonic()-start, workers=workers)


def pilot(workers=2):
    ids = load(EXPERIMENT/'protocol.json')['pilot_indices']
    a = run_stage('blind', workers, ids); b = run_stage('oracle', workers, ids)
    seconds = [sum(load(OUTPUT/stage/f'{i:04}.json')['seconds'] for stage in ['blind', 'oracle']) for i in ids]
    projected = float(2*np.mean(seconds)*2520/2)
    save(RECORD/'pilot-validation.json', dict(passed=projected <= 3600, indices=ids, per_frame_seconds=seconds,
        projected_total_two_worker_seconds=projected, cost_safety_factor=2., workers=workers,
        blind_truth_file_guard=True, oracle_truth_input='Explicit known count/times only; no true energies or charges', stages=[a, b]))
    print('ORACLE_PILOT', projected, projected <= 3600, flush=True)


def infer(stage, workers=2):
    verify('execution-freeze.json')
    assert load(RECORD/'execution-plan-validation.json')['passed']
    if stage == 'oracle': verify_predictions('blind')
    indices = load(EXPERIMENT/'execution-plan.json')['indices']
    result = run_stage(stage, workers, indices)
    files = [OUTPUT/stage/f'{i:04}.json' for i in indices]; assert len(files) == 1540 and all(p.exists() for p in files)
    save(RECORD/(stage+'-freeze.json'), dict(files={str(p.relative_to(OUTPUT)): digest(p) for p in files},
        stage=stage, true_E_never_supplied=True, true_t_supplied=(stage == 'oracle'), count=len(files)))
    save(RECORD/(stage+'-validation.json'), dict(passed=True, **result, frames=len(files)))
    print('ORACLE_STAGE_COMPLETE', stage, result, flush=True)


def gate():
    verify(); verify('execution-freeze.json'); verify_predictions('blind')
    assert not (EXPERIMENT/'gate.json').exists()
    cases = load(OUTPUT/'inputs/cases.json')
    ids = [i for i, c in enumerate(cases) if c['fold'] == 0 and c['kind'] == 'single']
    assert len(ids) == 280
    g = calibrate_gate([load(OUTPUT/'blind'/f'{i:04}.json')['candidates'] for i in ids])
    save(EXPERIMENT/'gate.json', g)
    freeze('selection-freeze.json', [EXPERIMENT/'freeze.json', EXPERIMENT/'execution-freeze.json', EXPERIMENT/'gate.json'],
        fold0_single_indices=ids, evaluation_not_started=True)
    print('ORACLE_GATE_FROZEN', list(map(len, g['scores'])), flush=True)


def evaluate():
    verify(); verify('execution-freeze.json'); verify('selection-freeze.json'); verify_predictions('blind'); verify_predictions('oracle')
    cases = load(OUTPUT/'inputs/cases.json'); g = load(EXPERIMENT/'gate.json')
    included = set(load(EXPERIMENT/'execution-plan.json')['indices'])
    with np.load(OUTPUT/'inputs/evaluation.labels.npz') as f: t, q, e = f['times_ns'], f['charge_pe'], f['energy_MeV']
    scores, predictions, search = [], [], []
    for i, c in enumerate(cases):
        if i not in included: continue
        blind = load(OUTPUT/'blind'/f'{i:04}.json'); oracle = load(OUTPUT/'oracle'/f'{i:04}.json')
        outputs = dict(blind_energy_gate=select(blind, g), **{m: oracle[m] for m in METHODS[1:]})
        for method, p in outputs.items():
            score = score_events(t[i], q[i], e[i], p['times_ns'], p['energy_MeV'], c['gap_ns'])
            tt = t[i, (q[i] >= 500.) & (e[i] > 0)]
            close = np.abs(tt[:, None]-np.asarray(p['times_ns'])[None, :]) <= 10.
            used_t = {m['target_index'] for m in score['matches']}; used_p = {m['prediction_index'] for m in score['matches']}
            merge = len({j for a in range(len(tt)) if a not in used_t for j in used_p if close[a, j]})
            split = len({a for a in used_t if any(close[a, j] for j in range(len(p['times_ns'])) if j not in used_p)})
            scores.append(dict(score, index=i, method=method, fold=c['fold'], kind=c['kind'], condition=c['condition'], gap_ns=c['gap_ns'],
                dataset=c['dataset'], temporal_split=split, temporal_merge=merge, true_count_supplied=method != 'blind_energy_gate',
                true_time_supplied=method in ['fixed_true_time', 'free_from_true_time'],
                not_converged=not p['converged'], boundary=p['boundary']))
            predictions.append(dict(p, index=i, method=method))
        if c['kind'] == 'pair':
            fixed, free, known = [oracle[m] for m in ['fixed_true_time', 'free_from_true_time', 'true_count_blind_start']]
            search.append(dict(index=i, fold=c['fold'], condition=c['condition'], gap_ns=c['gap_ns'],
                blind_search_minus_true_start_cost=known['cost']-free['cost'],
                fixed_minus_released_cost=fixed['cost']-free['cost'],
                missing_blind_candidate=oracle['missing_blind_candidate']))
    write_lines(OUTPUT/'predictions.jsonl', predictions); write_lines(OUTPUT/'scores.jsonl', scores)
    write_lines(OUTPUT/'search-diagnostics.jsonl', search)
    grouped = defaultdict(list)
    for r in scores: grouped[r['method'], r['fold'], r['kind'], r['condition'], r['gap_ns']].append(r)
    summaries = []
    for (method, fold, kind, condition, gap), rr in grouped.items():
        matches = [m for r in rr for m in r['matches']]
        errors = error_summary([m['true_energy_MeV'] for m in matches], [m['estimated_energy_MeV'] for m in matches]) if matches else None
        summaries.append(dict(method=method, fold=fold, kind=kind, condition=condition, gap_ns=gap, frames=len(rr),
            **{key: sum(r[key] for r in rr) for key in ['target', 'missed', 'extra', 'output_count', 'strict_5', 'strict_10', 'strict_20',
                'frame_success_10', 'temporal_split', 'temporal_merge', 'not_converged', 'boundary']}, matched_energy=errors))
    aggregate = []
    for fold in [0, 1]:
        for method in METHODS:
            for weak, group in [(True, 'weak'), (False, 'strong')]:
                rr = [r for r in summaries if r['method'] == method and r['fold'] == fold and r['kind'] == 'pair' and ('weak' in r['condition']) == weak]
                if not rr: continue
                aggregate.append(dict(method=method, fold=fold, group=group, frames=sum(r['frames'] for r in rr),
                    strict_pair_success=float(np.mean([r['frame_success_10']/r['frames'] for r in rr])),
                    successful_pairs=sum(r['frame_success_10'] for r in rr),
                    not_converged=sum(r['not_converged'] for r in rr)))
    save(OUTPUT/'report.json', dict(summary=summaries, aggregate=aggregate, methods=METHODS,
        interpretation='Development diagnosis, common observed weighting; truth-assisted rows are not blind recovery or rigorous limits',
        gate_counts=list(map(len, g['scores'])), inherited_model_unchanged=True, default_changed=False, new_ROOTs=0))
    print('ORACLE_EVALUATED', aggregate, flush=True)


def main():
    p = argparse.ArgumentParser(); p.add_argument('command', choices=['prepare', 'pilot', 'blind', 'gate', 'oracle', 'evaluate'])
    p.add_argument('--workers', type=int, default=2); args = p.parse_args()
    if args.command in ['blind', 'oracle']: infer(args.command, args.workers)
    elif args.command == 'pilot': pilot(args.workers)
    else: globals()[args.command]()


if __name__ == '__main__': main()
