"""Frozen, reproducible known-K/t deconvolution pilot; no Geant4 production."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import time

import numpy as np
import scipy

from .alpha_pair_deconvolution import AgeRecovery, EnergyProfiles, FiniteConvolution, calibrate_profiles, spe_kernel

NAME = 'alpha_pair_deconvolution_20261030'
ROOT = Path(os.environ.get('GAGG_DECONV_PROJECT', Path(__file__).resolve().parents[2]))
EXP = ROOT/'pileup/experiments'/NAME
ART = ROOT/'jupyter/artifacts'/NAME
REC = ROOT/'jupyter/run_records'/NAME
OLD = ROOT/'pileup/experiments/alpha_pair_energy_calibrated_20261030'
OLDART = ROOT/'jupyter/artifacts/alpha_pair_energy_calibrated_20261030'
ARMS = ('rl128', 'recovery128', 'rl512', 'recovery512')


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def read_json(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hashes(paths):
    return {str(p.relative_to(ROOT)): sha(p) for p in paths}


def prepare():
    if (EXP/'protocol.json').exists():
        raise FileExistsError('Protocol already frozen; do not overwrite an experiment')
    for p in (EXP, ART/'inputs', ART/'calibration', ART/'cases', REC):
        p.mkdir(parents=True, exist_ok=True)
    for name in ('sensor-config.json', 'energy-model.json'):
        shutil.copyfile(OLD/name, EXP/name)
    with np.load(OLD/'template.npz') as f:
        background = float(f['baseline_mV'])+float(f['background_mV'])
    protocol = dict(created=datetime.now(timezone.utc).isoformat(), commit='09732062ddc851946e7215a481e4e4a12484fbe4',
        status='New uncommitted research after published checkpoint; actual source snapshot below',
        scope='Pure alpha, supplied TRUE count and times; only E_dep unknown. Not blind or alpha/gamma PSD.',
        literature=dict(url='https://arxiv.org/abs/2010.14886', doi='10.48550/arXiv.2010.14886',
                        location='PDF pp19-24, Eqs5.8/5.12; preprint, formal journal version not verified',
                        adaptation='Our SPE, age-distribution mean field, stopping rule and energy profiles; NOT full author MC reproduction'),
        input='Complete 1024x2ns analog voltage, enough dynamic range; no clipping rejection',
        background_mV=background, latent_grid='2ns counts at bin centers1+2n; voltage centers0.25+2n',
        rl_iterations=[128, 512], rl_initial_counts_per_bin=1., primary='recovery128', comparator='new_mean_new_noise',
        recovery=dict(cells=14400, tau_ns=29., dt_ns=2., tail_tau=12, cross_iterations=6,
                      max_inverse_iterations=40, relative_L1_tolerance=.005,
                      limitations='Uniform illumination, at most one firing/cell/bin, mean-field prompt crosstalk; no dark-history correction'),
        calibration='280 fold0 isolated noisy voltages; 80 fold1 empty voltages for transformed background; no pair training',
        profiles='4 E_dep knots; normalized by inherited q0(E); 4ns mean smoothing, -64..1720ns with last80ns taper; constant normalized shape outside knots',
        objective='Fixed diagonal working weights: exposure/max(10ns-smoothed transformed counts,1). Correlations not modeled; not true Poisson likelihood.',
        evaluation='First8 per each of8 existing conditions, 64frames. Exposed historical holdout, NOT untouched confirmation; mothers reused across conditions.',
        endpoint='All injected E_dep within5/10/20%; no Q selection. Bias,sigma68,RMSE,boundaries,nonconvergence; K/t metrics not measured.',
        expansion_gate='Expand to480 only if recovery128 pair success count>=GLS and each isolated group loses at most1 of8, numerical controls pass, estimated cost<1800s. Pilot results remain pilot; no parameter tuning.',
        display_ids=[40, 120], workers=6, BLAS_threads=1, no_new_ROOT=True, default_changed=False)
    write_json(EXP/'protocol.json', protocol)
    originals = read_json(OLDART/'calibration/mean.json')
    records = [r for r in originals if r['kind'] == 'single']
    joint = ROOT/'jupyter/artifacts/alpha_pair_joint_likelihood_20261026/calibration.observations.npz'
    with np.load(joint) as f:
        cal = f['voltage_mV'][[r['source_row'] for r in records]]
    np.savez_compressed(ART/'calibration/observations.npz', voltage_mV=cal)
    write_json(ART/'calibration/labels.json', records)
    shutil.copyfile(OLDART/'calibration/empty.npz', ART/'calibration/empty.npz')
    labels = read_json(OLDART/'labels.json')
    ids = [start+j for start in (0, 40, 80, 120, 160, 240, 320, 400) for j in range(8)]
    np.savez_compressed(ART/'inputs/ids.npz', ids=np.asarray(ids))
    with np.load(OLDART/'inputs/observations.npz') as f:
        np.savez_compressed(ART/'inputs/observations.npz', **{k: f[k][ids] for k in f.files}, ids=np.asarray(ids))
    selected = [labels[i] for i in ids]
    write_json(ART/'labels.json', selected)
    calroots = sorted({root for r in records for root in r['mother_roots']})
    evalroots = sorted({root for r in selected for root in r['mother_roots']})
    assert not set(calroots)&set(evalroots) and len(records) == 280
    write_json(EXP/'mother-split.json', dict(calibration=calroots, evaluation=evalroots,
                                            calibration_count=len(calroots), evaluation_count=len(evalroots), disjoint=True))
    write_json(EXP/'pilot-ids.json', ids)
    write_json(EXP/'input-freeze.json', hashes([joint, OLDART/'calibration/mean.json', OLDART/'inputs/observations.npz',
        OLDART/'labels.json', OLDART/'calibration/empty.npz', ART/'calibration/observations.npz', ART/'calibration/labels.json',
        ART/'calibration/empty.npz', ART/'inputs/observations.npz', ART/'labels.json']))
    write_json(REC/'environment.json', dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform()))
    for name in ('codex_0007_SiPM_Energy_Calibration.ipynb', 'Research_Plan_20260921.ipynb'):
        shutil.copyfile(ROOT/'jupyter'/name, REC/(name+'.before'))
    for p in (ROOT/'docs/literature-use.md', ROOT/'README.md'):
        shutil.copyfile(p, REC/(p.name+'.before'))
    print(json.dumps(dict(prepared=True, frames=64, calibration=280, background_mV=background)), flush=True)


def read_guard():
    """Inference forbids data-file reads after loading its allowed frozen model."""
    def audit(event, args):
        if event == 'open' and isinstance(args[0], (str, bytes, os.PathLike)):
            path = Path(os.fsdecode(args[0])); mode = args[1]
            if (mode is None or 'r' in mode) and path.suffix not in ('.py', '.pyc', '.so'):
                raise PermissionError(f'Inference data read blocked: {path.name}')
    sys.addaudithook(audit)
    try:
        (ART/'labels.json').read_text()
    except PermissionError:
        return True
    raise AssertionError('Truth read guard did not fire')


def initialize(inference=False):
    global OP, RECOVERY, BACKGROUND, MODELS, GUARD
    configuration = read_json(EXP/'sensor-config.json')
    OP = FiniteConvolution(spe_kernel(configuration), 1024)
    RECOVERY = AgeRecovery(configuration['sensor'])
    BACKGROUND = read_json(EXP/'protocol.json')['background_mV']
    MODELS = {}
    if inference:
        for arm in ARMS:
            with np.load(EXP/(arm+'-model.npz')) as f:
                MODELS[arm] = EnergyProfiles(**{k: f[k] for k in f.files})
        GUARD = read_guard()


def transform(voltage):
    before = time.monotonic()
    unfolded = OP.unfold(voltage, BACKGROUND)
    arrays, diagnostics = {}, {}
    for iterations, value in unfolded.items():
        arrays['rl'+str(iterations)] = value
        corrected, info = RECOVERY.invert(value)
        arrays['recovery'+str(iterations)] = corrected
        diagnostics[str(iterations)] = info
    return arrays, diagnostics, time.monotonic()-before


def calibration_job(item):
    identifier, voltage = item
    arrays, diagnostics, seconds = transform(voltage)
    return identifier, arrays, diagnostics, seconds


def calibrate(workers):
    if (EXP/'model-freeze.json').exists():
        raise FileExistsError('Models already frozen')
    start = time.monotonic()
    with np.load(ART/'calibration/observations.npz') as f:
        voltage = f['voltage_mV']
    with np.load(ART/'calibration/empty.npz') as f:
        empty = f['voltage_mV']
    output = {}
    with ProcessPoolExecutor(max_workers=workers, initializer=initialize) as pool:
        for count, (i, arrays, diagnostics, seconds) in enumerate(pool.map(calibration_job, enumerate(np.concatenate((voltage, empty))))):
            output[i] = (arrays, diagnostics, seconds)
            if count % 40 == 0:
                print(f'calibration {count+1}/{len(voltage)+len(empty)} wall={time.monotonic()-start:.1f}s', flush=True)
    model = read_json(EXP/'energy-model.json'); records = read_json(ART/'calibration/labels.json')
    initialize()
    archive = {}
    for arm in ARMS:
        all_arrays = np.stack([output[i][0][arm] for i in range(len(output))])
        archive[arm] = all_arrays
        response = calibrate_profiles(records, all_arrays[:len(records)], all_arrays[len(records):], model, OP.exposure)
        np.savez_compressed(EXP/(arm+'-model.npz'), **asdict(response))
    np.savez_compressed(ART/'calibration/transformed.npz', **archive)
    write_json(ART/'calibration/diagnostics.json', [dict(id=i, recovery=output[i][1], seconds=output[i][2]) for i in range(len(output))])
    write_json(EXP/'model-freeze.json', hashes([EXP/(a+'-model.npz') for a in ARMS]+[EXP/'protocol.json', EXP/'sensor-config.json', EXP/'energy-model.json']))
    write_json(REC/'calibration-run.json', dict(frames=len(output), seconds=time.monotonic()-start, workers=workers,
        inverse_failures=sum(not d['converged'] for _, diag, _ in output.values() for d in diag.values())))
    print('CALIBRATION_COMPLETE', flush=True)


def inference_job(item):
    identifier, voltage, known_times = item
    arrays, diagnostics, seconds = transform(voltage)
    start = time.monotonic()
    fits = {arm: MODELS[arm].fit(arrays[arm], known_times) for arm in ARMS}
    return dict(id=int(identifier), fits=fits, recovery=diagnostics, transform_seconds=seconds,
                fit_seconds=time.monotonic()-start, known_K_t=True, truth_E_Q_read=False, worker_truth_read_guard=GUARD), arrays


def run(workers):
    if list((ART/'cases').glob('*.json')):
        raise FileExistsError('Results exist; preserve the completed run')
    for path, expected in read_json(EXP/'model-freeze.json').items():
        assert sha(ROOT/path) == expected, path
    with np.load(ART/'inputs/observations.npz') as f:
        items = [(i, v, t[:k]) for i, v, t, k in zip(f['ids'], f['voltage_mV'], f['known_times_ns'], f['known_counts'])]
    start = time.monotonic(); curves = {arm: [] for arm in ARMS}
    with ProcessPoolExecutor(max_workers=workers, initializer=initialize, initargs=(True,)) as pool:
        for count, (result, arrays) in enumerate(pool.map(inference_job, items)):
            write_json(ART/'cases'/f"{result['id']:04d}.json", result)
            for arm in ARMS: curves[arm].append(arrays[arm])
            if count % 8 == 0: print(f'evaluation {count+1}/{len(items)} wall={time.monotonic()-start:.1f}s', flush=True)
    np.savez_compressed(ART/'transformed.npz', **{k: np.stack(v) for k, v in curves.items()})
    write_json(REC/'evaluation-run.json', dict(frames=len(items), fits=len(items)*len(ARMS), workers=workers, seconds=time.monotonic()-start))
    print('EVALUATION_COMPLETE', flush=True)


def main():
    p = argparse.ArgumentParser(); p.add_argument('mode', choices=['prepare', 'calibrate', 'run']); p.add_argument('--workers', type=int, default=6)
    args = p.parse_args()
    if args.mode == 'prepare': prepare()
    elif args.mode == 'calibrate': calibrate(args.workers)
    else: run(args.workers)


if __name__ == '__main__': main()
