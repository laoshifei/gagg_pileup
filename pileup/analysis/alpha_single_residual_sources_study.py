"""Frozen mother/repeat residual experiment; run from the repository root."""
from __future__ import annotations

import argparse
import ast
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

from .alpha_pair_digitizer_optimized import simulate_avalanches
from .alpha_pair_energy_calibrated import CalibratedEnergyResponse, learned_metric
from .alpha_pair_recovery_energy import RecoveryEnergyResponse
from .alpha_pair_spatial_recovery import SpatialRecovery
from .sipm_response import SensorConfig, PhotonHits, AvalancheResult, map_cells, read_photon_hits
from .sipm_readout import ReadoutConfig, synthesize_readout
from .alpha_single_residual_sources import (variance_components, component_intervals,
    opposite_half_targets, prefix_features, select_ridge, ridge_predict, paired_loss_ratio)

NAME = 'alpha_single_residual_sources_20261102'
ROOT = Path(os.environ.get('GAGG_RESIDUAL_SOURCES_ROOT', Path(__file__).resolve().parents[2])).resolve()
EXP = ROOT/'pileup/experiments'/NAME
ART = ROOT/'jupyter/artifacts'/NAME
REC = ROOT/'jupyter/run_records'/NAME
PARENT = ROOT/'pileup/experiments/alpha_pair_residual_projection_20261102'
WORKERS, REPEATS = 8, 16
HORIZONS = (8., 40., 200., 1920.)
GAPS = (8., 40., 200.)


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.partial')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        while block := f.read(1024*1024):
            h.update(block)
    return h.hexdigest()


def snapshot():
    target = EXP/'source/pileup/analysis'
    target.mkdir(parents=True)
    (target.parent/'__init__.py').write_text('')
    (target/'__init__.py').write_text('')
    pending, copied = ['alpha_single_residual_sources_study', 'test_alpha_single_residual_sources'], set()
    while pending:
        name = pending.pop()
        if name in copied:
            continue
        p = ROOT/'pileup/analysis'/f'{name}.py'
        tree = ast.parse(p.read_text())
        shutil.copyfile(p, target/p.name); copied.add(name)
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.level == 1:
                pending.extend([node.module] if node.module else [a.name for a in node.names])
    for name in ('sipm_age_mean.cpp', 'sipm_age_mean_jacobian.cpp'):
        shutil.copyfile(ROOT/'pileup/analysis'/name, target/name)
    (EXP/'bin').mkdir()
    command = ['g++', '-O3', '-std=c++17', '-fPIC', '-shared',
        str((target/'sipm_age_mean.cpp').relative_to(ROOT)),
        str((target/'sipm_age_mean_jacobian.cpp').relative_to(ROOT)),
        '-o', str((EXP/'bin/age.so').relative_to(ROOT))]
    subprocess.run(command, cwd=ROOT, check=True)
    save(EXP/'build.json', dict(command=command, binary_sha256=sha(EXP/'bin/age.so'),
        compiler=subprocess.check_output(['g++', '--version'], text=True).splitlines()[0]))
    (EXP/'run.py').write_text(
        "from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\n"
        "os.environ['GAGG_RESIDUAL_SOURCES_ROOT']=str(f.parents[2])\n"
        "os.environ['OPENBLAS_NUM_THREADS']='1'\nos.environ['OMP_NUM_THREADS']='1'\n"
        "sys.path.insert(0,str(f/'source'))\n"
        "if __name__=='__main__':\n    from pileup.analysis.alpha_single_residual_sources_study import main\n    main()\n")


def extract(row):
    config = read(EXP/'sensor-config.json')
    source = ROOT/row['root']
    assert sha(source) == row['root_original_sha256']
    metadata = read(ROOT/row['metadata_file'])
    assert metadata['config']['sensor'] == config['sensor']
    assert metadata['config']['readout'] == config['readout']
    sensor = SensorConfig(**config['sensor'])
    hits = read_photon_hits(source)
    active = map_cells(hits, sensor) >= 0
    assert active.any()
    first = float(hits.time_ns[active].min())
    np.savez_compressed(ART/'photons'/f"{row['id']:04d}.npz", time_ns=hits.time_ns-first,
        x_mm=hits.x_mm, y_mm=hits.y_mm, energy_eV=hits.energy_eV)
    save(EXP/'original-runinfo'/f"{row['id']:04d}.json", hits.metadata)
    shutil.copyfile(ROOT/row['metadata_file'], EXP/'original-metadata'/f"{row['id']:04d}.json")
    return dict(id=row['id'], root=row['root'], first_active_ns=first,
        active_photons=int(active.sum()), original_window_ns=hits.window_ns,
        aligned_clock='First active-area photon at128ns; original RunInfo separately preserved')


def prepare():
    for p in (EXP, ART, REC):
        p.mkdir(parents=True, exist_ok=False)
    for name in ('photons', 'voltage', 'diagnostics', 'features', 'pilot', 'replay', 'figures'):
        (ART/name).mkdir()
    for name in ('original-runinfo', 'original-metadata'):
        (EXP/name).mkdir()
    for name in ('template.npz', 'energy-model.json', 'response-parameters.json', 'mean-correction.npz',
                 'noise-calibration.json', 'sensor-config.json', 'charge-map-spatial12.npz', 'spatial-profile.npz'):
        shutil.copyfile(PARENT/name, EXP/name)
    previous = ROOT/'pileup/experiments/alpha_pair_energy_calibrated_20261030'
    pool = read(previous/'mother-provenance.json')['development']
    oldsplit = read(previous/'mother-split.json')
    mothers = []
    for group, source in [('weak', 1.625), ('high', 8.)]:
        candidates = [r for r in pool if r['fold'] == 1 and r['source_energy_MeV'] == source]
        candidates.sort(key=lambda r: hashlib.sha256(('735000:'+r['root']).encode()).hexdigest())
        assert len(candidates) >= 32
        for j, row in enumerate(candidates[:32]):
            i = len(mothers)
            mothers.append(dict(row, id=i, group=group, role='train' if j < 16 else 'validation',
                seed_base=735000000+i*100, repeats=REPEATS, time_ns=128., window_ns=4096.))
    roots = {r['root'] for r in mothers}
    assert len(roots) == 64 and roots <= set(oldsplit['noise'])
    assert not roots & (set(oldsplit['mean']) | set(oldsplit['evaluation']))
    save(EXP/'mothers.json', mothers)
    save(EXP/'split.json', dict(train=[r['root'] for r in mothers if r['role']=='train'],
        validation=[r['root'] for r in mothers if r['role']=='validation'],
        original_fold=1, disjoint_from_frozen_mean_training=True, disjoint_from_recent_evaluation=True,
        exposure='All64 were historically used for noise calibration; new split is only a holdout for this residual predictor, not fresh transport confirmation'))
    save(EXP/'protocol.json', dict(created_utc=datetime.now(timezone.utc).isoformat(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip(),
        git_status=subprocess.check_output(['git','status','--short'], cwd=ROOT, text=True),
        purpose='Distinguish conditional mean error, persistent mother variation and within-mother digitization; test observable-prefix prediction of persistent strong isolated residual',
        selection='SHA256(735000:ROOT path) ordered first32 fold1 mothers at each source energy1.625/8MeV; first16train,next16validation. Source energy is not E_dep; no outcome selection',
        mothers=64, repeats=REPEATS, digitizations=64*REPEATS, workers=WORKERS, blas_threads=1,
        source_rate='No Poisson scan; exactly one stored optical mother per offline4096ns waveform. Historical RunInfo/macros archived verbatim',
        observation='First2048ns,1024 samples,2ns period,0.5ns aperture at each period start; onset128ns,phase0; no ADC quantization/clipping',
        physical_model='Unchanged95/351ns f=.3 mean anchor,29ns recovery proxy,14400cells,12x12 frozen spatial prior. Actual old optical photons replayed with original times/positions/energies',
        residual='Event-ancestor signal minus causal isolated mean at true E_dep; full analog minus background and mean retained separately. Ancestry is diagnostic truth',
        variance='r_mk=b+u_m+epsilon_mk, conditional on group/fixed phase. Between sample variance of repeat means minus within/R; retain negative estimates;4000 mother-cluster bootstrap, fixed repeats',
        conditioning='Subtract each mother model at its actual E_dep before grouping. Remaining E-dependent mean error/geometry can contribute to estimated between term; not a pure microscopic optical-variance identification',
        projections='Frozen GLS influence at nominal7.5MeV strong/.25MeV weak,128ns onset and gaps8/40/200ns. Same operator across mothers;100delta lnE units, not percent energy error or physical resolution',
        nominal_energy_MeV=dict(high=7.5, weak=.25), gaps_ns=list(GAPS),
        prediction='Train only16strong mothers; validation16strong mothers. One noisy observed prefix predicts mean ancestor-residual projection of the opposite8repeat half. This isolates persistent mother response; target still has finite-repeat noise',
        features='Known onset/K=1; prefix-only unit-voltage least-squares E fit,logE plus4 normalized residual block means. No true E/source energy/photons/ancestry features. Prefix8/40/200ns or full1920ns after onset',
        horizons_ns=list(HORIZONS), prediction_candidates=['training_mean','voltage_energy_only','voltage_energy_plus_shape'],
        regression='Standardize inside training fold; ridge penalties .01/.1/1/10 chosen by leave-one-mother-out MSE on16training mothers. Repeats stay grouped; no validation choice',
        comparisons='Gap-matched uncontaminated prefix primary; full isolated window secondary information control unavailable in an overlapping pair. Fixed narrow source groups do not test general energy-spectrum transfer',
        prediction_gate=dict(shape_mse_ratio_to_mean=.90, ratio_ci95_upper=1., shape_mse_ratio_to_energy=.95),
        gate_scope='Engineering gate for a subsequent pair-model study only; no default replacement. Pass requires each criterion at both8/40ns. Thresholds are project proposals, not literature guarantees',
        pilot=dict(ids=[0,32], repeats=4, projected_seconds_limit=900., safety_factor=1.5),
        replay=dict(ids=[1,33], repeats=2),
        no_new_ROOT=True, default_changed=False, literature='ProjectL91 and existing source/review hashes; no new literature claim'))
    for name in ('docs/literature-use.md', 'README.md', 'jupyter/codex_0007_SiPM_Energy_Calibration.ipynb',
                 'jupyter/Research_Plan_20260921.ipynb', 'jupyter/Run_Catalog.ipynb'):
        p = ROOT/name
        shutil.copyfile(p, REC/(p.name+'.before'))
    shutil.copyfile(ROOT/'jupyter/run_records/paper_search_strong_shape_20261102/evidence.json', EXP/'literature-evidence.json')
    snapshot()
    save(EXP/'design-freeze.json', {str(p.relative_to(ROOT)):sha(p) for p in EXP.rglob('*')
        if p.is_file() and 'bin' not in p.relative_to(EXP).parts})
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        rows = list(pool.map(extract, mothers))
    save(EXP/'extraction.json', rows)
    import scipy
    save(EXP/'environment.json', dict(python=sys.version, executable=sys.executable,
        numpy=np.__version__, scipy=scipy.__version__, workers=WORKERS, blas_threads=1))
    inputs = [ROOT/r[k] for r in mothers for k in ('root','metadata_file','truth_file')]
    inputs += list((ART/'photons').glob('*.npz'))
    inputs += list((EXP/'original-runinfo').glob('*.json')) + list((EXP/'original-metadata').glob('*.json'))
    save(EXP/'input-freeze.json', {str(p.relative_to(ROOT)):sha(p) for p in inputs})
    print('Prepared64 mothers,1024 planned digitizations;16train/16validation per energy group', flush=True)


def verify():
    counts = {}
    for name in ('design-freeze.json', 'input-freeze.json'):
        values = read(EXP/name)
        for p, expected in values.items():
            assert sha(ROOT/p) == expected, p
        counts[name] = len(values)
    return counts


def initialize():
    if not hasattr(np, 'trapezoid'):
        np.trapezoid = np.trapz
    global CONFIG, SENSOR, READOUT, MODEL, NOISE
    CONFIG = read(EXP/'sensor-config.json')
    SENSOR, READOUT = SensorConfig(**CONFIG['sensor']), ReadoutConfig(**CONFIG['readout'])
    anchor = CalibratedEnergyResponse.load(RecoveryEnergyResponse.load(EXP), EXP/'mean-correction.npz')
    with np.load(EXP/'charge-map-spatial12.npz') as f:
        mapping = {k:f[k].copy() for k in f.files}
    with np.load(EXP/'spatial-profile.npz') as f:
        profile = f['profile_12'].copy()
    MODEL = SpatialRecovery(anchor, CONFIG, EXP/'bin/age.so', mapping, profile)
    NOISE = read(EXP/'noise-calibration.json')


def digitize(job):
    row, repeats, directory = job
    start = time.monotonic()
    with np.load(ART/'photons'/f"{row['id']:04d}.npz") as f:
        times = f['time_ns']+128.
        keep = times < 4096.
        hits = PhotonHits(time_ns=times[keep], x_mm=f['x_mm'][keep], y_mm=f['y_mm'][keep],
                          energy_eV=f['energy_eV'][keep], window_ns=4096.)
    voltage, signal, charges = [], [], []
    for repeat in range(repeats):
        seed = row['seed_base']+repeat
        av = simulate_avalanches(hits, SENSOR, seed)
        voltage.append(synthesize_readout(av, READOUT, seed+10000000).arrays['analog_mV'][:1024])
        records = av.records
        owner = records['kind']==0
        for j in np.flatnonzero(records['parent_id'] >= 0):
            owner[j] = owner[records['parent_id'][j]]
        sub = records[owner]
        ancestor = AvalancheResult(sub, 4096., av.full_charge_C, np.array([]), {})
        signal.append(synthesize_readout(ancestor, READOUT, 0).arrays['signal_mV'][:1024])
        charges.append(float(sub['charge_C'].sum()/av.full_charge_C))
    voltage, signal = np.asarray(voltage), np.asarray(signal)
    assert np.isfinite(voltage).all() and np.isfinite(signal).all()
    assert np.max(abs(signal[:, :64])) < 1e-10
    if directory == 'production':
        vpath, dpath = ART/'voltage', ART/'diagnostics'
    else:
        vpath, dpath = ART/directory/'voltage', ART/directory/'diagnostics'
    for p in (vpath, dpath):
        p.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(vpath/f"{row['id']:04d}.npz", voltage_mV=voltage)
    np.savez_compressed(dpath/f"{row['id']:04d}.npz", signal_mV=signal, charge_pe=charges)
    result = dict(id=row['id'], group=row['group'], repeats=repeats, seconds=time.monotonic()-start,
        seed_base=row['seed_base'], readout_seed_offset=10000000, photons_after_window=int((~keep).sum()),
        finite=True, pre_event_signal_zero=True)
    save(REC/directory/f"{row['id']:04d}.json", result)
    return result


def generate(stage):
    verify()
    rows = read(EXP/'mothers.json')
    protocol = read(EXP/'protocol.json')
    count = REPEATS
    if stage in ('pilot','replay'):
        spec = protocol[stage]
        rows = [r for r in rows if r['id'] in spec['ids']]
        count = spec['repeats']
    if stage == 'production':
        assert read(EXP/'pilot-gate.json')['passed']
    paths = list((ART/('voltage' if stage=='production' else stage+'/voltage')).glob('*.npz'))
    assert not paths, 'Keep existing outputs; use another experiment for changed inputs'
    start = time.monotonic()
    with ProcessPoolExecutor(max_workers=WORKERS, initializer=initialize) as pool:
        results = []
        for row in pool.map(digitize, [(r,count,stage) for r in rows]):
            results.append(row)
            if len(results)%8 == 0:
                print(stage, len(results), '/', len(rows), flush=True)
    report = dict(rows=results, digitizations=len(rows)*count, wall_seconds=time.monotonic()-start)
    save(REC/(stage+'.json'), report)
    if stage == 'pilot':
        estimate = sum(r['seconds']/count*32*REPEATS for r in results)/WORKERS*1.5
        save(EXP/'pilot-gate.json', dict(projected_seconds=estimate, limit_seconds=900., passed=bool(estimate<900.)))
        assert estimate < 900., estimate
    if stage == 'production':
        save(EXP/'waveform-freeze.json', {str(p.relative_to(ROOT)):sha(p) for directory in ('voltage','diagnostics')
            for p in (ART/directory).glob('*.npz')})
    if stage in ('production','replay'):
        compare = 'pilot' if stage=='production' else 'replay'
        n = protocol[compare]['repeats']
        for i in protocol[compare]['ids']:
            for folder in ('voltage','diagnostics'):
                with np.load(ART/folder/f'{i:04d}.npz') as a, np.load(ART/compare/folder/f'{i:04d}.npz') as b:
                    assert a.files == b.files
                    for k in a.files:
                        np.testing.assert_array_equal(a[k][:n], b[k])
        save(EXP/(compare+'-equality.json'), dict(passed=True, exact_arrays=True, repeats_per_mother=n,
                                                ids=protocol[compare]['ids']))
    print(stage, report['digitizations'], 'digitizations complete', flush=True)


def feature_initialize():
    initialize()
    def audit(event, args):
        if event != 'open' or not isinstance(args[0], (str,bytes,os.PathLike)):
            return
        p = Path(os.fsdecode(args[0])).resolve()
        forbidden = p.suffix == '.root' or p.name in ('mothers.json','split.json')
        forbidden |= p.is_relative_to(ART/'diagnostics') or p.is_relative_to(ART/'photons')
        forbidden |= p.is_relative_to(EXP/'original-runinfo') or p.is_relative_to(EXP/'original-metadata')
        if forbidden:
            raise PermissionError('Voltage feature worker cannot read event truth: '+str(p))
    sys.addaudithook(audit)
    for p in (EXP/'mothers.json', ART/'diagnostics/0000.npz', ART/'photons/0000.npz'):
        try:
            p.open('rb')
        except PermissionError:
            pass
        else:
            raise AssertionError('Truth read guard failed')


def feature_one(i):
    start = time.monotonic()
    with np.load(ART/'voltage'/f'{i:04d}.npz') as f:
        voltage = f['voltage_mV']
    features, energy, convergence, boundary = [], [], [], []
    for horizon in HORIZONS:
        rows = [prefix_features(MODEL, y, horizon, noise_floor_mV=READOUT.noise_rms_mV) for y in voltage]
        features.append([r[0] for r in rows])
        energy.append([r[1]['energy_MeV'] for r in rows])
        convergence.append([r[1]['converged'] for r in rows])
        boundary.append([r[1]['boundary'] for r in rows])
        assert all(r[1]['aperture_end_ns']<=128.+horizon for r in rows)
    np.savez_compressed(ART/'features'/f'{i:04d}.npz', features=features, energy_MeV=energy,
                        converged=convergence, boundary=boundary)
    return dict(id=i, seconds=time.monotonic()-start, input_truth_guard=True,
                converged=int(np.sum(convergence)), boundary=int(np.sum(boundary)), fits=REPEATS*len(HORIZONS))


def features():
    verify()
    for p,d in read(EXP/'waveform-freeze.json').items():
        assert sha(ROOT/p)==d, p
    assert not list((ART/'features').glob('*.npz'))
    start = time.monotonic()
    with ProcessPoolExecutor(max_workers=WORKERS, initializer=feature_initialize) as pool:
        rows = []
        for r in pool.map(feature_one, range(64)):
            rows.append(r)
            if len(rows)%8==0:
                print('voltage features',len(rows),'/64',flush=True)
    save(REC/'features.json', dict(rows=rows, wall_seconds=time.monotonic()-start,
        fits=sum(r['fits'] for r in rows), true_energy_features=False, known_onset_and_K1=True))
    save(EXP/'feature-freeze.json', {str(p.relative_to(ROOT)):sha(p) for p in (ART/'features').glob('*.npz')})


def influence_operators():
    """Fixed nominal geometry, diagnostic truth coordinates only."""
    transform = MODEL.anchor.transform(np.eye(1024)).T
    result, checks = {}, []
    for group, energy in [('weak', .25), ('high', 7.5)]:
        dv = (MODEL.isolated([128.], [energy*np.exp(1e-4)])[0][0] -
              MODEL.isolated([128.], [energy*np.exp(-1e-4)])[0][0])/2e-4
        metric = learned_metric(MODEL.anchor, [128.], [energy], NOISE)
        operator = metric.apply(transform)
        jac = (operator@dv)[:, None]
        influence = np.linalg.pinv(jac, rcond=1e-12)@operator
        error = float(np.max(abs(influence@dv-1)))
        assert error < 1e-9
        checks.append(dict(name=group+'_single', tangent_identity_error=error))
        result[group] = [influence[0]]
    for gap in GAPS:
        times = [128., 128.+gap]
        energy = np.array([7.5, .25])
        mean, gradient = MODEL.mean_and_gradient(times, np.log(energy), MODEL.prepare(times), strength=1.)
        metric = learned_metric(MODEL.anchor, times, energy, NOISE)
        operator = metric.apply(transform)
        jac = operator@gradient.T
        influence = np.linalg.pinv(jac, rcond=1e-12)@operator
        error = float(np.max(abs(influence@gradient.T-np.eye(2))))
        assert error < 1e-9
        result['high'].append(influence[1])
        checks.append(dict(name=f'pair_{gap:g}', tangent_identity_error=error,
                           singular_values=np.linalg.svd(jac, compute_uv=False).tolist()))
    save(EXP/'projection-controls.json', dict(checks=checks, passed=True,
        fixed_reference=True, truth_E_used_for_diagnostic_residual_only=True))
    np.savez_compressed(ART/'influences.npz', **{k:np.array(v) for k,v in result.items()})
    return {k:np.array(v) for k,v in result.items()}


def analyze():
    verify(); initialize()
    assert not (EXP/'result-summary.json').exists()
    for name in ('waveform-freeze.json','feature-freeze.json'):
        for p,d in read(EXP/name).items():
            assert sha(ROOT/p)==d, p
    operators = influence_operators()
    mothers = read(EXP/'mothers.json')
    stat_rows, energy_rows, prediction_rows = [], [], []
    regression_models = {}
    features_all, signals_all, voltages_all, expected_all = [], [], [], []
    feature_meta = []
    for row in mothers:
        i = row['id']
        with np.load(ART/'voltage'/f'{i:04d}.npz') as f: v = f['voltage_mV'].copy()
        with np.load(ART/'diagnostics'/f'{i:04d}.npz') as f: s = f['signal_mV'].copy()
        with np.load(ART/'features'/f'{i:04d}.npz') as f:
            features_all.append(f['features'].copy())
            feature_meta.append({k:f[k].copy() for k in ('energy_MeV','converged','boundary')})
        expected_all.append(MODEL.isolated([128.], [row['edep_MeV']])[0][0])
        signals_all.append(s); voltages_all.append(v)
    observed, signals, expected, feats = map(np.asarray,
        (voltages_all, signals_all, expected_all, features_all))
    residual = signals-expected[:,None,:]
    full = observed-MODEL.background-expected[:,None,:]
    remainder = observed-MODEL.background-signals
    closure = float(np.max(abs(residual+remainder-full)))
    assert closure < 1e-9
    summary_arrays = dict(signal_residual_mV=residual, full_residual_mV=full,
        mean_prediction_mV=expected, signal_mV=signals, voltage_mV=observed,
        time_ns=MODEL.readout.adc_centers)
    for group in ('weak','high'):
        idx = np.array([r['id'] for r in mothers if r['group']==group])
        projection_names = ['isolated_logE'] + ([f'weak_logE_gap{int(g)}' for g in GAPS] if group=='high' else [])
        for kind, values in [('event_signal', residual), ('full_voltage', full)]:
            projected = 100*np.einsum('mrt,qt->mrq', values[idx], operators[group])
            summary_arrays[group+'_'+kind+'_projected_100logE'] = projected
            stats = variance_components(projected)
            intervals = component_intervals(projected, seed=735210+(group=='high')*10+(kind=='full_voltage'))
            for j, name in enumerate(projection_names):
                b, w = float(stats['between'][j]), float(stats['within'][j])
                stat_rows.append(dict(group=group, kind=kind, projection=name, mothers=32, repeats=REPEATS,
                    mean=float(stats['mean'][j]), mean_ci95=intervals['mean'][:,j].tolist(),
                    within_variance=w, within_ci95=intervals['within'][:,j].tolist(),
                    between_variance_unclipped=b, between_ci95=intervals['between'][:,j].tolist(),
                    mean16_variance=float(stats['variance_of_repeat_mean'][j]),
                    within_fraction_single=w/(w+b) if w+b>0 else None,
                    between_fraction_mean8=b/(b+w/8) if b>=0 and b+w/8>0 else None,
                    units='mean:100 delta ln E; variance:(100 delta ln E)^2',
                    moment_identity_error=stats['closure_error']))
        for kind, values in [('event_signal',residual),('full_voltage',full)]:
            stats = variance_components(values[idx])
            for key in ('mean','within','between','variance_of_repeat_mean'):
                summary_arrays[f'{group}_{kind}_{key}'] = stats[key]
        for role in ('train','validation'):
            subset = [r['id'] for r in mothers if r['group']==group and r['role']==role]
            for h,horizon in enumerate(HORIZONS):
                errors = np.array([100*(feature_meta[i]['energy_MeV'][h]/mothers[i]['edep_MeV']-1) for i in subset])
                energy_rows.append(dict(group=group, role=role, horizon_ns=horizon,
                    mothers=len(subset), repeat_fits=errors.size, bias_percent=float(errors.mean()),
                    rmse_percent=float(np.sqrt(np.mean(errors**2))), error_sd_percent=float(errors.std(ddof=1)),
                    within10_count=int(np.sum(abs(errors)<=10)),
                    boundary=int(sum(feature_meta[i]['boundary'][h].sum() for i in subset)),
                    nonconverged=int(sum((~feature_meta[i]['converged'][h]).sum() for i in subset)),
                    estimator='Single-prefix unit-voltage LS diagnostic,known onset/K=1; not default full GLS'))
    strong = [r['id'] for r in mothers if r['group']=='high']
    train = np.array([j for j,i in enumerate(strong) if mothers[i]['role']=='train'])
    val = np.array([j for j,i in enumerate(strong) if mothers[i]['role']=='validation'])
    strong_scores = summary_arrays['high_event_signal_projected_100logE'][:,:,1:]
    targets = opposite_half_targets(strong_scores)
    strong_features = feats[strong]
    for j,gap in enumerate(GAPS):
        for horizon in (gap,1920.):
            h = HORIZONS.index(horizon)
            x, y = strong_features[:,h], targets[:,:,j]
            baseline = np.full_like(y[val], y[train].mean())
            base_loss = ((baseline-y[val])**2).mean(1)
            pred = {'training_mean':baseline}
            key = f'gap{int(gap)}_prefix{int(horizon)}'
            regression_models[key] = {}
            for method, columns in [('voltage_energy_only',slice(0,1)), ('voltage_energy_plus_shape',slice(None))]:
                fit, cv = select_ridge(x[train,:,columns], y[train])
                pred[method] = ridge_predict(fit, x[val,:,columns])
                regression_models[key][method] = dict(fit=fit, selection=cv)
            losses = {method:((value-y[val])**2).mean(1) for method,value in pred.items()}
            shape_mean = paired_loss_ratio(losses['voltage_energy_plus_shape'], base_loss, seed=735310+j*10+h)
            shape_energy = paired_loss_ratio(losses['voltage_energy_plus_shape'],losses['voltage_energy_only'],seed=735410+j*10+h)
            energy_mean = paired_loss_ratio(losses['voltage_energy_only'],base_loss,seed=735510+j*10+h)
            passed = (shape_mean['ratio'] < .9 and shape_mean['ci95'][1] < 1. and shape_energy['ratio'] < .95)
            row = dict(gap_ns=gap, horizon_ns=horizon, accessible_before_weak=horizon==gap,
                training_mothers=len(train), validation_mothers=len(val), repeat_predictions=len(val)*REPEATS,
                shape_vs_mean=shape_mean, shape_vs_energy=shape_energy, energy_vs_mean=energy_mean,
                gate_passed=bool(passed), baseline_training_mean=float(y[train].mean()),
                validation_target_mean=float(y[val].mean()),
                opposite_half_target_noise_variance=float(strong_scores[val,:,j].var(1,ddof=1).mean()/8),
                baseline_bias=float((baseline-y[val]).mean()),
                energy_predictor_bias=float((pred['voltage_energy_only']-y[val]).mean()),
                shape_predictor_bias=float((pred['voltage_energy_plus_shape']-y[val]).mean()),
                units='Target and prediction:100 delta ln E; not nonlinear relative E error')
            prediction_rows.append(row)
            summary_arrays[key+'_target'] = y[val]
            for method,p in pred.items():
                summary_arrays[key+'_'+method] = p
            summary_arrays[key+'_validation_mother_ids'] = np.array(strong)[val]
    primary = [r for r in prediction_rows if r['accessible_before_weak'] and r['gap_ns'] in (8.,40.)]
    passed = all(r['gate_passed'] for r in primary)
    summary = dict(variance=stat_rows, isolated_energy=energy_rows, prediction=prediction_rows,
        decision=dict(primary_both_gaps_passed=passed,
            action='Eligible to design a separate pair validation; no demonstrated recovery yet' if passed else
                   'Do not promote this prefix shape predictor; inspect remaining mean/variance limits before a new pair model',
            no_blind_performance_measured=True, default_changed=False),
        exact_voltage_closure_mV=closure, mother_count=64, formal_digitizations=1024,
        limitations=['Historical noise-pool reuse; validation only for this new predictor',
            'Only two narrow source-energy groups; source energy is not deposited energy',
            'Known onset/K=1; opposite-half target is persistent mother response with finite averaging noise',
            'Full isolated-window features unavailable in an overlapping pair',
            'Fixed nominal logE linear projections are not actual pair-energy reconstruction',
            'No physical/mean parameter changes; covariance/predictor uncertainty not propagated into real energy fits'])
    save(EXP/'result-summary.json', summary)
    save(EXP/'regression-models.json', regression_models)
    np.savez_compressed(ART/'analysis.npz', **summary_arrays)
    save(EXP/'analysis-freeze.json', {str(p.relative_to(ROOT)):sha(p) for p in
        [EXP/'result-summary.json',EXP/'regression-models.json',EXP/'projection-controls.json',ART/'analysis.npz',ART/'influences.npz']})
    print(json.dumps(dict(primary_prediction_gate=passed,
        accessible=[dict(gap_ns=r['gap_ns'],shape_vs_mean=r['shape_vs_mean'],shape_vs_energy=r['shape_vs_energy'],passed=r['gate_passed']) for r in primary]),indent=2),flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['prepare','pilot','production','replay','features','analyze','verify'])
    stage = parser.parse_args().stage
    if stage == 'prepare': prepare()
    elif stage in ('pilot','production','replay'): generate(stage)
    elif stage == 'features': features()
    elif stage == 'analyze': analyze()
    else: print(json.dumps(verify()))


if __name__ == '__main__':
    main()
