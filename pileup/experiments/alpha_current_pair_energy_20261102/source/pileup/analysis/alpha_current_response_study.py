"""Frozen current-environment isolated controls before fresh-mother validation."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np

from . import alpha_mean_shape_study as frozen
from .alpha_current_response import paired_mother_interval, energy_metrics
from .alpha_transport_environment_study import load_root
from .sipm_response import read_photon_hits, SensorConfig, map_cells

NAME = 'alpha_current_response_20261102'
ROOT = Path(os.environ.get('GAGG_CURRENT_RESPONSE_ROOT', Path(__file__).resolve().parents[2])).resolve()
EXP = ROOT / 'pileup/experiments' / NAME
ART = ROOT / 'jupyter/artifacts' / NAME
REC = ROOT / 'jupyter/run_records' / NAME
PARENT = ROOT / 'pileup/experiments/alpha_mean_shape_20261102'
TRANSPORT = ROOT / 'pileup/experiments/alpha_transport_environment_20261102'
read, save, sha = frozen.read, frozen.save, frozen.sha
PHASES = (0., .5, 1., 1.5)


def initialize(fitting=False):
    frozen.ROOT, frozen.EXP, frozen.ART, frozen.REC = ROOT, EXP, ART, REC
    if fitting:
        frozen.fit_initialize()
    else:
        frozen.initialize(True)


def prepare():
    assert read(TRANSPORT / 'environment-status.json')['main_build_working']
    for folder in (EXP, ART, REC):
        folder.mkdir(parents=True, exist_ok=False)
    for part in ('original-runinfo', 'original-macros', 'bin'):
        (EXP / part).mkdir()
    for part in ('photons', 'voltage/single', 'truth/single', 'figures', 'fits'):
        (ART / part).mkdir(parents=True)
    models = ('template.npz', 'energy-model.json', 'response-parameters.json', 'mean-correction.npz',
              'noise-calibration.json', 'sensor-config.json', 'charge-map-spatial12.npz', 'spatial-profile.npz',
              'shape-correction.npz', 'nominal-influences.npz')
    for name in models:
        shutil.copyfile(PARENT / name, EXP / name)
    shutil.copyfile(PARENT / 'bin/age.so', EXP / 'bin/age.so')
    assert sha(EXP / 'bin/age.so') == read(PARENT / 'build.json')['binary_sha256']
    shutil.copytree(PARENT / 'source', EXP / 'source', ignore=shutil.ignore_patterns('__pycache__'))
    for name in ('alpha_current_response.py', 'alpha_current_response_study.py', 'test_alpha_current_response.py',
                 'alpha_transport_environment.py', 'alpha_transport_environment_study.py'):
        shutil.copyfile(ROOT / 'pileup/analysis' / name, EXP / 'source/pileup/analysis' / name)
    (EXP / 'source/analysis').mkdir()
    (EXP / 'source/analysis/__init__.py').write_text('')
    shutil.copyfile(ROOT / 'analysis/single_particle_sipm.py', EXP / 'source/analysis/single_particle_sipm.py')
    shutil.copytree(TRANSPORT / 'source/g4', EXP / 'source/g4')
    save(EXP / 'model-freeze.json', {str((EXP / name).relative_to(ROOT)): sha(EXP / name) for name in models})
    save(EXP / 'build.json', dict(reused_current_library=str((PARENT / 'bin/age.so').relative_to(ROOT)),
         binary_sha256=sha(EXP / 'bin/age.so'), compiler=read(PARENT / 'build.json')['compiler'],
         rebuild_command=['g++', '-O3', '-std=c++17', '-fPIC', '-shared',
            str((EXP / 'source/pileup/analysis/sipm_age_mean.cpp').relative_to(ROOT)),
            str((EXP / 'source/pileup/analysis/sipm_age_mean_jacobian.cpp').relative_to(ROOT)),
            '-o', str((EXP / 'bin/age.so').relative_to(ROOT))]))
    mothers = []
    native = sorted(read(TRANSPORT / 'manifest.json')['runs'], key=lambda r: r['id'])
    sensor = SensorConfig(**read(EXP / 'sensor-config.json')['sensor'])
    for environment in ('historical', 'current'):
        for row in native:
            key = 'reference_root' if environment == 'historical' else 'root'
            expected = row['reference_sha256'] if environment == 'historical' else row['root_sha256']
            root = ROOT / row[key]
            assert sha(root) == expected
            info, primary, photons, voltage, scalars = load_root(root)
            assert json.loads(info['sipm_config_json'])['sensor'] == read(EXP / 'sensor-config.json')['sensor']
            assert json.loads(info['sipm_config_json'])['readout'] == read(EXP / 'sensor-config.json')['readout']
            hits = read_photon_hits(root)
            active = map_cells(hits, sensor) >= 0
            first = float(hits.time_ns[active].min())
            i = len(mothers)
            # Match L93: retain outside-area photons; only the alignment uses active hits.
            np.savez_compressed(ART / 'photons' / f'{i:04d}.npz', time_ns=hits.time_ns-first,
                                x_mm=hits.x_mm, y_mm=hits.y_mm, energy_eV=hits.energy_eV)
            save(EXP / 'original-runinfo' / f'{i:04d}.json', info)
            macro = EXP / 'original-macros' / f'{i:04d}.mac'
            macro.write_text(info['macro_contents'])
            mothers.append(dict(id=i, pair=row['id'], group=row['group'], environment=environment,
                root=row[key], root_sha256=expected, seed=row['seed'], source_energy_MeV=row['source_energy_MeV'],
                edep_MeV=primary['edep_total_MeV'], first_active_ns=first,
                photons=len(hits.time_ns), active_photons=int(active.sum()),
                macro=str(macro.relative_to(ROOT)), historical_model_exposure=True,
                independent_new_mother=False))
    tasks = []
    for mother in mothers:
        for phase in PHASES:
            tasks.append(dict(id=len(tasks), kind='single', parents=[mother['id']], times_ns=[128+phase],
                phase_ns=phase, seed_base=739000000+mother['pair']*100, repeats=8,
                environment=mother['environment'], group=mother['group'], pair=mother['pair']))
    save(EXP / 'mothers.json', mothers)
    save(EXP / 'single-plan.json', tasks)
    save(EXP / 'protocol.json', dict(created_utc=datetime.now(timezone.utc).isoformat(),
         base_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
         git_status=subprocess.check_output(['git', 'status', '--short'], cwd=ROOT, text=True),
         scope='Engineering transfer controls only; known K=1 and first-active-photon time; no new G4 or fresh independent mothers',
         input='16 fixed historical/current same-seed pairs, four source energies x4; historical arm uses L93 training mothers',
         models='Copy unchanged L93 old spatial12 anchor and area-preserving mean correction; no retraining or threshold selection',
         primary_gate=dict(current_each_group_phase_abs_mean_energy_bias_percent=5.,
             corrected_to_old_individual_energy_rmse_ratio=1.1, corrected_phase_bias_span_pp=1.,
             max_nonconverged=0, native_readout_audit='Existing exact charge/current/voltage and half-dt convergence checks'),
         gate_scope='Inherited L93 engineering tolerances; n=4 mothers/group cannot establish generalization or statistical equivalence',
         secondary='Describe paired environment bias shifts and phase0 strong weak-direction projections; no small-sample projection gate',
         projection='Fixed nominal old-model weak-logE operators at gaps8/40/200ns; 100delta ln E, not nonlinear pair energy error',
         readout='Native stored avalanches -> Python readout, compare same arrays; separate from fresh stochastic offline avalanche simulation',
         observation='4096ns offline; first2048ns/1024samples;2nsperiod,.5nsaperture; firstactive at128+phase; keep all photons before sensor rejection',
         phases_ns=list(PHASES), repeats=8, formal_digitizations=1024, formal_fits=2304,
         seeds='739000000+pairID*100+repeat; common seed across environment/phase does not ensure identical avalanche streams; readout+10000000',
         uncertainty='4000 paired mother-cluster bootstrap draws with seed739900+group*10+methodID; four mothers/group, descriptive only',
         workers=8, BLAS_threads=1, pilot_tasks=[64,112], pilot_repeats=2, cost_limit_seconds=900.,
         default_changed=False, fresh_campaign='336strong+8pilot remains separate and unrun; engineering pass only permits planning current-env transfer validation'))
    shutil.copyfile(TRANSPORT / 'environment.json', EXP / 'transport-environment.json')
    shutil.copyfile(PARENT / 'literature-evidence.json', EXP / 'inherited-literature-evidence.json')
    for name in ('README.md', 'docs/literature-use.md', 'jupyter/codex_0007_SiPM_Energy_Calibration.ipynb',
                 'jupyter/Research_Plan_20260921.ipynb', 'jupyter/Run_Catalog.ipynb'):
        p = ROOT / name
        shutil.copyfile(p, REC / (p.name+'.before'))
    (EXP / 'run.py').write_text("from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\n"
        "os.environ['GAGG_CURRENT_RESPONSE_ROOT']=str(f.parents[2])\nos.environ['OPENBLAS_NUM_THREADS']='1'\nos.environ['OMP_NUM_THREADS']='1'\n"
        "sys.path.insert(0,str(f/'source'))\nif __name__=='__main__':\n    from pileup.analysis.alpha_current_response_study import main\n    main()\n")
    paths = [p for p in EXP.rglob('*') if p.is_file() and 'bin' not in p.relative_to(EXP).parts]
    save(EXP / 'design-freeze.json', {str(p.relative_to(ROOT)): sha(p) for p in paths})
    paths = [ROOT / m['root'] for m in mothers] + list((ART / 'photons').glob('*.npz'))
    save(EXP / 'input-freeze.json', {str(p.relative_to(ROOT)): sha(p) for p in paths})
    print('Frozen 16 environment pairs, 128 tasks, 1024 digits and 2304 known-K/t fits', flush=True)


def verify():
    counts = {}
    for name in ('design-freeze.json', 'input-freeze.json', 'model-freeze.json'):
        data = read(EXP / name)
        for path, digest in data.items():
            assert sha(ROOT / path) == digest, path
        counts[name] = len(data)
    assert sha(EXP / 'bin/age.so') == read(EXP / 'build.json')['binary_sha256']
    return counts


def native_one(mother):
    from analysis.single_particle_sipm import load_root as load_native, audit_response
    started = time.monotonic()
    out = audit_response(load_native(ROOT / mother['root']))
    return dict(id=mother['id'], environment=mother['environment'], pair=mother['pair'],
                seconds=time.monotonic()-started, **out)


def audit():
    verify()
    assert not (EXP / 'native-audit.json').exists()
    with ProcessPoolExecutor(max_workers=8) as pool:
        rows = list(pool.map(native_one, read(EXP / 'mothers.json')))
    save(EXP / 'native-audit.json', dict(passed=True, rows=rows, roots=32,
         new_digitizations=0, description='Two readout evaluations per stored native avalanche record; no resampled avalanche simulation'))
    print('Native avalanche/charge, Python readout and half-dt controls passed on 32 ROOTs', flush=True)


def generate(stage):
    verify()
    tasks = read(EXP / 'single-plan.json')
    if stage == 'pilot':
        tasks = [dict(tasks[i], repeats=2) for i in (64, 112)]
        output = 'pilot'
    else:
        assert read(EXP / 'pilot-gate.json')['passed'] and read(EXP / 'native-audit.json')['passed']
        output = 'formal'
    assert not (EXP / (stage+'-generation.json')).exists()
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=8, initializer=initialize) as pool:
        rows = []
        for row in pool.map(frozen.simulate, [(task, output) for task in tasks]):
            rows.append(row)
            if len(rows) % 16 == 0:
                print(stage, len(rows), '/', len(tasks), flush=True)
    save(EXP / (stage+'-generation.json'), dict(rows=rows, seconds=time.monotonic()-started,
                                               digitizations=sum(r['digitizations'] for r in rows)))
    if stage == 'pilot':
        estimate = max(r['seconds'] for r in rows)/2*1024/8*1.5
        save(EXP / 'pilot-gate.json', dict(passed=estimate < 900, projected_seconds=estimate))
        assert estimate < 900, estimate
    else:
        for i in (64, 112):
            for part in ('truth', 'voltage'):
                with np.load(ART / f'pilot/{part}/single/{i:04d}.npz') as a, np.load(ART / f'{part}/single/{i:04d}.npz') as b:
                    for key in a.files:
                        assert np.array_equal(a[key], b[key][:2]), (i, part, key)
        paths = list((ART / 'voltage/single').glob('*.npz'))+list((ART / 'truth/single').glob('*.npz'))
        save(EXP / 'digitization-freeze.json', {str(p.relative_to(ROOT)): sha(p) for p in paths})
    print(stage, 'completed', round(time.monotonic()-started, 2), 's', flush=True)


def fit():
    verify()
    for path, digest in read(EXP / 'digitization-freeze.json').items():
        assert sha(ROOT / path) == digest
    assert not list((ART / 'fits').glob('*.json'))
    tasks = read(EXP / 'single-plan.json')
    with ProcessPoolExecutor(max_workers=2, initializer=initialize, initargs=(True,)) as pool:
        pilot = list(pool.map(frozen.fit_one, [(i, tasks[i]['times_ns'][0], 'fit-pilot') for i in (64, 112)]))
    estimate = max(r['seconds'] for r in pilot)*len(tasks)/8*1.5
    save(EXP / 'fit-pilot-gate.json', dict(passed=estimate < 900, projected_seconds=estimate, rows=pilot))
    assert estimate < 900, estimate
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=8, initializer=initialize, initargs=(True,)) as pool:
        rows = []
        for row in pool.map(frozen.fit_one, [(t['id'], t['times_ns'][0], 'fits') for t in tasks]):
            rows.append(row)
            if len(rows) % 16 == 0:
                print('fit', len(rows), '/', len(tasks), flush=True)
    for i in (64, 112):
        a, b = read(ART / 'fit-pilot' / f'{i:04d}.json'), read(ART / 'fits' / f'{i:04d}.json')
        assert [{k: v for k, v in r.items() if k != 'seconds'} for r in a['rows']] == [
            {k: v for k, v in r.items() if k != 'seconds'} for r in b['rows']]
    save(EXP / 'fits.json', dict(fits=sum(r['fits'] for r in rows), pilot_fits=sum(r['fits'] for r in pilot),
                                rows=rows, seconds=time.monotonic()-started, pilot_equal=True))
    save(EXP / 'fit-output-freeze.json', {str(p.relative_to(ROOT)): sha(p) for p in (ART / 'fits').glob('*.json')})


def score():
    verify(); initialize()
    for name in ('digitization-freeze.json', 'fit-output-freeze.json'):
        for path, digest in read(EXP / name).items():
            assert sha(ROOT / path) == digest
    tasks, mothers = read(EXP / 'single-plan.json'), read(EXP / 'mothers.json')
    details, energy_rows, paired, projections, failures = [], [], [], [], []
    arrays = {'time_ns': frozen.OLD.readout.adc_centers}
    for task in tasks:
        mother = mothers[task['parents'][0]]
        fits = read(ART / 'fits' / f"{task['id']:04d}.json")
        assert fits['truth_guard'] and not fits['true_E_input'] and fits['known_K1_t']
        for fit in fits['rows']:
            details.append(dict(environment=mother['environment'], pair=mother['pair'], group=mother['group'],
                phase_ns=task['phase_ns'], method=fit['method'], repeat=fit['repeat'],
                relative_error_percent=100*(fit['energy_MeV'][0]/mother['edep_MeV']-1),
                converged=fit['converged'], boundary=fit['boundary']))
    for env in ('historical', 'current'):
        for group in range(4):
            for phase in PHASES:
                for method in ('old', 'corrected'):
                    for kind in ('individual', 'mean'):
                        chosen = [r for r in details if r['environment'] == env and r['group'] == group and
                                  r['phase_ns'] == phase and r['method'] == method and
                                  (r['repeat'] >= 0 if kind == 'individual' else r['repeat'] == -1)]
                        row = energy_metrics([r['relative_error_percent'] for r in chosen],
                                             [r['boundary'] for r in chosen], [r['converged'] for r in chosen])
                        energy_rows.append(dict(environment=env, group=group, phase_ns=phase, method=method,
                                                kind=kind, mothers=4, **row))
    for group in range(4):
        for phase in PHASES:
            def select(method, kind):
                return next(r for r in energy_rows if r['environment'] == 'current' and r['group'] == group
                            and r['phase_ns'] == phase and r['method'] == method and r['kind'] == kind)
            for method in ('old', 'corrected'):
                if abs(select(method, 'mean')['bias_percent']) > 5:
                    failures.append(f'group{group}/phase{phase}/{method}: mean energy bias >5%')
            if select('corrected', 'individual')['rmse_percent'] > 1.1*select('old', 'individual')['rmse_percent']:
                failures.append(f'group{group}/phase{phase}: corrected individual RMSE >1.1 baseline')
        span = np.ptp([r['bias_percent'] for r in energy_rows if r['environment'] == 'current' and
                       r['group'] == group and r['method'] == 'corrected' and r['kind'] == 'mean'])
        if span > 1:
            failures.append(f'group{group}: phase bias span >1pp')
        for method_id, method in enumerate(('old', 'corrected')):
            values = []
            for env in ('historical', 'current'):
                values.append(np.array([[r['relative_error_percent'] for r in details if r['environment'] == env
                    and r['pair'] == pair and r['phase_ns'] == 0 and r['method'] == method and r['repeat'] >= 0]
                    for pair in range(group*4, group*4+4)]))
            paired.append(dict(group=group, method=method, units='energy bias percentage points; current-historical',
                 **paired_mother_interval(*values, seed=739900+group*10+method_id)))
    if any(r['nonconverged'] for r in energy_rows if r['environment'] == 'current'):
        failures.append('Current-environment unconverged energy fit')
    with np.load(EXP / 'nominal-influences.npz') as data:
        operator = data['high'][1:]
    wave_rows = []
    for env in ('historical', 'current'):
        for group in range(4):
            residuals = {'old': [], 'corrected': []}; actual = []
            for mother in [r for r in mothers if r['environment'] == env and r['group'] == group]:
                task_id = mother['id']*4
                with np.load(ART / 'truth/single' / f'{task_id:04d}.npz') as data:
                    mean = data['components_mV'].mean(axis=0)[0]
                actual.append(mean)
                for method, model in [('old', frozen.OLD), ('corrected', frozen.NEW)]:
                    residuals[method].append(mean-model.isolated([128.], [mother['edep_MeV']])[0][0])
            arrays[f'{env}_{group}_actual_mV'] = np.array(actual)
            for method in residuals:
                residual = np.array(residuals[method])
                arrays[f'{env}_{group}_{method}_residual_mV'] = residual
                wave_rows.append(dict(environment=env, group=group, method=method,
                     mean_waveform_mse_mV2=float(np.mean(residual.mean(0)**2)), mother_mse_mV2=float(np.mean(residual**2))))
                if group == 3:
                    value = 100*residual@operator.T
                    arrays[f'{env}_{method}_projections_100logE'] = value
                    projections.extend(dict(environment=env, method=method, gap_ns=gap, mothers=4,
                       mean=float(value[:, j].mean()), rms=float(np.sqrt(np.mean(value[:, j]**2))))
                       for j, gap in enumerate((8., 40., 200.)))
    np.savez_compressed(ART / 'summary.npz', **arrays)
    save(EXP / 'fit-details.json', details)
    save(EXP / 'result-summary.json', dict(energy=energy_rows, paired_environment_shift=paired,
        waveform=wave_rows, projection=projections, formal_digitizations=1024, pilot_digitizations=4,
        formal_fits=2304, pilot_fits=36, native_audits=32, new_ROOTs=0, new_independent_mothers=0))
    save(EXP / 'gate.json', dict(passed=not failures, failures=failures,
         scope='Small exposed-sample engineering checks, not independent model confirmation',
         next='Separate current-environment fresh-mother transfer confirmation' if not failures else
              'Stop fresh campaign; retain failures and diagnose without retuning on these controls'))
    save(EXP / 'result-freeze.json', {str(p.relative_to(ROOT)): sha(p) for p in
          [EXP / 'result-summary.json', EXP / 'fit-details.json', EXP / 'gate.json', ART / 'summary.npz']})
    print(json.dumps(read(EXP / 'gate.json'), ensure_ascii=False, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['prepare', 'verify', 'audit', 'pilot', 'single', 'fit', 'score'])
    stage = parser.parse_args().stage
    if stage == 'prepare': prepare()
    elif stage == 'verify': print(json.dumps(verify()))
    elif stage == 'audit': audit()
    elif stage in ('pilot', 'single'): generate(stage)
    elif stage == 'fit': fit()
    elif stage == 'score': score()


if __name__ == '__main__':
    main()
