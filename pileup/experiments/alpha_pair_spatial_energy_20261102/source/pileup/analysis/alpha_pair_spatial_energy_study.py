"""Frozen paired pilot of uniform/spatial known-K/t energy inference."""
from __future__ import annotations

import argparse
import ast
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

from .alpha_pair_causal_controls import compatibility, read, save, sha
from .alpha_pair_causal_recovery import CausalRecovery
from .alpha_pair_digitizer_optimized import simulate_avalanches
from .alpha_pair_energy_calibrated import CalibratedEnergyResponse
from .alpha_pair_recovery_energy import RecoveryEnergyResponse
from .alpha_pair_spatial_recovery import SpatialRecovery
from .alpha_pair_spatial_energy import common_inputs, fit_energy
from .sipm_readout import ReadoutConfig, synthesize_readout
from .sipm_response import SensorConfig, PhotonHits, AvalancheResult

NAME = 'alpha_pair_spatial_energy_20261102'
ROOT = Path(os.environ.get('GAGG_SPATIAL_ENERGY_ROOT', Path(__file__).resolve().parents[2])).resolve()
EXP = ROOT / 'pileup/experiments' / NAME
ART = ROOT / 'jupyter/artifacts' / NAME
REC = ROOT / 'jupyter/run_records' / NAME
PREV = ROOT / 'pileup/experiments/alpha_pair_spatial_recovery_20261102'
PREVART = ROOT / 'jupyter/artifacts/alpha_pair_spatial_recovery_20261102'
METHODS = ('uniform', 'spatial12')
WORKERS = 8


def snapshot():
    source = EXP / 'source/pileup/analysis'
    source.mkdir(parents=True)
    (source.parent / '__init__.py').write_text('')
    (source / '__init__.py').write_text('')
    pending, copied = ['alpha_pair_spatial_energy_study'], set()
    while pending:
        name = pending.pop()
        if name in copied:
            continue
        path = ROOT / 'pileup/analysis' / (name + '.py')
        tree = ast.parse(path.read_text())
        shutil.copyfile(path, source / path.name)
        copied.add(name)
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.level == 1:
                pending.extend([node.module] if node.module else [n.name for n in node.names])
    for name in ('sipm_age_mean.cpp', 'sipm_age_mean_jacobian.cpp'):
        shutil.copyfile(ROOT / 'pileup/analysis' / name, source / name)
    (EXP / 'run.py').write_text("from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\n"
        "os.environ['GAGG_SPATIAL_ENERGY_ROOT']=str(f.parents[2])\n"
        "os.environ['OPENBLAS_NUM_THREADS']='1'\nos.environ['OMP_NUM_THREADS']='1'\n"
        "sys.path.insert(0,str(f/'source'))\n"
        "if __name__=='__main__':\n    from pileup.analysis.alpha_pair_spatial_energy_study import main\n    main()\n")
    (EXP / 'bin').mkdir()
    command = ['g++', '-O3', '-std=c++17', '-fPIC', '-shared',
        str((source / 'sipm_age_mean.cpp').relative_to(ROOT)),
        str((source / 'sipm_age_mean_jacobian.cpp').relative_to(ROOT)),
        '-o', str((EXP / 'bin/age.so').relative_to(ROOT))]
    subprocess.run(command, cwd=ROOT, check=True)
    save(EXP / 'build.json', dict(command=command, binary_sha256=sha(EXP / 'bin/age.so'),
        compiler=subprocess.check_output(['g++', '--version'], text=True).splitlines()[0]))


def prepare():
    compatibility()
    assert read(PREV / 'single-gate.json')['passed'] and read(PREV / 'shared-gate.json')['passed']
    for p in (EXP, ART, REC):
        p.mkdir(parents=True, exist_ok=False)
    for name in ('photons', 'waveforms', 'diagnostics'):
        (ART / name).mkdir()
    configs = ('template.npz', 'energy-model.json', 'response-parameters.json', 'mean-correction.npz',
        'noise-calibration.json', 'sensor-config.json', 'charge-map-1.0.npz',
        'charge-map-spatial12.npz', 'spatial-profile.npz', 'evaluation-mothers.json',
        'mother-split.json', 'profile-summary.json')
    for name in configs:
        shutil.copyfile(PREV / name, EXP / name)
    parents = read(EXP / 'evaluation-mothers.json')
    truth = {r['parent']: r['energy_dep_MeV'] for r in parents}
    pairs = read(PREV / 'pair-plan.json')
    singles = [r for r in read(PREV / 'single-plan.json') if r['phase_ns'] == 0]
    cases = []
    for row in singles:
        cid = len(cases)
        cases.append(dict(id=cid, kind='single', condition=row['condition'], block=row['block'],
            parents=[row['parent']], times_ns=[row['time_ns']], source_id=row['id'], gap_ns=None,
            seed_base=row['seed_base'], energy_dep_MeV=[truth[row['parent']]]))
        with np.load(PREVART / 'single/single' / f"{row['id']:04d}.npz") as f:
            np.savez_compressed(ART / 'waveforms' / f'{cid:04d}.npz', voltage_mV=f['analog_mV'])
    for row in pairs:
        cases.append(dict(id=len(cases), kind='pair', condition=row['condition'], block=row['block'],
            parents=row['parents'], times_ns=row['times_ns'], source_id=row['id'], gap_ns=row['gap_ns'],
            seed_base=row['seed_base'], energy_dep_MeV=[truth[p] for p in row['parents']]))
    for parent in sorted({p for c in cases if c['kind']=='pair' for p in c['parents']}):
        name = f'{parent:04d}.npz'
        shutil.copyfile(PREVART / 'photons' / name, ART / 'photons' / name)
    save(EXP / 'cases.json', cases)
    tasks = []
    for c in cases:
        for repeat in [-1] + (list(range(8)) if c['kind']=='single' else [0, 1]):
            for method in METHODS:
                tasks.append(dict(id=len(tasks), case=c['id'], kind=c['kind'], repeat=repeat, method=method,
                    times_ns=c['times_ns'], waveform=str((ART / 'waveforms' / f"{c['id']:04d}.npz").relative_to(ROOT))))
    save(EXP / 'fit-plan.json', tasks)
    save(EXP / 'protocol.json', dict(created_utc=datetime.now(timezone.utc).isoformat(),
        commit=subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip(),
        git_status=subprocess.check_output(['git','status','--short'], cwd=ROOT, text=True),
        scope='Fixed pilot, pure alpha, known true K and t; energy is unknown. Same historical mothers as the response check, not independent confirmation.',
        inputs='16 phase-zero single cases x8 saved noisy analog repeats. Recreate20 pair cases x8 repeats with original photons/seeds, adding complete dark/electronic-noise analog voltage; compare all saved ancestry arrays exactly.',
        primary='40 individual noisy pairs:5conditions x4mothers xfirst2repeats, fixed before output.20 eight-repeat mean pairs diagnose bias.128 individual singles and16 means are controls.',
        sample_scope='This study ends at this fixed pilot; no outcome-dependent sample expansion. Repeats are not independent optical mothers.',
        cost_rationale='Pre-freeze synthetic spatial fit took96.15s for3starts/31evaluations. Fixed pilot first; at most8workers with BLAS single-threaded.',
        model='Inherit unchanged uniform and primary12x12,1ns response,128-mother D4 prior, isolated anchor, own charge maps and all sensor parameters from preceding freeze.',
        weights='One common learned_metric from isolated NNLS voltage seed, times and inherited noise parameters. Fixed for both fits, not refitted or promoted to a calibrated spatial likelihood.',
        search='Identical log-energy bounds, NNLS start plus charge shares0.1/0.9; least_squares analytic gradient,max_nfev80 per start,ftol/xtol1e-8,gtol1e-5. Lowest objective selected including nonconverged candidates, all failures retained.',
        single_solver='Identical bounded scalar fit, xatol1e-10,maxiter150; interaction absent, hence models must be identical.',
        metrics='Individual E_dep bias/RMSE/SD, pair both-E-within10% AND converged, boundaries/nonconvergence, paired gain and cost differences. Known K/t gives no blind detection/merge/split efficiency.',
        uncertainty='10000 paired mother-block resamples; noisy endpoints additionally resample two repeats within each selected mother. Conditional on frozen prior/model; descriptive four-mother intervals only.',
        gates=dict(numerical_gradient_relative=3e-4, mean_component_abs_mV=1e-8, inverse_relative=.001,
                   single_group_mean_bias_percent=5., pilot_projected_seconds=3600.,
                   next_stage_short_strict_gain_pp=10., next_stage_short_rmse_ratio=.95,
                   next_stage_each_short_rmse_ratio=1.1, next_stage_short_mean_abs_bias_percent=10.),
        decision='Recommendation for later independent confirmation requires all next_stage gates on HW8/40; failure stops expansion. No default replacement in either outcome.',
        photon_truth_boundary='Inference jobs contain only waveform path, known times, method and row identifiers. No E, photon positions, ancestry, charges or response truth supplied; audited file opens reject label/diagnostic/photon inputs.',
        workers=WORKERS, BLAS_threads=1, new_ROOT=0, default_changed=False,
        literature='Project L88 and frozen code; paper-search configured library unavailable. No new external paper claim.'))
    import scipy
    save(EXP / 'environment.json', dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__))
    for path in [ROOT/'README.md', ROOT/'docs/literature-use.md', *[ROOT/'jupyter'/p for p in
        ['codex_0007_SiPM_Energy_Calibration.ipynb','Research_Plan_20260921.ipynb','Run_Catalog.ipynb']]]:
        shutil.copyfile(path, REC/(path.name+'.before'))
    snapshot()
    files = [p for p in EXP.rglob('*') if p.is_file() and 'bin' not in p.relative_to(EXP).parts]
    save(EXP/'freeze.json', {str(p.relative_to(ROOT)):sha(p) for p in files})
    save(EXP/'parent-input-freeze.json', {str(p.relative_to(ROOT)):sha(p) for p in [
        *list((ART/'photons').glob('*.npz')), *list((ART/'waveforms').glob('*.npz')),
        *list((PREVART/'pair/pair').glob('*.npz'))]})
    print('Frozen:36cases,408fits including means,40individual pairs; no truth E in fit jobs',flush=True)


def verify(inputs=False):
    names = ['freeze.json', 'parent-input-freeze.json'] + (['input-freeze.json'] if inputs else [])
    count = 0
    for name in names:
        for path, digest in read(EXP/name).items():
            assert sha(ROOT/path) == digest, path
            count += 1
    return count


def initialize():
    compatibility()
    global CONFIG, ANCHOR, MODELS, NOISE
    CONFIG = read(EXP/'sensor-config.json')
    ANCHOR = CalibratedEnergyResponse.load(RecoveryEnergyResponse.load(EXP), EXP/'mean-correction.npz')
    with np.load(EXP/'charge-map-1.0.npz') as f:
        mapping = {k:f[k].copy() for k in f.files}
    MODELS = {'uniform':CausalRecovery(ANCHOR,CONFIG,EXP/'bin/age.so',mapping)}
    with np.load(EXP/'spatial-profile.npz') as f:
        profile = f['profile_12'].copy()
    with np.load(EXP/'charge-map-spatial12.npz') as f:
        mapping = {k:f[k].copy() for k in f.files}
    MODELS['spatial12'] = SpatialRecovery(ANCHOR,CONFIG,EXP/'bin/age.so',mapping,profile)
    NOISE = read(EXP/'noise-calibration.json')


def make_pair(case):
    sensor, ro = SensorConfig(**CONFIG['sensor']), ReadoutConfig(**CONFIG['readout'])
    fields = {k:[] for k in ('time_ns','x_mm','y_mm','energy_eV')}
    owners = []
    for j,(parent,t0) in enumerate(zip(case['parents'],case['times_ns'])):
        with np.load(ART/'photons'/f'{parent:04d}.npz') as f:
            times = f['time_ns']+t0
            keep = times<4096.
            for key in fields:
                fields[key].append((times if key=='time_ns' else f[key])[keep])
        owners.append(np.full(int(keep.sum()),j))
    arrays = {k:np.concatenate(v) for k,v in fields.items()}
    owner = np.concatenate(owners)
    order = np.argsort(arrays['time_ns'],kind='stable')
    owner = owner[order]
    hits = PhotonHits(**{k:v[order] for k,v in arrays.items()},window_ns=4096.)
    voltage, signals, charges = [],[],[]
    for repeat in range(8):
        seed = case['seed_base']+repeat
        av = simulate_avalanches(hits,sensor,seed)
        voltage.append(synthesize_readout(av,ro,seed+10000000).arrays['analog_mV'])
        records = av.records
        ancestry = np.full(len(records),-1,int)
        direct = records['kind']==0
        ancestry[direct] = owner[records['photon_index'][direct]]
        for j in np.flatnonzero(records['parent_id']>=0):
            ancestry[j] = ancestry[records['parent_id'][j]]
        curves,q = [],[]
        for j in range(2):
            sub = records[ancestry==j]
            result = AvalancheResult(sub,4096.,av.full_charge_C,np.array([]),{})
            curves.append(synthesize_readout(result,ro,0).arrays['signal_mV'][:1024])
            q.append(float(sub['charge_C'].sum()/av.full_charge_C))
        signals.append(curves);charges.append(q)
    voltage,signals,charges = map(np.asarray,(voltage,signals,charges))
    with np.load(PREVART/'pair/pair'/f"{case['source_id']:04d}.npz") as f:
        assert np.array_equal(signals,f['on_shared']) and np.array_equal(charges,f['on_shared_charge'])
    assert np.isfinite(voltage).all()
    np.savez_compressed(ART/'waveforms'/f"{case['id']:04d}.npz",voltage_mV=voltage[:,:1024])
    np.savez_compressed(ART/'diagnostics'/f"{case['id']:04d}.npz",full_analog_mV=voltage,
                        event_signal_mV=signals,event_charge_pe=charges)
    return dict(case=case['id'],digitizations=8,ancestry_exact=True,full_samples=voltage.shape[1])


def generate():
    verify()
    assert not (EXP/'input-freeze.json').exists()
    cases = [c for c in read(EXP/'cases.json') if c['kind']=='pair']
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=WORKERS,initializer=initialize) as pool:
        rows = list(pool.map(make_pair,cases))
    save(REC/'digitization.json',dict(rows=rows,digitizations=160,seconds=time.monotonic()-started,
        ancestry_exact=True,new_ROOT=0,reused_single_waveforms=128))
    paths = [*list((ART/'waveforms').glob('*.npz')),*list((ART/'diagnostics').glob('*.npz'))]
    save(EXP/'input-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})
    print('160 complete pair digitizations;20 ancestry controls exact',flush=True)


def inverse_one(job):
    name,times,energy = job
    model = MODELS[name]
    waveform = model.components(times,energy).sum(0)+model.background
    common = common_inputs(MODELS['uniform'],waveform,times,NOISE)
    result = fit_energy(model,waveform,times,common)
    result.update(model=name, times_ns=times, true_synthetic_MeV=energy,
                  max_relative_error=float(np.max(abs(np.array(result['energy_MeV'])/energy-1))))
    return result


def numerical():
    verify();initialize()
    checks = []
    for name,model in MODELS.items():
        for times,energy in [([128.,136.],[7.,.25]),([128.,168.],[.25,7.])]:
            loge = np.log(energy);p = model.prepare(times)
            mean,grad = model.mean_and_gradient(times,loge,p,strength=1.)
            component_error = float(abs(mean-model.components(times,energy).sum(0)).max())
            direction = np.array([.7,-.3]);h=2e-5
            finite = (model.components(times,np.exp(loge+h*direction)).sum(0)-
                      model.components(times,np.exp(loge-h*direction)).sum(0))/(2*h)
            error = float(np.linalg.norm(grad.T@direction-finite)/np.linalg.norm(finite))
            checks.append(dict(model=name,times_ns=times,mean_component_abs_mV=component_error,gradient_relative=error))
    jobs = [(name,t,e) for name in METHODS for t,e in [
        ([128.],[.25]),([128.],[7.]),([128.,136.],[7.,.25]),([128.,168.],[.25,7.])]]
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=WORKERS,initializer=initialize) as pool:
        inverses = list(pool.map(inverse_one,jobs))
    gates = read(EXP/'protocol.json')['gates']
    passed = all(r['mean_component_abs_mV']<gates['mean_component_abs_mV'] and
                 r['gradient_relative']<gates['numerical_gradient_relative'] for r in checks)
    passed = passed and all(r['converged'] and r['max_relative_error']<gates['inverse_relative'] for r in inverses)
    save(REC/'numerical.json',dict(passed=bool(passed),checks=checks,inverses=inverses,seconds=time.monotonic()-started))
    print('Energy gradient and inverse controls',passed,flush=True)
    assert passed


def inference_initialize():
    initialize()
    global OPENED
    OPENED = set()
    forbidden = [ART/'photons',ART/'diagnostics',EXP/'cases.json',EXP/'evaluation-mothers.json',PREVART]
    def audit(event,args):
        if event!='open' or not isinstance(args[0],(str,bytes,os.PathLike)):
            return
        path = Path(os.fsdecode(args[0])).resolve()
        if any(path==p or p in path.parents for p in forbidden):
            raise PermissionError('Inference attempted to open labels/photons/ancestry: '+str(path))
        if ART in path.parents:
            OPENED.add(str(path.relative_to(ROOT)))
    sys.addaudithook(audit)


def fit_one(task):
    assert set(task)=={'id','case','kind','repeat','method','times_ns','waveform'}
    OPENED.clear()
    with np.load(ROOT/task['waveform']) as f:
        assert f.files==['voltage_mV']
        waveform = f['voltage_mV'].mean(0) if task['repeat']==-1 else f['voltage_mV'][task['repeat']]
    common = common_inputs(MODELS['uniform'],waveform,task['times_ns'],NOISE)
    result = fit_energy(MODELS[task['method']],waveform,task['times_ns'],common)
    result.update(task, read_inputs=sorted(OPENED), true_energy_input=False,
                  metric_signature=sha_array(np.r_[common['metric'].factor.ravel(),common['metric'].modes.ravel()]))
    return result


def sha_array(array):
    import hashlib
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def fit_stage(mode,output=None,ids=None):
    verify(True)
    assert read(REC/'numerical.json')['passed']
    tasks = read(EXP/'fit-plan.json')
    if mode=='single':
        tasks = [t for t in tasks if t['kind']=='single']
    elif mode=='pilot':
        tasks = [t for t in tasks if t['case']==17 and t['repeat']==0]
    elif mode=='pairs':
        assert read(EXP/'single-gate.json')['passed'] and read(REC/'pilot.json')['passed']
        tasks = [t for t in tasks if t['kind']=='pair']
    elif mode=='replay':
        assert output and ids
        tasks = [t for t in tasks if t['id'] in ids]
        assert {t['id'] for t in tasks}==set(ids)
    folder = Path(output) if output else ART/mode
    folder.mkdir(parents=True,exist_ok=False)
    started = time.monotonic();rows=[]
    with ProcessPoolExecutor(max_workers=WORKERS,initializer=inference_initialize) as pool:
        futures = [pool.submit(fit_one,t) for t in tasks]
        for future in as_completed(futures):
            row = future.result()
            save(folder/f"{row['id']:04d}.json",row);rows.append(row)
            print(mode,len(rows),'/',len(tasks),'elapsed',round(time.monotonic()-started,1),flush=True)
    record = dict(tasks=len(tasks),seconds=time.monotonic()-started,workers=WORKERS)
    if mode=='pilot':
        projected = float(sum(r['seconds'] for r in rows)*60/WORKERS*1.5)
        record.update(projected_seconds=projected,limit_seconds=3600.,passed=projected<3600.,rows=rows)
        save(REC/'pilot.json',record)
        assert record['passed']
    elif mode=='replay':
        comparisons=[]
        for row in rows:
            previous=read(ART/('single' if row['kind']=='single' else 'pairs')/f"{row['id']:04d}.json")
            relative=float(np.max(abs(np.array(row['energy_MeV'])/previous['energy_MeV']-1)))
            error=abs(row['cost']-previous['cost'])/max(abs(previous['cost']),1.)
            comparisons.append(dict(id=row['id'],energy_relative=relative,cost_relative=error,
                convergence_equal=row['converged']==previous['converged'],passed=relative<1e-10 and error<1e-10 and row['converged']==previous['converged']))
        record.update(comparisons=comparisons,passed=all(r['passed'] for r in comparisons))
        save(folder/'replay.json',record);assert record['passed']
    else:
        save(REC/(mode+'-run.json'),record)


def main():
    compatibility()
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['prepare','generate','numerical','single','pilot','pairs','replay'])
    parser.add_argument('--output');parser.add_argument('--ids',type=int,nargs='+')
    args=parser.parse_args()
    if args.mode in ('prepare','generate','numerical'):
        globals()[args.mode]()
    else:
        fit_stage(args.mode,args.output,args.ids)


if __name__=='__main__':
    main()
