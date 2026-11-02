"""Frozen residual, weighting and search diagnostics on 20 saved pair means.

All use true K/t. Projection, residual deletion and truth-start fits additionally
use true energy/ancestor signals. None is a blind reconstruction result.
"""
from __future__ import annotations

import argparse
import ast
from concurrent.futures import ProcessPoolExecutor, as_completed
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

from .alpha_pair_causal_recovery import CausalRecovery
from .alpha_pair_energy_calibrated import CalibratedEnergyResponse
from .alpha_pair_recovery_energy import RecoveryEnergyResponse
from .alpha_pair_spatial_recovery import SpatialRecovery
from .alpha_pair_spatial_energy import common_inputs, fit_energy
from .alpha_pair_residual_projection import (TERMS, WEIGHTS, DiagonalMetric,
    fixed_metrics, residual_terms, project, paired_bootstrap, symmetric_change)

NAME = 'alpha_pair_residual_projection_20261102'
ROOT = Path(os.environ.get('GAGG_RESIDUAL_ROOT', Path(__file__).resolve().parents[2])).resolve()
EXP = ROOT/'pileup/experiments'/NAME
ART = ROOT/'jupyter/artifacts'/NAME
REC = ROOT/'jupyter/run_records'/NAME
PARENT = ROOT/'pileup/experiments/alpha_pair_spatial_energy_20261102'
PARENTART = ROOT/'jupyter/artifacts/alpha_pair_spatial_energy_20261102'
SHAREDART = ROOT/'jupyter/artifacts/alpha_pair_spatial_recovery_20261102'
METHODS = ('uniform', 'spatial12')
WORKERS = 8


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def signature(metric):
    return hashlib.sha256(np.concatenate((metric.factor.ravel(), metric.modes.ravel())).tobytes()).hexdigest()


def snapshot():
    target = EXP/'source/pileup/analysis'
    target.mkdir(parents=True)
    (target.parent/'__init__.py').write_text('')
    (target/'__init__.py').write_text('')
    pending, copied = ['alpha_pair_residual_study'], set()
    while pending:
        name = pending.pop()
        if name in copied:
            continue
        source = ROOT/'pileup/analysis'/f'{name}.py'
        tree = ast.parse(source.read_text())
        shutil.copyfile(source, target/source.name)
        copied.add(name)
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.level == 1:
                pending.extend([node.module] if node.module else [n.name for n in node.names])
    for name in ('sipm_age_mean.cpp', 'sipm_age_mean_jacobian.cpp'):
        shutil.copyfile(ROOT/'pileup/analysis'/name, target/name)
    (EXP/'bin').mkdir()
    command = ['g++','-O3','-std=c++17','-fPIC','-shared',
        str((target/'sipm_age_mean.cpp').relative_to(ROOT)),
        str((target/'sipm_age_mean_jacobian.cpp').relative_to(ROOT)),
        '-o',str((EXP/'bin/age.so').relative_to(ROOT))]
    subprocess.run(command, cwd=ROOT, check=True)
    save(EXP/'build.json', dict(command=command, binary_sha256=sha(EXP/'bin/age.so'),
        compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0]))
    (EXP/'run.py').write_text(
        "from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\n"
        "os.environ['GAGG_RESIDUAL_ROOT']=str(f.parents[2])\n"
        "os.environ['OPENBLAS_NUM_THREADS']='1'\nos.environ['OMP_NUM_THREADS']='1'\n"
        "sys.path.insert(0,str(f/'source'))\n"
        "if __name__=='__main__':\n    from pileup.analysis.alpha_pair_residual_study import main\n    main()\n")


def prepare():
    for directory in (EXP, ART, REC):
        directory.mkdir(parents=True, exist_ok=False)
    for name in ('template.npz','energy-model.json','response-parameters.json','mean-correction.npz',
                 'noise-calibration.json','sensor-config.json','charge-map-1.0.npz',
                 'charge-map-spatial12.npz','spatial-profile.npz'):
        shutil.copyfile(PARENT/name, EXP/name)
    cases = [c for c in read(PARENT/'cases.json') if c['kind']=='pair']
    save(EXP/'cases.json', cases)
    tasks = []
    def add(case, method, weight='full', delete=None, start='truth'):
        tasks.append(dict(id=len(tasks),case=case['id'],method=method,weight=weight,
                          delete_term=delete,start=start))
    for case in cases:
        for method in METHODS:
            add(case, method)
    short = [c for c in cases if c['condition']=='high_weak' and c['gap_ns'] in (8.,40.)]
    for case in short:
        for weight in ('diagonal','unit_fir'):
            add(case,'spatial12',weight=weight,start='voltage')
        for term in TERMS:
            add(case,'spatial12',delete=term)
    assert len(tasks)==88
    save(EXP/'fit-plan.json',tasks)
    pilot = [t['id'] for t in tasks if t['case']==17 and (
        (t['weight']=='full' and t['delete_term'] is None) or
        t['delete_term']=='strong_isolated' or t['weight']=='diagonal')]
    assert len(pilot)==4
    replay = [t['id'] for t in tasks if t['case']==32 and (
        (t['weight']=='full' and t['delete_term'] is None) or
        t['delete_term']=='shared_correction' or t['weight']=='unit_fir')]
    assert len(replay)==4
    save(EXP/'protocol.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),
        commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        git_status=subprocess.check_output(['git','status','--short'],text=True),
        scope='Truth-assisted diagnostic,20 eight-repeat pair means,4 previously exposed mother pairs,5conditions. No new ROOT or digitization; no independent confirmation.',
        residual_terms=list(TERMS),
        decomposition='Istrong-Astrong; Iweak-Aweak; (sum S-sum I)-(mu-sum A); Y-background-sum S. Exact per-repeat voltage closure. Terms are correlated; shared/isolated streams need not have matching avalanche realizations.',
        projection='SVD of W*T*dmu/dlogE at true E,relative rank cutoff1e-12. Source contributions add in logE. Also bounded linear solve; exp transforms are applied only to total changes.',
        geometry='Profile weak derivative against strong. Weighted derivative geometry is not a calibrated Fisher/CR bound or universal physical resolution limit.',
        weights='Original full GLS; diagonal marginal variance; unit FIR weights. All apply the same inherited FIR transform. Parameters fixed from original mean voltage, not recomputed after deletions. Unit FIR is not raw-voltage OLS.',
        symmetric_split='Half-sum influence times residual change plus half-sum residual times influence change. Exact descriptive order-symmetric algebra, not unique causal attribution.',
        nonlinear='40 true-E-start full fits;16 spatial short-gap weight contrasts with original3 voltage starts;32 spatial short-gap oracle single-term deletions with one true-E start. Fixed logE bounds and max_nfev80/start inherited. No outcome-dependent fit expansion.',
        truth_boundary='Projection and deletion use true E and ancestor signals; oracle search uses true E only as start. Weight contrasts use voltage/known t for inference; truth is used only for scoring. Not blind K/t results.',
        bootstrap='10000 paired mother/repeat resamples of linear contributions,seed731100+condition index;fixed gradients and weights. Descriptive4mother intervals; optical mother variation remains in eight-repeat means.',
        numeric_cost_tolerance='max(1e-5,1e-7*abs(oldcost));within one weighting only,not statistical significance or evidence of global optimality.',
        pilot_ids=pilot,replay_ids=replay,pilot_projected_seconds_limit=3600,
        validation=dict(residual_closure_mV=1e-9,old_cost_relative=1e-7,
                        synthetic_projection=1e-9,small_step_absolute_log=1e-4),
        workers=WORKERS,BLAS_threads=1,new_ROOT=0,new_digitization=0,default_changed=False,
        literature='Reuse project L89 and model derivation. Configured paper library unavailable; no new external paper claim.',
        decision='Describe source/weight/search sensitivities. Do not replace defaults, refit physical parameters or launch new confirmation from these exposed mothers.'))
    import scipy
    save(EXP/'environment.json',dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__))
    for path in [ROOT/'README.md',ROOT/'docs/literature-use.md',*[ROOT/'jupyter'/name for name in
        ('codex_0007_SiPM_Energy_Calibration.ipynb','Research_Plan_20260921.ipynb','Run_Catalog.ipynb')]]:
        shutil.copyfile(path,REC/(path.name+'.before'))
    snapshot()
    paths = [p for p in EXP.rglob('*') if p.is_file() and 'bin' not in p.relative_to(EXP).parts]
    save(EXP/'freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})
    inputs = [PARENT/'cases.json',PARENT/'fit-plan.json',PARENT/'decision.json']
    for case in cases:
        inputs.extend([PARENTART/'waveforms'/f"{case['id']:04d}.npz",
                       PARENTART/'diagnostics'/f"{case['id']:04d}.npz",
                       SHAREDART/'pair/pair'/f"{case['source_id']:04d}.npz"])
        inputs.extend(PARENTART/'pairs'/f"{288+(case['id']-16)*6+j:04d}.json" for j in range(2))
    save(EXP/'input-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in inputs})
    print('Frozen20means,120linear projections,88nonlinear diagnostic fits',flush=True)


def verify():
    counts = {}
    for name in ('freeze.json','input-freeze.json'):
        values = read(EXP/name)
        for path,digest in values.items():
            if sha(ROOT/path)!=digest:
                raise AssertionError(f'Hash mismatch: {path}')
        counts[name]=len(values)
    return counts


def initialize():
    if not hasattr(np,'trapezoid'):
        np.trapezoid=np.trapz
    global ANCHOR, MODELS, NOISE, TRANSFORM, CASES
    ANCHOR=CalibratedEnergyResponse.load(RecoveryEnergyResponse.load(EXP),EXP/'mean-correction.npz')
    config=read(EXP/'sensor-config.json')
    with np.load(EXP/'charge-map-1.0.npz') as f:
        mapping={k:f[k].copy() for k in f.files}
    MODELS={'uniform':CausalRecovery(ANCHOR,config,EXP/'bin/age.so',mapping)}
    with np.load(EXP/'charge-map-spatial12.npz') as f:
        mapping={k:f[k].copy() for k in f.files}
    with np.load(EXP/'spatial-profile.npz') as f:
        profile=f['profile_12'].copy()
    MODELS['spatial12']=SpatialRecovery(ANCHOR,config,EXP/'bin/age.so',mapping,profile)
    NOISE=read(EXP/'noise-calibration.json')
    TRANSFORM=ANCHOR.transform(np.eye(ANCHOR.observed_samples)).T
    CASES={c['id']:c for c in read(EXP/'cases.json')}


def numerical():
    verify();initialize()
    rng=np.random.default_rng(731000)
    vector=rng.normal(size=ANCHOR.observed_samples)
    matrix_error=float(np.max(abs(TRANSFORM@vector-ANCHOR.transform(vector))))
    n=30; gradient=rng.normal(size=(2,n));delta=np.array([.02,-.03])
    residual=gradient.T@delta
    terms=np.stack([residual*.7,residual*.2,residual*.1,np.zeros(n)])[None]
    result,influence,jac=project(np.eye(n),DiagonalMetric(np.ones(n)),gradient,
                                 terms,np.array([1.,2.]),(.1,9.),1)
    span_error=float(np.max(abs(np.array(result['delta_log_energy'])-delta)))
    raw=rng.normal(size=n);orthogonal=raw-jac@(np.linalg.pinv(jac)@raw)
    orthogonal_error=float(np.max(abs(influence@orthogonal)))
    constant=np.ones((4,8,2,4,2))*np.arange(4)[None,None,None,:,None]
    samples,ci=paired_bootstrap(constant,seed=1,draws=100)
    bootstrap_error=float(np.max(abs(samples-constant[0,0])))
    check=dict(transform_matrix_error=matrix_error,synthetic_span_error=span_error,
               orthogonal_error=orthogonal_error,constant_bootstrap_error=bootstrap_error)
    assert max(check.values())<1e-9
    save(EXP/'numerical-checks.json',dict(passed=True,checks=check))
    print('Algebra and transform controls passed',flush=True)


def load_case(case):
    with np.load(PARENTART/'waveforms'/f"{case['id']:04d}.npz") as f:
        voltage=f['voltage_mV'].copy()
    with np.load(PARENTART/'diagnostics'/f"{case['id']:04d}.npz") as f:
        shared=f['event_signal_mV'].copy()
    with np.load(SHAREDART/'pair/pair'/f"{case['source_id']:04d}.npz") as f:
        isolated=f['on_isolated'].copy()
        assert np.array_equal(shared,f['on_shared'])
    assert voltage.shape==(8,1024)
    return voltage,shared,isolated


def projection_one(case):
    started=time.monotonic()
    voltage,shared,isolated=load_case(case)
    y=voltage.mean(0);times=case['times_ns'];energy=np.asarray(case['energy_dep_MeV'])
    weak=int(np.argmin(energy))
    common=common_inputs(MODELS['uniform'],y,times,NOISE)
    metrics=fixed_metrics(common['metric']);sig=signature(common['metric'])
    rows=[];arrays=dict(voltage_mean_mV=y,shared_mean_mV=shared.mean(0),isolated_mean_mV=isolated.mean(0))
    for index,(name,model) in enumerate(MODELS.items()):
        anchor=model.components(times,energy,interaction=False)
        mean,gradient=model.mean_and_gradient(times,np.log(energy),model.prepare(times),strength=1.)
        terms,closure=residual_terms(voltage,shared,isolated,anchor,mean,model.background,weak)
        previous=read(PARENTART/'pairs'/f"{288+(case['id']-16)*6+index:04d}.json")
        assert sig==previous['metric_signature']
        assert np.array_equal(common['starts'],np.asarray(previous['start_energy_MeV']))
        oldmean=model.components(times,previous['energy_MeV']).sum(0)
        wr=common['metric'].apply(ANCHOR.transform(oldmean-y+model.background))
        cost=float(wr@wr);cost_error=abs(cost-previous['cost'])/max(1.,abs(previous['cost']))
        assert cost_error<1e-7
        arrays[name+'_terms_mV']=terms
        arrays[name+'_mean_mV']=mean
        arrays[name+'_gradient_mV']=gradient
        for weighting,metric in metrics.items():
            result,influence,jac=project(TRANSFORM,metric,gradient,terms,energy,model.energy_bounds,weak)
            arrays[name+'_'+weighting+'_influence']=influence
            small_error=None
            if case['id'] in (16,17) and weighting=='full':
                step=np.array([.001,-.002])
                shifted=model.components(times,energy*np.exp(step)).sum(0)
                small_error=float(np.max(abs(influence@(shifted-mean)-step)))
                assert small_error<1e-4
            result.update(case=case['id'],condition=case['condition'],gap_ns=case['gap_ns'],
                block=case['block'],method=name,weight=weighting,weak_index=weak,
                true_energy_MeV=energy.tolist(),old_fit_energy_MeV=previous['energy_MeV'],
                old_fit_bias_percent=(100*(np.asarray(previous['energy_MeV'])/energy-1)).tolist(),
                old_fit_cost=previous['cost'],old_cost_relative_error=cost_error,
                metric_signature=sig,residual_closure_mV=closure,small_step_absolute_log_error=small_error)
            rows.append(result)
    differences={weight:symmetric_change(arrays['uniform_'+weight+'_influence'],
        arrays['spatial12_'+weight+'_influence'],arrays['uniform_terms_mV'].mean(0).sum(0),
        arrays['spatial12_terms_mV'].mean(0).sum(0)) for weight in WEIGHTS}
    (ART/'projections').mkdir(exist_ok=True)
    np.savez_compressed(ART/'projections'/f"{case['id']:04d}.npz",**arrays)
    save(ART/'projections'/f"{case['id']:04d}.json",dict(case=case['id'],rows=rows,
        symmetric_changes=differences,seconds=time.monotonic()-started))
    return dict(case=case['id'],projections=len(rows),seconds=time.monotonic()-started)


def pool_run(function,jobs,label):
    started=time.monotonic();rows=[]
    with ProcessPoolExecutor(max_workers=WORKERS,initializer=initialize) as pool:
        futures=[pool.submit(function,job) for job in jobs]
        for future in as_completed(futures):
            row=future.result();rows.append(row)
            with (REC/(label+'.jsonl')).open('a') as f:
                f.write(json.dumps(row,allow_nan=False)+'\n')
            print(f'{label}: {len(rows)}/{len(jobs)}',flush=True)
    save(REC/(label+'.json'),dict(rows=rows,seconds=time.monotonic()-started))
    return rows


def projections():
    verify()
    assert read(EXP/'numerical-checks.json')['passed']
    assert not (ART/'projections').exists()
    pool_run(projection_one,read(EXP/'cases.json'),'projections')
    paths=[p for p in (ART/'projections').glob('*') if p.is_file()]
    save(EXP/'projection-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})


def fit_one(job):
    task,batch=job
    output=ART/batch/f"{task['id']:04d}.json"
    if output.exists():
        raise FileExistsError(output)
    case=CASES[task['case']]
    with np.load(PARENTART/'waveforms'/f"{case['id']:04d}.npz") as f:
        y=f['voltage_mV'].mean(0)
    times=case['times_ns'];model=MODELS[task['method']]
    common=common_inputs(MODELS['uniform'],y,times,NOISE)
    metric_sig=signature(common['metric'])
    common['metric']=fixed_metrics(common['metric'])[task['weight']]
    if task['start']=='truth':
        common['starts']=np.asarray(case['energy_dep_MeV'],float)[None]
    if task['delete_term'] is not None:
        assert task['start']=='truth'
        with np.load(ART/'projections'/f"{case['id']:04d}.npz") as f:
            y=y-f[task['method']+'_terms_mV'][:,TERMS.index(task['delete_term'])].mean(0)
    result=fit_energy(model,y,times,common)
    energy=np.asarray(case['energy_dep_MeV']);weak=int(np.argmin(energy))
    error=100*(np.asarray(result['energy_MeV'])/energy-1)
    result.update(task,condition=case['condition'],gap_ns=case['gap_ns'],block=case['block'],
        batch=batch,true_energy_MeV=energy.tolist(),weak_index=weak,bias_percent=error.tolist(),
        weak_bias_percent=float(error[weak]),strict=bool(result['converged'] and np.max(abs(error))<=10),
        original_full_metric_signature=metric_sig,known_K_t=True,
        true_energy_in_fit=(task['start']=='truth'),ancestor_in_fit=(task['delete_term'] is not None))
    save(output,result)
    return dict(id=task['id'],seconds=result['seconds'],method=task['method'],
                starts=result['starts'],converged=result['converged'])


def pilot():
    verify()
    protocol=read(EXP/'protocol.json')
    tasks=[t for t in read(EXP/'fit-plan.json') if t['id'] in protocol['pilot_ids']]
    assert not (ART/'pilot').exists()
    rows=pool_run(fit_one,[(t,'pilot') for t in tasks],'pilot')
    costs={}
    for method in METHODS:
        costs[method]=max(r['seconds']/r['starts'] for r in rows if r['method']==method)
    all_tasks=read(EXP/'fit-plan.json')
    estimate=1.5*sum(costs[t['method']]*(3 if t['start']=='voltage' else 1) for t in all_tasks)/WORKERS
    gate=dict(projected_seconds=estimate,limit=protocol['pilot_projected_seconds_limit'],
              passed=bool(estimate<=protocol['pilot_projected_seconds_limit']),pilot_ids=protocol['pilot_ids'])
    save(EXP/'pilot-gate.json',gate)
    if not gate['passed']:
        raise RuntimeError('Fixed cost gate failed; no full run started')


def run(batch='fits'):
    verify()
    assert read(EXP/'pilot-gate.json')['passed']
    for path,digest in read(EXP/'projection-freeze.json').items():
        assert sha(ROOT/path)==digest,path
    tasks=read(EXP/'fit-plan.json')
    if batch=='replay':
        tasks=[t for t in tasks if t['id'] in read(EXP/'protocol.json')['replay_ids']]
    assert not (ART/batch).exists()
    pool_run(fit_one,[(t,batch) for t in tasks],batch)
    if batch=='fits':
        for taskid in read(EXP/'protocol.json')['pilot_ids']:
            a,b=[read(ART/p/f'{taskid:04d}.json') for p in ('pilot','fits')]
            for key in ('energy_MeV','cost','converged','boundary','start_energy_MeV'):
                assert a[key]==b[key],(taskid,key)
        print('All4pilot/formal fits exactly identical',flush=True)
    else:
        for task in tasks:
            a,b=[read(ART/p/f"{task['id']:04d}.json") for p in ('replay','fits')]
            for key in ('energy_MeV','cost','converged','boundary','start_energy_MeV'):
                assert a[key]==b[key],(task['id'],key)
        save(EXP/'replay-validation.json',dict(passed=True,ids=[t['id'] for t in tasks]))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=('prepare','numerical','project','pilot','run','replay','execute','verify'))
    stage=parser.parse_args().stage
    if stage=='execute':
        numerical();projections();pilot();run();run('replay');print('Diagnostic pipeline complete',flush=True)
    elif stage=='run':run()
    elif stage=='replay':run('replay')
    elif stage=='project':projections()
    elif stage=='verify':print(verify())
    else:globals()[stage]()


if __name__=='__main__':
    main()
