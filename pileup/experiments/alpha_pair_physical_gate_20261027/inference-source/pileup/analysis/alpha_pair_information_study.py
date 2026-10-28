"""Frozen known-K information analysis using exposed alpha-pair parameters."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

from .alpha_pair_information import converged_information

PROJECT = Path(os.environ.get('GAGG_INFORMATION_PROJECT', Path(__file__).resolve().parents[2])).resolve()
NAME = 'alpha_pair_information_20261027'
EXPERIMENT = PROJECT / 'pileup/experiments' / NAME
OUTPUT = PROJECT / 'jupyter/artifacts' / NAME
RECORD = PROJECT / 'jupyter/run_records' / NAME
PREVIOUS = PROJECT / 'pileup/experiments/alpha_pair_model_closure_20261026'
INPUT = PROJECT / 'jupyter/artifacts/alpha_pair_model_closure_20261026'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def serializable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def save(path, value):
    Path(path).write_text(json.dumps(value, default=serializable, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def rows(path):
    return [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]


def models():
    from .alpha_pair_joint_likelihood import ConservingResponse, JointNoise
    return (ConservingResponse.load(EXPERIMENT/'joint-response.npz'),
            JointNoise.load(EXPERIMENT/'joint_conserved_nll-noise.npz'))


def check_freeze():
    for name, expected in load(EXPERIMENT/'freeze.json')['files'].items():
        assert digest(EXPERIMENT/name) == expected, name
    for name, expected in load(RECORD/'input-provenance.json')['files'].items():
        assert digest(PROJECT/name) == expected, name


def prepare():
    assert not any(p.exists() for p in (EXPERIMENT, OUTPUT, RECORD))
    for p in (EXPERIMENT, OUTPUT, RECORD):
        p.mkdir(parents=True)
    previous_validation = load(PROJECT/'jupyter/run_records/alpha_pair_model_closure_20261026/final-validation.json')
    assert previous_validation['passed']
    notebook = PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    assert digest(notebook) == previous_validation['notebook_sha256']
    shutil.copyfile(notebook, RECORD/'notebook.before.ipynb')
    for name in ('joint-response.npz', 'joint_conserved_nll-noise.npz', 'settings.json'):
        shutil.copyfile(PREVIOUS/name, EXPERIMENT/name)
    cases = load(INPUT/'native.cases.json')
    with np.load(INPUT/'native-labels.truth.npz') as f:
        t, q = f['times_ns'], f['charge_pe']
    tasks = []
    for i, case in enumerate(cases):
        if case['kind'] == 'pair':
            assert np.all(q[i] >= 500.)
            tasks.append(dict(task_id=len(tasks), role='native_parameters', native_index=i,
                condition=case['condition'], gap_ns=case['gap_ns'], block=case['block'],
                replicate=case['replicate'], times_ns=t[i].tolist(), charge_pe=q[i].tolist()))
    assert len(tasks) == 360
    for condition in ('high_weak', 'weak_high'):
        for weak_charge in (500., 700., 1500., 3000.):
            for gap in (2., 4., 8., 12., 20., 40., 80., 120., 200., 400.):
                charges = [18000., weak_charge] if condition == 'high_weak' else [weak_charge, 18000.]
                tasks.append(dict(task_id=len(tasks), role='nominal_grid', native_index=None,
                    condition=condition, gap_ns=gap, times_ns=[160.37, 160.37+gap], charge_pe=charges,
                    weak_charge_pe=weak_charge, strong_charge_pe=18000.))
    for phase in (.13, .37, .73, 1.13):
        tasks.append(dict(task_id=len(tasks), role='phase_check', native_index=None,
            condition='high_weak', gap_ns=40., times_ns=[160.+phase, 200.+phase],
            charge_pe=[18000., 700.], phase_ns=phase))
    assert len(tasks) == 444
    save(EXPERIMENT/'tasks.truth.json', tasks)
    protocol = dict(created=datetime.now(timezone.utc).isoformat(),
        question='Separate selecting two outputs, locating both onsets, and estimating both charges; quantify local nuisance coupling',
        scope='Pure alpha;95/351ns fixed;2ns x1024 analog voltage; adequate dynamic range; existing blind K<=2 results are unchanged',
        data='760 exposed native and1520 matched evaluation frames from model closure;360 distinct numerical pair parameter points; no new waveform fits or ROOTs',
        truth='All information calculations condition on true Q/t and K=2. Only evaluation and diagnosis, never a blind input or operational interval',
        coordinates=['t1/(10ns)', 't2/(10ns)', 'ln(Q1/pe)', 'ln(Q2/pe)'],
        information='Fij=mu_i.T inv(C) mu_j +0.5 tr(inv(C) C_i inv(C) C_j), for the1021-dimensional FIR marginal',
        reference='Inverse expected Fisher information is a local regular/unbiased bound; finite, biased, selected, nonlinear and misspecified estimators need not attain it',
        conditioning=['All four parameters unknown', 'Both true times fixed: invert charge submatrix', 'Strong Q/t fixed: invert weak submatrix; strong stochastic noise remains'],
        derivatives=dict(initial_step=1e-3, primary_step=5e-4, adaptive_steps=[2.5e-4,1.25e-4],
            maximum_relative_matrix_change=.02, maximum_sigma_fractional_change=.02,
            scope='Step in scaled time and logQ; retain unresolved points and piecewise-template limitations'),
        inverse='No pseudoinverse/ridge; record nonpositive information or condition number>1e12 as unresolved',
        tasks=444, native_points=360, nominal_points=80, phase_points=4,
        nominal_grid='Qstrong18000pe;Qweak500/700/1500/3000pe;gaps2/4/8/12/20/40/80/120/200/400ns;both orders;first onset160.37ns',
        nominal_warning='Conditional model design coordinates, not new avalanche realizations or measured performance; single sampling phase',
        pilot_task_ids=[0,160,280],
        empirical='Reuse native240 and Gaussian480 truth-assisted pair fits; normalize known-K and true-start weak Q/t errors and fixed-time Q errors by local information',
        thresholds='Both trueQ>=500pe,eachQ within20%,eacht within min(10ns,gap/4),no extra output; no gate retuning',
        reference_regions='Chi-square2 95% radius and standard-normal95% interval are local Gaussian references, not calibrated confidence intervals',
        figures=['0004-77-1','0004-78-1','0004-79-1','0004-80-1','0004-80-2'],
        decision='Diagnostic only; no automatic replacement or new independent production',
        literature='Existing Liddick2012 notes: detection vs estimation; Fixsen2004 notes: conditional means/covariances. Fisher derivation is project mathematics, not a reported GAGG result.')
    save(EXPERIMENT/'protocol.json', protocol)
    inputs = ['native.cases.json','native-labels.truth.npz','labels.cases.json','labels.truth.npz',
              'scores.jsonl','diagnostic-summary.json','summary.json','report.json']
    provenance = {str((INPUT/name).relative_to(PROJECT)):digest(INPUT/name) for name in inputs}
    provenance.update({str((PREVIOUS/name).relative_to(PROJECT)):digest(PREVIOUS/name)
                       for name in ('freeze.json','report-freeze.json','joint-response.npz','joint_conserved_nll-noise.npz')})
    save(RECORD/'input-provenance.json', dict(files=provenance, inherited_validation=previous_validation,
        notebook_before_sha256=digest(notebook), native_Q='Realized avalanche-equivalent pe',
        matched_Q='Conditioning parameter equal to native numerical Q, not a new avalanche charge'))
    save(RECORD/'environment.json', dict(python=sys.version, executable=sys.executable,
        dependencies={name:importlib.metadata.version(name) for name in ('numpy','scipy','matplotlib','nbformat','pandas')},
        cpu_affinity=sorted(os.sched_getaffinity(0)), workers=4,
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip()))
    (RECORD/'git-status.before.txt').write_text(subprocess.check_output(
        ['git','--no-optional-locks','-c','core.fileMode=false','status','--short','--','pileup/analysis','README.md','docs/literature-use.md'],
        cwd=PROJECT,text=True))
    # Capture project modules actually needed to load and evaluate this model.
    response, noise = models()
    from .alpha_pair_information import model_distribution, coordinates
    model_distribution(response, noise, coordinates(t[0], q[0]))
    loaded = set()
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename:
            path = Path(filename).resolve()
            if path.is_relative_to(PROJECT) and path.suffix == '.py' and ('analysis' in path.relative_to(PROJECT).parts):
                loaded.add(path)
    loaded.add(PROJECT/'pileup/analysis/alpha_pair_information_study.py')
    loaded.add(PROJECT/'pileup/analysis/test_alpha_pair_information.py')
    for path in sorted(loaded):
        destination = EXPERIMENT/'inference-source'/path.relative_to(PROJECT)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, destination)
    for package in ('pileup', 'pileup/analysis', 'analysis'):
        destination = EXPERIMENT/'inference-source'/package/'__init__.py'
        destination.parent.mkdir(parents=True, exist_ok=True)
        if not destination.exists():
            destination.write_text('')
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
os.environ['GAGG_INFORMATION_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_information_study import main
    main()
''')
    files = {str(p.relative_to(EXPERIMENT)):digest(p) for p in (EXPERIMENT/'inference-source').rglob('*.py')}
    for name in ('joint-response.npz','joint_conserved_nll-noise.npz','settings.json','protocol.json','tasks.truth.json','run.py'):
        files[name] = digest(EXPERIMENT/name)
    save(EXPERIMENT/'freeze.json', dict(created=protocol['created'], files=files,
        source_capture='Project modules loaded by model evaluation, plus study and tests; no new native physics run',
        git_head=load(RECORD/'environment.json')['git_head']))
    print('INFORMATION_FROZEN',len(files),'files',len(tasks),'tasks',flush=True)


def initialize():
    global RESPONSE, NOISE
    RESPONSE, NOISE = models()


def calculate(task):
    started = time.monotonic()
    try:
        value = converged_information(RESPONSE, NOISE, task['times_ns'], task['charge_pe'])
        return dict(task, information=value, status='ok', seconds=time.monotonic()-started)
    except np.linalg.LinAlgError as error:
        return dict(task, information=None, status='unresolved', error=str(error), seconds=time.monotonic()-started)


def run(workers=4, pilot=False):
    assert 1 <= workers <= 8
    check_freeze()
    if not pilot:
        assert load(RECORD/'pilot-validation.json')['passed']
    tasks = load(EXPERIMENT/'tasks.truth.json')
    ids = load(EXPERIMENT/'protocol.json')['pilot_task_ids'] if pilot else list(range(len(tasks)))
    destination = OUTPUT/'information.truth.jsonl'
    existing = rows(destination) if destination.exists() else []
    done = {r['task_id'] for r in existing}
    assert len(done) == len(existing)
    pending = [tasks[i] for i in ids if i not in done]
    started = time.monotonic()
    with destination.open('a') as stream:
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn'), initializer=initialize) as pool:
            futures = [pool.submit(calculate, task) for task in pending]
            for j, future in enumerate(as_completed(futures), 1):
                stream.write(json.dumps(future.result(), default=serializable, allow_nan=False)+'\n')
                stream.flush()
                if j % 20 == 0:
                    print('INFORMATION',j,'/',len(pending),flush=True)
        os.fsync(stream.fileno())
    save(RECORD/('pilot-run.json' if pilot else 'information-run.json'), dict(passed=True,
        requested=len(ids), new_tasks=len(pending), workers=workers, seconds=time.monotonic()-started,
        output_sha256=digest(destination), truth_used=True, new_blind_fits=0, new_ROOTs=0))
    print('INFORMATION_COMPLETE',len(ids),flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare','run'))
    parser.add_argument('--pilot', action='store_true')
    parser.add_argument('--workers',type=int,default=4)
    args = parser.parse_args()
    if args.command == 'prepare':
        prepare()
    else:
        run(args.workers,args.pilot)


if __name__ == '__main__':
    main()
