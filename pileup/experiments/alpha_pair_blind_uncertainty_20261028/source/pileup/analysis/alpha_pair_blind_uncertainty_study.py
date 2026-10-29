"""Frozen calibration and independent confirmation of blind joint diagnostics."""
from __future__ import annotations

import argparse
import ast
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import io
import json
import multiprocessing as mp
import os
from pathlib import Path
import subprocess
import sys
import time
import unittest

import numpy as np

from .alpha_pair_blind_uncertainty import (
    DEFAULTS, domain_for, exterior_search, waveform_objective, in_domain,
    empirical_level, diagnostic_label,
)
from .alpha_pair_joint_likelihood import ConservingResponse, JointNoise, waveform_candidates, select

PROJECT = Path(os.environ.get('GAGG_UNCERTAINTY_PROJECT', Path(__file__).resolve().parents[2])).resolve()
NAME = 'alpha_pair_blind_uncertainty_20261028'
EXPERIMENT = PROJECT/'pileup/experiments'/NAME
OUTPUT = PROJECT/'jupyter/artifacts'/NAME
RECORD = PROJECT/'jupyter/run_records'/NAME
NOTEBOOK = PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
PRIOR = PROJECT/'pileup/experiments/alpha_pair_selection_20261027'


def now(): return datetime.now(timezone.utc).isoformat()
def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def load(path): return json.loads(Path(path).read_text())


def normalize(value):
    if isinstance(value, dict): return {str(k): normalize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [normalize(v) for v in value]
    if isinstance(value, np.ndarray): return normalize(value.tolist())
    if isinstance(value, np.generic): return normalize(value.item())
    return value


def write_new(path, payload):
    """Explicit-offset shared-folder write with a complete readback."""
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        offset = 0
        while offset < len(payload):
            block = payload[offset:offset+256*1024]
            n = os.pwrite(fd, block, offset)
            if n <= 0 or n > len(block): raise OSError(f'Unexpected write length: {n}')
            offset += n
        os.ftruncate(fd, len(payload)); os.fsync(fd)
    finally: os.close(fd)
    if path.read_bytes() != payload: raise OSError(f'Readback differs: {path}')


def save(path, value):
    write_new(path, (json.dumps(normalize(value), ensure_ascii=False, indent=2, allow_nan=False)+'\n').encode())


def save_array(path, **arrays):
    stream = io.BytesIO(); np.savez_compressed(stream, **arrays); write_new(path, stream.getvalue())


def rows(path):
    return [json.loads(s) for s in Path(path).read_text().splitlines() if s.strip()]


def design(split):
    npair, ncontrol = (32, 16) if split == 'calibration' else (64, 32)
    groups = [('pair40', [18000., 800.], [200.3, 240.3], npair),
              ('pair200', [18000., 800.], [200.3, 400.3], npair), ('empty', [], [], ncontrol)]
    groups += [(f'single{q}', [float(q)], [200.3], ncontrol) for q in (800, 3000, 8000, 18800)]
    cases = []
    for band, q, t, n in groups:
        for _ in range(n):
            index = len(cases)
            cases.append(dict(index=index, band=band, kind='pair' if len(q)==2 else 'single' if q else 'empty',
                count=len(q), charge_pe=q, times_ns=t, gap_ns=t[1]-t[0] if len(t)==2 else 0.,
                seed_words=[202610280701, 1 if split=='calibration' else 2, index, 997]))
    return cases


def source_closure(names):
    """Copy only actual Python dependencies; no traversal of old run trees."""
    pending = list(names); found = set()
    while pending:
        name = pending.pop()
        if name in found: continue
        path = PROJECT/'pileup/analysis'/f'{name}.py'
        if not path.is_file(): raise FileNotFoundError(path)
        found.add(name)
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.ImportFrom) and node.level == 1 and node.module:
                dependency = node.module.split('.')[0]
                if dependency not in found: pending.append(dependency)
    return sorted(found)


def prepare():
    previous = load(PROJECT/'jupyter/run_records/research_next_steps_20261028/validation.json')
    assert previous['status']=='passed' and digest(NOTEBOOK)==previous['notebook_sha256']
    assert digest(PROJECT/'docs/literature-use.md')==previous['docs_sha256']
    for folder in (EXPERIMENT, OUTPUT, RECORD): folder.mkdir(parents=True, exist_ok=False)
    for path, name in ((NOTEBOOK, 'notebook.before.ipynb'), (PROJECT/'docs/literature-use.md', 'literature.before.md'),
                       (PROJECT/'README.md', 'README.before.md')):
        write_new(RECORD/name, path.read_bytes()); print('BACKUP', name, flush=True)
    source_names = source_closure(['alpha_pair_blind_uncertainty_study', 'test_alpha_pair_blind_uncertainty',
                                  'alpha_pair_model_sampler', 'test_alpha_pair_model_sampler', 'alpha_pileup_fit'])
    for rel in ('pileup/__init__.py', 'pileup/analysis/__init__.py'):
        p = PROJECT/rel; write_new(EXPERIMENT/'source'/rel, p.read_bytes() if p.exists() else b'')
    for name in source_names:
        write_new(EXPERIMENT/'source/pileup/analysis'/f'{name}.py', (PROJECT/'pileup/analysis'/f'{name}.py').read_bytes())
    inputs = {}
    for name in ('joint-response.npz', 'noise.npz', 'original-gate.json'):
        source = PRIOR/name; inputs[str(source.relative_to(PROJECT))] = digest(source)
        write_new(EXPERIMENT/name, source.read_bytes())
    protocol = dict(created=now(), question='After the frozen blind count decision, does an empirically calibrated joint parameter cost region contain truth, and do finite exterior searches reveal ambiguity?',
        interpretation='New project diagnostic; no theorem of nominal confidence, no global optimization or non-identifiability claim.',
        reference='Unchanged waveform_candidates K2 cost, not the global maximum-likelihood fit; select uses original native gate alpha0.025 unchanged.',
        region='D(y):ordered t in observed K1-trigger local domain enlarged to contain original K2,charges50..50000pe. R_c={theta in D:J(y;theta)-J(y;theta_hat_K2)<=c}. The baseline candidate stays unchanged.',
        calibration='144 new frames:32 each40/200ns pairs;16 each empty,singles800/3000/8000/18800pe. Only selected TRUE pairs calibrate c=max(0,ceil(.95*n)-th ordered joint score); at least20 selected pairs and finite c required. Out-of-domain truths count as +inf.',
        confirmation='Freeze c before generating288 separate frames:64 each40/200ns pairs,32 each of the five controls. Fixed stratified design, not a universal input mixture.',
        selection_scope='Calibrated mixture is the selected-pair mixture induced by this design and gate. Report per-gap conditional containment and all-pair selected-and-contained separately; no uniform nuisance or nominal95% guarantee.',
        quality='Search eight closed exterior boxes at or beyond either ordered time error10ns or relative charge error20%, denominator possible true Q. Use Qhat/1.2,Qhat/.8. Explicit in-region exterior witness means ambiguity. No witness plus search stability/convergence is only a diagnostic label, never a proof.',
        settings=DEFAULTS, quality_changes_predictions=False, empirical_quantile=.95, minimum_selected_calibration_pairs=20,
        model='Same frozen conditional correlated Gaussian mean/covariance; FIR first3 residuals zero; fixed alpha95/351ns,1024x2ns float64 voltage,ample dynamic range.',
        truth='Truth only in generation and scoring/calibration after inference SHA256 freezes. Workers reject truth/cases/ROOT opens; no real K,Q,t or labels enter inference.',
        pilot_indices=[0,1,2,32,33,34,64,80,128], example_indices=[0,1,64,65], workers_limit=8, blas_threads=1,
        stopping='Analytic tests and complete pilot cost/domain checks precede main calibration. If empirical c is unbounded/insufficient, do not launch full confirmation; report calibration limitation. Numerical failures remain in scoring denominators.',
        scope_limits='Purealpha model-only K<=2; no native photons,new ROOT,training,free tau,alpha/gamma or rate scan. Q equivalentpe not depositedenergy. No inference improvement claimed.',
        new_ROOTs=0, new_training=0, default_changed=False)
    save(EXPERIMENT/'protocol.json', protocol)
    for split in ('calibration', 'test'): save(EXPERIMENT/(split+'-design.truth.json'), design(split))
    wrapper = '''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
os.environ['GAGG_UNCERTAINTY_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_blind_uncertainty_study import main
    main()
'''
    write_new(EXPERIMENT/'run.py', wrapper.encode())
    frozen = {str(p.relative_to(EXPERIMENT)): digest(p) for p in (EXPERIMENT/'source').rglob('*.py')}
    for n in ('run.py', 'protocol.json', 'calibration-design.truth.json', 'test-design.truth.json',
              'joint-response.npz', 'noise.npz', 'original-gate.json'): frozen[n] = digest(EXPERIMENT/n)
    save(EXPERIMENT/'freeze.json', dict(created=now(), files=frozen))
    save(RECORD/'environment.json', dict(created=now(), python=sys.version, executable=sys.executable,
        dependencies={n:importlib.metadata.version(n) for n in ('numpy','scipy','matplotlib','pandas','nbformat')},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'], cwd=PROJECT, text=True).strip(),
        scoped_changes='Before changes:docs/literature-use.md +20 lines L68; no staged paths. New study sources are uncommitted and copied exactly. Execute-bit changes ignored.',
        status_scan='Broad shared-folder status/listing interrupted without edits; relevant files checked by scoped diff.',
        inputs=inputs, notebook_before_sha256=digest(RECORD/'notebook.before.ipynb'),
        literature_before_sha256=digest(RECORD/'literature.before.md'), source_names=source_names))
    print('UNCERTAINTY_PREPARED', len(frozen), flush=True)


def check_freeze():
    for name, h in load(EXPERIMENT/'freeze.json')['files'].items():
        assert digest(EXPERIMENT/name)==h, name


def models(): return ConservingResponse.load(EXPERIMENT/'joint-response.npz'), JointNoise.load(EXPERIMENT/'noise.npz')


def tests():
    check_freeze(); suite = unittest.TestSuite()
    for name in ('test_alpha_pair_blind_uncertainty', 'test_alpha_pair_model_sampler'):
        suite.addTests(unittest.defaultTestLoader.loadTestsFromName('pileup.analysis.'+name))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    assert result.wasSuccessful()
    save(RECORD/'tests.json', dict(passed=True, tests=result.testsRun, created=now()))


def generate(split):
    from .alpha_pair_model_sampler import GaussianFactor, lift_innovations
    check_freeze(); assert load(RECORD/'tests.json')['passed']
    if split=='test':
        validation = load(RECORD/'calibration-freeze.json'); assert validation['passed']
        assert digest(EXPERIMENT/'calibration.json')==validation['sha256']
    cases = load(EXPERIMENT/(split+'-design.truth.json'))
    folder = OUTPUT/split; folder.mkdir(exist_ok=False)
    response, noise = models(); waves=[]; q=np.zeros((len(cases),2)); t=q.copy(); maxerror=0.
    for c in cases:
        curves=response.components(c['times_ns'], c['charge_pe'])
        factor=GaussianFactor.from_noise(noise, curves, c['charge_pe'])
        z=factor.draw(np.random.default_rng(np.random.SeedSequence(c['seed_words'])))
        residual=lift_innovations(z, noise.coefficients)
        maxerror=max(maxerror, float(np.max(np.abs(noise.transform(residual)-z))))
        assert np.isclose(noise.score(residual,curves,c['charge_pe']), np.sum(factor.whiten(z)**2)+factor.logdet, rtol=1e-10, atol=1e-8)
        waves.append(curves.sum(0)+response.baseline+noise.background_mV+residual)
        q[c['index'],:c['count']]=c['charge_pe']; t[c['index'],:c['count']]=c['times_ns']
    y=np.asarray(waves, dtype=np.float64)
    assert y.shape==(len(cases),1024) and np.isfinite(y).all() and maxerror<1e-9
    assert len({row.tobytes() for row in y})==len(y)
    save_array(folder/'observations.npz', voltage_mV=y)
    save_array(folder/'labels.truth.npz', charge_pe=q, times_ns=t, count=np.array([c['count'] for c in cases]))
    save(folder/'labels.cases.json', cases)
    save(RECORD/(split+'-generation.json'), dict(passed=True,created=now(),frames=len(cases),max_fir_error=maxerror,
        files={p.name:digest(p) for p in folder.iterdir()},source_freeze_sha256=digest(EXPERIMENT/'freeze.json'),
        calibration_sha256=digest(EXPERIMENT/'calibration.json') if split=='test' else None))
    print('GENERATED', split, len(y), flush=True)


def block_labels(event, args):
    if event=='open':
        path=str(args[0])
        if '.truth.' in path or path.endswith(('.cases.json','.root')) or '/truth/' in path:
            raise PermissionError('Inference attempted label access: '+path)


def initialize(split):
    global RESPONSE, NOISE, OBS, GATE, SETTINGS
    RESPONSE, NOISE=models(); sys.addaudithook(block_labels)
    try:
        with (EXPERIMENT/'calibration-design.truth.json').open():
            raise AssertionError('Label-access audit failed')
    except PermissionError:
        pass
    GATE=load(EXPERIMENT/'original-gate.json'); SETTINGS=load(EXPERIMENT/'protocol.json')['settings']
    with np.load(OUTPUT/split/'observations.npz') as f: OBS=f['voltage_mV']


def infer_one(index):
    start=time.monotonic(); y=OBS[index]
    candidates=waveform_candidates(RESPONSE,NOISE,y); selected=select(candidates,GATE)
    objective=waveform_objective(RESPONSE,NOISE,y)
    error=max(abs(objective(c['times_ns'],c['charge_pe'])-c['cost'])/max(1.,abs(c['cost'])) for c in candidates)
    assert error<1e-9
    domain=domain_for(RESPONSE,candidates[1],candidates[2]) if len(candidates)==3 else None
    exterior=None
    if len(selected['times_ns'])==2:
        exterior=exterior_search(objective,candidates[1],selected,domain,SETTINGS)
        assert exterior['cost_replay_error']<1e-9
        for box in exterior['verified']:
            if box['best'] is not None:
                f=box['best']; assert np.isfinite(f['cost'])
                assert domain['time'][0]-1e-8<=f['times_ns'][0]<=f['times_ns'][1]+1e-8<=domain['time'][1]+1e-8
                assert all(domain['charge'][0]-1e-7<=q<=domain['charge'][1]+1e-7 for q in f['charge_pe'])
    return dict(index=index,candidates=candidates,selected=selected,domain=domain,exterior=exterior,
                seconds=time.monotonic()-start,cost_replay_error=error,truth_inputs=[])


def guarded_infer_one(index):
    try:
        return infer_one(index)
    except (ArithmeticError, ValueError, np.linalg.LinAlgError) as exc:
        # Explicit numerical failures stay in the frame denominator. Security,
        # I/O, assertion and programming errors stop the run for investigation.
        return dict(index=index,candidates=[],selected=dict(times_ns=[],charge_pe=[],converged=False),
            domain=None,exterior=None,seconds=0.,cost_replay_error=0.,truth_inputs=[],
            numerical_failure=type(exc).__name__+': '+str(exc))


def infer(split, workers, pilot=False):
    check_freeze(); assert 1<=workers<=8
    check=load(RECORD/(split+'-generation.json')); assert check['passed']
    assert digest(OUTPUT/split/'observations.npz')==check['files']['observations.npz']
    if split=='test': assert digest(EXPERIMENT/'calibration.json')==check['calibration_sha256']
    prior=[]; tasks=list(range(check['frames']))
    if pilot:
        assert split=='calibration'; tasks=load(EXPERIMENT/'protocol.json')['pilot_indices']
    elif split=='calibration':
        v=load(RECORD/'pilot-validation.json'); assert v['passed']
        assert digest(OUTPUT/split/'pilot.jsonl')==v['sha256']
        prior=rows(OUTPUT/split/'pilot.jsonl'); done={r['index'] for r in prior}; tasks=[i for i in tasks if i not in done]
    target=OUTPUT/split/('pilot.jsonl' if pilot else 'predictions.jsonl')
    start=time.monotonic(); collected=list(prior)
    with target.open('x') as stream, ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize,initargs=(split,)) as pool:
        for row in prior: stream.write(json.dumps(normalize(row),allow_nan=False)+'\n')
        futures={pool.submit(guarded_infer_one,i):i for i in tasks}
        for future in as_completed(futures):
            row=future.result(); collected.append(row)
            stream.write(json.dumps(normalize(row),allow_nan=False)+'\n'); stream.flush()
            print('UNCERTAINTY',split,'pilot' if pilot else 'main',len(collected),len(tasks)+len(prior),
                  'index',row['index'],'K',len(row['selected']['times_ns']),'s',round(row['seconds'],2),flush=True)
    assert len(collected)==len(tasks)+len(prior) and len({r['index'] for r in collected})==len(collected)
    save(RECORD/('pilot-validation.json' if pilot else split+'-inference.json'),dict(passed=True,created=now(),
        frames=len(collected),selected_pairs=sum(len(r['selected']['times_ns'])==2 for r in collected),
        stable_exterior=sum(bool(r['exterior'] and r['exterior']['stable']) for r in collected),
        seconds=time.monotonic()-start,workers=workers,inherited_pilot=len(prior),sha256=digest(target),
        max_cost_replay_error=max(r['cost_replay_error'] for r in collected),truth_inputs=[],labels_blocked=True))
    print('INFERENCE_COMPLETE',split,len(collected),flush=True)


def score_rows(split):
    """Only called after primary predictions have frozen; no feedback to fits."""
    from .alpha_pileup_fit import score
    check=load(RECORD/(split+'-inference.json')); assert check['passed']
    assert digest(OUTPUT/split/'predictions.jsonl')==check['sha256']
    predictions={r['index']:r for r in rows(OUTPUT/split/'predictions.jsonl')}
    cases=load(OUTPUT/split/'labels.cases.json')
    with np.load(OUTPUT/split/'observations.npz') as f: waves=f['voltage_mV']
    with np.load(OUTPUT/split/'labels.truth.npz') as f: q=f['charge_pe']; t=f['times_ns']
    response,noise=models(); scored=[]
    for c in cases:
        i=c['index']; r=predictions[i]; n=c['count']; fitted=r['selected']
        truth=dict(first_active_ns=t[i,:n],emission_ns=t[i,:n],charge_pe=q[i,:n])
        evaluation=score(fitted,truth,[0.,2048.],tolerance_ns=10.)
        strict=bool(evaluation['recovered_count']==n and evaluation['unmatched_predictions']==0)
        joint=None; admissible=None
        if n==2:
            admissible=bool(r['domain'] is not None and in_domain(t[i],q[i],r['domain']))
            if admissible:
                objective=waveform_objective(response,noise,waves[i])
                joint=float(objective(t[i],q[i])-r['candidates'][2]['cost'])
        scored.append(dict(index=i,band=c['band'],kind=c['kind'],true_count=n,
            predicted_count=len(fitted['times_ns']),strict_success=strict,evaluation=normalize(evaluation),
            truth_in_domain=admissible,joint_delta=joint,baseline_converged=fitted.get('converged',True),
            numerical_failure=r.get('numerical_failure')))
    return scored, predictions


def calibrate():
    check_freeze(); assert not (OUTPUT/'test').exists()
    scored, predictions=score_rows('calibration')
    subset=[r for r in scored if r['true_count']==2 and r['predicted_count']==2]
    values=[r['joint_delta'] if r['truth_in_domain'] else float('inf') for r in subset]
    calibration=empirical_level(values)
    calibration.update(created=now(),basis='Selected true calibration pairs only; fixed stratified mixture',
        selected_per_gap={band:sum(r['band']==band for r in subset) for band in ('pair40','pair200')},
        out_of_domain=sum(not r['truth_in_domain'] for r in subset),test_generated=False,
        prediction_sha256=digest(OUTPUT/'calibration/predictions.jsonl'),
        truth_used_for_calibration_after_predictions=True,quality_rule='Fixed diagnostic_label; no threshold optimization on event accuracy')
    for r in scored: r['quality_label']=diagnostic_label(predictions[r['index']],calibration['level'])
    save(OUTPUT/'calibration/scores.truth.json',scored)
    save(EXPERIMENT/'calibration.json',calibration)
    save(RECORD/'calibration-freeze.json',dict(passed=calibration['level'] is not None,created=now(),
        sha256=digest(EXPERIMENT/'calibration.json'),calibration_scores_sha256=digest(OUTPUT/'calibration/scores.truth.json'),
        no_confirmation_generated=True,prediction_sha256=calibration['prediction_sha256']))
    print('CALIBRATION_FROZEN',json.dumps(calibration,ensure_ascii=False),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['prepare','tests','generate-calibration','pilot','infer-calibration','calibrate','generate-test','infer-test'])
    parser.add_argument('--workers',type=int,default=2); args=parser.parse_args()
    if args.command=='prepare': prepare()
    elif args.command=='tests': tests()
    elif args.command.startswith('generate-'): generate(args.command.split('-',1)[1])
    elif args.command=='pilot': infer('calibration',args.workers,pilot=True)
    elif args.command=='calibrate': calibrate()
    else: infer(args.command.split('-',1)[1],args.workers)


if __name__=='__main__': main()
