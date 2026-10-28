"""Separate event selection from charge inference with frozen Gaussian controls.

The physical fitter is unchanged. Fresh single-pulse calibration precedes a
separate confirmation sample; old interval curves are diagnostic data only.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import importlib.metadata
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import unittest

import numpy as np
from scipy.stats import binomtest

from .alpha_pair_information_study import load, save, digest, rows, serializable
from .alpha_pair_stability_study import copy_chunks
from .alpha_pair_joint_likelihood import ConservingResponse, JointNoise, select, evidence, calibrate_gate

PROJECT = Path(os.environ.get('GAGG_SELECTION_PROJECT', Path(__file__).resolve().parents[2])).resolve()
NAME = 'alpha_pair_selection_20261027'
EXPERIMENT = PROJECT / 'pileup/experiments' / NAME
OUTPUT = PROJECT / 'jupyter/artifacts' / NAME
RECORD = PROJECT / 'jupyter/run_records' / NAME
OLD_NAME = 'alpha_pair_contour_coverage_20261027'
OLD_EXP = PROJECT / 'pileup/experiments' / OLD_NAME
OLD_OUT = PROJECT / 'jupyter/artifacts' / OLD_NAME
OLD_REC = PROJECT / 'jupyter/run_records' / OLD_NAME
NOTEBOOK = PROJECT / 'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
ALPHAS = (.025, .05, .10)
METHODS = tuple(f'{origin}_{int(1000*a):03d}' for origin in ('original', 'matched') for a in ALPHAS)


def now():
    return datetime.now(timezone.utc).isoformat()


def fraction(k, n):
    if not n:
        return dict(numerator=int(k), denominator=0, fraction=None, binomial95=None)
    ci = binomtest(int(k), int(n)).proportion_ci(.95, method='exact')
    return dict(numerator=int(k), denominator=int(n), fraction=k/n, binomial95=[ci.low, ci.high])


def inclusion_partition(included, selected):
    included, selected = np.asarray(included, bool), np.asarray(selected, bool)
    if included.shape != selected.shape or included.ndim != 1:
        raise ValueError('Aligned one-dimensional frame masks required')
    n = len(included)
    return dict(all=fraction(included.sum(), n), selected=fraction((included & selected).sum(), selected.sum()),
                unselected=fraction((included & ~selected).sum(), (~selected).sum()),
                selected_and_included=fraction((included & selected).sum(), n),
                selected_fraction=fraction(selected.sum(), n))


def design(split):
    # These points are fixed before either new split is generated.
    groups = [(str(q), [float(q)], [200.3]) for q in (800, 3000, 8000, 18800)]
    if split == 'test':
        groups = [('empty', [], [])] + groups + [(f'pair{g}', [18000., 800.], [200.3, 200.3+g]) for g in (40, 200)]
    cases = []
    for band, q, t in groups:
        for _ in range(64):
            i = len(cases)
            cases.append(dict(index=i, split=split, kind='pair' if len(q)==2 else 'single' if q else 'empty',
                band=band, condition='high_weak' if len(q)==2 else band, gap_ns=t[1]-t[0] if len(t)==2 else 0.,
                charge_pe=q, times_ns=t, count=len(q), seed_words=[202610272001, 1 if split=='calibration' else 2, i, 995]))
    return cases


def prepare():
    previous = load(OLD_REC/'final-validation.json')
    assert previous['passed'] and digest(NOTEBOOK) == previous['notebook_sha256']
    for name, h in previous['source_sha256'].items():
        assert digest(PROJECT/name) == h, name
    for name, h in previous['documentation_sha256'].items():
        assert digest(PROJECT/name) == h, name
    for folder in (EXPERIMENT, OUTPUT, RECORD):
        folder.mkdir(parents=True, exist_ok=False)
    for source, name in ((NOTEBOOK, 'notebook.before.ipynb'), (PROJECT/'README.md', 'README.before.md'),
                         (PROJECT/'docs/literature-use.md', 'literature.before.md')):
        copy_chunks(source, RECORD/name)
    for source, name in (('joint-response.npz','joint-response.npz'), ('joint_conserved_nll-noise.npz','noise.npz'),
                         ('coverage-baseline-gate.json','original-gate.json')):
        copy_chunks(OLD_EXP/source, EXPERIMENT/name)
    shutil.copytree(OLD_EXP/'inference-source', EXPERIMENT/'source')
    for name in ('alpha_pair_selection_study.py', 'test_alpha_pair_selection.py', 'test_alpha_pair_model_sampler.py', 'alpha_pair_risk_report.py'):
        copy_chunks(PROJECT/'pileup/analysis'/name, EXPERIMENT/'source/pileup/analysis'/name)
    protocol = dict(created=now(), question='Does event-count selection bias recovered weak charge and interval inclusion; can matched-null calibration improve detection without unacceptable extras?',
        scope='Pure alpha,fixed95/351ns,full1024x2ns float64 voltage,ample dynamic range; unchanged waveform_candidates and select; local K<=2.',
        design='256 calibration singles:64 each800/3000/8000/18800pe.448 independent confirmation:64 each empty,the four singles,pair40ns,pair200ns. Pairs18000/800pe,t1=200.3ns.',
        distribution='Frozen response and candidate-dependent JointNoise Gaussian in FIR space,conditional first-three residual samples zero; no native photons or new ROOT.',
        calibration='Same original observed fitted-Q bins1500/5000/12500pe; only calibration single predictions used. Matched calibration freezes before confirmation generation.',
        methods=list(METHODS), alphas=list(ALPHAS), primary_comparison=['original_025','matched_025'],
        comparison='Other four thresholds are prespecified tradeoffs, not selected using test truth. No default update. Rank alpha is a decision parameter, not a charge confidence level or a guaranteed FPR.',
        old_interval_data='All32 previous pairs and their original interval sets are reused only for exploratory selection decomposition; no fresh interval coverage measurement.',
        denominators='Report P(C),P(S),P(C|S),P(C and S),P(C|not S) separately; undefined conditional rates stay null, not zero.',
        recovery='Both Q errors<=20%,times<=min(10ns,gap/4),no unmatched predictions. Also ordinary10ns matched Q/time errors,missed,extra,split,merge.',
        truth='Simulation/scoring only; primary workers block truth,cases,ROOT. No true count/times/charges in fits or selection. All candidates fitted before labels are read.',
        limitations='Fixed parameter points; empirical fitted-Q bin calibration is not uniform nuisance-parameter coverage; finite optimizer; no energy calibration,free tau,alpha-gamma,K>=3 or rate scan.',
        pilot_indices=[0,64,128,192], example_test_indices=[320,321,384,385], workers_limit=8, blas_threads=1,
        new_ROOTs=0, new_training=0, fitter_changed=False, new_interval_profiles=0)
    save(EXPERIMENT/'protocol.json', protocol)
    for split in ('calibration','test'):
        save(EXPERIMENT/(split+'-design.truth.json'), design(split))
    wrapper = '''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for key in ('SELECTION','CONTOUR','STABILITY','EVIDENCE','GATE','RISK','INFORMATION'):os.environ['GAGG_'+key+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_selection_study import main
    main()
'''
    (EXPERIMENT/'run.py').write_text(wrapper)
    files = {str(p.relative_to(EXPERIMENT)):digest(p) for p in (EXPERIMENT/'source').rglob('*.py')}
    for n in ('run.py','protocol.json','calibration-design.truth.json','test-design.truth.json','joint-response.npz','noise.npz','original-gate.json'):
        files[n] = digest(EXPERIMENT/n)
    save(EXPERIMENT/'freeze.json', dict(created=now(), files=files))
    old = {}
    for name in ('coverage/predictions.jsonl','coverage/observations.npz','coverage/labels.truth.npz','coverage/labels.cases.json',
                 'coverage-intervals.json','coverage-frames.json','report.json'):
        old[name] = digest(OLD_OUT/name)
    save(RECORD/'input-provenance.json', dict(old_results=old, previous_validation=previous,
        previous_validation_sha256=digest(OLD_REC/'final-validation.json'), previous_freeze_sha256=digest(OLD_EXP/'coverage-freeze.json')))
    status=subprocess.check_output(['git','-c','core.fileMode=false','status','--short','--untracked-files=no','--','README.md','docs/literature-use.md','pileup/analysis'],cwd=PROJECT,text=True,timeout=60)
    save(RECORD/'environment.json', dict(created=now(), python=sys.version, executable=sys.executable,
        dependencies={n:importlib.metadata.version(n) for n in ('numpy','scipy','matplotlib','pandas','nbformat')},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),scoped_git_status=status,
        uncommitted='Existing changes retained; actual source snapshot,protocol and old documentation backups saved. No commit or push.'))
    print('SELECTION_PREPARED',len(files),flush=True)


def check_freeze():
    for name,h in load(EXPERIMENT/'freeze.json')['files'].items():
        assert digest(EXPERIMENT/name)==h,name


def models():
    return ConservingResponse.load(EXPERIMENT/'joint-response.npz'), JointNoise.load(EXPERIMENT/'noise.npz')


def tests():
    check_freeze()
    suite=unittest.TestSuite()
    for n in ('test_alpha_pair_selection','test_alpha_pair_model_sampler'):
        suite.addTests(unittest.defaultTestLoader.loadTestsFromName('pileup.analysis.'+n))
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    assert result.wasSuccessful()
    save(RECORD/'tests.json',dict(passed=True,tests=result.testsRun,created=now()))


def generate(split):
    from .alpha_pair_model_sampler import GaussianFactor, lift_innovations
    check_freeze();assert load(RECORD/'tests.json')['passed']
    if split=='test':
        gatecheck=load(RECORD/'gate-validation.json');assert gatecheck['passed']
        assert digest(EXPERIMENT/'gates.json')==gatecheck['sha256']
    cases=load(EXPERIMENT/(split+'-design.truth.json'))
    folder=OUTPUT/split;folder.mkdir(exist_ok=False)
    response,noise=models();observations=[];q=np.zeros((len(cases),2));t=q.copy();errors=[]
    for c in cases:
        curves=response.components(c['times_ns'],c['charge_pe']);factor=GaussianFactor.from_noise(noise,curves,c['charge_pe'])
        z=factor.draw(np.random.default_rng(np.random.SeedSequence(c['seed_words'])))
        residual=lift_innovations(z,noise.coefficients)
        y=curves.sum(0)+response.baseline+noise.background_mV+residual
        error=float(np.max(np.abs(noise.transform(residual)-z)));assert error<1e-9;errors.append(error)
        cost=noise.score(residual,curves,c['charge_pe']);expected=float(np.sum(factor.whiten(z)**2)+factor.logdet)
        assert np.isclose(cost,expected,rtol=1e-10,atol=1e-8)
        observations.append(y);q[c['index'],:c['count']]=c['charge_pe'];t[c['index'],:c['count']]=c['times_ns']
    y=np.asarray(observations,dtype=np.float64)
    assert y.shape==(len(cases),1024) and np.isfinite(y).all()
    assert len({v.tobytes() for v in y})==len(y)
    np.savez_compressed(folder/'observations.npz',voltage_mV=y)
    np.savez_compressed(folder/'labels.truth.npz',charge_pe=q,times_ns=t,count=np.array([c['count'] for c in cases]))
    save(folder/'labels.cases.json',cases)
    save(RECORD/(split+'-generation.json'),dict(passed=True,created=now(),frames=len(cases),max_fir_error=max(errors),
        files={p.name:digest(p) for p in folder.iterdir()},freeze_sha256=digest(EXPERIMENT/'freeze.json'),
        gates_sha256=digest(EXPERIMENT/'gates.json') if split=='test' else None))
    print('GENERATED',split,len(y),flush=True)


def initialize(split):
    global RESPONSE,NOISE,OBS
    from .alpha_pair_charge_allocation import forbid_labels
    RESPONSE,NOISE=models()
    sys.addaudithook(forbid_labels)
    with np.load(OUTPUT/split/'observations.npz') as f:OBS=f['voltage_mV']


def infer_one(index):
    from .alpha_pair_joint_likelihood import waveform_candidates
    began=time.monotonic();y=OBS[index]
    fits=waveform_candidates(RESPONSE,NOISE,y)
    error=0.
    for fit in fits:
        c=RESPONSE.components(fit['times_ns'],fit['charge_pe'])
        cost=NOISE.score(y-RESPONSE.baseline-NOISE.background_mV-c.sum(0),c,fit['charge_pe'])
        error=max(error,abs(cost-fit['cost'])/max(1.,abs(cost)))
    assert error<1e-9
    return dict(index=index,candidates=fits,seconds=time.monotonic()-began,cost_replay_error=error,truth_inputs=[])


def infer(split,workers,pilot=False):
    check_freeze();assert 1<=workers<=8
    check=load(RECORD/(split+'-generation.json'));assert check['passed']
    for name,h in check['files'].items():assert digest(OUTPUT/split/name)==h
    if split=='test':assert digest(EXPERIMENT/'gates.json')==check['gates_sha256']
    prior=[];tasks=list(range(check['frames']))
    if pilot:
        assert split=='calibration';tasks=load(EXPERIMENT/'protocol.json')['pilot_indices']
    elif split=='calibration':
        validation=load(RECORD/'pilot-validation.json');assert validation['passed']
        assert digest(OUTPUT/split/'pilot.jsonl')==validation['sha256']
        prior=rows(OUTPUT/split/'pilot.jsonl');done={r['index'] for r in prior};tasks=[i for i in tasks if i not in done]
    target=OUTPUT/split/('pilot.jsonl' if pilot else 'predictions.jsonl');began=time.monotonic();result=list(prior)
    with target.open('x') as stream,ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize,initargs=(split,)) as pool:
        for row in prior:stream.write(json.dumps(row,default=serializable,allow_nan=False)+'\n')
        futures=[pool.submit(infer_one,i) for i in tasks]
        for future in as_completed(futures):
            row=future.result();result.append(row);stream.write(json.dumps(row,default=serializable,allow_nan=False)+'\n');stream.flush()
            if pilot or len(result)%16==0 or len(result)==check['frames']:
                print('SELECTION_INFER',split,len(result),len(prior)+len(tasks),row['index'],round(row['seconds'],2),flush=True)
    expected=len(tasks)+len(prior);assert len(result)==expected and len({r['index'] for r in result})==expected
    save(RECORD/('pilot-validation.json' if pilot else split+'-inference.json'),dict(passed=True,created=now(),frames=expected,
        sha256=digest(target),seconds=time.monotonic()-began,workers=workers,inherited_pilot=len(prior),truth_inputs=[],label_access_blocked=True,
        max_cost_replay_error=max(r['cost_replay_error'] for r in result)))
    print('SELECTION_INFER_COMPLETE',split,expected,flush=True)


def fit_gates():
    check_freeze();v=load(RECORD/'calibration-inference.json');assert v['passed'] and v['frames']==256
    assert digest(OUTPUT/'calibration/predictions.jsonl')==v['sha256']
    assert not (OUTPUT/'test').exists(),'Confirmation must not be generated before gate freeze'
    assert not (EXPERIMENT/'gates.json').exists()
    rr=rows(OUTPUT/'calibration/predictions.jsonl');matched=calibrate_gate([r['candidates'] for r in rr])
    old=load(EXPERIMENT/'original-gate.json');gates={}
    for origin,gate in (('original',old),('matched',matched)):
        for a in ALPHAS:gates[f'{origin}_{int(1000*a):03d}']=dict(gate,alpha=a)
    assert tuple(gates)==METHODS
    save(EXPERIMENT/'gates.json',gates)
    save(RECORD/'gate-validation.json',dict(passed=True,created=now(),sha256=digest(EXPERIMENT/'gates.json'),
        calibration_prediction_sha256=v['sha256'],fitted_bins=[len(s) for s in matched['scores']],test_generated=False,
        no_test_truth_used=True,original_gate_unchanged=True,methods=list(METHODS)))
    print('GATES_FROZEN',[len(s) for s in matched['scores']],flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=('prepare','tests','generate-calibration','pilot','infer-calibration','gates','generate-test','infer-test'))
    p.add_argument('--workers',type=int,default=6);args=p.parse_args()
    if args.command=='prepare':prepare()
    elif args.command=='tests':tests()
    elif args.command.startswith('generate-'):generate(args.command.split('-',1)[1])
    elif args.command=='pilot':infer('calibration',args.workers,pilot=True)
    elif args.command=='gates':fit_gates()
    else:infer(args.command.split('-',1)[1],args.workers)


if __name__=='__main__':main()
