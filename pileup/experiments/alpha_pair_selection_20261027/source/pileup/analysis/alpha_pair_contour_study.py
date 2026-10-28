"""Frozen contour refinement before an independent model-coverage pilot."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
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

from .alpha_pair_information_study import load,save,digest,rows,serializable
from .alpha_pair_stability_study import copy_chunks
from .alpha_pair_evidence_study import paths
from .alpha_pair_profile_contours import converge_contours
from .alpha_pair_profile_stability import point_cost

PROJECT=Path(os.environ.get('GAGG_CONTOUR_PROJECT',Path(__file__).resolve().parents[2])).resolve()
NAME='alpha_pair_contour_coverage_20261027'
EXPERIMENT=PROJECT/'pileup/experiments'/NAME
OUTPUT=PROJECT/'jupyter/artifacts'/NAME
RECORD=PROJECT/'jupyter/run_records'/NAME
PREVIOUS=PROJECT/'pileup/experiments/alpha_pair_profile_stability_20261027'
INPUT=PROJECT/'jupyter/artifacts/alpha_pair_profile_stability_20261027'
OLD_RECORD=PROJECT/'jupyter/run_records/alpha_pair_profile_stability_20261027'
NOTEBOOK=PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
SETTINGS=dict(levels=[1.,4.,9.],bracket_fraction=.01,endpoint_fraction=.01,reference_tolerance=.01,
    max_rounds=12,local=dict(starts=6,maxiter=250),grid=dict(starts=4,maxiter=150),verify=dict(starts=8,maxiter=300))
TASKS=[dict(domain=d,index=i) for d,i in (('continuous',160),('continuous',190),('native',440),('native',480),('continuous',200))]


def prepare():
    prior=load(OLD_RECORD/'final-validation.json');assert prior['passed']
    assert digest(NOTEBOOK)==prior['notebook_sha256']
    for name,value in prior['source_sha256'].items():assert digest(PROJECT/name)==value,name
    for name,value in prior['documentation_sha256'].items():assert digest(PROJECT/name)==value,name
    for folder in (EXPERIMENT,OUTPUT,RECORD):folder.mkdir(parents=True,exist_ok=False)
    for source,target in ((NOTEBOOK,'notebook.before.ipynb'),(PROJECT/'README.md','README.before.md'),(PROJECT/'docs/literature-use.md','literature.before.md')):
        copy_chunks(source,RECORD/target)
    for name in ('joint-response.npz','joint_conserved_nll-noise.npz'):copy_chunks(PREVIOUS/name,EXPERIMENT/name)
    protocol=dict(created=datetime.now(timezone.utc).isoformat(),question='Resolve two changing contour boundaries before repeated-sampling coverage checks',
        scope='Pure alpha,fixed95/351ns,1024x2ns full analog voltage,knownK2,ordered laterQ2. Same mean,covariance,logdet,charge50-50000pe and observed-trigger local time domain.',
        settings=SETTINGS,tasks=TASKS,
        selection='Two previously failed Gaussian middle cases160/190 plus first native40/native200/Gaussianlong controls440/480/200. Selected numerical diagnostics, not representative recovery statistics.',
        method='Refit sampled crossing endpoints and minimum neighbours on unchanged grid; bisect crossings and minimum cells until logQ bracket width<=log1p(.01); verify boundary nodes with eight observed starts,300 iterations. Maximum12 rounds.',
        stop='After stronger verification, every level has same component count and edge flags,endpoint factor<=1.01,reference gain<=.01,and all crossing brackets<=1%. Failure retained rather than forcing success.',
        diagnosis='Compare same-grid refit gains with newly evaluated midpoint cost versus old logQ-linear interpolation; does not prove a unique physical cause or a global optimum.',
        coverage='Only after all five diagnostics pass, freeze a separate model-only repeated-sampling protocol and source before generating fresh noise. No adaptive threshold chosen using that evaluation.',
        workers_limit=8,blas_threads=1,truth='No truth/cases/ROOT access in primary workers; scoring only after hashes freeze.',
        new_ROOTs=0,new_training=0,seed=None,seed_note='PhaseA deterministic old-data refinement; later independent-noise design gets separate seed freeze.',
        figures=['0004-108-1','0004-109-1','0004-110-1','0004-111-1'])
    save(EXPERIMENT/'protocol.json',protocol)
    inherited=load(OLD_RECORD/'input-provenance.json')
    inputs=dict(inherited['inputs'])
    for name in ('profiles.jsonl','report.json','profile-scores.json','truth-probes.jsonl'):
        inputs[str((INPUT/name).relative_to(PROJECT))]=digest(INPUT/name)
    save(RECORD/'input-provenance.json',dict(inputs=inputs,prior_validation=prior,inherited_sha256=digest(OLD_RECORD/'input-provenance.json'),
        prior_freezes={n:digest(PREVIOUS/n) for n in ('freeze.json','report-freeze.json','delivery-freeze.json')},notebook_sha256=digest(NOTEBOOK),mother_groups=inherited['mother_groups']))
    status=subprocess.check_output(['git','-c','core.fileMode=false','status','--short','--untracked-files=no','--','README.md','docs/literature-use.md','pileup/analysis'],cwd=PROJECT,text=True,timeout=60)
    save(RECORD/'environment.json',dict(created=protocol['created'],python=sys.version,executable=sys.executable,
        dependencies={k:importlib.metadata.version(k) for k in ('numpy','scipy','matplotlib','pandas','nbformat')},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),scoped_git_status=status,
        uncommitted='Earlier studies retained; exact source snapshots and documentation/notebook backups capture state. No commit or push.',workers_limit=8,blas_threads=1))
    shutil.copytree(PREVIOUS/'inference-source',EXPERIMENT/'inference-source')
    for name in ('alpha_pair_profile_contours.py','test_alpha_pair_profile_contours.py','alpha_pair_contour_study.py'):
        copy_chunks(PROJECT/'pileup/analysis'/name,EXPERIMENT/'inference-source/pileup/analysis'/name)
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for key in ('CONTOUR','STABILITY','EVIDENCE','GATE','RISK','INFORMATION'):os.environ['GAGG_'+key+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_contour_study import main
    main()
''')
    files={str(p.relative_to(EXPERIMENT)):digest(p) for p in (EXPERIMENT/'inference-source').rglob('*.py')}
    for n in ('run.py','protocol.json','joint-response.npz','joint_conserved_nll-noise.npz'):files[n]=digest(EXPERIMENT/n)
    save(EXPERIMENT/'freeze.json',dict(files=files,created=protocol['created']))
    print('CONTOUR_PREPARED',len(files),flush=True)


def check_freeze():
    for name,value in load(EXPERIMENT/'freeze.json')['files'].items():assert digest(EXPERIMENT/name)==value,name


def models():
    from .alpha_pair_joint_likelihood import ConservingResponse,JointNoise
    return ConservingResponse.load(EXPERIMENT/'joint-response.npz'),JointNoise.load(EXPERIMENT/'joint_conserved_nll-noise.npz')


def initialize():
    global RESPONSE,NOISE,OBS,CACHE,OLD,CONFIG
    from .alpha_pair_charge_allocation import forbid_labels
    RESPONSE,NOISE=models();sys.addaudithook(forbid_labels)
    OBS={};CACHE={}
    for domain in ('native','continuous'):
        with np.load(paths(domain)['observation']) as f:OBS[domain]=f['voltage_mV']
        CACHE[domain]={r['index']:r for r in rows(paths(domain)['candidates'])}
    OLD={(r['domain'],r['index']):r['stages']['contour'] for r in rows(INPUT/'profiles.jsonl')}
    CONFIG=load(EXPERIMENT/'protocol.json')['settings']


def task(t):
    d,i=t['domain'],t['index']
    def progress(r):
        print('CONTOUR_PROGRESS',d,i,r['stage'],r['points'],'unresolved',r['quality']['unresolved_cells'],flush=True)
    r=converge_contours(RESPONSE,NOISE,OBS[d][i],CACHE[d][i],OLD[d,i],CONFIG,progress=progress)
    error=0.
    for c in r['stages'].values():
        assert c['reference_cost']==min(f['cost'] for f in c['fits'])
        for f in c['fits']:
            actual=point_cost(RESPONSE,NOISE,OBS[d][i],f)
            error=max(error,abs(actual-f['cost'])/max(1.,abs(actual)))
            assert c['bounds_ns'][0]-1e-8<=f['times_ns'][0]<=f['times_ns'][1]<=c['bounds_ns'][1]+1e-8
    assert error<1e-9
    return dict(r,**t,cost_replay_error=error)


def tests():
    check_freeze();suite=unittest.TestSuite()
    for n in ('test_alpha_pair_profile_stability','test_alpha_pair_profile_contours'):
        suite.addTests(unittest.defaultTestLoader.loadTestsFromName('pileup.analysis.'+n))
    result=unittest.TextTestRunner(verbosity=2).run(suite);assert result.wasSuccessful()
    save(RECORD/'numerical-tests.json',dict(passed=True,tests=result.testsRun,core_sha256=digest(EXPERIMENT/'inference-source/pileup/analysis/alpha_pair_profile_contours.py')))


def run(workers):
    check_freeze();assert load(RECORD/'numerical-tests.json')['passed'] and 1<=workers<=8
    for name,value in load(RECORD/'input-provenance.json')['inputs'].items():assert digest(PROJECT/name)==value,name
    began=time.monotonic();result=[];target=OUTPUT/'refined-profiles.jsonl'
    with target.open('x') as stream,ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize) as pool:
        futures=[pool.submit(task,t) for t in load(EXPERIMENT/'protocol.json')['tasks']]
        for f in as_completed(futures):
            r=f.result();result.append(r);stream.write(json.dumps(r,default=serializable,allow_nan=False)+'\n');stream.flush()
            print('CONTOUR_DONE',len(result),r['domain'],r['index'],r['stable'],round(r['seconds'],1),flush=True)
    save(RECORD/'refinement-validation.json',dict(passed=True,frames=len(result),stable=sum(r['stable'] for r in result),
        seconds=time.monotonic()-began,workers=workers,sha256=digest(target),maximum_cost_replay_error=max(r['cost_replay_error'] for r in result),
        label_files_blocked_in_workers=True,primary_truth_inputs=[],coverage_gate_passed=all(r['stable'] for r in result)))
    print('CONTOUR_REFINEMENT_COMPLETE',len(result),sum(r['stable'] for r in result),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('prepare','tests','run'));p.add_argument('--workers',type=int,default=5)
    args=p.parse_args()
    if args.command=='prepare':prepare()
    elif args.command=='tests':tests()
    else:run(args.workers)


if __name__=='__main__':main()
