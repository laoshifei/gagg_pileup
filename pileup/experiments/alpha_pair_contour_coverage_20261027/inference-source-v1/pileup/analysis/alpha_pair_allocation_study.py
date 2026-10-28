"""Frozen paired development comparison of component and allocation searches."""
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

import numpy as np

from .alpha_pair_information_study import load,save,digest,rows
from .alpha_pair_charge_allocation import VARIANTS,GAPS,SHARES

PROJECT=Path(os.environ.get('GAGG_ALLOCATION_PROJECT',Path(__file__).resolve().parents[2])).resolve()
NAME='alpha_pair_charge_allocation_20261027'
EXPERIMENT=PROJECT/'pileup/experiments'/NAME
OUTPUT=PROJECT/'jupyter/artifacts'/NAME
RECORD=PROJECT/'jupyter/run_records'/NAME
CLOSURE=PROJECT/'jupyter/artifacts/alpha_pair_model_closure_20261026'
POOL=PROJECT/'jupyter/artifacts/alpha_pair_weak_search_20261026'
MODEL=PROJECT/'pileup/experiments/alpha_pair_model_closure_20261026'


def models():
    from .alpha_pair_joint_likelihood import ConservingResponse,JointNoise
    return ConservingResponse.load(EXPERIMENT/'joint-response.npz'),JointNoise.load(EXPERIMENT/'joint_conserved_nll-noise.npz')


def check_freeze():
    for name,expected in load(EXPERIMENT/'freeze.json')['files'].items():assert digest(EXPERIMENT/name)==expected,name
    provenance=load(RECORD/'input-provenance.json')
    for name,expected in provenance['files'].items():assert digest(PROJECT/name)==expected,name
    for name,expected in provenance['prepared_files'].items():assert digest(OUTPUT/name)==expected,name


def prepare():
    assert not any(p.exists() for p in (EXPERIMENT,OUTPUT,RECORD))
    previous=load(PROJECT/'jupyter/run_records/alpha_pair_information_20261027/final-validation.json')
    assert previous['passed']
    notebook=PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    assert digest(notebook)==previous['notebook_sha256']
    for folder in (EXPERIMENT,OUTPUT,RECORD):folder.mkdir(parents=True)
    shutil.copyfile(notebook,RECORD/'notebook.before.ipynb')
    for name in ('joint-response.npz','joint_conserved_nll-noise.npz','settings.json','baseline-gate.json'):
        shutil.copyfile(MODEL/name,EXPERIMENT/name)
    native_cases=load(CLOSURE/'native.cases.json');matched_cases=load(CLOSURE/'labels.cases.json');pool_cases=load(POOL/'labels.cases.json')
    cases=[dict(c,domain='native',role='gate',native_index=-1,block=-1,draw=-1) for c in pool_cases[:280]]
    cases+=native_cases
    cases+=matched_cases[280:1040]
    assert len(cases)==1800 and all(c['role']=='gate' for c in cases[:280])
    assert all(c['domain']=='matched' and c['draw']==0 for c in cases[1040:])
    arrays={}
    for key,filename in [('voltage_mV','observations.npz'),('charge_pe','labels.truth.npz'),('times_ns','labels.truth.npz')]:
        with np.load(POOL/filename) as f:gate=f[key][:280]
        with np.load(CLOSURE/('native-'+filename)) as f:native=f[key]
        with np.load(CLOSURE/filename) as f:matched=f[key][280:1040]
        arrays[key]=np.concatenate([gate,native,matched])
    np.savez_compressed(OUTPUT/'observations.npz',voltage_mV=arrays['voltage_mV'])
    np.savez_compressed(OUTPUT/'labels.truth.npz',charge_pe=arrays['charge_pe'],times_ns=arrays['times_ns'])
    save(OUTPUT/'labels.cases.json',cases)
    old_pool={r['index']:r['candidates'] for r in rows(POOL/'predictions.baseline.jsonl')}
    old_native={r['index']:r['candidates'] for r in rows(CLOSURE/'native.predictions.jsonl')}
    old_matched={r['index']:r['candidates'] for r in rows(CLOSURE/'matched.predictions.jsonl')}
    baseline=[old_pool[i] for i in range(280)]+[old_native[i] for i in range(760)]+[old_matched[i] for i in range(280,1040)]
    with (OUTPUT/'predictions.baseline.jsonl').open('w') as stream:
        for i,candidates in enumerate(baseline):stream.write(json.dumps(dict(index=i,candidates=candidates))+'\n')
    groups={}
    for i,c in enumerate(cases):groups.setdefault((c['role'],c['domain'],c['kind'],c['condition'],c['gap_ns']),i)
    assert len(groups)==32
    source_groups=load(PROJECT/'jupyter/run_records/alpha_pair_model_closure_20261026/input-provenance.json')['mother_groups']
    active=[set(source_groups[k]) for k in ('training','gate','latest')]
    assert [len(s) for s in active]==[280,280,320]
    assert all(not active[i]&active[j] for i in range(3) for j in range(i))
    protocol=dict(created=datetime.now(timezone.utc).isoformat(),
        scope='Pure alpha; fixed95/351ns;2ns x1024 full analog; sufficient dynamic range; local unknown K<=2',
        question='Does rotating/rescaling the charge coordinates or profiling total charge improve finite blind search and strict event recovery?',
        input='280 independent native gate singles +760 exposed native evaluation +760 matched Gaussian draw0 evaluation; development/diagnostic only; no new optical mothers or noise draws',
        variants=list(VARIANTS),
        unchanged='Exact same conditional mean, candidate covariance, logdet, K0/K1 candidates, individual Q50..50000pe and time bounds; no prior or true labels',
        coordinates='t1/10,t2/10,100*ln(Q1+Q2),10*r; Q1=L(S)+(U(S)-L(S))*r; L=max(50,S-50000),U=min(50000,S-50)',
        coordinate_warning='Rotation/rescaling changes finite optimizer behavior and finite-difference displacements, not the likelihood or available information; no Jacobian added to an MLE objective',
        local='One70-iteration refinement of the original K2 candidate, time radius24ns, in each coordinate system; keep baseline if better',
        profile=dict(gaps_ns=list(GAPS),shares_first=list(SHARES),
            total='Numerical bounded scalar minimization of ln(total Q) at fixed times and share, max24 iterations,xatol1e-4; retain initial value too; full admissible charge domain',
            times='Observed single-candidate onset plus/minus each gap, chronological order; no true times',
            selection='Best share per gap/direction, best4 gaps per direction; same starts in both coordinates',
            budget='Eight20-iteration screens then two70-iteration polishes, plus the local-only result; same iteration limits, report actual evaluation/time differences',
            limits='Finite scalar/local search, no global-minimum claim; profile cost is shared preprocessing, additional to cached baseline'),
        gates='Primary: each method recalibrates its K2 evidence on the SAME280 native singles; same fitted-singleQ bins1500/5000/12500pe and plus-one/tie rank<=.025. Also report original-gate sensitivity without selecting by it.',
        frames=1800,pilot_indices=list(groups.values()),
        evaluation='Both trueQ>=500pe,eachQ within20%,each time within min(10ns,gap/4),no extra output; report all counts and matched errors; no failed-fit exclusion',
        selection=dict(native_extra_max=.05,native_single_success_min=.95,native_single_absolute_median_Q_bias_max=.05,
            native_high_weak40_gain_min=.10,native_high_weak200_no_drop=True,native_weak_high8_max_drop=.05,
            native_mid_mid_each_min=.95,native_six_weak_groups_mean_no_drop=True,
            matched='Diagnostic stress check only, not an independent confirmation or gate training source',
            ranking='If native constraints pass, largest native high_weak40+200 successes, then six-group weak mean, then listed variant order',
            decision='Provisional development candidate only; require continuous-Q/mismatch checks and fresh-parent confirmation before any default replacement'),
        bootstrap=dict(replicates=2000,seed=2026102703,design='40 parent blocks, first20/last20 stratified, common weights over methods/conditions/domains; descriptive paired development intervals; no independent-frame coverage claim'),
        examples='Preselected native high_weak replicate0 at40 and200ns; full0-2048ns; observed, conditional truth mean, blind selected fits and weak fitted components labeled',
        figures=['0004-82-1','0004-83-1','0004-83-2','0004-84-1','0004-85-1'],
        literature='Reuse Fixsen2004 conditional mean/covariance and Liddick2012 detection/precision distinction, AI notes pending review; present reparameterization/profile procedure is a project proposal')
    save(EXPERIMENT/'protocol.json',protocol)
    input_paths=[]
    for folder,names in [(POOL,['observations.npz','labels.truth.npz','labels.cases.json','predictions.baseline.jsonl']),
                         (CLOSURE,['native-observations.npz','native-labels.truth.npz','native.cases.json','native.predictions.jsonl','observations.npz','labels.truth.npz','labels.cases.json','matched.predictions.jsonl','scores.jsonl'])]:
        input_paths += [folder/name for name in names]
    save(RECORD/'input-provenance.json',dict(files={str(p.relative_to(PROJECT)):digest(p) for p in input_paths},
        prepared_files={p.name:digest(p) for p in OUTPUT.iterdir() if p.is_file()},mother_groups=source_groups,
        active_mothers_disjoint=True,notebook_before_sha256=digest(notebook),previous_validation=previous,
        native_Q='realized avalanche-equivalent pe',matched_Q='same numerical conditioning parameter, not new avalanche charge'))
    save(RECORD/'environment.json',dict(python=sys.version,executable=sys.executable,
        dependencies={n:importlib.metadata.version(n) for n in ('numpy','scipy','matplotlib','pandas','nbformat','ipython','uproot','awkward')},
        workers=4,cpu_affinity=sorted(os.sched_getaffinity(0)),blas_threads=1,
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        uncommitted='Prior information work retained; new allocation core, tests and study captured below; no detector/source change'))
    (RECORD/'git-status.before.txt').write_text(subprocess.check_output(['git','--no-optional-locks','-c','core.fileMode=false','status','--short','--','pileup/analysis','README.md','docs/literature-use.md'],cwd=PROJECT,text=True))
    response,noise=models()
    response.components([100.,140.],[18000.,700.])
    from . import alpha_pileup_fit
    loaded=set()
    for module in tuple(sys.modules.values()):
        filename=getattr(module,'__file__',None)
        if filename:
            p=Path(filename).resolve()
            if p.is_relative_to(PROJECT) and p.suffix=='.py' and 'analysis' in p.relative_to(PROJECT).parts:loaded.add(p)
    loaded.add(PROJECT/'pileup/analysis/test_alpha_pair_charge_allocation.py')
    loaded.add(PROJECT/'pileup/analysis/alpha_pair_allocation_study.py')
    for p in sorted(loaded):
        dest=EXPERIMENT/'inference-source'/p.relative_to(PROJECT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    for package in ('pileup','pileup/analysis','analysis'):
        dest=EXPERIMENT/'inference-source'/package/'__init__.py';dest.parent.mkdir(parents=True,exist_ok=True)
        if not dest.exists():dest.write_text('')
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for name in ('ALLOCATION','INFORMATION'):
    os.environ['GAGG_'+name+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_allocation_study import main
    main()
''')
    files={str(p.relative_to(EXPERIMENT)):digest(p) for p in (EXPERIMENT/'inference-source').rglob('*.py')}
    for name in ('joint-response.npz','joint_conserved_nll-noise.npz','settings.json','baseline-gate.json','protocol.json','run.py'):files[name]=digest(EXPERIMENT/name)
    save(EXPERIMENT/'freeze.json',dict(created=protocol['created'],files=files,source_capture='Loaded inference dependencies, core, study and tests; existing observations only'))
    print('ALLOCATION_PREPARED',1800,'frames',len(files),'frozen files',flush=True)


def initialize():
    global RESPONSE,NOISE
    RESPONSE,NOISE=models()
    from .alpha_pair_charge_allocation import forbid_labels
    sys.addaudithook(forbid_labels)


def fit_one(task):
    from .alpha_pair_charge_allocation import extend_allocation
    i,voltage,baseline=task;started=time.monotonic()
    result=extend_allocation(RESPONSE,NOISE,voltage,baseline)
    error=0.
    for variant in result['variants'].values():
        candidates=variant['candidates'];assert candidates[:2]==baseline[:2]
        if len(baseline)<3:continue
        best=candidates[2];assert best['cost']<=baseline[2]['cost']+1e-10
        q=np.asarray(best['charge_pe']);t=np.asarray(best['times_ns'])
        assert q.min()>=50.-1e-7 and q.max()<=50000.+1e-7 and np.isfinite(q).all()
        assert t.min()>=0. and t.max()<=RESPONSE.sample_grid[-1]
        curves=RESPONSE.components(t,q)
        replay=NOISE.score(np.asarray(voltage)-RESPONSE.baseline-NOISE.background_mV-curves.sum(0),curves,q)
        relative=abs(replay-best['cost'])/max(1.,abs(replay));error=max(error,relative)
        assert relative<1e-7
    return dict(index=i,seconds=time.monotonic()-started,maximum_cost_replay_error=error,**result)


def infer(workers=4,pilot=False):
    assert 1<=workers<=8
    check_freeze()
    if not pilot:assert load(RECORD/'pilot-validation.json')['passed']
    with np.load(OUTPUT/'observations.npz') as f:y=f['voltage_mV']
    baseline={r['index']:r['candidates'] for r in rows(OUTPUT/'predictions.baseline.jsonl')}
    indices=load(EXPERIMENT/'protocol.json')['pilot_indices'] if pilot else list(range(len(y)))
    destination=OUTPUT/'predictions.allocation.jsonl'
    previous=rows(destination) if destination.exists() else [];done={r['index'] for r in previous};assert len(done)==len(previous)
    pending=[i for i in indices if i not in done];started=time.monotonic()
    with destination.open('a') as stream:
        with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize) as pool:
            futures=[pool.submit(fit_one,(i,y[i],baseline[i])) for i in pending]
            for n,future in enumerate(as_completed(futures),1):
                stream.write(json.dumps(future.result(),allow_nan=False)+'\n');stream.flush()
                if n%20==0:os.fsync(stream.fileno());print('ALLOCATION',n,'/',len(pending),flush=True)
        os.fsync(stream.fileno())
    save(RECORD/('pilot-inference.json' if pilot else 'inference.json'),dict(passed=True,requested=len(indices),new_frames=len(pending),
        seconds=time.monotonic()-started,workers=workers,prediction_sha256=digest(destination),freeze_sha256=digest(EXPERIMENT/'freeze.json'),
        truth_inputs=[],worker_labels_blocked=True,new_ROOTs=0))
    print('ALLOCATION_COMPLETE',len(indices),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('prepare','infer'))
    parser.add_argument('--pilot',action='store_true');parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    if args.command=='prepare':prepare()
    else:infer(args.workers,args.pilot)


if __name__=='__main__':main()
