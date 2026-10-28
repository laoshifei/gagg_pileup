"""Frozen diagnostic experiment: detection tradeoffs and later-Q profiles."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
import importlib.metadata
import json
import multiprocessing
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

from .alpha_pair_information_study import load,save,digest,rows
from .alpha_pair_evidence_profile import profile,fit_at_charge,seeds_for,rank_grid

PROJECT=Path(os.environ.get('GAGG_EVIDENCE_PROJECT',Path(__file__).resolve().parents[2])).resolve()
NAME='alpha_pair_evidence_profile_20261027'
EXPERIMENT=PROJECT/'pileup/experiments'/NAME
OUTPUT=PROJECT/'jupyter/artifacts'/NAME
RECORD=PROJECT/'jupyter/run_records'/NAME
PREVIOUS=PROJECT/'pileup/experiments/alpha_pair_physical_gate_20261027'
INPUT=PROJECT/'jupyter/artifacts/alpha_pair_physical_gate_20261027'
OLD_RECORD=PROJECT/'jupyter/run_records/alpha_pair_physical_gate_20261027'
DATA=PROJECT/'jupyter/artifacts/alpha_pair_risk_training_v2_20261027'


def paths(domain):
    return dict(observation=DATA/('native-observations.npz' if domain=='native' else 'evaluation/observations.npz'),
        truth=DATA/('native-labels.truth.npz' if domain=='native' else 'evaluation/labels.truth.npz'),
        cases=DATA/('native-labels.cases.json' if domain=='native' else 'evaluation/labels.cases.json'),
        candidates=INPUT/(domain+'-candidates.jsonl'))


def models():
    from .alpha_pair_joint_likelihood import ConservingResponse,JointNoise
    return ConservingResponse.load(EXPERIMENT/'joint-response.npz'),JointNoise.load(EXPERIMENT/'joint_conserved_nll-noise.npz')


def check_freeze():
    for name,value in load(EXPERIMENT/'freeze.json')['files'].items():assert digest(EXPERIMENT/name)==value,name


def prepare():
    prior=load(OLD_RECORD/'final-validation.json');assert prior['passed']
    notebook=PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    assert digest(notebook)==prior['notebook_sha256']
    for name,value in prior['source_sha256'].items():assert digest(PROJECT/name)==value,name
    for folder in (EXPERIMENT,OUTPUT,RECORD):folder.mkdir(parents=True,exist_ok=False)
    for source,target in ((notebook,'notebook.before.ipynb'),(PROJECT/'README.md','README.before.md'),(PROJECT/'docs/literature-use.md','literature.before.md')):
        shutil.copyfile(source,RECORD/target)
    for name in ('joint-response.npz','joint_conserved_nll-noise.npz','settings.json','gate-calibrations.json'):
        shutil.copyfile(PREVIOUS/name,EXPERIMENT/name)
    tasks=[]
    for domain in ('native','continuous'):
        cases=load(paths(domain)['cases'])
        ids=[i for i,c in enumerate(cases) if c['kind']=='pair' and c['condition']=='high_weak' and
            (c['gap_ns'] in (40.,200.) if domain=='native' else c['band'] in ('middle','long'))]
        assert len(ids)==80
        tasks.extend(dict(domain=domain,index=i) for i in ids)
    assert [t['index'] for t in tasks[:80]]==list(range(440,520))
    assert [t['index'] for t in tasks[80:]]==list(range(160,240))
    save(EXPERIMENT/'tasks.json',tasks)
    gates=load(EXPERIMENT/'gate-calibrations.json')
    protocol=dict(created=datetime.now(timezone.utc).isoformat(),question='Quantify the cost of relaxing count decisions and finite conditional later-charge uncertainty',
        scope='Pure alpha; fixed95/351ns; full1024x2ns analog voltage; ample dynamic range; original unknownK<=2 count study, separate forcedK2 profiles',
        data='280 prior native calibration,760 exposed native and1000 exposed continuous Gaussian evaluation; no new ROOTs, draws or training',
        threshold=dict(methods=['baseline','learned_gate','local_physical','mlp_physical'],alpha_grid=rank_grid(gates),
            rule='Only alpha changes; fixed caches, per-method calibration, presence, Q>=500 and physical positive-evidence guards remain',
            selected_alpha=.025,extra_limit=.05,interpretation='Full descriptive sweep; best under observed control limits is retrospective development envelope, not new validated default'),
        profile=dict(tasks=160,grid_count=33,base_charge_range_pe=[50.,50000.],inserted='500, originalQ2 and .5/.8/1/1.2/2 x hybridQ2 clipped to fit domain',
            time_domain='Both ordered times in [max(0,pilot_t-96),min(2046,pilot_t+512)] ns; pilot from observed K1. Cached out-of-domain candidates flagged, not silently assumed covered.',
            nuisance='For each fixed later Q2 optimize Q1 and both times; neither true times nor true Q enter primary profiles',
            optimization='Two best raw starts from cached fits, total-charge adjustment, warm continuation and0/4/8/16/32/64/128/256/512ns gaps;50 L-BFGS-B iterations; forward and reverse passes',
            envelope='Retain initial values and fixed-nuisance slice in numerical envelope; same candidate-dependent covariance and logdet. No claim of global convergence.',
            levels=[1.,4.,9.],interpretation='Finite compatibility contours, not68/95/99% intervals or posterior probabilities. Keep disjoint segments and truncation flags.',
            probe='After primary profiles are frozen, fix true later Q2 only in a separately labelled scoring probe, refit nuisance from observed seeds; no probe influences profiles, count gates or default fits',
            interpolation='LogQ piecewise-linear contour endpoints approximate; actual trueQ probe avoids using interpolation to measure truth inclusion',
            decision='No default change. Report reference-level truth inclusion, widths, boundary truncation, finite-search errors, and gate-selection conditioning separately.'),
        pilot=[tasks[i] for i in (0,40,80,120)],workers=8,blas_threads=1,
        metrics='Original strict pair Qwithin20%, twithinmin(10ns,gap/4), no extras; missing/extra/split/merge and matched bias/sigma68/time errors retained',
        figures=['0004-98-1','0004-98-2','0004-99-1','0004-100-1','0004-101-1'],
        examples=[tasks[i] for i in (0,40,80,120)],
        literature='Existing local Liddick2012 separate detection and parameter accuracy; Baker/Cousins1984 distinguish point estimation, interval coverage and goodness of fit. No direct Poisson histogram likelihood transplanted to correlated voltage.')
    save(EXPERIMENT/'protocol.json',protocol)
    input_files={str(p.relative_to(PROJECT)):digest(p) for domain in ('native','continuous') for p in paths(domain).values()}
    for name in ('scores.jsonl','report.json'):
        input_files[str((INPUT/name).relative_to(PROJECT))]=digest(INPUT/name)
    save(RECORD/'input-provenance.json',dict(inputs=input_files,prior_validation=prior,
        inherited_provenance_sha256=digest(OLD_RECORD/'input-provenance.json'),
        prior_freezes={name:digest(PREVIOUS/name) for name in ('freeze.json','report-freeze.json','delivery-freeze.json','delivery-v2-freeze.json')},
        notebook_sha256=digest(notebook),mother_groups=load(OLD_RECORD/'input-provenance.json')['mother_groups'],new_ROOTs=0,new_training=0,new_draws=0))
    save(RECORD/'environment.json',dict(python=sys.version,executable=sys.executable,workers=8,blas_threads=1,
        dependencies={k:importlib.metadata.version(k) for k in ('numpy','scipy','matplotlib','pandas','nbformat')},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        uncommitted='Prior information/allocation/risk/physical-gate work preserved; this experiment adds separately frozen diagnostic code'))
    shutil.copytree(PREVIOUS/'inference-source',EXPERIMENT/'inference-source')
    for name in ('alpha_pair_evidence_profile.py','alpha_pair_evidence_study.py','test_alpha_pair_evidence_profile.py'):
        shutil.copyfile(PROJECT/'pileup/analysis'/name,EXPERIMENT/'inference-source/pileup/analysis'/name)
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for name in ('EVIDENCE','GATE','RISK','INFORMATION'):os.environ['GAGG_'+name+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_evidence_study import main
    main()
''')
    files={str(p.relative_to(EXPERIMENT)):digest(p) for p in (EXPERIMENT/'inference-source').rglob('*.py')}
    for name in ('run.py','tasks.json','protocol.json','joint-response.npz','joint_conserved_nll-noise.npz','settings.json','gate-calibrations.json'):files[name]=digest(EXPERIMENT/name)
    save(EXPERIMENT/'freeze.json',dict(files=files,created=protocol['created']))
    print('EVIDENCE_PROFILE_PREPARED',len(files),len(protocol['threshold']['alpha_grid']),flush=True)


def initialize():
    global RESPONSE,NOISE,OBS,CACHE,PROFILES
    from .alpha_pair_charge_allocation import forbid_labels
    RESPONSE,NOISE=models();sys.addaudithook(forbid_labels)
    OBS={};CACHE={};PROFILES={}
    for domain in ('native','continuous'):
        with np.load(paths(domain)['observation']) as f:OBS[domain]=f['voltage_mV']
        CACHE[domain]={r['index']:r for r in rows(paths(domain)['candidates'])}
    if (OUTPUT/'profiles.jsonl').exists():PROFILES={(r['domain'],r['index']):r for r in rows(OUTPUT/'profiles.jsonl')}


def run_task(task):
    domain=task['domain'];i=task['index'];started=time.monotonic();row=CACHE[domain][i]
    result=profile(RESPONSE,NOISE,OBS[domain][i],row)
    maximum=0.
    from .alpha_pair_physical_gate import replay_cost
    for fit in result['fits']:
        actual=replay_cost(RESPONSE,NOISE,OBS[domain][i],fit)
        err=abs(actual-fit['cost'])/max(1.,abs(actual));maximum=max(maximum,err);assert err<1e-9
        assert result['bounds_ns'][0]-1e-8<=fit['times_ns'][0]<=fit['times_ns'][1]<=result['bounds_ns'][1]+1e-8
    return dict(result,**task,seconds=time.monotonic()-started,cost_replay_error=maximum)


def run_probe(task):
    domain=task['domain'];i=task['index'];q2=task['true_later_Q_pe'];row=CACHE[domain][i];curve=PROFILES[domain,i]
    best=min(curve['fits'],key=lambda r:r['cost'])
    nearest=sorted(curve['fits'],key=lambda r:abs(np.log(r['charge_pe'][1]/q2)))[:2]
    seeds=seeds_for(row['sets']['baseline'][1],row['sets']['baseline'][2],row['sets']['mlp'][2],q2,curve['bounds_ns'],best)+nearest
    started=time.monotonic();fit=fit_at_charge(RESPONSE,NOISE,OBS[domain][i],q2,seeds,curve['bounds_ns'])
    return dict(**task,fit=fit,delta_cost=fit['cost']-curve['reference_cost'],seconds=time.monotonic()-started,
        truth_assistance='Fixed true later Q only for scoring; true times and true Q1 not provided; no result returned to primary profile/gates')


def calculate(command,workers):
    check_freeze();protocol=load(EXPERIMENT/'protocol.json');provenance=load(RECORD/'input-provenance.json')
    for name,value in provenance['inputs'].items():assert digest(PROJECT/name)==value,name
    assert 1<=workers<=8
    tasks=protocol['pilot'] if command=='pilot' else load(EXPERIMENT/'tasks.json')
    pilot=command=='pilot';probe=command=='probe';prior=[]
    if not pilot:
        assert load(RECORD/'pilot-validation.json')['passed'] and load(RECORD/'numerical-tests.json')['passed']
    if command=='infer':
        prior=rows(OUTPUT/'pilot-profiles.jsonl');done={(r['domain'],r['index']) for r in prior}
        tasks=[r for r in tasks if (r['domain'],r['index']) not in done]
    if probe:
        assert load(RECORD/'inference.json')['passed']
        assert digest(OUTPUT/'profiles.jsonl')==load(RECORD/'inference.json')['sha256']
        truths={}
        for domain in ('native','continuous'):
            with np.load(paths(domain)['truth']) as f:truths[domain]=f['charge_pe']
        tasks=[dict(r,true_later_Q_pe=float(truths[r['domain']][r['index'],1])) for r in tasks]
    target=OUTPUT/('truth-probes.jsonl' if probe else 'pilot-profiles.jsonl' if pilot else 'profiles.jsonl')
    started=time.monotonic();results=list(prior)
    with target.open('x') as stream:
        for row in prior:stream.write(json.dumps(row,allow_nan=False)+'\n')
        with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn'),initializer=initialize) as pool:
            futures=[pool.submit(run_probe if probe else run_task,task) for task in tasks]
            for future in as_completed(futures):
                row=future.result();results.append(row);stream.write(json.dumps(row,allow_nan=False)+'\n');stream.flush()
                print('EVIDENCE_'+command.upper(),len(results),len(tasks)+len(prior),row['domain'],row['index'],round(row['seconds'],2),flush=True)
    assert len(results)==(4 if pilot else 160)
    if not probe:assert max(r['cost_replay_error'] for r in results)<1e-9
    save(RECORD/('probe-validation.json' if probe else 'pilot-validation.json' if pilot else 'inference.json'),
        dict(passed=True,frames=len(results),workers=workers,seconds=time.monotonic()-started,sha256=digest(target),
            source_sha256=digest(__file__),primary_profiles_sha256=digest(OUTPUT/'profiles.jsonl') if probe else None,
            primary_truth_inputs=[],probe_truth_input='Only fixed true later Q; scoring diagnostic' if probe else None,
            label_files_blocked_in_workers=True,maximum_cost_replay_error=None if probe else max(r['cost_replay_error'] for r in results),
            inherited_pilot_frames=len(prior),failed_trials=None if probe else sum(r['nonconverged_trials'] for r in results)))
    print('EVIDENCE_'+command.upper()+'_COMPLETE',len(results),round(time.monotonic()-started,2),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=('prepare','pilot','infer','probe'))
    parser.add_argument('--workers',type=int,default=8);args=parser.parse_args()
    if args.command=='prepare':prepare()
    else:calculate(args.command,args.workers)


if __name__=='__main__':main()
