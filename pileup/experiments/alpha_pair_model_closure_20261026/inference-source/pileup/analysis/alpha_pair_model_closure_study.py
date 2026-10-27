"""Matched-Gaussian versus exposed native alpha pairs, with frozen inference."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from importlib.metadata import version
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

from .alpha_pair_resolution_data import load, save, digest

PROJECT = Path(os.environ.get('GAGG_MODEL_PROJECT', Path(__file__).resolve().parents[2])).resolve()
NAME = 'alpha_pair_model_closure_20261026'
EXPERIMENT = PROJECT / 'pileup/experiments' / NAME
OUTPUT = PROJECT / 'jupyter/artifacts' / NAME
RECORD = PROJECT / 'jupyter/run_records' / NAME
PREVIOUS = PROJECT / 'pileup/experiments/alpha_pair_weak_search_20261026'
POOL = PROJECT / 'jupyter/artifacts/alpha_pair_weak_search_20261026'
CONDITIONS = ('mid_mid', 'high_weak', 'weak_high')
ENTROPY = (20261026, 7101)


def rows(path):
    return [json.loads(line) for line in Path(path).read_text().splitlines()]


def check_freeze():
    for name, sha in load(EXPERIMENT / 'freeze.json')['files'].items():
        assert digest(EXPERIMENT / name) == sha, name


def models():
    from .alpha_pair_joint_likelihood import ConservingResponse, JointNoise
    return (ConservingResponse.load(EXPERIMENT / 'joint-response.npz'),
            JointNoise.load(EXPERIMENT / 'joint_conserved_nll-noise.npz'))


def prepare():
    assert not any(p.exists() for p in (EXPERIMENT, OUTPUT, RECORD)), 'Use an independent destination'
    cases = load(POOL / 'labels.cases.json'); blocks = load(POOL / 'mother-blocks.truth.json')
    indices = [i for i, c in enumerate(cases) if c['dataset_role'] == 'latest'
               and (c['kind'] != 'pair' or c['condition'] in CONDITIONS)]
    assert len(indices) == 760
    with np.load(POOL / 'observations.npz') as f: y = f['voltage_mV']
    with np.load(POOL / 'labels.truth.npz') as f: q = f['charge_pe']; t = f['times_ns']
    baseline = {r['index']: r['candidates'] for r in rows(POOL / 'predictions.baseline.jsonl')}
    native_cases = [dict(cases[i], pool_index=i, native_index=j, block=blocks[i], domain='native', role='evaluation', draw=-1)
                    for j, i in enumerate(indices)]
    matched_cases = [dict(cases[i], pool_index=i, native_index=-1, block=-1, domain='matched', role='gate', draw=-1,
                          seed_entropy=[*ENTROPY, 0, i]) for i in range(280)]
    for draw in range(2):
        matched_cases += [dict(c, domain='matched', draw=draw, seed_entropy=[*ENTROPY, 1, draw, c['pool_index']])
                          for c in native_cases]
    assert len(matched_cases) == 1800
    groups = {}
    for i, c in enumerate(matched_cases):
        groups.setdefault((c['role'], c['draw'], c['kind'], c['condition'], c['gap_ns']), i)
    pilot = list(groups.values()); assert len(pilot) == 32
    diagnostics = [i for i, c in enumerate(matched_cases) if c['role'] == 'evaluation' and c['kind'] == 'pair'
                   and c['condition'] in ('high_weak', 'weak_high')]
    assert len(diagnostics) == 480
    boundary_groups = [('pair', 'mid_mid', 8.), ('pair', 'high_weak', 8.), ('pair', 'high_weak', 40.),
                       ('pair', 'high_weak', 200.), ('pair', 'weak_high', 8.), ('empty', 'all', None),
                       ('single', 'weak', None), ('single', 'high', None)]
    boundary_indices = [next(i for i, c in enumerate(matched_cases) if c['role'] == 'evaluation' and c['draw'] == 0
                             and (c['kind'], c['condition'], c['gap_ns']) == group) for group in boundary_groups]
    for folder in (EXPERIMENT, OUTPUT, RECORD): folder.mkdir(parents=True)
    shutil.copyfile(PROJECT / 'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb', RECORD / 'notebook.before.ipynb')
    for name in ('joint-response.npz', 'joint_conserved_nll-noise.npz', 'settings.json', 'baseline-gate.json'):
        shutil.copyfile(PREVIOUS / name, EXPERIMENT / name)
    np.savez_compressed(OUTPUT / 'native-observations.npz', voltage_mV=y[indices])
    np.savez_compressed(OUTPUT / 'native-labels.truth.npz', charge_pe=q[indices], times_ns=t[indices])
    source_ids = [c['pool_index'] for c in matched_cases]
    np.savez_compressed(OUTPUT / 'labels.truth.npz', charge_pe=q[source_ids], times_ns=t[source_ids])
    save(OUTPUT / 'native.cases.json', native_cases); save(OUTPUT / 'labels.cases.json', matched_cases)
    with (OUTPUT / 'native.predictions.jsonl').open('w') as stream:
        for j, i in enumerate(indices): stream.write(json.dumps(dict(index=j, pool_index=i, candidates=baseline[i])) + '\n')
    mapping = {i: j for j, i in enumerate(indices)}
    native_diagnostics = [dict(r, index=mapping[r['index']], pool_index=r['index'])
                          for r in rows(POOL / 'diagnostics.truth.jsonl') if r['kind'] == 'parameters' and r['index'] in mapping]
    assert len(native_diagnostics) == 240
    with (OUTPUT / 'native-diagnostics.truth.jsonl').open('w') as stream:
        for r in native_diagnostics: stream.write(json.dumps(r) + '\n')
    protocol = dict(created=datetime.now().astimezone().isoformat(),
        question='Separate finite-search/statistical limitations within the fitted Gaussian surrogate from native response/noise mismatch',
        scope='Pure alpha; frozen95/351ns; analog voltage; sufficient dynamic range; 2ns x1024; unknown local K<=2',
        native='760 previously exposed latest frames: mid_mid/high_weak/weak_high x8/40/200ns x40 pairs, plus80empty and4x80native singles; all development/diagnostic, no new optical parents',
        matched='Same numerical true Q/t parameters; two independent Gaussian draws per760 native frame;280 separate gate-mother parameter points each with one independent draw',
        meaning_of_Q='Native Q is realized avalanche-equivalent pe; matched Q is a conditioning parameter set to that numerical value, not a newly simulated avalanche charge or deposited energy',
        generation='Exact C=B-VV^T in the1021-dimensional FIR marginal; zero residual prefix of3 voltage samples; inverse fixed FIR; no clipping or waveform rescaling',
        seed_entropy=list(ENTROPY), rng='NumPy PCG64 / SeedSequence; role, draw and source index appended to entropy',
        unchanged='Frozen mean, covariance, logdet, waveform_candidates, fit bounds/budgets and Q>=500 selection guard; no retraining or extra weak-start search',
        gates='Primary: unchanged native280-single gate on both domains; secondary: matched280-single gate on both domains. Calibration only, no evaluation-label tuning',
        matched_inference_frames=1800, pilot_indices=pilot, diagnostic_indices=diagnostics,
        diagnostic_pilot=[next(i for i in diagnostics if matched_cases[i]['draw'] == draw and matched_cases[i]['condition'] == 'high_weak'
                               and matched_cases[i]['gap_ns'] == 40.) for draw in range(2)],
        diagnostics='Native240 cached and matched480: known K from existing candidate, true Q/t initialization radius48ns maxiter200, fixed true times with7 charge starts maxiter100; explicit truth assistance, no theoretical bound',
        boundary_indices=boundary_indices, boundary_prefixes_mV=[[.3, -.3, .3], [-.3, .3, -.3]],
        preflight='Five sampler tests;4096draws at4 prespecified parameter points; projected covariance, mean, whitening, dense algebra and FIR round-trip;32 inference pilots and2 diagnostics before full inference',
        success='Both trueQ>=500pe; eachQ error<=20%; each t error<=min(10ns,gap/4); no extra output; also retain misses/splits/merges/extra outputs and matched Q/t errors',
        bootstrap=dict(replicates=2000, seed=276107101,
            design='40 native mother blocks, first20/last20 stratified, common weights over conditions/gaps/domains; average2 matched noise draws per parameter point; conditional descriptive intervals'),
        residual_check='Truth-parameter whitening, quadratic/1021, empirical skewness/kurtosis/lag correlations; Gaussian comparison uses the frozen model, not fitted residual degrees of freedom',
        decision='Diagnostic only; no model or threshold selected for deployment, no automatic newROOT production; noisy Gaussian agreement is not a universal resolution limit',
        figures=['0004-72-1','0004-72-2','0004-73-1','0004-74-1','0004-75-1'],
        examples='Prespecified latest high_weak replicate0 at40 and200ns; matched draw0; full0-2048ns, conditional weak means and native avalanche components labelled separately',
        literature='Fixsen2004 existing reviewed-source note: conditional means/covariances, TES5000repeats perenergy, AI pending human review; exact Gaussian sampling and closure comparison are project implementations')
    save(EXPERIMENT / 'protocol.json', protocol)
    parent_provenance = load(PROJECT / 'jupyter/run_records/alpha_pair_weak_search_20261026/input-provenance.json')
    sources = [POOL / s for s in ('observations.npz','labels.truth.npz','labels.cases.json','mother-blocks.truth.json',
                                  'predictions.baseline.jsonl','diagnostics.truth.jsonl')]
    save(RECORD / 'input-provenance.json', dict(files={str(p.relative_to(PROJECT)): digest(p) for p in sources},
        prepared_files={p.name: digest(p) for p in OUTPUT.iterdir() if p.is_file()},
        mother_groups=parent_provenance['mother_groups'], native_pool_indices=indices,
        parent_freeze_sha256=digest(PREVIOUS / 'freeze.json'), new_optical_parents=0))
    save(RECORD / 'environment.json', dict(python=sys.version, executable=sys.executable,
        dependencies={n: version(n) for n in ('numpy','scipy','matplotlib','nbformat','pandas')},
        cpu_affinity=sorted(os.sched_getaffinity(0)), max_workers=8, blas_threads=1,
        git_head=subprocess.check_output(['git','rev-parse','HEAD'], cwd=PROJECT,text=True).strip(),
        uncommitted='Previous weak-search work retained; new sampler/test/study source snapshotted for this experiment',
        no_G4_or_detector_parameter_change=True, no_new_roots=True))
    shutil.copytree(PREVIOUS / 'inference-source', EXPERIMENT / 'inference-source')
    for filename in ('alpha_pair_model_sampler.py','test_alpha_pair_model_sampler.py','alpha_pair_model_closure_study.py'):
        shutil.copyfile(PROJECT / 'pileup/analysis' / filename, EXPERIMENT / 'inference-source/pileup/analysis' / filename)
    (EXPERIMENT / 'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for name in ('MODEL','SEARCH','JOINT','CONFIRM','NULL','CONDITIONAL','RESOLUTION','PAIR','COMPONENT'):
    os.environ['GAGG_'+name+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_model_closure_study import main
    main()
''')
    files = [p for p in (EXPERIMENT / 'inference-source').rglob('*') if p.is_file()]
    files += [EXPERIMENT / n for n in ('joint-response.npz','joint_conserved_nll-noise.npz','settings.json','baseline-gate.json','protocol.json','run.py')]
    save(EXPERIMENT / 'freeze.json', dict(created=datetime.now().astimezone().isoformat(),
        files={str(p.relative_to(EXPERIMENT)): digest(p) for p in files},
        native_physics='No new G4 invocation; native inputs retain their previous source/macro/RunInfo records'))
    print('MODEL_CLOSURE_PREPARED', len(native_cases), len(matched_cases), len(files), flush=True)


def generate():
    from .alpha_pair_model_sampler import GaussianFactor, conditional_mean, lift_innovations
    check_freeze(); assert not (OUTPUT / 'observations.npz').exists()
    response, noise = models(); cases = load(OUTPUT / 'labels.cases.json'); protocol = load(EXPERIMENT / 'protocol.json')
    with np.load(OUTPUT / 'labels.truth.npz') as f: q=f['charge_pe']; t=f['times_ns']
    waveform = np.empty((len(cases), 1024)); roundtrip = 0.; minimum_schur = 1.
    for i, c in enumerate(cases):
        valid = q[i] > 0.; curves = response.components(t[i,valid], q[i,valid])
        factor = GaussianFactor.from_noise(noise, curves, q[i,valid])
        minimum_schur = min(minimum_schur, 1.-float(factor.eigenvalues.max(initial=0.)))
        rng = np.random.default_rng(np.random.SeedSequence(c['seed_entropy']))
        z = factor.draw(rng); residual = lift_innovations(z, noise.coefficients)
        roundtrip = max(roundtrip, float(np.max(np.abs(noise.transform(residual)-z))))
        waveform[i] = conditional_mean(response, noise, t[i,valid], q[i,valid]) + residual
    assert np.isfinite(waveform).all() and roundtrip < 1e-9
    np.savez_compressed(OUTPUT / 'observations.npz', voltage_mV=waveform)
    boundary=[]; metadata=[]
    for index in protocol['boundary_indices']:
        for prefix in protocol['boundary_prefixes_mV']:
            delta = lift_innovations(np.zeros(1021), noise.coefficients, prefix)
            boundary.append(waveform[index] + delta)
            metadata.append(dict(main_index=index, prefix_mV=prefix, transform_max_error_mV=float(np.max(abs(noise.transform(delta)))),
                max_difference_after100ns_mV=float(np.max(abs(delta[response.sample_grid>=100.])))))
    np.savez_compressed(OUTPUT / 'boundary-observations.npz', voltage_mV=np.asarray(boundary))
    save(OUTPUT / 'boundary.cases.json', metadata)
    seeds=[tuple(c['seed_entropy']) for c in cases];assert len(set(seeds))==len(seeds)
    save(RECORD / 'generation.json', dict(passed=True, waveform_sha256=digest(OUTPUT / 'observations.npz'),
        boundary_waveform_sha256=digest(OUTPUT / 'boundary-observations.npz'), seed_streams_unique=True,
        matched_frames=len(cases), maximum_FIR_roundtrip_error_mV=roundtrip, minimum_schur_eigenvalue=minimum_schur,
        FIR_roots=np.abs(np.roots(noise.coefficients)).tolist(), zero_residual_prefix_samples=3,
        boundary_after100ns_max_mV=max(r['max_difference_after100ns_mV'] for r in metadata), no_clipping=True))
    print('MATCHED_GAUSSIAN_GENERATED', len(cases), 'FIR_ERROR', roundtrip, flush=True)


def preflight():
    from .alpha_pair_model_sampler import GaussianFactor, lift_innovations
    check_freeze(); response, noise=models(); cases=load(OUTPUT/'labels.cases.json')
    with np.load(OUTPUT/'labels.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
    selected=[next(i for i,c in enumerate(cases) if c['role']=='evaluation' and c['draw']==0 and
                   (c['kind'],c['condition'],c['gap_ns'])==group) for group in
              [('empty','all',None),('single','weak',None),('pair','high_weak',40.),('pair','mid_mid',8.)]]
    results=[]
    for j,i in enumerate(selected):
        ids=q[i]>0.;curves=response.components(t[i,ids],q[i,ids]);band,V=noise.covariance(curves,q[i,ids])
        factor=GaussianFactor.from_covariance(band,V);rng=np.random.default_rng(np.random.SeedSequence([*ENTROPY,99,j]))
        draws=factor.draw(rng,4096);white=factor.whiten(draws);n=factor.size
        probes=rng.normal(size=(n,16));probes/=np.linalg.norm(probes,axis=0)
        action=band[0,:,None]*probes
        for k in range(1,len(band)):
            action[k:]+=band[k,:n-k,None]*probes[:-k];action[:-k]+=band[k,:n-k,None]*probes[k:]
        action-=V@(V.T@probes)
        predicted=np.sum(probes*action,axis=0);projected=draws@probes
        relative_variance=projected.var(axis=0,ddof=1)/predicted-1.
        mean_z=projected.mean(0)/np.sqrt(predicted)
        expected_diagonal=band[0]-np.sum(V*V,axis=1)
        diagonal_rmse=float(np.sqrt(np.mean((draws.var(0,ddof=1)/expected_diagonal-1.)**2)))
        replay=noise.score(draws[0],noise.transform(curves),q[i,ids],transformed=True,details=True)
        assert abs(replay['quadratic']-float(white[0]@white[0]))<1e-8
        assert abs(replay['logdet']-factor.logdet)<1e-8
        assert abs(np.mean(white**2)-1.)<.005 and np.max(abs(relative_variance))<.12
        assert np.max(abs(mean_z))<.08 and diagonal_rmse<.04
        error=float(np.max(abs(noise.transform(lift_innovations(draws[0],noise.coefficients))-draws[0])))
        assert error<1e-9
        results.append(dict(index=i,kind=cases[i]['kind'],condition=cases[i]['condition'],gap_ns=cases[i]['gap_ns'],
            draws=4096,dimension=n,mean_white_square=float(np.mean(white**2)),max_projection_variance_error=float(np.max(abs(relative_variance))),
            max_projection_mean_z=float(np.max(abs(mean_z))),diagonal_relative_variance_rmse=diagonal_rmse,
            maximum_roundtrip_error_mV=error))
    save(RECORD/'sampler-preflight.json',dict(passed=True,results=results,diagnostic_only=True))
    print('SAMPLER_PREFLIGHT_PASSED',len(results),'x4096 draws',flush=True)


def initialize(blind=True):
    global RESPONSE,NOISE
    RESPONSE,NOISE=models()
    if blind:
        from .alpha_pair_conditional_study import forbid_labels
        sys.addaudithook(forbid_labels)


def fit_one(task):
    from .alpha_pair_joint_likelihood import waveform_candidates
    index, voltage=task;started=time.monotonic()
    return dict(index=index,candidates=waveform_candidates(RESPONSE,NOISE,voltage),seconds=time.monotonic()-started)


def infer(workers=4,pilot=False,boundary=False):
    check_freeze();assert load(RECORD/'sampler-preflight.json')['passed'];assert 1<=workers<=8
    filename='boundary-observations.npz' if boundary else 'observations.npz'
    with np.load(OUTPUT/filename) as f:y=f['voltage_mV']
    indices=load(EXPERIMENT/'protocol.json')['pilot_indices'] if pilot and not boundary else list(range(len(y)))
    path=OUTPUT/('boundary.predictions.jsonl' if boundary else 'matched.predictions.jsonl')
    existing=rows(path) if path.exists() else [];done={r['index'] for r in existing};assert len(done)==len(existing)
    started=time.monotonic();tasks=[i for i in indices if i not in done]
    with path.open('a') as stream:
        with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize) as pool:
            futures=[pool.submit(fit_one,(i,y[i])) for i in tasks]
            for n,future in enumerate(as_completed(futures),1):
                stream.write(json.dumps(future.result(),allow_nan=False)+'\n');stream.flush()
                if n%50==0:print('CLOSURE_INFER',n,'/',len(tasks),flush=True)
        os.fsync(stream.fileno())
    stage='boundary' if boundary else 'pilot' if pilot else 'matched'
    save(RECORD/(stage+'-inference.json'),dict(passed=True,frames=len(indices),new_frames=len(tasks),workers=workers,
        seconds=time.monotonic()-started,worker_labels_blocked=True,truth_inputs=[],
        prediction_sha256=digest(path),waveform_sha256=digest(OUTPUT/filename),freeze_sha256=digest(EXPERIMENT/'freeze.json')))
    print('CLOSURE_INFERENCE_COMPLETE',stage,len(indices),flush=True)


def diagnose_one(task):
    from .alpha_pair_profile_diagnostics import truth_diagnosis
    i,y,t,q,baseline=task;started=time.monotonic()
    return dict(index=i,seconds=time.monotonic()-started,**truth_diagnosis(RESPONSE,NOISE,y,t,q,baseline))


def diagnostics(workers=4,pilot=False):
    check_freeze();protocol=load(EXPERIMENT/'protocol.json')
    indices=protocol['diagnostic_pilot'] if pilot else protocol['diagnostic_indices']
    with np.load(OUTPUT/'observations.npz') as f:y=f['voltage_mV']
    with np.load(OUTPUT/'labels.truth.npz') as f:q=f['charge_pe'];t=f['times_ns']
    predictions={r['index']:r['candidates'] for r in rows(OUTPUT/'matched.predictions.jsonl')}
    assert all(i in predictions for i in indices)
    path=OUTPUT/'matched-diagnostics.truth.jsonl';existing=rows(path) if path.exists() else []
    done={r['index'] for r in existing};assert len(done)==len(existing);started=time.monotonic()
    tasks=[i for i in indices if i not in done]
    with path.open('a') as stream:
        with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize,initargs=(False,)) as pool:
            futures=[pool.submit(diagnose_one,(i,y[i],t[i],q[i],predictions[i])) for i in tasks]
            for n,future in enumerate(as_completed(futures),1):
                stream.write(json.dumps(future.result(),allow_nan=False)+'\n');stream.flush()
                if n%20==0:print('CLOSURE_DIAGNOSTIC',n,'/',len(tasks),flush=True)
        os.fsync(stream.fileno())
    save(RECORD/('diagnostic-pilot.json' if pilot else 'diagnostic-run.json'),dict(passed=True,frames=len(indices),new_frames=len(tasks),
        workers=workers,seconds=time.monotonic()-started,truth_used=True,output_sha256=digest(path)))
    print('CLOSURE_DIAGNOSTICS_COMPLETE',len(indices),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('prepare','generate','preflight','infer','diagnostics'))
    parser.add_argument('--workers',type=int,default=4);parser.add_argument('--pilot',action='store_true');parser.add_argument('--boundary',action='store_true')
    args=parser.parse_args()
    if args.command=='infer':infer(args.workers,args.pilot,args.boundary)
    elif args.command=='diagnostics':diagnostics(args.workers,args.pilot)
    else:globals()[args.command]()


if __name__=='__main__':main()
