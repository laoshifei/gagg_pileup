"""Frozen development experiment for a causal energy-response surrogate."""
from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import hashlib
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
from scipy.optimize import minimize_scalar
from scipy.signal import lfilter

from .alpha_pair_recovery_energy import RecoveryEnergyResponse, coupled_charge_inverse, refine_energy
from .sipm_energy_calibration import score_events, error_summary

PROJECT = Path(os.environ.get('GAGG_RECOVERY_ENERGY_PROJECT', Path(__file__).resolve().parents[2])).resolve()
NAME = 'alpha_pair_recovery_energy_20261029'
EXPERIMENT = PROJECT/'pileup/experiments'/NAME
OUTPUT = PROJECT/'jupyter/artifacts'/NAME
RECORD = PROJECT/'jupyter/run_records'/NAME
GAPS = [8., 40., 200.]
VARIANTS = ['static_nll', 'coupled_charge_only', 'isolated_energy_refit', 'recovery_energy_refit']


def load(path): return json.loads(Path(path).read_text())
def lines(path): return [json.loads(s) for s in Path(path).read_text().splitlines()]
def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name+'.partial')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n'); temp.replace(path)
def write_lines(path, rows):
    with Path(path).open('w') as stream:
        for row in rows: stream.write(json.dumps(row, ensure_ascii=False, allow_nan=False)+'\n')
def freeze(name, files, **metadata):
    save(EXPERIMENT/name, dict(created=datetime.now().astimezone().isoformat(),
        files={str(p.relative_to(EXPERIMENT)): digest(p) for p in files}, **metadata))
def verify(name='freeze.json'):
    for path, sha in load(EXPERIMENT/name)['files'].items():
        assert digest(EXPERIMENT/path) == sha, path


def prepare():
    from .alpha_pair_joint_likelihood import ConservingResponse, JointNoise
    for folder in [EXPERIMENT, OUTPUT, RECORD]: folder.mkdir(parents=True, exist_ok=False)
    train = PROJECT/'jupyter/artifacts/alpha_pair_joint_likelihood_20261026'
    history = PROJECT/'jupyter/artifacts/sipm_energy_calibration_v2_20261028'
    guard = PROJECT/'pileup/experiments/sipm_energy_guardband_20261029'
    old = PROJECT/'pileup/experiments/alpha_pair_joint_confirmation_20261026'
    inputs = OUTPUT/'inputs'; inputs.mkdir()
    sources = {'development.truth.npz': train/'calibration-ancestry.truth.npz',
               'development.cases.json': train/'calibration.cases.json',
               'development.energy.json': history/'calibration-records.json',
               'evaluation.labels.npz': history/'evaluation.labels.npz',
               'evaluation.cases.json': history/'evaluation.cases.json',
               'observations.npz': PROJECT/'jupyter/artifacts/alpha_pair_joint_confirmation_20261026/observations.npz'}
    for target, source in sources.items(): shutil.copyfile(source, inputs/target)
    cases = load(inputs/'development.cases.json')
    records = load(inputs/'development.energy.json')
    roots = [{r['root'] for r in records if r['fold'] == f} for f in [0, 1]]
    testroots = {r['root'] for r in load(history/'evaluation-parent-labels.json')}
    assert len(roots[0]) == len(roots[1]) == 280 and len(testroots) == 320
    assert not roots[0]&roots[1] and not (roots[0]|roots[1])&testroots
    original = ConservingResponse.load(old/'joint-response.npz')
    noise = JointNoise.load(old/'joint_conserved_nll-noise.npz')
    np.savez_compressed(EXPERIMENT/'template.npz', base_grid=original.base.grid, base_kernel=original.base.kernel,
        base_area=original._base_area, correction_grid=original.grid, mean_coefficients=original.coefficients,
        gain=original.base.area_per_pe, baseline_mV=original.baseline, ar=noise.coefficients,
        noise_lags=noise.noise_lags, noise_parameters=noise.parameters, lag1=noise.lag1,
        cross=noise.cross, background_mV=noise.background_mV)
    for name in ['energy-model.json', 'sensor-config.json']: shutil.copyfile(guard/name, EXPERIMENT/name)
    sensor = load(EXPERIMENT/'sensor-config.json')['sensor']
    save(EXPERIMENT/'initial-response-parameters.json', dict(strength=0., recovery_ns=sensor['recovery_ns'],
        cells=int(sensor['cells_x']*sensor['cells_y']),
        role='Uncalibrated zero-neighbor baseline; physical detector unchanged'))
    response = RecoveryEnergyResponse.load(EXPERIMENT, strength=0.)
    checks = []
    for times, q in [([128.], [700.]), ([128., 168.], [17000., 700.]), ([128., 328.], [2700., 17000.])]:
        predicted = response.legacy_components(times, q)
        expected = original.components(times, q)
        np.testing.assert_allclose(predicted, expected, rtol=1e-12, atol=1e-10)
        # Verify the copied banded weighting against the old implementation
        # with charge-total conditioning deliberately disabled.
        from dataclasses import replace
        band, _ = replace(noise, kappa=0.).covariance(expected, q)
        from scipy.linalg import cholesky_banded
        np.testing.assert_allclose(response.fixed_metric(dict(times_ns=times, charge_pe=q)),
                                   cholesky_banded(band, lower=True), rtol=1e-11, atol=1e-11)
        checks.append(dict(times_ns=times, charge_pe=q, compatible=True))
    starts_source = PROJECT/'jupyter/artifacts/sipm_energy_guardband_20261029/historical-energy.predictions.jsonl'
    source_rows = [r for r in lines(starts_source) if r['method']=='joint_conserved_nll']
    assert len(source_rows) == 1600 and {r['index'] for r in source_rows} == set(range(1600))
    write_lines(inputs/'blind-starts.jsonl', [dict(index=r['index'], times_ns=r['times_ns'], charge_pe=r['charge_pe']) for r in source_rows])
    write_lines(inputs/'static.predictions.jsonl', [dict(index=r['index'], times_ns=r['times_ns'], charge_pe=r['charge_pe'],
        energy_MeV=r['energy_MeV'], method='static_nll') for r in source_rows])
    evaluation_cases = load(inputs/'evaluation.cases.json')
    pilot_indices = []
    for kind, condition, gap in [('empty','empty',None), ('single','weak',None), ('single','high',None),
                                  ('pair','high_weak',8.), ('pair','weak_high',40.), ('pair','high_high',8.)]:
        eligible = [i for i,c in enumerate(evaluation_cases) if c['kind']==kind and
                    (kind=='empty' or c['condition']==condition) and c['gap_ns']==gap]
        assert eligible, (kind, condition, gap)
        pilot_indices.append(eligible[0])
    protocol = dict(created=datetime.now().astimezone().isoformat(), git_baseline=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        scope='Development prototype and retrospective replay; not a fresh independent transport confirmation',
        default_algorithm_changed=False, new_ROOTs=0, particles='Pure central alpha, fixed95/351ns, local K<=2, full analog range',
        response='Frozen isolated conserving mean at q0(E); fixed electronics AR3; causal 29ns neighbor-charge proxy; multiply innovations by exp(-strength*neighbor); inverse AR3 to voltage',
        hypothesis='Project empirical surrogate, not Rosado Eq19 or microscopic state inference',
        physical_config_changed=False, calibrate='One strength in [0,8], 900 fold0 pair frames; equal condition/gap/component strata; fractional Q residual at known development E/t',
        train_mothers=280, heldout_development_mothers=280, historical_evaluation_mothers=320,
        model_validation='900 fold1 pairs plus280 singles; no strength selection or retuning on heldout/evaluation data',
        waveform_refit='E/t jointly fitted; original BLIND selected K retained; radius24ns,60 L-BFGS-B iterations',
        weights='Fixed old banded noise at initial blind Q/t; kappa=0 (unconditioned total); same weights for both new waveform fits; no new likelihood or CI claim',
        methods=VARIANTS, endpoints='E within10%, time<=min(10ns,gap/4), no extras for pair success; all Qtrue>=500pe; report counts, bias, sigma68, split/merge',
        energy_bounds='Frozen guard-calibration domain; clamping only for starts, optimizer boundary hits retained explicitly',
        gate='No new event-count gate; original blind K is a fixed design condition, not true K',
        pilot_indices=pilot_indices, runtime_gate='All invariants and truth-read guards pass; projected1600-frame two-worker runtime <=1800s',
        native_expansion_gate='Separate future experiment only after heldout weak charge RMSE improves>=10% and blind energy controls are assessed',
        full_samples=3072, observed_samples=1024, sample_spacing_ns=2., first_sample_ns=.25,
        max_workers=2, no_new_training_network=True,
        restrictions=['Existing optical mothers already exposed','No gamma/free tau/high rate/K>=3','No blind count improvement claimed','No confidence interval from this fixed weighting metric'])
    save(EXPERIMENT/'protocol.json', protocol)
    save(EXPERIMENT/'input-hashes.json', {str(p.relative_to(OUTPUT)): digest(p) for p in inputs.iterdir()})
    save(RECORD/'provenance.json', dict(sources={name:dict(path=str(p.relative_to(PROJECT)),sha256=digest(p)) for name,p in sources.items()},
        blind_start_source=str(starts_source.relative_to(PROJECT)), blind_start_sha256=digest(starts_source),
        train_roots=sorted(roots[0]), heldout_roots=sorted(roots[1]), evaluation_roots=sorted(testroots), mother_disjoint=True))
    save(RECORD/'compatibility.json', dict(passed=True, old_mean_and_unconditioned_metric_reproduced=checks))
    save(RECORD/'environment.json', dict(python=sys.version, executable=sys.executable,
        packages={name:version(name) for name in ['numpy','scipy','matplotlib','nbformat']}, threads=1, max_workers=2))
    src=EXPERIMENT/'source/pileup/analysis'; src.mkdir(parents=True)
    (src.parent/'__init__.py').write_text(''); (src/'__init__.py').write_text('')
    for name in ['alpha_pair_recovery_energy.py','alpha_pair_recovery_energy_study.py','sipm_energy_calibration.py']:
        shutil.copyfile(PROJECT/'pileup/analysis'/name, src/name)
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
os.environ['GAGG_RECOVERY_ENERGY_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_recovery_energy_study import main
    main()
''')
    freeze('freeze.json', [p for p in EXPERIMENT.rglob('*') if p.is_file()], calibration_started=False)
    for name in ['codex_0007_SiPM_Energy_Calibration.ipynb','Research_Plan_20260921.ipynb']:
        shutil.copyfile(PROJECT/'jupyter'/name, RECORD/(name+'.before'))
    print('RECOVERY_PREPARED',len(cases),'development frames;1600 historical replay frames',flush=True)


def calibrate():
    verify(); assert not (EXPERIMENT/'response-parameters.json').exists()
    response = RecoveryEnergyResponse.load(EXPERIMENT, strength=0.)
    cases=load(OUTPUT/'inputs/development.cases.json')
    records={r['root']:r for r in load(OUTPUT/'inputs/development.energy.json')}
    with np.load(OUTPUT/'inputs/development.truth.npz') as f:
        true_q=f['charge_pe']; true_t=f['times_ns']; components=f['components_mV']
    indices=[i for i,c in enumerate(cases) if c['kind']=='pair' and c['gap_ns'] in GAPS]
    assert len(indices)==1800
    def build(ids):
        zz=[];hh=[];q0=[];keys=[]
        for i in ids:
            c=cases[i]; energies=[records[p]['edep_MeV'] for p in c['mother_roots']]
            z,h,q=response.drivers(true_t[i],energies)
            zz.append(z.astype('f4'));hh.append(h.astype('f4'));q0.append(q)
            keys.append((c['condition'],c['gap_ns']))
        return np.asarray(zz),np.asarray(hh),np.asarray(q0),keys
    train_ids=[i for i in indices if cases[i]['fold']==0]
    z,h,q0,keys=build(train_ids)
    groups=defaultdict(list)
    for j,key in enumerate(keys):groups[key].append(j)
    weights=np.zeros(len(z))
    for ids in groups.values():weights[ids]=1/(len(groups)*len(ids))
    target=true_q[train_ids]
    def loss(k):
        predicted=np.sum(z*np.exp(-k*h),axis=-1)/response._innovation_per_pe
        return float(np.sum(weights[:,None]*((predicted-target)/q0)**2)/2)
    initial_training_loss=loss(0.)
    opt=minimize_scalar(loss,bounds=(0.,8.),method='bounded',options=dict(xatol=1e-5,maxiter=60))
    assert opt.success and np.isfinite(opt.fun)
    strength=float(opt.x if opt.fun<loss(0.) else 0.)
    params=load(EXPERIMENT/'initial-response-parameters.json')
    params.update(strength=strength, role='Calibrated once on fold0 true component charges',
        training_frames=len(train_ids), objective='Stratum-balanced mean fractional Q MSE', boundary=bool(strength<1e-4 or strength>7.9999))
    save(EXPERIMENT/'response-parameters.json',params)
    freeze('method-freeze.json',[EXPERIMENT/'freeze.json',EXPERIMENT/'response-parameters.json',EXPERIMENT/'template.npz',EXPERIMENT/'energy-model.json'],
           inference_started=False, strength_fit_before_heldout_evaluation=True)
    results=[]
    for fold in [0,1]:
        ids=[i for i in indices if cases[i]['fold']==fold]
        if fold==1:z,h,q0,keys=build(ids)
        for k,method in [(0.,'isolated'),(strength,'recovery')]:
            predicted_z=z*np.exp(-k*h)
            qpred=np.sum(predicted_z,axis=-1)/response._innovation_per_pe
            target_z=lfilter(response.ar,[1.],components[ids],axis=-1)
            relative_innovation_mse=np.sum((predicted_z[:,:,:1024]-target_z)**2,axis=-1)/np.maximum(np.sum(target_z**2,axis=-1),1e-12)
            relative=(qpred-true_q[ids])/true_q[ids]
            for j,i in enumerate(ids):
                for component in [0,1]:
                    results.append(dict(index=i,fold=fold,condition=cases[i]['condition'],gap_ns=cases[i]['gap_ns'],
                        component=component,method=method,true_charge_pe=float(true_q[i,component]),predicted_charge_pe=float(qpred[j,component]),
                        relative_charge_error=float(relative[j,component]),relative_innovation_mse=float(relative_innovation_mse[j,component])))
        print('RECOVERY_MEAN_VALIDATION',fold,len(ids),flush=True)
    write_lines(OUTPUT/'mean-charge-diagnostics.jsonl',results)
    summary=[]
    for fold in [0,1]:
        for group in ['weak','strong']:
            for method in ['isolated','recovery']:
                rows=[r for r in results if r['fold']==fold and r['method']==method and ('weak' in r['condition'])==(group=='weak')]
                summary.append(dict(fold=fold,group=group,method=method,n=len(rows),
                    relative_rmse=float(np.sqrt(np.mean([r['relative_charge_error']**2 for r in rows]))),
                    relative_innovation_rmse=float(np.sqrt(np.mean([r['relative_innovation_mse'] for r in rows]))),
                    median_relative=float(np.median([r['relative_charge_error'] for r in rows]))))
    save(OUTPUT/'calibration-report.json',dict(strength=strength, training_loss_before=initial_training_loss,
        training_optimum_loss=float(opt.fun), summary=summary, selected_on_fold0=True,
        heldout_not_used_for_tuning=True, parameter_at_boundary=params['boundary']))
    print('RECOVERY_CALIBRATED',strength,summary,flush=True)


def initialize():
    global RESPONSE
    RESPONSE=RecoveryEnergyResponse.load(EXPERIMENT)
    def guard(event,args):
        if event!='open' or not args:return
        value=str(args[0]).lower()
        if any(word in value for word in ['.truth','.labels','/truth/','/labels/','.root','development.energy','development.cases','evaluation.cases']):
            raise PermissionError('Inference cannot read truth: '+value)
    sys.addaudithook(guard)
    for path in ['forbidden.truth.npz','forbidden.labels.npz','forbidden.root']:
        try:open(path,'rb')
        except PermissionError:pass
        else:raise AssertionError('Truth-read guard failed')


def one(task):
    i,voltage,initial=task;start=time.monotonic()
    result={}
    result['coupled_charge_only']=coupled_charge_inverse(RESPONSE,initial['times_ns'],initial['charge_pe'])
    result['isolated_energy_refit']=refine_energy(RESPONSE,voltage,initial,0.)
    result['recovery_energy_refit']=refine_energy(RESPONSE,voltage,initial,RESPONSE.strength)
    for row in result.values():assert len(row['times_ns'])==len(initial['times_ns'])
    return dict(index=i,methods=result,seconds=time.monotonic()-start)


def infer(workers=2,pilot=False):
    verify();verify('method-freeze.json');assert 1<=workers<=2
    if not pilot:assert load(RECORD/'pilot-validation.json')['passed']
    with np.load(OUTPUT/'inputs/observations.npz') as f:y=f['voltage_mV']
    starts={r['index']:r for r in lines(OUTPUT/'inputs/blind-starts.jsonl')}
    assert y.shape==(1600,1024)
    indices=load(EXPERIMENT/'protocol.json')['pilot_indices'] if pilot else list(range(1600))
    destination=OUTPUT/'inference';destination.mkdir(exist_ok=True)
    pending=[i for i in indices if not (destination/f'{i:04}.json').exists()]
    started=time.monotonic()
    with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize) as pool:
        futures=[pool.submit(one,(i,y[i],starts[i])) for i in pending]
        for count,future in enumerate(as_completed(futures),1):
            row=future.result();save(destination/f"{row['index']:04}.json",row)
            if count%40==0:print('RECOVERY_INFER',count,'/',len(pending),flush=True)
    if pilot:
        results=[load(destination/f'{i:04}.json') for i in indices]
        # Wall-clock upper estimate uses the slowest complete 3-variant frame.
        maximum=max(r['seconds'] for r in results);projected=maximum*1600/2
        save(RECORD/'pilot-validation.json',dict(passed=projected<=1800,frames=len(indices),max_seconds=maximum,
            projected_two_worker_seconds=projected,workers=workers,truth_guard=True))
        print('RECOVERY_PILOT',projected<=1800,'projected seconds',projected,flush=True)
    else:
        assert len(list(destination.glob('*.json')))==1600
        files=[destination/f'{i:04}.json' for i in range(1600)]
        save(RECORD/'predictions-freeze.json',dict(files={str(p.relative_to(OUTPUT)):digest(p) for p in files},
            truth_read_by_workers=False, methods=VARIANTS[1:], stage='Frozen before joining evaluation labels'))
        save(RECORD/'inference-validation.json',dict(passed=True,frames=1600,new_frames=len(pending),
            seconds=time.monotonic()-started,workers=workers,truth_guard=True,original_blind_count_retained=True))
        print('RECOVERY_INFERENCE_COMPLETE',len(files),flush=True)


def evaluate():
    verify('method-freeze.json')
    for file,sha in load(RECORD/'predictions-freeze.json')['files'].items():assert digest(OUTPUT/file)==sha
    cases=load(OUTPUT/'inputs/evaluation.cases.json')
    with np.load(OUTPUT/'inputs/evaluation.labels.npz') as f:q=f['charge_pe'];t=f['times_ns'];e=f['energy_MeV'];block=f['block']
    rows=[dict(r,method='static_nll') for r in lines(OUTPUT/'inputs/static.predictions.jsonl')]
    for i in range(1600):
        result=load(OUTPUT/'inference'/f'{i:04}.json')
        rows += [dict(r,index=i,method=m) for m,r in result['methods'].items()]
    scores=[]
    for row in rows:
        i=row['index'];c=cases[i]
        score=score_events(t[i],q[i],e[i],row['times_ns'],[np.nan if x is None else x for x in row['energy_MeV']],c['gap_ns'])
        scores.append(dict(score,index=i,method=row['method'],kind=c['kind'],condition=c['condition'],gap_ns=c['gap_ns'],
            block=int(block[i]),split=int(score['target']==1 and score['output_count']>1),
            merge=int(score['target']>=2 and score['output_count']==1)))
    write_lines(OUTPUT/'predictions.jsonl',rows);write_lines(OUTPUT/'scores.jsonl',scores)
    summaries=[]
    grouped=defaultdict(list)
    for row in scores:grouped[row['method'],row['kind'],row['condition'],row['gap_ns']].append(row)
    for (method,kind,condition,gap),ss in grouped.items():
        matched=[m for r in ss for m in r['matches'] if m['relative_energy_error'] is not None]
        summaries.append(dict(method=method,kind=kind,condition=condition,gap_ns=gap,frames=len(ss),
            **{name:sum(r[name] for r in ss) for name in ['target','missed','extra','output_count','strict_10','frame_success_10','split','merge','energy_out_of_domain']},
            matched_energy=error_summary([m['true_energy_MeV'] for m in matched],[m['estimated_energy_MeV'] for m in matched]) if matched else None))
    aggregate={}
    for method in VARIANTS:
        rr=[r for r in summaries if r['method']==method]
        aggregate[method]={group:float(np.mean([r['frame_success_10']/r['frames'] for r in rr if r['kind']=='pair' and ('weak' in r['condition'])==(group=='weak')])) for group in ['weak','strong']}
        aggregate[method]['boundary_frames']=sum(r.get('boundary',False) for r in rows if r['method']==method)
        aggregate[method]['not_converged']=sum(not r.get('converged',True) for r in rows if r['method']==method)
    save(OUTPUT/'report.json',dict(calibration=load(OUTPUT/'calibration-report.json'),summary=summaries,aggregate=aggregate,
        interpretation='Retrospective, fixed blind-selected count. No new independent transport or event-count selection.',
        inference=load(RECORD/'inference-validation.json')))
    print('RECOVERY_EVALUATED',aggregate,flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['prepare','calibrate','pilot','infer','evaluate'])
    parser.add_argument('--workers',type=int,default=2);args=parser.parse_args()
    if args.command in ['pilot','infer']:infer(args.workers,args.command=='pilot')
    else:globals()[args.command]()


if __name__=='__main__':main()
