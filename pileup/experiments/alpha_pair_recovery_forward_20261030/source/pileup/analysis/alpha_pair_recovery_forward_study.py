"""Native-voltage diagnostic, frozen before examining its fitted energies."""
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import argparse
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

from .alpha_pair_deconvolution import AgeRecovery, EnergyProfiles, FiniteConvolution, spe_kernel
from .alpha_pair_recovery_forward import CompiledAgeRecovery, NativeForward, METHODS, calibrate_voltage_residual
from .alpha_pair_recovery_energy import RecoveryEnergyResponse
from .alpha_pair_energy_calibrated import CalibratedEnergyResponse, learned_metric

NAME='alpha_pair_recovery_forward_20261030'
ROOT=Path(os.environ.get('GAGG_FORWARD_PROJECT',Path(__file__).resolve().parents[2]))
EXP=ROOT/'pileup/experiments'/NAME
ART=ROOT/'jupyter/artifacts'/NAME
REC=ROOT/'jupyter/run_records'/NAME
PREV=ROOT/'pileup/experiments/alpha_pair_deconvolution_20261030'
PREVART=ROOT/'jupyter/artifacts/alpha_pair_deconvolution_20261030'
OLD=ROOT/'pileup/experiments/alpha_pair_energy_calibrated_20261030'


def read(p): return json.loads(Path(p).read_text())
def write(p,v): Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def hashes(paths):return {str(p.relative_to(ROOT)):sha(p) for p in paths}


def build():
    (EXP/'bin').mkdir(parents=True,exist_ok=True)
    source=EXP/'source/pileup/analysis/sipm_age_mean.cpp'; target=EXP/'bin/libsipm_age_mean.so'
    if not target.exists() or not (EXP/'build.json').exists() or read(EXP/'build.json')['source_sha256'] != sha(source):
        command=['g++','-O3','-std=c++17','-fPIC','-shared',str(source),'-o',str(target)]
        subprocess.run(command,check=True)
        write(EXP/'build.json',dict(command=command,source_sha256=sha(source),binary_sha256=sha(target),
            compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0],fast_math=False))


def prepare():
    if (EXP/'protocol.json').exists():raise FileExistsError('Frozen experiment exists')
    for p in (EXP,ART/'inputs',ART/'cases',REC):p.mkdir(parents=True,exist_ok=True)
    for name in ('sensor-config.json','energy-model.json','recovery128-model.npz','mother-split.json'):
        shutil.copyfile(PREV/name,EXP/name)
    for name in ('template.npz','response-parameters.json','mean-correction.npz','noise-calibration.json'):
        shutil.copyfile(OLD/name,EXP/name)
    shutil.copyfile(PREVART/'inputs/observations.npz',ART/'inputs/observations.npz')
    shutil.copyfile(PREVART/'labels.json',ART/'labels.json')
    protocol=dict(created=datetime.now(timezone.utc).isoformat(),commit='09732062ddc851946e7215a481e4e4a12484fbe4',
        prior_local_stage='alpha_pair_deconvolution_20261030, not committed; actual snapshots/hashes preserved',
        background_mV=read(PREV/'protocol.json')['background_mV'],
        scope='Pure alpha, TRUE K/t supplied, E_dep fitted; same64 historical diagnostic frames,32 optical mothers; not untouched validation',
        motivation='Prior transformed-domain recovery pilot4/32 pairs versus voltageGLS9/32; test refolding and independently calibrated native mean',
        methods=list(METHODS),primary='native_gls_calibrated',references=['recovery128','new_mean_new_noise'],
        physical_model='Identical age distribution and SPE kernel to preceding pilot; no revised physical constants. C++ acceleration numerically compared to Python.',
        voltage_model='SPE*Age(frozen recovered-light profiles plus empty transformed mean), minus the empty response, plus fixed voltage background. Centering avoids counting algorithmic empty bias twice.',
        mean_calibration='280fold0 isolated noisy voltage residuals at true calibration E/t;4energy knots, normalized by inheritedq0(E),4ns smoothing,-64..1720ns,last80ns taper; additive correction then nonnegative total signal. No pair/evaluation labels.',
        metrics='OLS native voltage; GLS frozen previous calibrated FIR covariance at SAME voltage-derived recovery128 energy seed for both GLS arms. No covariance refit, no dependence on candidate E.',
        seed='Recompute RL128/recovery and transformed-domain energy fit from voltage plus knownK/t; do not read previous fitted answers',
        search='Shared initial+five charge shares for pairs, one initial for singles;100least-squares evaluations per start; same inherited E bounds',
        endpoint='All injected E_dep within5/10/20%; preserve every boundary/nonconverged case; biases,sigma68,RMSE,gained/lost. No unknownK/t or PSD metrics.',
        plan='64frames only; no new optical mothers or parameter tuning from fitted outputs. No full480 expansion in this diagnostic.',
        cost_gate_seconds=1800,pilot_ids=[40,160],pilot_use='Cost/numerical guard only; fixed parameters unchanged',
        display_ids=[40,120],workers=6,BLAS_threads=1,no_new_ROOT=True,default_changed=False)
    write(EXP/'protocol.json',protocol)
    for name in ('codex_0007_SiPM_Energy_Calibration.ipynb','Research_Plan_20260921.ipynb'):
        shutil.copyfile(ROOT/'jupyter'/name,REC/(name+'.before'))
    for p in (ROOT/'README.md',ROOT/'docs/literature-use.md'):shutil.copyfile(p,REC/(p.name+'.before'))
    write(EXP/'input-freeze.json',hashes([ART/'inputs/observations.npz',ART/'labels.json',PREVART/'calibration/labels.json',
         PREVART/'calibration/observations.npz',PREV/'model-freeze.json',OLD/'model-freeze.json']))
    build();print('PREPARED',flush=True)


def guard_reads():
    def audit(event,args):
        if event=='open' and isinstance(args[0],(str,bytes,os.PathLike)):
            path=Path(os.fsdecode(args[0]));mode=args[1]
            if (mode is None or 'r' in mode) and path.suffix not in ('.py','.pyc','.so'):
                raise PermissionError('Inference data read blocked: '+path.name)
    sys.addaudithook(audit)
    try:(ART/'labels.json').read_text()
    except PermissionError:return True
    raise AssertionError('Read guard failed')


def initialize(inference=False):
    global MODEL,OP,PROFILE,AGE,BACKGROUND,LEGACY,NOISE,GUARD
    cfg=read(EXP/'sensor-config.json');BACKGROUND=read(EXP/'protocol.json')['background_mV']
    with np.load(EXP/'recovery128-model.npz') as f:PROFILE=EnergyProfiles(**{k:f[k] for k in f.files})
    OP=FiniteConvolution(spe_kernel(cfg),1024);AGE=CompiledAgeRecovery(cfg['sensor'],EXP/'bin/libsipm_age_mean.so')
    correction=None
    if (EXP/'voltage-correction.npz').exists():
        with np.load(EXP/'voltage-correction.npz') as f:correction={k:f[k] for k in f.files}
    MODEL=NativeForward(PROFILE,OP,AGE,BACKGROUND,correction)
    base=RecoveryEnergyResponse.load(EXP);LEGACY=CalibratedEnergyResponse.load(base,EXP/'mean-correction.npz')
    NOISE=read(EXP/'noise-calibration.json')
    GUARD=guard_reads() if inference else False


def calibrate():
    if (EXP/'model-freeze.json').exists():raise FileExistsError('Calibration already frozen')
    initialize();start=time.monotonic()
    records=read(PREVART/'calibration/labels.json')
    with np.load(PREVART/'calibration/observations.npz') as f:voltage=f['voltage_mV']
    correction=calibrate_voltage_residual(MODEL,records,voltage)
    np.savez_compressed(EXP/'voltage-correction.npz',**correction)
    write(REC/'calibration-run.json',dict(frames=len(records),seconds=time.monotonic()-start,energy_knots_MeV=correction['energy_knots'].tolist()))
    write(EXP/'model-freeze.json',hashes([EXP/name for name in ('protocol.json','sensor-config.json','energy-model.json',
        'recovery128-model.npz','template.npz','response-parameters.json','mean-correction.npz','noise-calibration.json','voltage-correction.npz')]))
    print('CALIBRATED',flush=True)


def controls():
    initialize();start=time.monotonic();cfg=read(EXP/'sensor-config.json')['sensor'];rng=np.random.default_rng(473)
    comparisons=[]
    for dt in (1.,2.):
        slow=AgeRecovery(cfg,dt);fast=CompiledAgeRecovery(cfg,EXP/'bin/libsipm_age_mean.so',dt)
        for x in [np.zeros(128),rng.gamma(2.,50.,1024),np.r_[np.zeros(40),np.ones(30)*500.,np.zeros(954)]]:
            expected=slow.forward(x);actual=fast.forward(x)
            np.testing.assert_allclose(actual,expected,rtol=5e-13,atol=5e-11)
            comparisons.append(float(abs(actual-expected).max()))
    np.testing.assert_allclose(MODEL.signal([],[]),np.zeros(1024),atol=1e-12)
    np.testing.assert_allclose(MODEL.signal([180.,220.],[7.,.3],corrected=True),MODEL.signal([220.,180.],[.3,7.],corrected=True),atol=1e-9)
    inverse=[]
    for t,e,start_e in [([180.],[.25],[.4]),([180.],[7.5],[6.]),([180.,220.],[7.5,.25],[6.,.4])]:
        for corrected in (False,True):
            y=BACKGROUND+MODEL.signal(t,e,corrected=corrected)
            metric=learned_metric(LEGACY,t,start_e,NOISE)
            method='native_gls_calibrated' if corrected else 'native_gls_centered'
            fit=MODEL.fit(y,t,start_e,method=method,metric=metric,transform=LEGACY.transform)
            error=float(np.max(np.abs(np.asarray(fit['energy_MeV'])/e-1)))
            assert fit['converged'] and error<1e-5,(method,error)
            inverse.append(dict(times_ns=t,energies_MeV=e,method=method,max_relative_error=error))
    write(REC/'numerical-controls.json',dict(passed=True,cpp_python_max_abs_error=max(comparisons),comparisons=len(comparisons),
        centered_empty_zero=True,component_permutation=True,noiseless_inverse=inverse,seconds=time.monotonic()-start))
    print('CONTROLS_PASSED',flush=True)


def job(item):
    identifier,voltage,times=item;start=time.monotonic()
    released=OP.unfold(voltage,BACKGROUND,(128,))[128];light,diag=AGE.invert(released)
    seed=PROFILE.fit(light,times)
    metric=learned_metric(LEGACY,times,seed['energy_MeV'],NOISE)
    fits={m:MODEL.fit(voltage,times,seed['energy_MeV'],method=m,metric=metric,transform=LEGACY.transform) for m in METHODS}
    return dict(id=int(identifier),fits=fits,seed=seed,recovery=diag,seconds=time.monotonic()-start,
                known_K_t=True,true_E_Q_input=False,worker_truth_read_guard=GUARD)


def run(workers):
    build()
    assert read(REC/'numerical-controls.json')['passed']
    for name,h in read(EXP/'model-freeze.json').items():assert sha(ROOT/name)==h,name
    with np.load(ART/'inputs/observations.npz') as f:
        items=[(i,v,t[:k]) for i,v,t,k in zip(f['ids'],f['voltage_mV'],f['known_times_ns'],f['known_counts'])]
    if list((ART/'cases').glob('*.json')):raise FileExistsError('Preserve existing results')
    pilot_ids=read(EXP/'protocol.json')['pilot_ids'];pilot=[x for x in items if int(x[0]) in pilot_ids]
    start=time.monotonic()
    with ProcessPoolExecutor(max_workers=workers,initializer=initialize,initargs=(True,)) as pool:
        pilot_results=list(pool.map(job,pilot))
        estimate=max(r['seconds'] for r in pilot_results)*len(items)/workers*1.5
        write(REC/'pilot-cost.json',dict(estimated_seconds=estimate,gate_seconds=1800,frames=len(pilot_results),
                                       seconds=[r['seconds'] for r in pilot_results],passed=estimate<1800))
        for r in pilot_results:write(ART/'cases'/f"{r['id']:04d}.json",r)
        print(f'PILOT_COST estimate={estimate:.1f}s',flush=True)
        if estimate>=1800:raise RuntimeError('Cost gate failed; pilot preserved, no main run launched')
        rest=[x for x in items if int(x[0]) not in pilot_ids]
        for count,result in enumerate(pool.map(job,rest)):
            write(ART/'cases'/f"{result['id']:04d}.json",result)
            if count%8==0:print(f'evaluation {count+3}/{len(items)}, wall={time.monotonic()-start:.1f}s',flush=True)
    write(REC/'evaluation-run.json',dict(frames=len(items),fits=len(items)*len(METHODS),seconds=time.monotonic()-start,workers=workers))
    print('COMPLETE',flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['prepare','calibrate','controls','run']);parser.add_argument('--workers',type=int,default=6)
    args=parser.parse_args()
    if args.mode=='prepare':prepare()
    elif args.mode=='calibrate':calibrate()
    elif args.mode=='controls':controls()
    else:run(args.workers)


if __name__=='__main__':main()
