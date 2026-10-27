"""Ancestry-checked joint residual diagnostics; development truth only."""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import replace
import shutil

import numpy as np
from scipy.signal import lfilter

from .alpha_pair_joint_data import PROJECT, EXPERIMENT, OUTPUT, RECORD, SOURCES
from .alpha_pair_resolution_data import load, save, digest


def poles(config):
    """Configured electronics poles, including the separate RC low pass."""
    r = config['readout']; coefficients = np.array([1.])
    for tau in (r['rise_ns'], r['decay_ns'], r['lowpass_ns']):
        if tau > 0:
            coefficients = np.convolve(coefficients, [1., -np.exp(-r['sample_period_ns']/tau)])
    return coefficients


def innovations(values, coefficients):
    # Fixed transform, marginalize the first three transient samples.
    return lfilter(coefficients, [1.], values, axis=-1)[..., len(coefficients)-1:]


def ancestry():
    from .sipm_response import AvalancheResult
    from .sipm_readout import ReadoutConfig, synthesize_readout
    target = OUTPUT/'calibration-ancestry.truth.npz'
    if target.exists():
        raise FileExistsError(target)
    bootstrap = EXPERIMENT/'preliminary-data-source'; bootstrap.mkdir()
    shutil.copyfile(PROJECT/'pileup/analysis/alpha_pair_joint_data.py', bootstrap/'alpha_pair_joint_data.py')
    records = {}
    for dataset, _ in SOURCES:
        for r in load(PROJECT/'jupyter/artifacts'/dataset/'native-records.json')['test_single']:
            records[dataset, r['seed']] = r
    cases = load(OUTPUT/'calibration.cases.json')
    with np.load(OUTPUT/'calibration.truth.npz') as f:
        arrays = {k:f[k] for k in f.files}
    checks = []; reference = None
    for i, c in enumerate(cases):
        if c['kind'] != 'single':
            continue
        r = records[c['dataset'], c['seed']]
        meta_path = PROJECT/r['metadata']; meta = load(meta_path); cfg = meta['config']
        if reference is None:
            reference = cfg
        assert cfg == reference, 'Historical digitizers must agree'
        readout = replace(ReadoutConfig(**cfg['readout']), baseline_mV=0., noise_rms_mV=0.)
        with np.load(PROJECT/r['truth']) as f:
            at=f['avalanche_time_ns']; aq=f['avalanche_charge_C']; owner=f['avalanche_owner']; signal=f['signal_mV']
        rec = np.empty(len(at), dtype=[('time_ns','f8'),('charge_C','f8')])
        rec['time_ns']=at; rec['charge_C']=aq
        q0=cfg['sensor']['gain_electrons']*1.602176634e-19
        wave=lambda block:synthesize_readout(AvalancheResult(block,meta['duration_ns'],q0,np.array([]),{}),readout,0).arrays['signal_mV']
        primary=wave(rec[owner>=0]); dark=wave(rec[owner<0])
        error=float(np.max(np.abs(primary+dark-signal)))
        assert error<1e-8, (r['seed'],error)
        assert np.isclose(aq[owner>=0].sum()/q0, arrays['charge_pe'][i,0], rtol=1e-10)
        start=c['crop_first_index']; crop=slice(start,start+1024)
        arrays['components_mV'][i,0]=primary[crop]
        checks.append(dict(index=i,seed=r['seed'],metadata=r['metadata'],metadata_sha256=digest(meta_path),
            closure_error_mV=error,crop_dark_peak_mV=float(np.max(dark[crop])),
            crop_dark_sum_mV=float(np.sum(dark[crop]))))
        if len(checks)%80==0:print('ANCESTRY_CHECKED',len(checks),flush=True)
    assert len(checks)==560
    np.savez_compressed(target,**arrays)
    assert reference==load(EXPERIMENT/'sensor-config.json'), 'Current copy disagrees with original metadata'
    save(EXPERIMENT/'historical-sensor-config.json',reference)
    save(RECORD/'ancestry-correction.json',dict(reason='Native signal_mV includes dark ancestry; preliminary pooled single components were therefore not primary-only.',
        preliminary_kept=True,preliminary_used_for_calibration=False,corrected_file=str(target.relative_to(PROJECT)),
        corrected_sha256=digest(target),checks=checks,readout_config_verified_against_all_native_metadata=True,
        ar_coefficients=poles(reference).tolist(),source_sha256=digest(__file__)))
    print('ANCESTRY_CORRECTED',len(checks),max(c['closure_error_mV'] for c in checks),flush=True)


def diagnose():
    from .alpha_pair_conditional import ConditionalResponse, calibrate
    with np.load(OUTPUT/'calibration.observations.npz') as f:y=f['voltage_mV'].astype(float)
    with np.load(OUTPUT/'calibration-ancestry.truth.npz') as f:
        q=f['charge_pe'];t=f['times_ns'];components=f['components_mV'].astype(float);fold=f['fold']
    cases=load(OUTPUT/'calibration.cases.json');train=fold==0
    response, details=calibrate(load(EXPERIMENT/'settings.json'),components[train],t[train],q[train],y[train])
    response.save(EXPERIMENT/'joint-response.npz')
    save(RECORD/'mean-calibration.json',dict(details,fold=0,frames=int(train.sum()),
        input_sha256=digest(OUTPUT/'calibration-ancestry.truth.npz'),source_sha256=digest(__file__)))
    coefficients=poles(load(EXPERIMENT/'historical-sensor-config.json'))
    groups=defaultdict(list)
    for i,c in enumerate(cases):groups[c['fold'],c['kind'],c['condition'],c['gap_ns']].append(i)
    background=y-response.baseline-components.sum(1)
    ibackground=innovations(background,coefficients)
    summaries=[];models={}
    for name,model in [('old',ConditionalResponse.load(EXPERIMENT/'response.npz')),('joint',response)]:
        mean=np.zeros_like(components)
        for i in range(len(y)):
            ids=np.flatnonzero(q[i]>=500.)
            mean[i,ids]=model.components(t[i,ids],q[i,ids])
        residual=components-mean;ir=innovations(residual,coefficients)
        imean=innovations(mean,coefficients)
        np.savez_compressed(OUTPUT/f'calibration-{name}-residuals.truth.npz',
            mean_mV=mean.astype(np.float32),innovations=ir.astype(np.float32),
            innovations_mean=imean.astype(np.float32))
        for key,ids in groups.items():
            if key[1]=='empty':continue
            rr=residual[ids];ii=ir[ids]
            # Raw-grid trace sums require no interpolation; fixed random onset phase is retained.
            # Center within the experimental group for empirical covariance, report mean error separately.
            r=rr-rr.mean(0);a=ii-ii.mean(0)
            var=np.mean(a*a,axis=(0,2));cross=float(np.mean(a[:,0]*a[:,1]))
            active=np.abs(imean[ids])>0.05
            lag=[]
            for k in range(9):
                u=a[...,k:];v=a[...,:a.shape[-1]-k]
                lag.append(float(np.sum(u*v)/max(np.sqrt(np.sum(u*u)*np.sum(v*v)),1e-30)))
            row=dict(model=name,fold=key[0],kind=key[1],condition=key[2],gap_ns=key[3],frames=len(ids),
                component_raw_mse=np.mean(rr*rr,axis=(0,2)).tolist(),component_raw_bias_norm=np.sqrt(np.mean(rr.mean(0)**2,axis=1)).tolist(),
                component_innovation_variance=var.tolist(),same_time_cross_covariance=cross,
                cross_correlation=cross/max(np.sqrt(np.prod(var)),1e-30),
                cross_variance_fraction=2*cross/max(float(var.sum()),1e-30),innovation_lag_correlation=lag,
                mean_innovation_negative_fraction=float(np.mean(imean[ids][active]<0)) if active.any() else 0.)
            summaries.append(row)
        models[name]=dict(component_raw_rmse=float(np.sqrt(np.mean(residual[fold==1]**2))),
            component_innovation_rmse=float(np.sqrt(np.mean(ir[fold==1]**2))))
    # Calibrate stationary dark/electronics covariance only from training frames.
    z=ibackground[train];lags=[float(np.mean(z[:,k:]*z[:,:z.shape[1]-k])) for k in range(9)]
    save(OUTPUT/'joint-diagnosis.json',dict(summary=summaries,models=models,
        noise_innovation_lags=lags,noise_innovation_mean=float(z.mean()),
        coefficients=coefficients.tolist(),fold0_training_frames=int(train.sum()),
        note='Empirical group centering on original sample grids. Reused optical mothers across gaps: descriptive moments, not independent uncertainty estimates.'))
    np.savez_compressed(OUTPUT/'background-innovations.truth.npz',innovations=ibackground.astype(np.float32))
    print('JOINT_DIAGNOSIS',models,'BACKGROUND_LAGS',lags,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['ancestry','diagnose'])
    args=parser.parse_args()
    ancestry() if args.command=='ancestry' else diagnose()
