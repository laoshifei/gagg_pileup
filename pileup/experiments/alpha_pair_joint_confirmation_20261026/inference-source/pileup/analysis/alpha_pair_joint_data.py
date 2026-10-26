"""Mother-disjoint development data for joint-pulse likelihood calibration."""
from __future__ import annotations

import argparse
from datetime import datetime
from importlib.metadata import version
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

from .alpha_pair_resolution_data import load,save,digest

PROJECT=Path(os.environ.get('GAGG_JOINT_PROJECT',Path(__file__).resolve().parents[2])).resolve()
NAME='alpha_pair_joint_likelihood_20261026'
EXPERIMENT=PROJECT/'pileup/experiments'/NAME
OUTPUT=PROJECT/'jupyter/artifacts'/NAME
RECORD=PROJECT/'jupyter/run_records'/NAME
PREVIOUS=PROJECT/'pileup/experiments/alpha_pair_null_gate_20261023'
SOURCES=[('alpha_pair_conditional_20261022','test'),('alpha_pair_residual_20261023','test'),
         ('alpha_pair_null_gate_20261023/development','labels')]


def prepare():
    for p in (EXPERIMENT,OUTPUT,RECORD):p.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb',RECORD/'notebook.before.ipynb')
    for name in ('response.npz','settings.json','calibration.json','model-null-calibration.json'):
        shutil.copyfile(PREVIOUS/name,EXPERIMENT/name)
    config=load(PROJECT/'g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json')
    save(EXPERIMENT/'sensor-config.json',config)
    plan=dict(created=datetime.now().astimezone().isoformat(),scope='Pure alpha; fixed95/351ns; local unknown K<=2',
        data_role='Previously exposed development only; no fresh independent test yet',
        split='Native mother replicate//2 modulo2; retain only frames whose mothers belong to one fold',
        training_fold=0,gate_fold=1,gate='Only held-out native singles; no training mothers',
        development_validation='Previous null-gate test: all ten combinations at8/40/200ns and all controls',
        time_definition='First active-SiPM photon; Q is realized avalanche-equivalent pe, not energy',
        diagnoses=['Per-component conditional-mean residuals','Same-time and lagged cross covariance',
                   'Readout-pole innovations from configured2/20/2ns poles','Finite-record charge conservation'],
        input='Full analog voltage; no ADC clipping or quantization selection',
        native_single_components='Recover genuine signal_mV from native truth export; old packed zeros are placeholders',
        selection=dict(control_extra_fraction_max=.05,single_success_min=.95,median_abs_charge_bias_max=.05,
            weak_average_must_exceed_both_legacy_blind_methods=True,strong_max_drop_percentage_points=2.),
        literature='Fixsen2004 conditional mean/covariance idea only; TES simulation, AI note pending review; implementation is project-specific')
    save(EXPERIMENT/'data-plan.json',plan)
    waves=[];charges=[];times=[];components=[];all_cases=[];source_files={};native_files={};fold_mothers=[set(),set()]
    for dataset,stem in SOURCES:
        folder=PROJECT/'jupyter/artifacts'/dataset
        native=load(folder/'native-records.json')['test_single']
        parents=load(folder/'physical-parents.json')['records'] if stem=='labels' else native
        cases=load(folder/(stem+'.cases.json'))
        wave=folder/('observations.npz' if stem=='labels' else stem+'.waveforms.npz')
        truth=folder/(stem+'.truth.npz')
        with np.load(wave) as f:y=f['voltage_mV']
        with np.load(truth) as f:q=f['charge_pe'];t=f['times_ns'];comp=f['components_mV']
        assert len(y)==len(cases)==len(q)
        native_by_seed={r['seed']:r for r in native}
        for i,c in enumerate(cases):
            if c['kind']=='single':mother_records=[native_by_seed[c['seed']]]
            else:mother_records=[parents[j] for j in c['parents']]
            folds={r['replicate']//2%2 for r in mother_records}
            if len(folds)>1:continue
            fold=next(iter(folds)) if folds else c['replicate']//2%2
            mother_keys=[r['root'] for r in mother_records]
            fold_mothers[fold].update(mother_keys)
            curves=comp[i].copy()
            if c['kind']=='single':
                r=mother_records[0];path=PROJECT/r['truth'];start=c['crop_first_index']
                with np.load(path) as f:curves[0]=f['signal_mV'][start:start+1024]
                native_files[r['truth']]=dict(sha256=digest(path),root=r['root'],root_sha256=r['root_original_sha256'],seed=r['seed'])
            waves.append(y[i]);charges.append(q[i]);times.append(t[i]);components.append(curves)
            all_cases.append(dict(c,dataset=dataset,source_index=i,fold=fold,mother_roots=mother_keys))
        source_files.update({str(p.relative_to(PROJECT)):digest(p) for p in
            (wave,truth,folder/(stem+'.cases.json'),folder/'native-records.json')})
        print('JOINT_DATA_SOURCE',dataset,len(all_cases),flush=True)
    assert not fold_mothers[0]&fold_mothers[1]
    assert len(fold_mothers[0])==len(fold_mothers[1])==280
    y=np.asarray(waves,np.float32);q=np.asarray(charges);t=np.asarray(times);comp=np.asarray(components,np.float32)
    np.savez_compressed(OUTPUT/'calibration.observations.npz',voltage_mV=y)
    np.savez_compressed(OUTPUT/'calibration.truth.npz',charge_pe=q,times_ns=t,components_mV=comp,
                        fold=np.asarray([c['fold'] for c in all_cases]))
    save(OUTPUT/'calibration.cases.json',all_cases)
    validation=PROJECT/'jupyter/artifacts/alpha_pair_null_gate_20261023/test'
    valid_cases=load(validation/'labels.cases.json')
    indices=[i for i,c in enumerate(valid_cases) if c['kind']!='pair' or c['gap_ns'] in (8.,40.,200.)]
    valid_mothers={r['root'] for r in load(validation/'native-records.json')['test_single']}
    assert not (fold_mothers[0]|fold_mothers[1])&valid_mothers and len(indices)==1600
    save(EXPERIMENT/'validation-indices.json',indices)
    save(RECORD/'data-provenance.json',dict(input_files=source_files,native_single_files=native_files,
        fold_mothers=[sorted(s) for s in fold_mothers],validation_mothers=sorted(valid_mothers),
        mothers_disjoint=True,frames=len(all_cases),
        source_code_sha256=digest(Path(__file__)),
        array_hashes={p.name:digest(p) for p in OUTPUT.glob('calibration.*')}))
    save(RECORD/'environment.json',dict(python=sys.version,executable=sys.executable,
        dependencies={name:version(name) for name in ('numpy','scipy','matplotlib','nbformat','pandas')},
        workers=4,blas_threads=1,git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        preexisting_changes='Previous weak-diagnostic and empirical-covariance work is uncommitted and preserved; README, literature-use and eight modules plus two experiment directories'))
    print('JOINT_DATA_PREPARED',len(y),'frames',len(fold_mothers[0]),len(fold_mothers[1]),flush=True)


if __name__=='__main__':
    argparse.ArgumentParser(description=__doc__).parse_args();prepare()
