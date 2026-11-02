"""Verify frozen conditional analysis, grouped identities and append-only notes."""
from pathlib import Path
import json
import os
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent
os.environ['GAGG_SHAPE_STABILITY_ROOT']=str(HERE.parents[2])
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
sys.path.insert(0,str(HERE/'source'))
import numpy as np
import nbformat
from pileup.analysis.alpha_mean_shape_stability_study import ROOT,EXP,ART,REC,read,save,sha,verify,load_arrays
from pileup.analysis.alpha_mean_shape_stability import bootstrap_mothers


def main():
    checks=verify()
    for name in ('results/primary/output-freeze.json','report-source-freeze.json','report-output-freeze.json'):
        manifest=read(EXP/name)
        for p,d in manifest.items():assert sha(ROOT/p)==d,p
        checks[name]=len(manifest)
    tests=subprocess.run([sys.executable,'-B','-m','unittest','pileup.analysis.test_alpha_mean_shape_stability'],
                         cwd=HERE/'source',capture_output=True,text=True)
    assert tests.returncode==0,tests.stdout+tests.stderr
    (REC/'unit-tests.log').write_text(tests.stdout+tests.stderr)
    summary=read(EXP/'results/primary/summary.json');parts=read(EXP/'results/primary/partitions.json')
    boots=read(EXP/'results/primary/bootstrap.json');data=load_arrays(ART/'inputs.npz')
    mothers=read(EXP/'mothers.json');folds=np.array(read(EXP/'folds.json')['assignments'])
    assert len(mothers)==128 and len({r['root'] for r in mothers})==128
    assert all(r['role']=='train' for r in mothers)
    assert data['signal_repeats_mV'].shape==(128,8,1024)
    np.testing.assert_array_equal(data['signal_repeats_mV'].mean(1),data['means_mV'])
    assert folds.shape==(10,128) and len(parts)==10 and len(boots)==1000
    for partition in folds:
        for group in range(4):np.testing.assert_array_equal(np.bincount(partition[data['groups']==group]),[8]*4)
    for i,row in enumerate(boots):
        sample=bootstrap_mothers(data['groups'],folds[0],737100000+i)
        for fold in range(4):
            train=sample[folds[0,sample]!=fold];test=sample[folds[0,sample]==fold]
            assert not set(train)&set(test)
            model=row['models'][fold]
            assert model['unique_training_mothers']==len(set(train)) and model['unique_test_mothers']==len(set(test))
    for i,row in enumerate(parts):
        a=load_arrays(ART/'results/primary'/f'partition{i:02d}.npz')
        high=data['groups']==3
        for method,key in [('old','old_projection_100logE'),('corrected','corrected_projection_100logE')]:
            residual=a[method+'_residual_mV']
            np.testing.assert_allclose(100*residual[high]@data['operator'].T,a[key],rtol=1e-13,atol=1e-11)
            for j,r in enumerate(row['projection']):
                assert abs(r[method+'_mean']-a[key][:,j].mean())<1e-10
                assert abs(r[method+'_mse']-np.mean(a[key][:,j]**2))<1e-8
        assert np.max(abs(a['profiles'].sum(2)))<1e-12
    for r in summary['variance']:
        assert abs(r['observed_mean_variance']-(r['between_variance_unclipped']+r['within_variance']/8))<1e-9
    assert not summary['development_stability_passed']
    assert summary['failures']==['gap8:insufficient development stability','gap40:insufficient development stability']
    assert summary['minimum_scale']==1 and summary['new_digitizations']==0 and summary['energy_fits']==0
    assert read(EXP/'numerical.json')['full_recipe_exactly_reproduces_L93']
    replay=read(EXP/'replay.json');assert replay['passed'] and replay['exact_arrays']==12
    exact=0
    for i in (0,9):
        a=load_arrays(ART/'results/primary'/f'partition{i:02d}.npz');b=load_arrays(ART/'replay'/f'partition{i:02d}.npz')
        for k in a:np.testing.assert_array_equal(a[k],b[k]);exact+=1
    for r in read(REC/'pilot.json')['rows']:assert r['result']==boots[r['id']]
    for row in read(REC/'replay.json')['rows'][2:]:assert row['result']==boots[row['id']]
    algebra=read(EXP/'pooled-mean-control.json')
    assert abs(algebra['corrected_pooled_mean'])<1e-14 and abs(algebra['corrected_to_old_mse']-16/9)<1e-14
    notebook_checks={}
    for filename,count,added in [('codex_0007_SiPM_Energy_Calibration.ipynb',185,8),
                                 ('Research_Plan_20260921.ipynb',103,1),('Run_Catalog.ipynb',64,1)]:
        nb=read(ROOT/'jupyter'/filename);before=read(REC/(filename+'.before'))
        assert len(before['cells'])==count and len(nb['cells'])==count+added
        assert nb['cells'][:count]==before['cells']
        assert {k:v for k,v in nb.items() if k!='cells'}=={k:v for k,v in before.items() if k!='cells'}
        nbformat.validate(nbformat.from_dict(nb));images=0
        for cell in nb['cells'][count:]:
            source=cell['source'] if isinstance(cell['source'],str) else ''.join(cell['source'])
            if cell['cell_type']=='code':
                assert cell['execution_count'] is not None
                assert not any(o['output_type']=='error' for o in cell['outputs'])
                images+=sum('image/png' in o.get('data',{}) for o in cell['outputs'])
            else:assert not re.search(r'\]\([^\n)]*\.(?:md|pdf|ipynb)(?:[#:][^)]*)?\)',source,re.I)
        assert images==(3 if added==8 else 0)
        notebook_checks[filename]=dict(preserved_cells=count,added_cells=added,new_images=images)
    for name in ('README.md','docs/literature-use.md'):
        p=ROOT/name;assert p.read_text().startswith((REC/(p.name+'.before')).read_text())
    for p in (EXP/'tables').glob('*.csv'):assert b'\r' not in p.read_bytes()
    assert read(REC/'primary.json')['worker_confirmation_read_guard']
    result=dict(passed=True,hash_checks=checks,unit_tests=4,partitions=10,refit_bootstraps=1000,
        all_bootstrap_identities_isolated=True,all_partition_metrics_recomputed=True,
        exact_replay_arrays=exact,selected_summary_replays=5,pilot_equal=True,full_L93_recipe_equal=True,
        pooled_mean_counterexample_verified=True,formal_calibrations=4040,total_study_calibrations=4093,
        new_digitizations=0,new_ROOT=0,energy_fits=0,notebooks=notebook_checks)
    save(EXP/'validation.json',result);print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':main()
