"""Read-only verification of diagnostic algebra, provenance and report."""
import base64
import json
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import re
import sys
import nbformat
import numpy as np
import awkward as ak
import uproot
from uproot.source.file import MemmapSource

here=Path(__file__).resolve().parent;root=here.parents[2]
os.environ['GAGG_CURRENT_ENERGY_DIAGNOSIS_ROOT']=str(root)
sys.path.insert(0,str(here/'source'))
from pileup.analysis import alpha_current_energy_diagnosis_study as study
from pileup.analysis.alpha_current_energy_diagnosis import local_operator,repeat_covariance
read,sha=study.read,study.sha
checks=dict(frozen_inputs=study.verify())
print('Frozen design verified',file=sys.stderr,flush=True)
for name in ('original-input-freeze.json','provenance-freeze.json','diagnostic-output-freeze.json',
             'result-freeze.json','report-source-freeze.json','report-output-freeze.json'):
    for path,digest in read(here/name).items():assert sha(root/path)==digest,(name,path)
assert read(here/'unit-tests.json')['passed']==3
assert read(here/'reproduction-check.json')['passed']
assert read(here/'replay.json')['passed'] and read(here/'replay.json')['cases']==2
assert read(here/'run-completion.json')['cases']==40
plan=read(here/'case-plan.json');rows=read(here/'details.json');summary=read(here/'result-summary.json')
assert len(plan)==40 and len(rows)==160 and sum(r['reference_fit'] is not None for r in rows)==112
assert len({p for c in plan for p in c['parents']})==16 and len(summary['blocks'])==8
assert [c['block'] for c in plan[::5]]==summary['blocks']
for case in plan:
    cid=case['id']
    with np.load(study.ART/'diagnostics'/f'{cid:04d}.npz') as data:
        for row in (r for r in rows if r['case']==cid):
            prefix=row['method']+'_'+row['kind']
            inverse,geometry=local_operator(data[prefix+'_jacobian'])
            assert geometry==row['geometry']
            np.testing.assert_array_equal(inverse[1],data[prefix+'_weak_influence'])
            projected=data[prefix+'_projected_repeats']
            selected=projected[0] if row['kind']=='individual' else projected.mean(0)
            np.testing.assert_array_equal(selected,row['projection_terms'])
            np.testing.assert_allclose(selected.sum(0),row['projection_total'],rtol=1e-10,atol=1e-9)
            assert repeat_covariance(projected[:,:,1])==row['repeat_covariance']
            assert row['truth_assisted'] and row['closure_mV']<1e-8 and row['projection_closure']<1e-8
            if row['reference_fit']:
                f=row['reference_fit'];assert f['converged'] and f['cost_replay_relative_error']<1e-8
                assert f['cost_improvement']>=-1e-8
index=np.random.default_rng(742100000).integers(8,size=(10000,8))
for group in summary['groups']:
    selected=[r for r in rows if all(r[k]==group[k] for k in ('condition','gap_ns','kind','method'))]
    selected.sort(key=lambda r:summary['blocks'].index(r['block']))
    terms=np.array([r['projection_terms'] for r in selected])[:,:,1]
    np.testing.assert_array_equal(terms.mean(0),group['weak_terms_mean'])
    np.testing.assert_array_equal(np.percentile(terms[index].mean(1),[2.5,97.5],axis=0).T,group['weak_terms_ci95'])
for row in read(here/'input-mothers.json'):
    info=read(here/'original-runinfo'/f"{row['id']:04d}.json")
    assert info['macro_contents']==(root/row['diagnostic_macro']).read_text()
    with uproot.open(root/row['root'],handler=MemmapSource,num_workers=1,
                     decompression_executor=uproot.TrivialExecutor(),
                     interpretation_executor=uproot.TrivialExecutor()) as f:
        assert ak.to_list(f['RunInfo'].arrays(library='ak'))[0]==info
        edep=float(ak.to_numpy(f['PrimaryPhoton']['edep_total_MeV'].array(library='ak'))[0])
        assert edep==row['edep_MeV'] and info['random_seed']==row['seed']
        assert info['n_events']==1 and info['t_length_ps']==20000000
print('160 projections,bootstrap and16ROOT records verified',file=sys.stderr,flush=True)
counts={}
for name,before_count,after_count in [('codex_0007_SiPM_Energy_Calibration.ipynb',230,240),
            ('Research_Plan_20260921.ipynb',108,109),('Run_Catalog.ipynb',69,70)]:
    before=read(study.REC/(name+'.before'));after=read(root/'jupyter'/name)
    assert len(after['cells'])==after_count and before['cells']==after['cells'][:before_count]
    assert {k:v for k,v in before.items() if k!='cells'}=={k:v for k,v in after.items() if k!='cells'}
    nbformat.validate(nbformat.from_dict(after))
    for cell in after['cells'][before_count:]:
        if cell['cell_type']=='markdown':
            assert not re.search(r'\]\([^)]*\.(?:md|pdf|ipynb)(?:[#?)])',''.join(cell['source']),re.I)
        else:
            assert cell['execution_count'] is not None and any('image/png' in x.get('data',{}) for x in cell['outputs'])
            number=cell['id'].rsplit('-',1)[1]
            images=list((study.ART/'figures').glob(f'fig_0007_{number}_1_*.png'));assert len(images)==1
            output=next(x['data']['image/png'] for x in cell['outputs'] if 'image/png' in x.get('data',{}))
            assert base64.b64decode(''.join(output) if isinstance(output,list) else output)==images[0].read_bytes()
    counts[name]=after_count
checks.update(passed=True,cases=40,independent_mother_pairs=8,reused_ROOTs=16,diagnostic_rows=160,
    reference_cost_checks=112,diagnostic_replay_cases=2,reproduction_diagnostic_cases=1,new_fits=0,new_ROOTs=0,new_digitizations=0,
    notebooks=counts,diagnostic_only=True,default_changed=False)
print(json.dumps(checks,ensure_ascii=False,indent=2))
