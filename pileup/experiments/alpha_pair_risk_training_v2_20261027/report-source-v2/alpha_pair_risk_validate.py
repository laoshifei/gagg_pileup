"""Final source, split, scoring and executed-notebook checks."""
from __future__ import annotations

import ast
import base64
from datetime import datetime
import hashlib
import json
import re
import subprocess
import sys

import nbformat
import numpy as np

from .alpha_pair_risk_study import PROJECT,EXPERIMENT,OUTPUT,RECORD,NAME,plan_case
from .alpha_pair_information_study import load,save,digest,rows


def validate():
    required=('numerical-tests','pilot-validation','preflight-head-amendment','train-generation','validation-generation',
              'evaluation-generation','native-features','continuous-baseline','native-inference','continuous-inference',
              'science-validation','figure-delivery','figure-visual-review','notebook-delivery','document-delivery','dependency-repair')
    for name in required:assert load(RECORD/(name+'.json'))['passed'],name
    manifests={}
    for name in ('freeze','evaluation-freeze','dependency-freeze','training-freeze','report-freeze','report-v2-freeze'):
        manifest=load(EXPERIMENT/(name+'.json'));manifests[name]=manifest
        for path,expected in manifest['files'].items():assert digest(EXPERIMENT/path)==expected,path
    assert manifests['evaluation-freeze']['training_source_freeze_sha256']==digest(EXPERIMENT/'freeze.json')
    assert manifests['dependency-freeze']['evaluation_freeze_sha256']==digest(EXPERIMENT/'evaluation-freeze.json')
    frozen=manifests['training-freeze'];launch=load(RECORD/'evaluation-launch.json')
    assert launch['training_freeze_sha256']==digest(EXPERIMENT/'training-freeze.json')
    assert datetime.fromisoformat(frozen['frozen_at'])<datetime.fromisoformat(launch['launched_at'])
    assert not frozen['continuous_evaluation_generated'] and launch['models_chosen_before_evaluation']
    old=PROJECT/'pileup/experiments/alpha_pair_risk_training_20261027'
    for name in ('freeze','evaluation-freeze'):
        for path,expected in load(old/(name+'.json'))['files'].items():assert digest(old/path)==expected,path
    source_hashes={}
    groups={
        'inference-source/pileup/analysis':('alpha_pair_risk_learning.py','alpha_pair_risk_study.py','test_alpha_pair_risk_learning.py'),
        'evaluation-source/pileup/analysis':('alpha_pair_risk_learning.py','alpha_pair_risk_study.py','alpha_pair_risk_inference.py'),
        'report-source-v2':('alpha_pair_risk_report.py','alpha_pair_risk_notebook.py','alpha_pair_risk_documents.py','alpha_pair_risk_validate.py')}
    for folder,names in groups.items():
        for name in names:
            path=PROJECT/'pileup/analysis'/name;ast.parse(path.read_text())
            assert digest(path)==digest(EXPERIMENT/folder/name)
            source_hashes[str(path.relative_to(PROJECT))]=digest(path)
    provenance=load(RECORD/'input-provenance.json')
    for name,expected in provenance['input_files'].items():assert digest(PROJECT/name)==expected,name
    for name,expected in provenance['native_files'].items():assert digest(OUTPUT/name)==expected,name
    for role in ('train-pilot','train','validation','evaluation'):
        generation=load(RECORD/(role+'-generation.json'))
        for name,expected in generation['files'].items():assert digest(OUTPUT/name)==expected,name
    seeds={};counts={}
    for role,size in (('train',9600),('validation',960),('evaluation',1000)):
        cases=load(OUTPUT/role/'labels.cases.json');assert len(cases)==size
        seeds[role]={c['seed'] for c in cases};assert len(seeds[role])==size
        with np.load(OUTPUT/role/'labels.truth.npz') as data:
            q=data['charge_pe'];t=data['times_ns'];k=data['count']
        with np.load(OUTPUT/role/'features.npz') as data:
            assert data['features'].shape==(size,1408) and np.isfinite(data['features']).all()
        for i,case in enumerate(cases):
            assert case==plan_case(role,i),(role,i)
            assert k[i]==case['count']
            np.testing.assert_array_equal(q[i,:k[i]],case['charge_pe'])
            np.testing.assert_array_equal(t[i,:k[i]],case['times_ns'])
        counts[role]={str(n):int((k==n).sum()) for n in (0,1,2)}
    assert not(seeds['train']&seeds['validation'] or seeds['train']&seeds['evaluation'] or seeds['validation']&seeds['evaluation'])
    # Audit the actual frozen worker initializer in a fresh process; never disable the hook.
    program='''import sys,json
sys.path.insert(0,sys.argv[1])
from pileup.analysis.alpha_pair_risk_inference import initialize
initialize()
blocked=[]
for p in ('forbidden.truth.npz','forbidden.cases.json','forbidden.root','/tmp/truth/forbidden.npy','native-records.json'):
    try:open(p,'rb')
    except PermissionError:blocked.append(p)
    else:raise AssertionError(p)
print(json.dumps(blocked))
'''
    blocked=json.loads(subprocess.check_output([sys.executable,'-c',program,str(EXPERIMENT/'evaluation-source')],cwd=PROJECT,text=True))
    assert len(blocked)==5
    for domain in ('native','continuous'):
        record=load(RECORD/(domain+'-inference.json'))
        assert record['truth_inputs']==[] and record['label_access_blocked']
        assert record['prediction_sha256']==digest(OUTPUT/(domain+'-predictions.jsonl'))
        assert record['network_outputs_sha256']==digest(OUTPUT/(domain+'-network-outputs.npz'))
        assert len(rows(OUTPUT/(domain+'-predictions.jsonl')))==(1040 if domain=='native' else 1000)
    science=load(RECORD/'science-validation.json')
    for name,expected in science['output_hashes'].items():assert digest(OUTPUT/name)==expected,name
    assert len(rows(OUTPUT/'scores.jsonl'))==10560 and len(rows(OUTPUT/'forced-K2.jsonl'))==5760
    document=load(RECORD/'document-delivery.json')
    for name,expected in document['files'].items():assert digest(PROJECT/name)==expected,name
    previous=load(PROJECT/'jupyter/run_records/alpha_pair_charge_allocation_20261027/final-validation.json')
    for name,expected in previous['source_sha256'].items():assert digest(PROJECT/name)==expected,name
    for name,backup in (('README.md','README.before.md'),('docs/literature-use.md','literature.before.md')):
        assert digest(RECORD/backup)==previous['documentation_sha256'][name]
        assert (PROJECT/name).read_text().startswith((RECORD/backup).read_text().rstrip()+'\n\n')
    older=load(PROJECT/'jupyter/run_records/alpha_pair_information_20261027/final-validation.json')
    for name,expected in older['source_sha256'].items():assert digest(PROJECT/name)==expected,name
    path=PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    notebook=nbformat.read(path,as_version=4);nbformat.validate(notebook)
    before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4);nbformat.validate(before)
    assert notebook.cells[:431]==before.cells and len(before.cells)==431
    new=notebook.cells[431:];delivery=load(RECORD/'notebook-delivery.json')
    assert len(new)==delivery['new_cells'] and digest(path)==delivery['notebook_sha256']
    assert all(c.id.startswith('risk-20261027-') for c in new)
    assert len({c.id for c in notebook.cells})==len(notebook.cells)
    anchors=re.findall(r'id="codex-figure-([^\"]+)"','\n'.join(c.source for c in notebook.cells if c.cell_type=='markdown'))
    assert len(anchors)==len(set(anchors))
    codes=[c for c in new if c.cell_type=='code'];assert len(codes)==delivery['executed_cells']
    assert all(c.execution_count and c.outputs and not any(o.output_type=='error' for o in c.outputs) for c in codes)
    figures=load(RECORD/'figure-delivery.json')['figures'];review=load(RECORD/'figure-visual-review.json')
    assert figures==review['figures']
    for name,expected in figures.items():assert digest(OUTPUT/'figures'/(name+'.png'))==expected
    embedded=[hashlib.sha256(base64.b64decode(o.data['image/png'])).hexdigest() for c in codes for o in c.outputs if 'image/png' in o.get('data',{})]
    assert len(embedded)==6 and set(embedded)==set(figures.values())
    for c in new:
        if c.cell_type=='markdown':assert not re.search(r'\]\([^)]*\.(?:pdf|ipynb|md)(?:[#?)][^)]*)?\)',c.source,re.I)
    ignored=[str(path.relative_to(PROJECT)),f'jupyter/artifacts/{NAME}/report.json',f'pileup/experiments/{NAME}/models/relative/model.npz']
    for name in ignored:assert subprocess.run(['git','check-ignore','-q',name],cwd=PROJECT).returncode==0,name
    subprocess.run(['git','-c','core.fileMode=false','diff','--check','--','README.md','docs/literature-use.md','pileup/analysis'],cwd=PROJECT,check=True)
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip()==load(RECORD/'environment.json')['git_head']
    save(RECORD/'final-validation.json',dict(passed=True,numerical_tests=6,training_frames=9600,validation_frames=960,
        new_Gaussian_evaluation=1000,exposed_native_evaluation=760,calibration_frames=280,primary_scores=10560,forced_K2_scores=5760,
        frozen_files={k:len(v['files']) for k,v in manifests.items()},source_sha256=source_hashes,dataset_counts=counts,
        parameter_noise_seeds_disjoint=True,weights_frozen_before_evaluation=True,truth_open_probes_blocked=blocked,
        original_431_cells_preserved=True,new_cells=len(new),executed_new_cells=len(codes),figures=figures,
        notebook_sha256=digest(path),prior_information_allocation_preserved=True,new_physical_mothers=0,new_ROOTs=0,
        documentation_sha256=document['files'],no_default_change=True,no_git_commit_or_push=True,
        limitations=['one training initialization','one calibrated conditional Gaussian model','exposed native transfer only','local K<=2; fixed tau; pure alpha']))
    print('RISK_FINAL_VALIDATION_PASSED',len(new),len(codes),flush=True)


if __name__=='__main__':validate()
