"""Validate delivered sources, paired outputs and append-only local notebooks."""
from pathlib import Path
import ast
import hashlib
import json
import re
import subprocess
import importlib.util

import nbformat
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
EXP=Path(__file__).resolve().parent
NAME=EXP.name
ART=ROOT/'jupyter/artifacts'/NAME
REC=ROOT/'jupyter/run_records'/NAME


def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    counts={};namespace=read(EXP/'namespace-revision.json')
    for name in ('freeze.json','input-freeze.json','projection-freeze.json','report-freeze.json','report-format-freeze.json'):
        records=read(EXP/name)
        for path,digest in records.items():
            actual=namespace['historical_report_hash_path_mapping'].get(path,path) if name=='report-freeze.json' else path
            assert sha(ROOT/actual)==digest,(name,path)
        counts[name]=len(records)
    for record in read(EXP/'report-format-revision.json')['tables']:
        assert sha(ROOT/record['path'])==record['after_sha256'],record['path']
        assert record['data_identical']
    for path,record in namespace['root_modules'].items():
        source=(ROOT/path).read_text();snapshot=(ROOT/record['snapshot']).read_text()
        assert sha(ROOT/path)==record['sha256'],path
        if path.endswith('_study.py'):
            source=source.replace('alpha_pair_residual_projection_study','alpha_pair_residual_study')
        else:
            source=source.replace('import alpha_pair_residual_projection_study as s','import alpha_pair_residual_study as s')
        assert source==snapshot,path
    core='pileup/analysis/alpha_pair_residual_projection.py'
    assert (ROOT/core).read_bytes()==(EXP/'source'/core).read_bytes()
    for old in ('alpha_pair_residual_study.py','alpha_pair_residual_report.py'):
        path='pileup/analysis/'+old
        assert (ROOT/path).read_bytes()==subprocess.check_output(['git','show','402f792:'+path],cwd=ROOT)
    for parent in ('alpha_pair_spatial_energy_20261102','alpha_pair_spatial_recovery_20261102'):
        delivery=read(ROOT/'pileup/experiments'/parent/'delivery-freeze.json')
        for group in ('source_configuration_reports','local_artifacts'):
            for path,digest in delivery[group].items():assert sha(ROOT/path)==digest,path
            counts[parent+'/'+group]=len(delivery[group])
    tasks=read(EXP/'fit-plan.json')
    fits=[read(ART/'fits'/f"{t['id']:04d}.json") for t in tasks]
    assert len(fits)==88 and len(list((ART/'fits').glob('*.json')))==88
    for task,row in zip(tasks,fits):
        for key,value in task.items():assert row[key]==value,(task['id'],key)
        assert row['true_energy_in_fit']==(task['start']=='truth')
        assert row['ancestor_in_fit']==(task['delete_term'] is not None)
        assert np.isfinite(row['energy_MeV']).all() and np.isfinite(row['cost'])
    protocol=read(EXP/'protocol.json')
    for batch,key in (('pilot','pilot_ids'),('replay','replay_ids')):
        for tid in protocol[key]:
            old,new=[read(ART/p/f'{tid:04d}.json') for p in ('fits',batch)]
            for field in ('energy_MeV','cost','converged','boundary','start_energy_MeV'):
                assert old[field]==new[field],(batch,tid,field)
    rows=[r for p in (ART/'projections').glob('*.json') for r in read(p)['rows']]
    assert len(rows)==120
    numerical={key:max(r[key] for r in rows if r[key] is not None) for key in
        ('residual_closure_mV','old_cost_relative_error','additivity_error','profile_error','small_step_absolute_log_error')}
    assert numerical['residual_closure_mV']<1e-9 and numerical['old_cost_relative_error']<1e-7
    assert numerical['small_step_absolute_log_error']<1e-4
    spec=importlib.util.spec_from_file_location('projection_checks',EXP/'source'/core)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    random=np.random.default_rng(731200).normal(size=(4,8))
    paired=np.stack((random,-random),axis=-1)
    draws,_=module.paired_bootstrap(paired,seed=731201,draws=1000)
    numerical['paired_bootstrap_cancellation_error']=float(np.max(abs(draws.sum(-1))))
    assert numerical['paired_bootstrap_cancellation_error']<1e-12
    notebooks=[]
    for name in ('codex_0007_SiPM_Energy_Calibration.ipynb','Research_Plan_20260921.ipynb','Run_Catalog.ipynb'):
        old=read(REC/(name+'.before'));new=read(ROOT/'jupyter'/name)
        assert new['cells'][:len(old['cells'])]==old['cells'],name
        assert new['metadata']==old['metadata'],name
        nbformat.validate(nbformat.from_dict(new))
        for cell in new['cells'][len(old['cells']):]:
            if cell['cell_type']=='markdown':
                text=''.join(cell['source'])
                assert not re.search(r'\]\([^)]*\.(?:ipynb|md|pdf)(?:[)#?])',text,re.I),name
        notebooks.append(dict(name=name,old_cells=len(old['cells']),new_cells=len(new['cells']),prefix_exact=True))
    nb=read(ROOT/'jupyter/codex_0007_SiPM_Energy_Calibration.ipynb')
    markdown='\n'.join(''.join(c['source']) for c in nb['cells'] if c['cell_type']=='markdown')
    for section in (94,95,96,97):assert markdown.count(f'图 0007-{section}-1：')==1,section
    figures=list((ART/'figures').glob('*.png'))
    assert len(figures)==4 and len(list((ART/'figures').glob('*.svg')))==4
    for path in figures:assert path.stat().st_size>10000
    for path in [*EXP.rglob('*.py'),*[ROOT/'pileup/analysis'/f'alpha_pair_residual_projection{s}.py' for s in ('','_study','_report')]]:
        ast.parse(path.read_text())
        assert not (path.stat().st_mode&0o111),str(path)
    for path in ['AGENTS.md','jupyter/Research_Plan_20260921.ipynb','jupyter/Run_Catalog.ipynb']:
        assert subprocess.run(['git','check-ignore','-q',path],cwd=ROOT).returncode==0,path
    staged=subprocess.check_output(['git','diff','--cached','--name-only'],cwd=ROOT,text=True).splitlines()
    assert not any(p=='AGENTS.md' or p.startswith('jupyter/') or p.endswith('.ipynb') for p in staged)
    subprocess.run(['git','diff','--check'],cwd=ROOT,check=True)
    result=dict(passed=True,hash_counts=counts,numerical=numerical,formal_fits=88,pilot_fits=4,replay_fits=4,
        projections=120,notebooks=notebooks,figures=4,default_changed=False,new_ROOT=0,
        namespace_revision_checked=True,old_tracked_modules_restored=True,
        limits='Historically exposed4mother pairs;known K/t. Truth used in projection/oracle fits;no independent or blind validation.')
    (EXP/'validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False))


if __name__=='__main__':main()
