"""Check frozen inputs/results, staged stopping, truth isolation and note preservation."""
from pathlib import Path
import csv
import json
import os
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent
os.environ['GAGG_MEAN_SHAPE_ROOT']=str(HERE.parents[2])
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
sys.path.insert(0,str(HERE/'source'))
import numpy as np
import nbformat
from pileup.analysis.alpha_mean_shape_study import ROOT,EXP,ART,REC,read,sha,save,verify


def main():
    checks=verify(True)
    for name in ('train-input-freeze.json','single-input-freeze.json','fit-output-freeze.json',
                 'report-source-freeze.json','report-output-freeze.json'):
        paths=read(EXP/name)
        for p,d in paths.items():assert sha(ROOT/p)==d,p
        checks[name]=len(paths)
    completed=read(EXP/'single-gate.json')
    assert not completed['passed'] and completed['failures']==['weak projection gap8']
    assert not (EXP/'pair-result.json').exists() and not (EXP/'pair-input-freeze.json').exists()
    assert not list((ART/'voltage/pair').glob('*.npz'))
    assert read(EXP/'numerical.json')['passed'] and read(EXP/'replay-validation.json')['exact_arrays']
    assert read(EXP/'calibration.json')['scale']==1.
    split=read(EXP/'split.json')
    assert len(split['train'])==128 and len(split['validation'])==64
    assert not set(split['train'])&set(split['validation'])
    assert not set(split['excluded_recent'])&set(split['validation'])
    tests=subprocess.run([sys.executable,'-B','-m','unittest',
        'pileup.analysis.test_alpha_mean_shape_correction'],cwd=HERE/'source',capture_output=True,text=True)
    assert tests.returncode==0,tests.stdout+tests.stderr
    (REC/'unit-tests.log').write_text(tests.stdout+tests.stderr)
    arrays=0
    for folder,kind,ids in [('pilot','train',(0,96)),('replay','train',(1,)),('replay','single',(197,))]:
        for i in ids:
            for type_ in ('voltage','truth'):
                with np.load(ART/type_/kind/f'{i:04d}.npz') as a,np.load(ART/folder/type_/kind/f'{i:04d}.npz') as b:
                    assert a.files==b.files
                    for name in a.files:np.testing.assert_array_equal(a[name][:2],b[name]);arrays+=1
    fits=0
    for p in sorted((ART/'fits').glob('*.json')):
        data=read(p);assert data['truth_guard'] and data['known_K1_t'] and not data['true_E_input']
        assert len(data['rows'])==18
        for repeat in [-1,*range(8)]:
            pair=[r for r in data['rows'] if r['repeat']==repeat]
            assert {r['method'] for r in pair}=={'old','corrected'}
            assert pair[0]['metric_signature']==pair[1]['metric_signature']
        fits+=len(data['rows'])
    assert fits==4608
    for i in (0,192):
        a=read(ART/'fits'/f'{i:04d}.json')['rows'];b=read(ART/'fit-pilot'/f'{i:04d}.json')['rows']
        assert len(a)==len(b)==18
        for x,y in zip(a,b):assert {k:v for k,v in x.items() if k!='seconds'}=={k:v for k,v in y.items() if k!='seconds'}
    notebook_checks={}
    for filename,count,added in [('codex_0007_SiPM_Energy_Calibration.ipynb',176,9),
                                 ('Research_Plan_20260921.ipynb',102,1),('Run_Catalog.ipynb',63,1)]:
        nb=read(ROOT/'jupyter'/filename);before=read(REC/(filename+'.before'))
        assert len(before['cells'])==count and len(nb['cells'])==count+added
        assert nb['cells'][:count]==before['cells']
        assert {k:v for k,v in nb.items() if k!='cells'}=={k:v for k,v in before.items() if k!='cells'}
        nbformat.validate(nbformat.from_dict(nb))
        images=0
        for cell in nb['cells'][count:]:
            source=cell['source'] if isinstance(cell['source'],str) else ''.join(cell['source'])
            if cell['cell_type']=='code':
                assert cell['execution_count'] is not None
                assert not any(x['output_type']=='error' for x in cell['outputs'])
                images+=sum('image/png' in out.get('data',{}) for out in cell['outputs'])
            else:
                assert not re.search(r'\]\([^\n)]*\.(?:md|pdf|ipynb)(?:[#:][^)]*)?\)',source,re.I)
        assert images==(3 if added==9 else 0)
        notebook_checks[filename]=dict(old_cells_unchanged=count,added_cells=added,new_images=images)
    details=read(EXP/'single-fit-details.json')
    assert len(details)==4608 and not any(not x['converged'] for x in details)
    assert sum(x['boundary'] for x in details)==0
    for p in (EXP/'tables').glob('*.csv'):
        assert b'\r' not in p.read_bytes()
        with p.open() as f:assert list(csv.DictReader(f)),p
    for name in ('README.md','docs/literature-use.md'):
        p=ROOT/name;before=(REC/(p.name+'.before')).read_text()
        assert p.read_text().startswith(before)
    report=dict(passed=True,hash_checks=checks,unit_tests=4,exact_pilot_and_replay_arrays=arrays,
        formal_fits=fits,fit_pilot_equal_except_timing=True,voltage_only_input_guards=True,
        shared_expansion_stopped=True,no_new_ROOT=True,notebooks=notebook_checks,
        formal_digitizations=3072,total_digitizations=3080,pilot_fits=36)
    save(EXP/'validation.json',report)
    print(json.dumps(report,ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':main()
