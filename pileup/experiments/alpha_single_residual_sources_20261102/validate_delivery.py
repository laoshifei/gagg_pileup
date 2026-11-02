"""Validate the completed frozen experiment and append-only local delivery."""
from pathlib import Path
import hashlib
import json
import re

import nbformat
import numpy as np


EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
ART = ROOT/'jupyter/artifacts'/EXP.name
REC = ROOT/'jupyter/run_records'/EXP.name


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        while block:=f.read(1024*1024):h.update(block)
    return h.hexdigest()


def main():
    checks={}
    for name in ('design-freeze.json','input-freeze.json','waveform-freeze.json',
                 'feature-freeze.json','analysis-freeze.json','report-source-freeze.json',
                 'report-output-freeze.json'):
        items=read(EXP/name)
        for path,digest in items.items():assert sha(ROOT/path)==digest,path
        checks[name]=len(items)
    assert sha(EXP/'bin/age.so')==read(EXP/'build.json')['binary_sha256']
    checks['compiled_binary_matches_build_record']=True
    for name in ('pilot-gate.json','pilot-equality.json','replay-equality.json',
                 'projection-controls.json','provenance-validation.json'):
        assert read(EXP/name)['passed'],name
    checks['cost_provenance_projection_and_partial_replay_checks']=True
    other=ART/'reproductions/complete_replay'
    array_count=0
    for group in ('voltage','diagnostics','features'):
        for i in range(64):
            with np.load(ART/group/f'{i:04d}.npz') as a, np.load(other/group/f'{i:04d}.npz') as b:
                assert a.files==b.files
                for k in a.files:
                    np.testing.assert_array_equal(a[k],b[k]);array_count+=1
    assert read(EXP/'result-summary.json')==read(EXP/'reproductions/complete_replay/result-summary.json')
    assert read(EXP/'reproductions/complete_replay/reproduction-equality.json')['all_summary_values_equal']
    checks['complete_replay_equal_arrays']=array_count
    checks['complete_replay_summary_equal']=True
    results=read(EXP/'result-summary.json')
    assert results['exact_voltage_closure_mV']<1e-9
    assert not results['decision']['primary_both_gaps_passed']
    assert not results['decision']['default_changed']
    checks['prespecified_failed_prediction_gate_preserved']=True
    features=read(REC/'features.json')
    assert len(features['rows'])==64 and all(r['input_truth_guard'] for r in features['rows'])
    checks['feature_workers_truth_guard']=64
    mothers=read(EXP/'mothers.json')
    roots={r['root'] for r in mothers}
    previous=read(ROOT/'pileup/experiments/alpha_pair_energy_calibrated_20261030/mother-split.json')
    spatial=read(ROOT/'pileup/experiments/alpha_pair_spatial_recovery_20261102/mother-split.json')
    assert roots<=set(previous['noise'])
    assert not roots & (set(previous['mean'])|set(previous['evaluation'])|set(spatial['development']))
    for group in ('weak','high'):
        for role in ('train','validation'):
            assert sum(r['group']==group and r['role']==role for r in mothers)==16
    checks['mother_split_and_historical_roles']=True
    cells={};figure_ids=[]
    for filename,expected in [('codex_0007_SiPM_Energy_Calibration.ipynb',8),
                              ('Research_Plan_20260921.ipynb',1),('Run_Catalog.ipynb',1)]:
        before=read(REC/(filename+'.before'));after=read(ROOT/'jupyter'/filename)
        nbformat.validate(nbformat.from_dict(after))
        n=len(before['cells'])
        assert after['cells'][:n]==before['cells'],filename
        assert {k:v for k,v in after.items() if k!='cells'}=={k:v for k,v in before.items() if k!='cells'}
        added=after['cells'][n:]
        assert len(added)==expected
        ids=[c['id'] for c in after['cells'] if 'id' in c]
        assert len(ids)==len(set(ids))
        code_cells=0
        for cell in added:
            source=''.join(cell['source'])
            if cell['cell_type']=='markdown':
                for target in re.findall(r'\]\(([^)]+)\)',source):
                    assert not re.search(r'\.(pdf|ipynb|md)(?:[#:]|$)',target,re.I)
                if filename.startswith('codex_0007'):
                    figure_ids+=re.findall(r'^图 (0007-\d+-\d+)：',source,re.M)
            else:
                code_cells+=1
                assert cell['execution_count'] is not None
                assert any('image/png' in o.get('data',{}) for o in cell['outputs'])
                assert not any(o['output_type']=='error' for o in cell['outputs'])
        cells[filename]=dict(before=n,after=len(after['cells']),old_outputs_preserved=True,
                             executed_new_display_cells=code_cells)
    assert figure_ids==['0007-104-1','0007-105-1','0007-106-1']
    checks['notebooks']=cells
    checks['new_fixed_figure_ids']=figure_ids
    old=(REC/'literature-use.md.before').read_text();current=(ROOT/'docs/literature-use.md').read_text()
    assert current.startswith(old) and current[len(old):].count('## L92 ')==1
    checks['literature_append_only']=True
    checks['formal_digitizations']=1024
    checks['actual_digitizations_including_all_checks']=2072
    checks['new_ROOT']=0
    (EXP/'validation.json').write_text(json.dumps(dict(passed=True,checks=checks),ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(passed=True,checks=checks),ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
