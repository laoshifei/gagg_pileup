"""Read-only verification of the fixed energy pilot and local report."""
from pathlib import Path
import json
import os
import re
import sys
import nbformat
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
os.environ['GAGG_CURRENT_PAIR_ENERGY_ROOT']=str(ROOT)
sys.path.insert(0,str(HERE/'source'))
from pileup.analysis import alpha_current_pair_energy_study as study
from pileup.analysis.alpha_current_pair_energy import selected_blocks,paired_statistics
from pileup.analysis.alpha_transport_environment_study import load_root

read,sha=study.read,study.sha
checks=dict(frozen_inputs=study.verify())
for name in ('upstream-input-freeze.json','resume-source-freeze.json','resume-retained-freeze.json',
             'single-fit-freeze.json','pair-fit-freeze.json','result-freeze.json','report-source-freeze.json','report-output-freeze.json'):
    for path,digest in read(HERE/name).items():assert sha(ROOT/path)==digest,(name,path)
assert read(HERE/'numerical.json')['passed'] and read(HERE/'numerical.json')['fits']==4
assert read(HERE/'single-control.json')['exact_fits_except_timing']==32
assert read(HERE/'replay.json')['passed'] and read(HERE/'replay.json')['fits']==6
assert read(HERE/'cost-preflight.json')['pilot_fits']==6
assert read(HERE/'unit-tests.json')['passed']==2
assert read(HERE/'reproduction-check.json')['passed']
assert read(HERE/'pair-fits.json')['retained_outputs']==29 and read(HERE/'pair-fits.json')['resumed_fits']==83
tasks=read(HERE/'fit-plan.json');cases=read(HERE/'cases.json')
assert len(tasks)==144 and len(cases)==48
assert read(HERE/'protocol.json')['cohort']['blocks']==selected_blocks()
for task in tasks:
    assert set(task)=={'id','case','kind','times_ns','repeat','method','voltage','voltage_sha256'}
    result=read(study.ART/'fits'/f"{task['id']:04d}.json")
    assert result['truth_guard'] and not result['true_E_input']
    assert result['known_K']==len(task['times_ns'])
    assert result['voltage_sha256']==sha(ROOT/task['voltage'])
for c in cases:
    assert sha(ROOT/c['source_voltage'])==sha(ROOT/c['voltage'])==c['voltage_sha256']
    with np.load(ROOT/c['voltage']) as data:assert data['voltage_mV'].shape==(8,1024)
    for repeat in c['repeats']:
        chosen=[t for t in tasks if t['case']==c['id'] and t['repeat']==repeat]
        assert len(chosen)==2
        a,b=[read(study.ART/'fits'/f"{t['id']:04d}.json") for t in chosen]
        for field in ('metric_signature','voltage_array_sha256','initial_MeV','start_energy_MeV','times_ns'):
            assert a[field]==b[field],(c['id'],repeat,field)
        assert a['starts']==b['starts']==(3 if c['kind']=='pair' else 1)

parents=read(HERE/'input-mothers.json');assert len(parents)==20
for row in parents:
    assert sha(ROOT/row['root'])==row['root_sha256']
    assert sha(ROOT/row['snapshot_macro'])==row['macro_sha256']
    info,primary,photons,voltage,scalars=load_root(ROOT/row['root'])
    assert info==read(HERE/'original-runinfo'/f"{row['id']:04d}.json")
    assert info['macro_contents']==(ROOT/row['snapshot_macro']).read_text()
    assert info['n_events']==1 and info['t_length_ps']==20000000
    assert info['random_seed']==row['seed'] and primary['edep_total_MeV']==row['edep_MeV']

details=read(HERE/'fit-details.json');summary=read(HERE/'result-summary.json')
assert len(details)==112
for row in details:
    assert np.isclose(row['weak_error_percent'],100*(row['weak_fit_MeV']/row['weak_true_MeV']-1))
    assert row['strict']==(max(abs(row['weak_error_percent']),abs(row['strong_error_percent']))<=10 and row['converged'])
blocks=selected_blocks();index=np.random.default_rng(741900000).integers(8,size=(10000,8))
for pooled in summary['short_pooled']:
    arrays={}
    for method in ('old','corrected'):
        rows=[[next(r for r in details if r['block']==block and r['method']==method and r['kind']==pooled['kind']
                    and r['condition']=='high_weak' and r['gap_ns']==gap) for gap in (8.,40.)] for block in blocks]
        arrays[method]=np.array([[r['weak_error_percent'] for r in col] for col in rows])
        arrays[method+'_success']=np.array([[r['strict'] for r in col] for col in rows])
    calculated=paired_statistics(arrays['old'],arrays['corrected'],arrays['old_success'],arrays['corrected_success'],index)
    assert calculated==pooled['comparison']

notebooks={}
for name,old_count,new_count in [('codex_0007_SiPM_Energy_Calibration.ipynb',220,230),
                              ('Research_Plan_20260921.ipynb',107,108),('Run_Catalog.ipynb',68,69)]:
    before=read(study.REC/(name+'.before'));after=read(ROOT/'jupyter'/name)
    assert len(after['cells'])==new_count and before['cells']==after['cells'][:old_count]
    assert {k:v for k,v in before.items() if k!='cells'}=={k:v for k,v in after.items() if k!='cells'}
    nbformat.validate(nbformat.from_dict(after))
    for cell in after['cells'][old_count:]:
        if cell['cell_type']=='markdown':
            assert not re.search(r'\]\([^)]*\.(?:md|pdf|ipynb)(?:[#?)])',''.join(cell['source']),re.I)
        else:
            assert cell['execution_count'] is not None and any('image/png' in o.get('data',{}) for o in cell['outputs'])
    notebooks[name]=new_count
checks.update(passed=True,reused_ROOTs=20,new_ROOTs=0,new_digitizations=0,independent_mother_pairs=8,
              formal_pair_fits=112,formal_single_fits=32,cost_pilot_fits=6,numerical_fits=4,replay_fits=6,
              total_fits=160,notebooks=notebooks,energy_gate_passed=read(HERE/'gate.json')['passed'])
print(json.dumps(checks,ensure_ascii=False,indent=2))
