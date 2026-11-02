"""Read-only audit of frozen stages, native provenance, statistics and reports."""
from pathlib import Path
import json
import os
import re
import sys

import nbformat
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
os.environ['GAGG_CURRENT_SHARED_ROOT']=str(ROOT)
sys.path.insert(0,str(HERE/'source'))
from pileup.analysis import alpha_current_shared_study as study
from pileup.analysis.alpha_transport_environment_study import load_root
from pileup.analysis.alpha_current_shared import make_plans

read,sha=study.read,study.sha
checks=dict(design=study.verify())
assert read(HERE/'unit-tests.json')['passed']==2
reproduction=read(HERE/'reproduction-check.json')
assert reproduction['passed'] and reproduction['existing_target_refused'] and reproduction['new_ROOTs']==0
shared=(HERE/'pair-result.json').exists()
names=['native-output-freeze.json','single-output-freeze.json','single-fit-freeze.json',
       'single-result-freeze.json','report-source-freeze.json','report-output-freeze.json']
if shared:names+=['prediction-freeze.json','pair-on-freeze.json','pair-off-freeze.json','pair-result-freeze.json']
for name in names:
    for path,digest in read(HERE/name).items():assert sha(ROOT/path)==digest,(name,path)

native,single,pairs=make_plans()
assert native==read(HERE/'native-plan.json') and single==read(HERE/'single-plan.json')
assert pairs==read(HERE/'pair-plan.json')
manifest=read(HERE/'manifest.json');mothers=read(HERE/'mothers.json')
assert not manifest['failed_attempts'] and len(manifest['runs'])==len(mothers)==136
assert [r['id'] for r in mothers]==list(range(136))
assert len({r['root'] for r in mothers})==136 and len({r['seed'] for r in mothers})==136
spacing=np.diff(sorted(r['started_epoch'] for r in mothers));assert spacing.min()>=1
for row in mothers:
    assert all(row[k]==v for k,v in native[row['id']].items())
    assert row['status']=='completed'
    assert re.fullmatch(r'g4/data/\d{8}_\d{2}h\d{2}m\d{2}s\.root',row['root'])
    info,primary,photons,voltage,scalars=load_root(ROOT/row['root'])
    assert info==read(HERE/'original-runinfo'/f"{row['id']:04d}.json")
    assert info['macro_contents']==(ROOT/row['macro']).read_text()
    assert sha(ROOT/row['macro'])==row['macro_sha256']
    assert info['random_seed']==row['seed'] and info['sipm_avalanche_base_seed']==row['sipm_seed']
    assert info['sipm_readout_base_seed']==row['readout_seed']
    assert info['t_length_ps']==20000000 and info['n_events']==1 and primary['primary_time_ps']==5000000
    assert np.isclose(primary['primary_energy_MeV'],row['source_energy_MeV'],rtol=2e-7)
    assert primary['edep_total_MeV']==row['edep_MeV']
    assert json.loads(info['sipm_config_json'])==read(HERE/'source/g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json')
    assert abs(scalars['charge_balance_error_C'])<1e-18 and all(np.isfinite(a).all() for a in voltage.values())
    with np.load(study.ART/'photons'/f"{row['id']:04d}.npz") as cache:
        assert np.array_equal(cache['time_ns'],photons['time_ps'].astype(float)/1000-row['first_active_ns'])
        for key in ('x_mm','y_mm','energy_eV'):assert np.array_equal(cache[key],photons[key])

pilot=read(HERE/'pilot-gate.json');replay=read(HERE/'replay.json')
assert pilot['passed'] and pilot['efficacy_outcomes_not_used'] and len(pilot['readout_audits'])==8
assert pilot['digitizations']==256 and pilot['fits']==72 and replay['passed']
assert replay['single_energy_fits']==36
assert read(HERE/'single-generation.json')['digitizations']==4096
assert read(HERE/'single-fits.json')['fits']==9216
count=0
for task in single:
    if task['role']!='formal':continue
    fit=read(study.ART/'fits'/f"{task['id']:04d}.json")
    assert fit['truth_guard'] and not fit['true_E_input'] and fit['known_K1_t']
    assert len(fit['rows'])==18
    for a,b in zip(fit['rows'][::2],fit['rows'][1::2]):
        assert a['method']=='old' and b['method']=='corrected' and a['repeat']==b['repeat']
        assert a['metric_signature']==b['metric_signature'] and a['initial_MeV']==b['initial_MeV']
        assert a['known_time_ns']==b['known_time_ns']
    count+=len(fit['rows'])
assert count==9216
result=read(HERE/'single-result.json')
with np.load(study.ART/'single-summary.npz') as data:
    for row in result['waveform']:
        prefix=f"g{row['group']}_p{row['phase_ns']}_"
        old=data[prefix+'old_residual_mV'];new=data[prefix+'corrected_residual_mV']
        assert old.shape==new.shape==(32,1024)
        assert np.isclose(row['group_mean_mse_ratio'],np.mean(new.mean(0)**2)/np.mean(old.mean(0)**2),rtol=1e-14)

digitizations=4096+256+replay['digitizations']
if shared:
    assert read(HERE/'single-gate.json')['passed'] and read(HERE/'pair-pilot-gate.json')['passed']
    assert all(r['max_direct_difference_mV']<1e-10 for r in read(HERE/'pair-pilot-gate.json')['prediction_rows'])
    assert read(HERE/'pair-on-generation.json')['digitizations']==15360
    assert read(HERE/'pair-off-generation.json')['digitizations']==3840
    assert len(read(HERE/'predictions.json')['rows'])==640
    details=read(HERE/'pair-details.json');assert len(details['on'])==1280 and len(details['off'])==160
    for row in read(HERE/'pair-result.json')['groups']:
        selected={m:[r for r in details['on'] if r['method']==m and all(r[k]==row[k] for k in ('condition','gap_ns','phase_ns'))]
                  for m in ('old','corrected')}
        assert all([r['block'] for r in selected[m]]==list(range(32)) for m in selected)
        for key,ratio in row['mse_ratios'].items():
            assert np.isclose(ratio,np.mean([r[key] for r in selected['corrected']])/np.mean([r[key] for r in selected['old']]),rtol=1e-14)
    assert replay['exact_arrays']==18 and replay['digitizations']==22
    digitizations+=19200+36
else:
    assert not read(HERE/'single-gate.json')['passed']
    assert not (HERE/'pair-pilot-gate.json').exists()
    assert not list((study.ART/'predictions').glob('*.npz')) and not list((study.ART/'voltage/pair').glob('*.npz'))
    assert replay['exact_arrays']==6 and replay['digitizations']==4

notebooks={}
for name,old_count,new_count in [('codex_0007_SiPM_Energy_Calibration.ipynb',210,220 if shared else 218),
                              ('Research_Plan_20260921.ipynb',106,107),('Run_Catalog.ipynb',67,68)]:
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

checks.update(passed=True,new_ROOTs=136,formal_mothers=128,pilot_mothers=8,shared_executed=shared,
              minimum_launch_spacing_seconds=float(spacing.min()),total_offline_digitizations=digitizations,
              formal_single_fits=9216,pilot_single_fits=72,replay_single_fits=36,
              total_single_fits=9324,notebooks=notebooks,fresh_rerun_inputs_verified=True)
print(json.dumps(checks,ensure_ascii=False,indent=2))
