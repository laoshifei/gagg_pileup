"""Read-only verification of both current-environment stages and fresh outputs."""
from pathlib import Path
import json
import os
import re
import sys

import nbformat
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
os.environ['GAGG_CURRENT_CONFIRMATION_ROOT'] = str(ROOT)
os.environ['GAGG_CURRENT_RESPONSE_ROOT'] = str(ROOT)
sys.path.insert(0, str(HERE / 'source'))
from pileup.analysis import alpha_current_confirmation_study as study
from pileup.analysis import alpha_current_response_study as control
from pileup.analysis.alpha_transport_environment_study import load_root

read, sha = study.read, study.sha
checks = dict(confirmation=study.verify(), controls=control.verify())
for experiment, names in [(HERE, ['native-output-freeze.json', 'digitization-freeze.json',
                                 'fit-output-freeze.json', 'result-freeze.json', 'report-source-freeze.json',
                                 'report-output-freeze.json']),
                          (control.EXP, ['digitization-freeze.json', 'fit-output-freeze.json',
                                         'result-freeze.json', 'analysis-environment-freeze.json'])]:
    for name in names:
        for path, digest in read(experiment / name).items():
            assert sha(ROOT / path) == digest, (name, path)

manifest = read(HERE / 'manifest.json')
assert not manifest['failed_attempts'] and len(manifest['runs']) == 344
assert len({r['root'] for r in manifest['runs']}) == 344
assert len({r['seed'] for r in manifest['runs']}) == 344
spacing = np.diff(sorted(r['started_epoch'] for r in manifest['runs']))
assert spacing.min() >= 1
assert {r['id'] for r in manifest['runs'] if r['role'] == 'pilot'} == set(range(8))
assert {r['id'] for r in manifest['runs'] if r['role'] == 'formal'} == set(range(8,344))
for row in manifest['runs']:
    assert row['status'] == 'completed'
    assert re.fullmatch(r'g4/data/\d{8}_\d{2}h\d{2}m\d{2}s\.root', row['root'])
    info, primary, photons, voltage, scalars = load_root(ROOT / row['root'])
    assert info['macro_contents'] == (ROOT / row['macro']).read_text()
    assert row['macro_sha256'] == sha(ROOT / row['macro'])
    assert info['random_seed'] == row['seed'] == 738000000+row['id']
    assert info['sipm_avalanche_base_seed'] == row['sipm_seed'] == row['seed']+100000
    assert info['sipm_readout_base_seed'] == row['readout_seed'] == row['seed']+200000
    assert info['t_length_ps'] == 20000000 and info['n_events'] == 1
    assert primary['primary_energy_MeV'] == 8 and primary['primary_time_ps'] == 5000000
    assert primary['edep_total_MeV'] == row['edep_MeV']
    assert np.isfinite(voltage['analog_mV']).all()
    with np.load(study.ART / 'photons' / f"{row['id']:04d}.npz") as cache:
        assert np.array_equal(cache['time_ns'], photons['time_ps'].astype(float)/1000-row['first_active_ns'])
        for key in ('x_mm','y_mm','energy_eV'):
            assert np.array_equal(cache[key], photons[key])

assert read(control.EXP / 'native-audit.json')['passed']
assert read(control.EXP / 'gate.json')['passed']
assert read(HERE / 'pilot-gate.json')['efficacy_outcomes_not_used_for_gate']
assert read(HERE / 'replay.json')['passed']
fit_counts = {}
for art, ids, label in [(control.ART, range(128), 'controls'), (study.ART, range(8,344), 'confirmation')]:
    count = 0
    for i in ids:
        data = read(art / 'fits' / f'{i:04d}.json')
        assert data['truth_guard'] and not data['true_E_input'] and data['known_K1_t']
        assert len(data['rows']) == 18
        for a,b in zip(data['rows'][::2],data['rows'][1::2]):
            assert a['method']=='old' and b['method']=='corrected' and a['repeat']==b['repeat']
            assert a['metric_signature']==b['metric_signature'] and a['initial_MeV']==b['initial_MeV']
            assert a['known_time_ns']==b['known_time_ns']
        count += len(data['rows'])
    fit_counts[label] = count
assert fit_counts == {'controls':2304,'confirmation':6048}

summary = read(HERE / 'result-summary.json')
with np.load(study.ART / 'summary.npz') as data:
    assert data['mother_ids'].tolist() == list(range(8,344))
    for j, row in enumerate(summary['projection']):
        old, new = data['old_projection_100logE'][:,j], data['corrected_projection_100logE'][:,j]
        assert np.isclose(row['mse_ratio'], np.mean(new**2)/np.mean(old**2), rtol=1e-14)
        assert np.isclose(row['old_mean'], old.mean(), rtol=1e-14)
        assert np.isclose(row['corrected_mean'], new.mean(), rtol=1e-14)

notebooks = {}
for name, old_count, new_count in [('codex_0007_SiPM_Energy_Calibration.ipynb',199,210),
                                  ('Research_Plan_20260921.ipynb',105,106),('Run_Catalog.ipynb',66,67)]:
    before = read(control.REC / (name+'.before'))
    after = read(ROOT / 'jupyter' / name)
    assert len(after['cells']) == new_count
    assert before['cells'] == after['cells'][:old_count]
    assert {k:v for k,v in before.items() if k!='cells'} == {k:v for k,v in after.items() if k!='cells'}
    nbformat.validate(nbformat.from_dict(after))
    for cell in after['cells'][old_count:]:
        if cell['cell_type']=='markdown':
            assert not re.search(r'\]\([^)]*\.(?:md|pdf|ipynb)(?:[#?)])',''.join(cell['source']),re.I)
        else:
            assert cell['execution_count'] is not None
            assert any('image/png' in output.get('data',{}) for output in cell['outputs'])
    notebooks[name] = new_count

checks.update(passed=True,new_ROOTs=344,formal_mothers=336,pilot_mothers=8,
              minimum_launch_spacing_seconds=float(spacing.min()),fit_counts=fit_counts,
              formal_digitizations=1024+2688,pilot_digitizations=4+64,replay_digitizations=4,
              formal_fits=2304+6048,pilot_fits=36+36,replay_fits=18,notebooks=notebooks)
print(json.dumps(checks,ensure_ascii=False,indent=2))
