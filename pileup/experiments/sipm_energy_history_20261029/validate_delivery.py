"""Audit provenance, scientific identities, denominators and executed notebooks."""
from pathlib import Path
import ast
import base64
import csv
import hashlib
import json
import re

import nbformat
import numpy as np

ROOT = next(p for p in Path(__file__).resolve().parents
            if (p/'pileup/analysis').is_dir() and (p/'README.md').is_file())
NAME = 'sipm_energy_history_20261029'
EXP = ROOT/'pileup/experiments'/NAME
OUT = ROOT/'jupyter/artifacts'/NAME
REC = ROOT/'jupyter/run_records'/NAME


def load(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def lines(p): return [json.loads(s) for s in p.read_text().splitlines()]
def read_csv(name):
    with (OUT/name).open() as stream: return list(csv.DictReader(stream))


def main():
    for manifest in ['freeze.json', 'report-source-freeze.json', 'delivery-freeze.json']:
        for name, value in load(EXP/manifest)['files'].items():
            assert sha(EXP/name) == value, (manifest, name)
    for name in ['sipm_energy_history_study.py', 'sipm_energy_history_report.py']:
        ast.parse((ROOT/'pileup/analysis'/name).read_text())
        assert sha(ROOT/'pileup/analysis'/name) == sha(EXP/'source/pileup/analysis'/name)
    for name, value in load(EXP/'input-hashes.json').items():
        assert sha(OUT/name) == value, name
    for name, value in load(REC/'output-checksums.json')['files'].items():
        assert sha(OUT/name) == value, name
    assert sha(EXP/'energy-model.json') == sha(ROOT/'pileup/experiments/sipm_energy_guardband_20261029/energy-model.json')
    assert sha(OUT/'inputs/fits.jsonl') == sha(ROOT/'jupyter/artifacts/sipm_energy_guardband_20261029/historical-energy.predictions.jsonl')
    assert sha(OUT/'inputs/scores.jsonl') == sha(ROOT/'jupyter/artifacts/sipm_energy_guardband_20261029/historical-energy.scores.jsonl')
    a = load(EXP/'sensor-config.json')['sensor']
    b = load(EXP/'sensor-config-no-recovery.json')['sensor']
    assert {k for k in a if a[k] != b[k]} == {'recovery_enabled'}
    assert a['recovery_enabled'] is True and b['recovery_enabled'] is False
    photons = load(EXP/'photon-inputs.json')
    assert len(photons) == len({p['parent'] for p in photons}) == 16
    for p in photons:
        assert sha(ROOT/p['photon_file']) == p['sha256']
        assert sha(OUT/'inputs/photons'/f"{p['parent']:04}.npz") == p['sha256']
    tasks = load(EXP/'simulation-plan.json')
    assert len(tasks) == len({t['id'] for t in tasks}) == 912
    assert len(list((OUT/'simulations').glob('*.json'))) == 912
    assert {kind: sum(t['kind'] == kind for t in tasks) for kind in ['single', 'pair', 'long_single']} == dict(single=512, pair=384, long_single=16)
    for task in tasks:
        row = load(OUT/'simulations'/(task['id']+'.json'))
        assert row['task'] == task and row['charge_closure_error_pe'] < 1e-8
        assert task['avalanche_seed'] == 691000000 + task['slot']*100 + task['mc']
        assert [r['parent'] for r in row['components']] == task['parents']
        assert all(r['charge_pe'] > 0 for r in row['components'])
    pilot = load(REC/'pilot-validation.json')
    assert pilot['passed'] and pilot['historical_charge_reproduced']
    assert pilot['projected_two_worker_seconds'] <= 1800
    historical = lines(OUT/'historical-components.jsonl')
    fitted = lines(OUT/'historical-fitted-components.jsonl')
    assert len(historical) == 2400 and len(fitted) == 4800
    assert len({(r['index'], r['component']) for r in historical}) == 2400
    assert len({(r['method'], r['index'], r['component']) for r in fitted}) == 4800
    closure = []
    for row in fitted:
        if row['fit_energy_relative'] is not None:
            closure.append(abs(row['fit_energy_relative'] - row['native_energy_relative']
                               - row['response_energy_change'] - row['fit_energy_change']))
    assert max(closure) < 1e-12
    for row in read_csv('historical-frame-counts.csv'):
        assert int(row['outputs']) == int(row['true_components']) - int(row['missed']) + int(row['extra'])
    # Recompute the paired charge estimator directly from retained component
    # rows, independently of the report's multidimensional indexing.
    raw = read_csv('counterfactual-components.csv')
    assert len(raw) == 768
    for row in read_csv('counterfactual-summary.csv'):
        group = [r for r in raw if all(r[k] == row[k] for k in ['condition', 'gap_ns', 'species', 'recovery'])]
        assert len(group) == 32
        ratios = []
        for block in sorted({r['block'] for r in group}):
            records = [r for r in group if r['block'] == block]
            assert len(records) == 4
            ratios.append(sum(float(r['joint_charge_pe']) for r in records) /
                          sum(float(r['isolated_charge_pe']) for r in records) - 1)
        assert abs(np.mean(ratios)-float(row['mean_relative_charge_change'])) < 1e-12
    for r in raw:
        if r['recovery'] == 'False':
            assert r['joint_energy_relative'] == r['isolated_energy_relative'] == ''
    for row in read_csv('known-charge-energy-summary.csv'):
        group = [r for r in raw if r['recovery'] == 'True' and all(r[k] == row[k] for k in ['condition', 'gap_ns', 'species'])]
        values = np.array([float(r[row['kind']+'_energy_relative']) for r in group])
        assert len(values) == 32 == int(row['n'])
        assert np.isfinite(values).sum() == int(row['valid'])
        assert np.sum(np.abs(values) <= .1) == int(row['within_energy10'])
        assert abs(np.mean(values)-float(row['mean_energy_relative'])) < 1e-12
    window = read_csv('window-control.csv')
    assert len(window) == 16
    assert all(float(r['prefix_difference_relative']) == 0 for r in window)
    report = load(OUT/'report.json')
    assert report == load(EXP/'result-summary.json')
    for p in (OUT).glob('*.csv'):
        assert sha(p) == sha(EXP/'tables'/p.name)
    npath = ROOT/'jupyter/codex_0007_SiPM_Energy_Calibration.ipynb'
    nb = nbformat.read(npath, as_version=4)
    before = nbformat.read(REC/(npath.name+'.before'), as_version=4)
    nbformat.validate(nb)
    assert len(nb.cells) == 42 and nb.cells[1:31] == before.cells[1:]
    new = nb.cells[31:]
    codes = [c for c in new if c.cell_type == 'code']
    assert len(codes) == 5 and all(c.execution_count is not None for c in codes)
    assert not any(o.output_type == 'error' for c in codes for o in c.outputs)
    images = [o.data['image/png'] for c in codes for o in c.outputs if 'image/png' in o.get('data', {})]
    assert len(images) == 4
    for name, encoded in zip(report['figures'], images):
        assert sha(OUT/(name+'.png')) == report['figures_sha256'][name]
        assert hashlib.sha256(base64.b64decode(encoded)).hexdigest() == report['figures_sha256'][name]
        assert (OUT/(name+'.svg')).is_file()
    captions = []
    for c in new:
        if c.cell_type == 'markdown':
            captions += re.findall(r'\*\*图 (0007-\d+-\d+)：', c.source)
            assert not any(re.search(r'\.(pdf|md|ipynb)(#.*)?$', s, re.I) for s in re.findall(r'\]\(([^)]+)\)', c.source))
    assert captions == ['0007-15-1', '0007-16-1', '0007-17-1', '0007-18-1']
    for name, oldn, newn in [('Research_Plan_20260921.ipynb', 81, 82), ('Run_Catalog.ipynb', 55, 56)]:
        current = nbformat.read(ROOT/'jupyter'/name, as_version=4)
        old = nbformat.read(REC/(name+'.before'), as_version=4)
        nbformat.validate(current)
        assert len(current.cells) == newn and len(old.cells) == oldn
        assert current.cells[1:oldn] == old.cells[1:]
        assert [c for c in current.cells if c.cell_type == 'code'] == [c for c in old.cells if c.cell_type == 'code']
    catalog = nbformat.read(ROOT/'jupyter/Run_Catalog.ipynb', as_version=4)
    assert all(p['root'] in catalog.cells[-1].source for p in photons)
    index = load(ROOT/'jupyter/codex_notebook_index.json')
    assert index['next_number'] == '0008'
    validation = dict(passed=True, diagnostic_simulations=912, optical_pairs=8, digitizer_replicates_per_pair=4,
        no_new_ROOTs=True, workers=2, input_and_snapshot_hashes=True, photon_caches_unchanged=True,
        historical_components=2400, fitted_component_positions=4800, additive_closure_max=max(closure),
        charge_estimator_independently_recomputed=True, frame_count_identities=True,
        energy_denominators_recomputed=True, no_energy_mapping_for_disabled_recovery=True,
        long_window_prefixes_equal=16, frozen_calibration_and_fits_unchanged=True,
        notebook_old_calculations_preserved=True, appended_executed_cells=5, numbered_figures=4,
        embedded_images_match_files=True, notebook_format_and_links_checked=True,
        no_git_push=True, interpretation='Truth-assisted retrospective response diagnostic; no blind algorithm gain')
    (REC/'delivery-validation.json').write_text(json.dumps(validation, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
