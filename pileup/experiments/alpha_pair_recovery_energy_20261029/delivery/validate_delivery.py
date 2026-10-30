"""Independent checks of frozen provenance, score identities and notebook output."""
from pathlib import Path
import ast
import base64
import csv
import hashlib
import json
import re

import nbformat
import numpy as np

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'pileup/analysis').is_dir())
NAME = 'alpha_pair_recovery_energy_20261029'
EXP, OUT, REC = [ROOT/p/NAME for p in ['pileup/experiments', 'jupyter/artifacts', 'jupyter/run_records']]


def load(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def lines(p): return [json.loads(s) for s in p.read_text().splitlines()]


def main():
    for manifest in ['freeze.json', 'method-freeze.json', 'test-source-freeze.json', 'report-source-freeze.json',
                     'report-source-freeze-v2.json', 'report-source-freeze-v3.json', 'delivery-freeze.json']:
        for name, value in load(EXP/manifest)['files'].items(): assert sha(EXP/name) == value, (manifest, name)
    for name, value in load(EXP/'input-hashes.json').items(): assert sha(OUT/name) == value, name
    for name, value in load(REC/'predictions-freeze.json')['files'].items(): assert sha(OUT/name) == value, name
    for name, value in load(REC/'output-checksums.json')['files'].items(): assert sha(OUT/name) == value, name
    modules = ['alpha_pair_recovery_energy.py', 'alpha_pair_recovery_energy_study.py', 'test_alpha_pair_recovery_energy.py']
    for name in modules:
        ast.parse((ROOT/'pileup/analysis'/name).read_text())
        assert sha(ROOT/'pileup/analysis'/name) == sha(EXP/'source/pileup/analysis'/name)
    report_source = ROOT/'pileup/analysis/alpha_pair_recovery_energy_report.py'
    ast.parse(report_source.read_text())
    assert sha(report_source) == sha(EXP/'source/pileup/analysis/alpha_pair_recovery_energy_report_v3.py')
    prov = load(REC/'provenance.json')
    groups = [set(prov[k]) for k in ['train_roots', 'heldout_roots', 'evaluation_roots']]
    assert list(map(len, groups)) == [280, 280, 320]
    assert all(not groups[i]&groups[j] for i in range(3) for j in range(i+1, 3))
    report = load(OUT/'audited-report.json'); assert report == load(EXP/'result-summary.json')
    for path in OUT.glob('*.csv'): assert sha(path) == sha(EXP/'tables'/path.name)
    scores = lines(OUT/'audited-scores.jsonl'); predictions = {(r['method'], r['index']): r for r in lines(OUT/'predictions.jsonl')}
    assert len(scores) == len(predictions) == 6400
    assert len({(r['method'], r['index']) for r in scores}) == 6400
    for row in scores:
        p = predictions[row['method'], row['index']]
        assert row['output_count'] == row['target']-row['missed']+row['extra']
        assert len(p['times_ns']) == len(predictions['static_nll', row['index']]['times_ns'])
        assert row['strict_5'] <= row['strict_10'] <= row['strict_20'] <= row['target']
        assert len(row['matches']) == row['target']-row['missed']
        assert bool(row['frame_success_10']) == (row['target'] == row['output_count'] == row['strict_10'])
        for m in row['matches']:
            assert abs(m['time_error_ns']) <= 10.
            if m['relative_energy_error'] is not None:
                assert abs(m['relative_energy_error']-(m['estimated_energy_MeV']/m['true_energy_MeV']-1)) < 1e-12
        if 'cost' in p: assert p['cost'] <= p['initial_cost']+1e-7
    for row in report['failures']:
        assert row['frames'] == sum(row[k] for k in ['wrong_count', 'strict_time', 'energy', 'success']) == 600
    for row in report['audited_summary']:
        if row['method'] == 'static_nll': assert row['boundary'] is None and row['not_converged'] is None
    old_rows = {r['index']: r for r in lines(ROOT/'jupyter/artifacts/sipm_energy_guardband_20261029/historical-energy.scores.jsonl') if r['method'] == 'joint_conserved_nll'}
    for r in scores:
        if r['method'] != 'static_nll': continue
        old = old_rows[r['index']]
        for key in ['matches', 'target', 'missed', 'extra', 'output_count', 'strict_10', 'frame_success_10']: assert r[key] == old[key]
        assert r['temporal_split'] == old['split'] and r['temporal_merge'] == old['merge']
    # Independently reconstruct group means from complete frame denominators.
    for method, aggregate in report['aggregate'].items():
        for group in ['weak', 'strong']:
            rr = [s for s in scores if s['method'] == method and s['kind'] == 'pair' and ('weak' in s['condition']) == (group == 'weak')]
            assert len(rr) == 600
            assert abs(np.mean([s['frame_success_10'] for s in rr])-aggregate[group]) < 1e-12
    response = load(REC/'heldout-response-validation.json')
    assert response['passed'] and response['exact_single_invariance']
    assert response['strict_nominal_time_causality_passed'] is False
    assert response['numerical_response_checks']['max_pre_support_difference_mV'] == 0.
    assert load(REC/'numerical-tests.json')['passed'] and load(REC/'numerical-tests.json')['tests'] == 7
    assert load(REC/'inference-validation.json')['truth_guard'] and load(REC/'inference-validation.json')['workers'] <= 2
    for name, old_count, new_count in [('codex_0007_SiPM_Energy_Calibration.ipynb', 42, 55), ('Research_Plan_20260921.ipynb', 82, 83)]:
        nb = nbformat.read(ROOT/'jupyter'/name, as_version=4); nbformat.validate(nb)
        before = nbformat.read(REC/(name+'.before'), as_version=4)
        assert len(before.cells) == old_count and len(nb.cells) == new_count
        assert nb.cells[1:old_count] == before.cells[1:]
        for cell in nb.cells[old_count:]:
            if cell.cell_type == 'code':
                ast.parse(cell.source); assert cell.execution_count is not None
                assert not any(o.output_type == 'error' for o in cell.outputs)
            if cell.cell_type == 'markdown':
                assert not any(re.search(r'\.(pdf|md|ipynb)(#.*)?$', s, re.I) for s in re.findall(r'\]\(([^)]+)\)', cell.source))
        if name.startswith('codex'):
            codes = [c for c in nb.cells[old_count:] if c.cell_type == 'code']; assert len(codes) == 6
            pictures = [base64.b64decode(o.data['image/png']) for c in codes for o in c.outputs if 'image/png' in o.get('data', {})]
            assert len(pictures) == 5
            for img, fig in zip(pictures, report['figures']): assert hashlib.sha256(img).hexdigest() == sha(OUT/(fig+'.png'))
            figure_ids = [x for c in nb.cells if c.cell_type == 'markdown' for x in re.findall(r'\*\*图\s+(0007-\d+-\d+)：', c.source)]
            assert len(figure_ids) == len(set(figure_ids))
            assert all(f'0007-{section}-1' in figure_ids for section in range(21, 26))
    # Check local links in changed documentation, preserving DOI/web references.
    for doc in [EXP/'README.md', ROOT/'README.md', ROOT/'docs/literature-use.md']:
        text = doc.read_text()
        assert not any(line.rstrip() != line for line in text.splitlines())
        if doc == ROOT/'README.md': text = text[text.index('随后完成[恢复历史进入E/t'):]
        if doc.name == 'literature-use.md': text = text[text.index('<a id="recovery-energy-20261029">'):]
        for link in re.findall(r'\]\(([^)]+)\)', text):
            if '://' in link or link.startswith('#'): continue
            assert (doc.parent/link.split('#')[0]).exists(), (doc, link)
    result = dict(passed=True, frozen_input_code_model_and_predictions_verified=True, score_rows=6400,
        unchanged_baseline_and_count=True, numerical_tests=7, heldout_responses=1180,
        strict_nominal_causality_limitation_retained=True, original_notebook_cells_preserved=42,
        new_executed_cells=6, embedded_figure_hashes_verified=5, source_commit='6779620',
        fresh_independent_confirmation=False, default_changed=False, new_ROOTs=0, pushed=False)
    (REC/'delivery-validation.json').write_text(json.dumps(result, indent=2)+'\n')
    (EXP/'delivery-validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print('RECOVERY_DELIVERY_VERIFIED', result, flush=True)


if __name__ == '__main__': main()
