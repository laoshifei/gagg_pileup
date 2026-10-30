"""Audit the declared information supplied to each diagnostic and final outputs."""
from pathlib import Path
import ast
import base64
import hashlib
import json
import re

import nbformat
import numpy as np

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'pileup/analysis').is_dir())
NAME = 'alpha_pair_energy_oracle_20261029'
EXP, OUT, REC = [ROOT/p/NAME for p in ['pileup/experiments', 'jupyter/artifacts', 'jupyter/run_records']]


def load(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def lines(p): return [json.loads(s) for s in p.read_text().splitlines()]


def main():
    for manifest in ['freeze.json', 'execution-freeze.json', 'selection-freeze.json', 'report-freeze.json', 'reference-freeze.json', 'delivery-freeze.json']:
        for name, value in load(EXP/manifest)['files'].items(): assert sha(EXP/name) == value, (manifest, name)
    for name, value in load(EXP/'input-hashes.json').items(): assert sha(OUT/name) == value, name
    for stage in ['blind', 'oracle']:
        frozen = load(REC/(stage+'-freeze.json')); assert frozen['count'] == 1540
        for name, value in frozen['files'].items(): assert sha(OUT/name) == value, name
    for name, value in load(REC/'output-checksums.json')['files'].items(): assert sha(OUT/name) == value, name
    for name in ['alpha_pair_energy_oracle.py', 'alpha_pair_energy_oracle_report.py', 'alpha_pair_energy_oracle_reference.py', 'test_alpha_pair_energy_oracle.py']:
        assert sha(ROOT/'pileup/analysis'/name) == sha(EXP/'source/pileup/analysis'/name)
        ast.parse((ROOT/'pileup/analysis'/name).read_text())
    assert sha(ROOT/'pileup/analysis/alpha_pair_energy_oracle_study.py') == sha(EXP/'source/pileup/analysis/alpha_pair_energy_oracle_study_v2.py')
    prov = load(REC/'provenance.json')
    assert len(prov['fold0_roots']) == len(prov['fold1_roots']) == 280
    assert not set(prov['fold0_roots']) & set(prov['fold1_roots'])
    old = ROOT/'pileup/experiments/alpha_pair_recovery_energy_20261029'
    for name, value in prov['inherited_model'].items(): assert sha(EXP/name) == sha(old/name) == value
    assert load(REC/'pilot-validation.json')['passed'] is False
    assert load(REC/'execution-plan-validation.json')['passed'] is True
    plan = load(EXP/'execution-plan.json'); ids = set(plan['indices']); assert len(ids) == 1540
    cases = load(OUT/'inputs/cases.json')
    assert ids == {i for i, c in enumerate(cases) if c['fold'] == 1 or c['kind'] == 'single'}
    assert sum(cases[i]['fold'] == 1 for i in ids) == 1260
    expected_gate = [i for i in ids if cases[i]['fold'] == 0 and cases[i]['kind'] == 'single']
    assert sorted(load(EXP/'selection-freeze.json')['fold0_single_indices']) == sorted(expected_gate)
    assert sum(map(len, load(EXP/'gate.json')['scores'])) == 280
    report = load(OUT/'audited-report.json'); assert report == load(EXP/'result-summary.json')
    scores = lines(OUT/'scores.jsonl'); predictions = {(r['method'], r['index']): r for r in lines(OUT/'predictions.jsonl')}
    assert len(scores) == len(predictions) == 6160
    methods = ['blind_energy_gate', 'true_count_blind_start', 'fixed_true_time', 'free_from_true_time']
    assert set(predictions) == {(m, i) for m in methods for i in ids}
    times = load(OUT/'inputs/oracle-times.truth.json')
    for i in ids:
        blind = load(OUT/'blind'/f'{i:04}.json'); oracle = load(OUT/'oracle'/f'{i:04}.json')
        k = len(times[i]); expected = blind['candidates'][k] if len(blind['candidates']) > k else blind['candidates'][0]
        assert oracle['true_count_blind_start'] == expected
        assert oracle['supplied_true_count'] == k
        np.testing.assert_array_equal(oracle['fixed_true_time']['times_ns'], sorted(times[i]))
        assert oracle['free_from_true_time']['cost'] <= oracle['fixed_true_time']['cost']+1e-7
        for method in methods:
            row = predictions[method, i]
            assert row['cost'] <= row['initial_cost']+1e-7
            assert np.isfinite(row['energy_MeV']).all() and np.isfinite(row['times_ns']).all()
    for row in scores:
        assert row['output_count'] == row['target']-row['missed']+row['extra']
        assert row['strict_5'] <= row['strict_10'] <= row['strict_20'] <= row['target']
        assert bool(row['frame_success_10']) == (row['target'] == row['output_count'] == row['strict_10'])
        assert row['true_count_supplied'] == (row['method'] != 'blind_energy_gate')
        assert row['true_time_supplied'] == (row['method'] in ['fixed_true_time', 'free_from_true_time'])
    for aggregate in report['aggregate']:
        rr = [s for s in scores if s['method'] == aggregate['method'] and s['fold'] == 1 and s['kind'] == 'pair' and
              ('weak' in s['condition']) == (aggregate['group'] == 'weak')]
        assert len(rr) == aggregate['frames']
        assert sum(s['frame_success_10'] for s in rr) == aggregate['successful_pairs']
        assert abs(np.mean([s['frame_success_10'] for s in rr])-aggregate['strict_pair_success']) < 1e-12
    blocks = load(REC/'mother-block-design.json'); assert blocks['count'] == 110
    assert sorted(j for group in blocks['strata'] for j in group) == list(range(110))
    root_sets = list(map(set, blocks['blocks'])); assert sum(map(len, root_sets)) == 280
    assert len(set.union(*root_sets)) == 280
    for i, b in blocks['frame_blocks'].items(): assert set(cases[int(i)]['mother_roots']) <= root_sets[b]
    for row in report['paired_intervals']:
        assert row['frames'] == row['both_success']+row['gained']+row['lost']+row['both_failed']
        assert abs((row['gained']-row['lost'])/row['frames']-row['difference']) < 1e-12
    identity = load(REC/'input-identity-validation.json')
    assert identity['passed'] and identity['isolated_roots'] == 280 and identity['energy_identity_checks'] == 1800
    reference = report['isolated_reference']
    assert reference == load(OUT/'isolated-reference.json')
    assert reference['pairs'] == 900 and reference['single_mothers'] == 280 and reference['true_energy_identity_verified']
    for row in reference['totals']:
        assert row['pairs'] == (300 if row['group'] == 'weak' else 600)
        assert row['shared_fixed_time_pair_energy10']-row['isolated_both_energy10'] == row['gained_after_pairing']-row['lost_after_pairing']
        agg = next(r for r in report['aggregate'] if r['method'] == 'fixed_true_time' and r['group'] == row['group'])
        assert row['shared_fixed_time_pair_energy10'] == agg['successful_pairs']
    numerical = load(REC/'numerical-validation.json')
    assert numerical['passed'] and numerical['tests'] == 5 and len(numerical['real_template_noiseless_checks']) == 12
    assert max(r['max_relative_error'] for r in numerical['real_template_noiseless_checks']) < 1e-4
    for path in OUT.glob('*.csv'): assert sha(path) == sha(EXP/'tables'/path.name)
    for name, old_count, new_count in [('codex_0007_SiPM_Energy_Calibration.ipynb', 55, 68), ('Research_Plan_20260921.ipynb', 83, 84)]:
        nb = nbformat.read(ROOT/'jupyter'/name, as_version=4); before = nbformat.read(REC/(name+'.before'), as_version=4)
        nbformat.validate(nb); assert len(nb.cells) == new_count and len(before.cells) == old_count
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
            assert len(pictures) == 4
            for data, filename in zip(pictures, report['figures']): assert hashlib.sha256(data).hexdigest() == sha(OUT/(filename+'.png'))
            numbers = [s for c in nb.cells if c.cell_type == 'markdown' for s in re.findall(r'\*\*图\s+(0007-\d+-\d+)：', c.source)]
            assert len(numbers) == len(set(numbers))
            assert all(f'0007-{i}-1' in numbers for i in [29, 30, 31, 32])
    result = dict(passed=True, scored_rows=6160, frames=1540, blind_frames=1540, oracle_frames=1540,
        inference_truth_separation_checked=True, same_frozen_blind_candidates_used_for_known_K=True,
        inherited_model_unchanged=True, numerical_tests=5, real_template_noiseless_cases=12,
        mother_components=110, validation_roots=280, old_notebook_cells_preserved=55,
        new_executed_cells=6, new_figures=4, new_ROOTs=0, default_changed=False,
        original_cost_gate_failure_preserved=True, independent_confirmation=False, committed=False, pushed=False)
    for p in [REC/'delivery-validation.json', EXP/'delivery-validation.json']: p.write_text(json.dumps(result, indent=2)+'\n')
    print('ORACLE_DELIVERY_VERIFIED', result, flush=True)


if __name__ == '__main__': main()
