"""Read-only verification of the saved migration controls and current build."""
from pathlib import Path
import json
import os
import re
import sys

import nbformat
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
os.environ['GAGG_TRANSPORT_ENV_ROOT'] = str(ROOT)
sys.path.insert(0, str(HERE / 'analysis-source'))
from pileup.analysis.alpha_transport_environment_study import load_root, read, sha, verify
from pileup.analysis.alpha_transport_environment import macro_semantics

result = verify()
frozen_files = 0
for name in ['output-freeze.json', 'main-replay-source-freeze.json',
             'main-smoke-source-freeze.json', 'report-source-freeze.json', 'report-output-freeze.json']:
    for path, digest in read(HERE / name).items():
        assert sha(ROOT / path) == digest, (name, path)
        frozen_files += 1
build = read(HERE / 'main-build.json')
assert build['ctest_passed'] and build['all_libraries_resolved']
assert sha(ROOT / build['binary']) == build['binary_sha256']
for path, digest in build['source_sha256'].items():
    assert sha(ROOT / path) == digest, path
assert read(HERE / 'build.json')['ctest_passed']

manifest = read(HERE / 'manifest.json')
replay = read(HERE / 'main-replay.json')
smoke = read(HERE / 'main-smoke.json')
assert not manifest['failed_attempts'] and len(manifest['runs']) == 16
assert replay['passed'] and len(replay['rows']) == 4 and smoke['passed']
all_runs = manifest['runs'] + replay['rows'] + [smoke]
assert len({row['root'] for row in all_runs}) == 21
starts = sorted(row['started_epoch'] for row in all_runs)
assert min(np.diff(starts)) >= 1.0
for row in all_runs:
    assert re.fullmatch(r'g4/data/\d{8}_\d{2}h\d{2}m\d{2}s\.root', row['root'])
    assert sha(ROOT / row['root']) == row['root_sha256']
    info, primary, photons, voltage, scalars = load_root(ROOT / row['root'])
    assert info['macro_contents'] == (ROOT / row['macro']).read_text()
    assert info['random_seed'] == row['seed'] and info['t_length_ps'] == 20000000
    assert np.isclose(primary['primary_energy_MeV'], row['source_energy_MeV'], rtol=2e-7, atol=1e-9)
    assert primary['edep_total_MeV'] == row['edep_MeV']
    assert primary['primary_time_ps'] == 5000000 and len(photons['time_ps']) > 0
    assert all(np.isfinite(array).all() for array in voltage.values())
    assert abs(scalars['charge_balance_error_C']) < 1e-18
    assert (ROOT / row['log']).is_file()

for row in manifest['runs']:
    assert sha(ROOT / row['reference_root']) == row['reference_sha256']
    old = read(HERE / 'reference-runinfo' / f"{row['id']:02d}.json")['run_info']
    current = read(HERE / 'runinfo' / f"{row['id']:02d}.json")
    assert macro_semantics(old['macro_contents']) == macro_semantics(current['macro_contents'])
    assert old['random_engine_state'] == current['random_engine_state']
    assert json.loads(old['sipm_config_json']) == json.loads(current['sipm_config_json'])
    for key in ['random_seed', 'sipm_avalanche_base_seed', 'sipm_readout_base_seed',
                't_length_ps', 'gagg_max_step_um', 'mylar_front_max_step_um', 'mylar_side_max_step_um']:
        assert old[key] == current[key], key

exact_arrays = 0
for row in replay['rows']:
    old = load_root(ROOT / row['reference_root'])
    current = load_root(ROOT / row['root'])
    assert old[0]['random_engine_state'] == current[0]['random_engine_state']
    assert old[1] == current[1] and old[4] == current[4]
    for column in (2, 3):
        for key in old[column]:
            assert np.array_equal(old[column][key], current[column][key]), (row['id'], key)
            exact_arrays += 1
assert exact_arrays == 28

notebook_counts = {}
for name, count, added in [('codex_0007_SiPM_Energy_Calibration.ipynb', 193, 6),
                           ('Research_Plan_20260921.ipynb', 104, 1), ('Run_Catalog.ipynb', 65, 1)]:
    before = read(ROOT / 'jupyter/run_records' / HERE.name / (name + '.before'))
    after = read(ROOT / 'jupyter' / name)
    assert before['cells'] == after['cells'][:count] and len(after['cells']) == count + added
    assert {k: v for k, v in before.items() if k != 'cells'} == {
        k: v for k, v in after.items() if k != 'cells'}
    nbformat.validate(nbformat.from_dict(after))
    for cell in after['cells'][count:]:
        source = ''.join(cell['source'])
        if cell['cell_type'] == 'markdown':
            assert not re.search(r'\]\([^)]*\.(?:md|ipynb|pdf)(?:[#?)])', source, re.I)
    notebook_counts[name] = len(after['cells'])
    if name.startswith('codex_'):
        code = [c for c in after['cells'][count:] if c['cell_type'] == 'code']
        assert len(code) == 2
        for cell in code:
            assert cell['execution_count'] is not None
            assert any('image/png' in output.get('data', {}) for output in cell['outputs'])
        text = '\n'.join(''.join(c['source']) for c in after['cells'] if c['cell_type'] == 'markdown')
        assert text.count('图 0007-120-1：') == 1 and text.count('图 0007-121-1：') == 1

status = read(HERE / 'environment-status.json')
assert status['native_root_count'] == 21 and status['independent_new_mothers'] == 0
assert not status['old_environment_exact_matching'] and not status['fresh_confirmation_run']
result.update(passed=True, extra_frozen_files=frozen_files, native_roots=21,
              preserved_reference_roots=16, current_exact_arrays=exact_arrays,
              minimum_launch_spacing_seconds=float(min(np.diff(starts))),
              main_source_files=len(build['source_sha256']), notebooks=notebook_counts)
print(json.dumps(result, ensure_ascii=False, indent=2))
