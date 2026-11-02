"""Check the documented single-alpha entry point with the current installation."""
from pathlib import Path
import os
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
os.environ['GAGG_TRANSPORT_ENV_ROOT'] = str(ROOT)
sys.path.insert(0, str(HERE / 'analysis-source'))
from pileup.analysis.alpha_transport_environment_study import load_root, read, save, sha

REC = ROOT / 'jupyter/run_records' / HERE.name
assert not (HERE / 'main-smoke.json').exists(), 'Keep the saved smoke result; do not overwrite it'
macro = ROOT / 'g4/macros/run_single_alpha_sipm.mac'
macro_text = macro.read_text()
copy = HERE / 'macros/main_single_alpha.mac'
assert not copy.exists()
copy.write_text(macro_text)
script = '''set -e
source /opt/root/root63410/bin/thisroot.sh
source /opt/geant4/geant4-v1132/bin/geant4.sh
exec ./gagg ../macros/run_single_alpha_sipm.mac
'''
plan = dict(command=['bash', '--noprofile', '--norc', '-c', script], cwd='g4/build',
            source_energy_MeV=5.486, seed=261010201, native_sipm_seed=261010301,
            native_readout_seed=261010401, window_ns=20000, primary_time_ns=5000,
            macro=str(copy.relative_to(ROOT)), macro_sha256=sha(copy),
            binary_sha256=sha(ROOT / 'g4/build/gagg'), new_independent_mothers=0,
            purpose='Operational check of the existing README single-alpha macro; not a new physics sample')
assert plan['binary_sha256'] == read(HERE / 'main-build.json')['binary_sha256']
save(HERE / 'main-smoke-plan.json', plan)
save(HERE / 'main-smoke-source-freeze.json', {
    str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__), copy, HERE / 'main-smoke-plan.json']})
before = set((ROOT / 'g4/data').glob('*.root'))
start = time.time()
log = REC / 'logs/main_single_alpha.log'
with log.open('w') as output:
    result = subprocess.run(plan['command'], cwd=ROOT / plan['cwd'], stdout=output,
                            stderr=subprocess.STDOUT)
new = set((ROOT / 'g4/data').glob('*.root')) - before
save(REC / 'main-smoke-attempt.json', dict(returncode=result.returncode,
     started_epoch=start, seconds=time.time()-start, new_roots=[str(p.relative_to(ROOT)) for p in new],
     log=str(log.relative_to(ROOT))))
assert result.returncode == 0 and len(new) == 1
root = new.pop()
info, primary, photons, voltage, scalars = load_root(root)
assert macro.read_text() == macro_text == info['macro_contents']
assert info['random_seed'] == plan['seed'] and info['t_length_ps'] == 20000000
assert info['sipm_avalanche_base_seed'] == plan['native_sipm_seed']
assert info['sipm_readout_base_seed'] == plan['native_readout_seed']
assert abs(primary['primary_energy_MeV'] - 5.486) < 1e-6
assert primary['primary_time_ps'] == 5000000 and primary['edep_total_MeV'] > 0
assert len(voltage['analog_mV']) > 0 and len(photons['time_ps']) > 0
assert abs(scalars['charge_balance_error_C']) < 1e-18
save(HERE / 'runinfo/main_single_alpha.json', info)
save(HERE / 'main-smoke.json', dict(passed=True, root=str(root.relative_to(ROOT)),
     root_sha256=sha(root), root_bytes=root.stat().st_size, macro=plan['macro'],
     seed=plan['seed'], source_energy_MeV=plan['source_energy_MeV'],
     edep_MeV=primary['edep_total_MeV'], photons=len(photons['time_ps']),
     samples=len(voltage['analog_mV']), started_epoch=start, seconds=time.time()-start,
     log=str(log.relative_to(ROOT)), native_roots=1, new_independent_mothers=0))
print('Documented single-alpha command passed:', root.relative_to(ROOT), flush=True)
