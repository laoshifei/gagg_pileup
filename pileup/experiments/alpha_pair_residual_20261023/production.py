"""Resume the frozen independent campaign; stop on any failed stage."""
from datetime import datetime
from pathlib import Path
import json
import os
import subprocess
import sys
import time

FOLDER = Path(__file__).resolve().parent
PROJECT = FOLDER.parents[2]
RECORD = PROJECT / 'jupyter/run_records/alpha_pair_residual_20261023'
OUTPUT = PROJECT / 'jupyter/artifacts/alpha_pair_residual_20261023'


def run_stage(name, script, *arguments):
    started = time.monotonic()
    command = [sys.executable, str(FOLDER / script), *arguments]
    with (RECORD / (name + '-process.log')).open('a') as log:
        log.write('\nSTART ' + datetime.now().astimezone().isoformat() + '\n')
        log.flush()
        process = subprocess.Popen(command, cwd=PROJECT, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            log.write(line)
            log.flush()
            print(name + ': ' + line.rstrip(), flush=True)
        code = process.wait()
    with (RECORD / 'production-stages.jsonl').open('a') as log:
        log.write(json.dumps(dict(stage=name, command=command, returncode=code,
            seconds=time.monotonic()-started, finished=datetime.now().astimezone().isoformat()))+'\n')
    if code:
        raise RuntimeError(f'{name} failed with code {code}; completed outputs are preserved')


def main():
    assert (FOLDER / 'freeze.json').exists(), 'Freeze selection before generating any new test parent'
    for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[key] = '1'
    manifest = json.loads((FOLDER / 'manifest.json').read_text())
    if len(manifest['runs']) != len(manifest['planned_runs']):
        # Equivalent photon extraction overlaps ROOT I/O. It never reads fits
        # or labels for algorithm selection and stays below the eight-core cap.
        with (RECORD / 'prefill-process.log').open('a') as log:
            prefill = subprocess.Popen([sys.executable, str(FOLDER / 'prefill_photons.py')],
                                       cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT)
            try:
                run_stage('native', 'run.py', 'native', '--workers', '4')
                if prefill.wait():
                    raise RuntimeError('Photon prefill failed; inspect prefill-process.log')
            finally:
                if prefill.poll() is None:
                    prefill.terminate()
                    prefill.wait()
    if not (RECORD / 'physical-validation.json').exists():
        run_stage('physical', 'run.py', 'physical', '--workers', '4')
    if not (RECORD / 'test-baseline-inference.json').exists():
        run_stage('baselines', 'run.py', 'seeds')
    run_stage('test-inference', 'run.py', 'infer', '--dataset', 'test', '--workers', '4')
    run_stage('test-report', 'report.py', 'test')
    run_stage('figures', 'report.py', 'figures')
    run_stage('science-validation', 'validate.py')
    print('PRODUCTION_AND_VALIDATION_COMPLETE; notebook delivery is a separate reviewed step', flush=True)


if __name__ == '__main__':
    main()
