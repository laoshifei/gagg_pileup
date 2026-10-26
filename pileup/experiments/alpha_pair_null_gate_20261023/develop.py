"""Resume weak-pulse development with durable stage records."""
from datetime import datetime
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

FOLDER = Path(__file__).resolve().parent
PROJECT = FOLDER.parents[2]
RECORD = PROJECT/'jupyter/run_records'/FOLDER.name


def run(name, *arguments):
    started = time.monotonic()
    source = PROJECT/'pileup/analysis/alpha_pair_null_study.py'
    snapshot = RECORD/(name+'-source.py')
    if not snapshot.exists(): shutil.copyfile(source, snapshot)
    command = [sys.executable, str(FOLDER/'run.py'), *arguments]
    with (RECORD/(name+'-process.log')).open('a') as log:
        log.write('\nSTART '+datetime.now().astimezone().isoformat()+'\n'); log.flush()
        process = subprocess.Popen(command, cwd=PROJECT, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            log.write(line); log.flush(); print(name+': '+line.rstrip(), flush=True)
        code = process.wait()
    with (RECORD/'development-stages.jsonl').open('a') as stream:
        stream.write(json.dumps(dict(stage=name, command=command, returncode=code,
            seconds=time.monotonic()-started, source_sha256=hashlib.sha256(snapshot.read_bytes()).hexdigest(),
            finished=datetime.now().astimezone().isoformat()))+'\n')
    if code: raise RuntimeError(name+' failed; preserve outputs and inspect the log')


if __name__ == '__main__':
    for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[key] = '1'
    manifest = json.loads((FOLDER/'development/manifest.json').read_text())
    if len(manifest['runs']) != len(manifest['planned_runs']):
        run('development-native', 'native', '--stage', 'development')
    if not (RECORD/'development/physical-validation.json').exists():
        run('development-physical', 'physical', '--stage', 'development')
    run('development-inference', 'infer', '--stage', 'development')
    run('development-selection', 'select')
