"""Run the predeclared independent test; never choose parameters from it."""
from datetime import datetime
from pathlib import Path
import json
import os
import subprocess
import sys
import time

FOLDER=Path(__file__).resolve().parent
PROJECT=FOLDER.parents[2]
RECORD=PROJECT/'jupyter/run_records'/FOLDER.name


def run(name,script,*arguments):
    started=time.monotonic()
    command=[sys.executable,str(FOLDER/script),*arguments]
    with (RECORD/(name+'-process.log')).open('a') as log:
        log.write('\nSTART '+datetime.now().astimezone().isoformat()+'\n');log.flush()
        process=subprocess.Popen(command,cwd=PROJECT,stdout=log,stderr=subprocess.STDOUT)
        code=process.wait()
    with (RECORD/'production-stages.jsonl').open('a') as stream:
        stream.write(json.dumps(dict(stage=name,command=command,returncode=code,
            seconds=time.monotonic()-started,finished=datetime.now().astimezone().isoformat()))+'\n')
    print('STAGE_COMPLETE',name,code,flush=True)
    if code:raise RuntimeError(name+' failed; inspect durable stage logs')


if __name__=='__main__':
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    assert (FOLDER/'freeze.json').exists()
    if not (FOLDER/'test/manifest.json').exists():run('test-prepare','run.py','create-test')
    manifest=json.loads((FOLDER/'test/manifest.json').read_text())
    if len(manifest['runs'])!=len(manifest['planned_runs']):run('test-native','run.py','native','--stage','test')
    if not (RECORD/'test/physical-validation.json').exists():run('test-physical','run.py','physical','--stage','test')
    run('test-model-inference','run.py','infer','--stage','test','--variant','model')
    run('test-observed-inference','run.py','infer','--stage','test','--variant','observed')
    run('test-report','report.py','test')
    run('test-figures','report.py','figures')
    print('INDEPENDENT_TEST_COMPLETE; scientific validation and notebook delivery follow',flush=True)
