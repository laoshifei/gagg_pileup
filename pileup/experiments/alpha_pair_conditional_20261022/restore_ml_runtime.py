"""Restore the previously used CPU wheels offline, verifying saved hashes."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys

PROJECT=Path(__file__).resolve().parents[3]
RECORD=PROJECT/'jupyter/run_records/alpha_pair_conditional_20261022'


def main():
    previous=PROJECT/'jupyter/run_records/alpha_correlated_global_20261020/ml-dependencies.json'
    manifest=json.loads(previous.read_text())
    wheels=Path('/tmp/gagg-conditional-wheels-20261022');wheels.mkdir(exist_ok=True)
    copied={}
    for item in manifest['wheels']:
        source=Path(item['cache']);name=Path(item['wheel']).name
        assert hashlib.sha256(source.read_bytes()).hexdigest()==item['sha256'],name
        target=wheels/name
        if not target.exists():shutil.copyfile(source,target)
        assert hashlib.sha256(target.read_bytes()).hexdigest()==item['sha256'],name
        copied[name]=item['sha256']
    target=Path('/tmp/gagg-correlated-ml-deps-20261020')
    command=[sys.executable,'-m','pip','install','--no-index','--no-cache-dir','--no-compile',
             '--find-links',str(wheels),'--target',str(target),'torch==2.14.1+cpu']
    with (RECORD/'ml-restore.log').open('w') as stream:
        result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
    record=dict(command=command,returncode=result.returncode,network=False,wheels=copied,
        previous_manifest=str(previous.relative_to(PROJECT)),target=str(target),
        reason='Old optional /tmp dependency directory was absent; data generation completed before the import failure')
    (RECORD/'ml-restore.json').write_text(json.dumps(record,indent=2)+'\n')
    assert result.returncode==0,'Inspect ml-restore.log'
    print('OFFLINE_RUNTIME_RESTORED',target,flush=True)


if __name__=='__main__':main()
