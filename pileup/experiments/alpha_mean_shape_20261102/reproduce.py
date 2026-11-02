"""Repeat the frozen stages in a fresh local workspace; never overwrite inputs."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys


def read(p):return json.loads(p.read_text())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--label',required=True)
    parser.add_argument('--prepare-only',action='store_true',help='Validate the isolated input mirror without new digitizations')
    args=parser.parse_args()
    if not re.fullmatch(r'[a-zA-Z0-9_-]+',args.label):parser.error('Use letters,digits,underscore or dash')
    original=Path(__file__).resolve().parent;root=original.parents[2]
    workspace=root/'jupyter/run_records'/original.name/'reproductions'/args.label/'workspace'
    workspace.mkdir(parents=True,exist_ok=False)
    relative=original.relative_to(root);experiment=workspace/relative
    design=read(original/'design-freeze.json');inputs=read(original/'input-freeze.json')
    # Small frozen source/configuration files are copied; immutable historical
    # data use individual symlinks. No source output directory is linked.
    for name,digest in design.items():
        source=root/name;target=workspace/name;target.parent.mkdir(parents=True,exist_ok=True)
        assert hashlib.sha256(source.read_bytes()).hexdigest()==digest,name
        shutil.copyfile(source,target)
    for name,digest in inputs.items():
        source=root/name;target=workspace/name;target.parent.mkdir(parents=True,exist_ok=True)
        if not target.exists():target.symlink_to(source)
    for name in ('design-freeze.json','input-freeze.json'):
        shutil.copyfile(original/name,experiment/name)
    (experiment/'bin').mkdir(exist_ok=True)
    if (original/'bin/age.so').exists():shutil.copyfile(original/'bin/age.so',experiment/'bin/age.so')
    else:subprocess.run(read(original/'build.json')['command'],cwd=workspace,check=True)
    for folder in ('photons','voltage/train','voltage/single','voltage/pair','truth/train','truth/single',
                   'truth/pair','fits','figures','pilot','replay'):
        (workspace/'jupyter/artifacts'/original.name/folder).mkdir(parents=True,exist_ok=True)
    (workspace/'jupyter/run_records'/original.name).mkdir(parents=True,exist_ok=True)
    runner=[sys.executable,'-B',str(experiment/'run.py')]
    def run(stage):subprocess.run(runner+[stage],cwd=workspace,check=True)
    run('verify')
    marker=dict(workspace=str(workspace.relative_to(root)),input_mirror_verified=True,
                prepare_only=args.prepare_only,original_results_untouched=True)
    if not args.prepare_only:
        for stage in ('pilot','train','calibrate','numerical','single','fit','single_report'):
            run(stage)
        if read(experiment/'single-gate.json')['passed']:
            for stage in ('pair_pilot','pair','pair_report'):run(stage)
        run('replay')
        names=['calibration.json','numerical.json','single-result.json','single-gate.json']
        if (original/'pair-result.json').exists():names+=['pair-result.json','pair-gate.json']
        equal={name:read(original/name)==read(experiment/name) for name in names}
        marker['summary_exact_equality']=equal
        assert all(equal.values()),equal
    (workspace.parent/'reproduction-status.json').write_text(json.dumps(marker,indent=2)+'\n')
    print(json.dumps(marker,indent=2),flush=True)


if __name__=='__main__':main()
