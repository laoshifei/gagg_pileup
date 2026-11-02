"""Prepare an empty project root for a numerical rerun from frozen files."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target-root',required=True,type=Path)
    target=parser.parse_args().target_root.resolve()
    if target.exists():parser.error('Target root must not exist; existing outputs are never overwritten')
    here=Path(__file__).resolve().parent;root=here.parents[2]
    names=('design-freeze.json','input-freeze.json','upstream-input-freeze.json','resume-source-freeze.json')
    files={}
    for name in names:files.update(json.loads((here/name).read_text()))
    for name,expected in files.items():
        p=Path(name)
        if p.is_absolute() or '..' in p.parts:raise ValueError(name)
        assert sha(root/p)==expected,name
    build=json.loads((here/'build.json').read_text())
    assert sha(here/'bin/age.so')==build['age_binary_sha256']
    target.mkdir(parents=True,exist_ok=False)
    for name in files:
        dest=target/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(root/name,dest)
    experiment=target/here.relative_to(root)
    for name in names:shutil.copyfile(here/name,experiment/name)
    (experiment/'bin').mkdir();shutil.copyfile(here/'bin/age.so',experiment/'bin/age.so')
    for part in ('fits','numerical','replay','figures'):
        (target/'jupyter/artifacts'/here.name/part).mkdir(parents=True,exist_ok=True)
    (target/'jupyter/run_records'/here.name).mkdir(parents=True,exist_ok=True)
    record=dict(passed=True,target_root=str(target),frozen_files_copied=len(files),
                new_simulations=0,fits_submitted=0,scope='run.py numerical,single,pair,score,replay;report requires original local Notebook context')
    (experiment/'reproduction.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
