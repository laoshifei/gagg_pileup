"""Copy frozen numerical inputs to a new root without overwriting results."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

here=Path(__file__).resolve().parent
root=here.parents[2]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target-root',required=True,type=Path)
    target=parser.parse_args().target_root.resolve()
    if target.exists():parser.error('Target must not exist')
    files={}
    for name in ('design-freeze.json','input-freeze.json'):
        files.update(json.loads((here/name).read_text()))
    for name,digest in files.items():
        path=Path(name)
        assert not path.is_absolute() and '..' not in path.parts
        assert hashlib.sha256((root/path).read_bytes()).hexdigest()==digest
    target.mkdir(parents=True)
    for name in files:
        p=target/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(root/name,p)
    experiment=target/here.relative_to(root)
    for name in ('design-freeze.json','input-freeze.json'):shutil.copyfile(here/name,experiment/name)
    (experiment/'bin').mkdir();shutil.copyfile(here/'bin/age.so',experiment/'bin/age.so')
    for part in ('diagnostics','replay','figures'):
        (target/'jupyter/artifacts'/here.name/part).mkdir(parents=True,exist_ok=True)
    (target/'jupyter/run_records'/here.name).mkdir(parents=True,exist_ok=True)
    print(json.dumps(dict(target_root=str(target),copied_files=len(files),
        scope='run.py verify,pilot,run,score,replay;no original ROOT needed for numerical diagnosis;report requires existing Notebook context')))


if __name__=='__main__':main()
