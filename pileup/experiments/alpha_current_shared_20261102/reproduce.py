"""Copy frozen inputs to a new project root for a fresh numerical rerun.

Existing roots are refused. No production, cleanup, Git operation or Notebook
write is performed. The copied run.py uses the original source snapshot.
"""
from pathlib import Path
import argparse
from datetime import datetime,timezone
import hashlib
import json
import shutil


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target-root',type=Path,required=True)
    args=parser.parse_args()
    here=Path(__file__).resolve().parent;root=here.parents[2]
    target=args.target_root.resolve()
    if target.exists():parser.error('Target root must not exist; existing outputs are never overwritten')
    freeze=json.loads((here/'design-freeze.json').read_text())
    relative=here.relative_to(root)
    for name,expected in freeze.items():
        path=Path(name)
        if path.is_absolute() or '..' in path.parts or not path.is_relative_to(relative):
            raise ValueError(f'Unexpected frozen path: {name}')
        assert digest(root/path)==expected,name
    build=json.loads((here/'build.json').read_text())
    for name,key in [('gagg','native_binary_sha256'),('age.so','age_binary_sha256')]:
        assert digest(here/'bin'/name)==build[key],name
    target.mkdir(parents=True,exist_ok=False)
    for name in freeze:
        dest=target/name;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(root/name,dest)
    experiment=target/relative
    shutil.copyfile(here/'design-freeze.json',experiment/'design-freeze.json')
    (experiment/'bin').mkdir();(experiment/'build').mkdir();(experiment/'macros').mkdir()
    (experiment/'original-runinfo').mkdir()
    for name in ('gagg','age.so'):shutil.copyfile(here/'bin'/name,experiment/'bin'/name)
    (experiment/'bin/gagg').chmod(0o755)
    for part in ('g4/data',*[f'jupyter/artifacts/{here.name}/{p}' for p in
                 ('photons','voltage/single','truth/single','voltage/pair','truth/pair','fits','predictions','figures')],
                 f'jupyter/run_records/{here.name}/logs',f'jupyter/run_records/{here.name}/attempts'):
        (target/part).mkdir(parents=True,exist_ok=True)
    planned=json.loads((experiment/'native-plan.json').read_text())
    (experiment/'manifest.json').write_text(json.dumps(dict(planned_runs=planned,runs=[],failed_attempts=[]),indent=2)+'\n')
    record=dict(created_utc=datetime.now(timezone.utc).isoformat(),original_snapshot=str(here),
                original_design_sha256=digest(here/'design-freeze.json'),fresh_root=str(target),
                production_submitted=False,old_outputs_copied=False,same_seeds=True,
                scope='Numerical stages through replay; reporting requires the local Notebook context')
    (experiment/'reproduction.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))
    print('From the new root, run the copied run.py stages in protocol order, respecting each gate.')


if __name__=='__main__':main()
