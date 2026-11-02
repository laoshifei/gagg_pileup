"""Run the frozen analysis under a new output label and compare all summaries."""
from pathlib import Path
import argparse
import json
import subprocess
import sys


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--label',required=True);args=p.parse_args()
    if not args.label.replace('_','').isalnum() or args.label=='primary':p.error('Use a new alphanumeric/underscore label')
    exp=Path(__file__).resolve().parent;root=exp.parents[2]
    assert not (exp/'results'/args.label).exists()
    subprocess.run([sys.executable,'-B',str(exp/'run.py'),'analyze','--label',args.label],cwd=root,check=True)
    checks={}
    for name in ('summary.json','partitions.json','bootstrap.json'):
        a=json.loads((exp/'results/primary'/name).read_text());b=json.loads((exp/'results'/args.label/name).read_text())
        checks[name]=a==b
    assert all(checks.values()),checks
    (exp/'results'/args.label/'reproduction-equality.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps(checks,indent=2),flush=True)


if __name__=='__main__':main()
