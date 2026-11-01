"""Replay selected frozen physical controls into a new local directory."""
from pathlib import Path
import argparse,importlib,os,sys,json
from concurrent.futures import ProcessPoolExecutor
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
folder=Path(__file__).resolve().parent
is_single=folder.name=='alpha_pair_causal_controls_20261101'
os.environ['GAGG_CAUSAL_CONTROLS_ROOT' if is_single else 'GAGG_CAUSAL_SHARED_ROOT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'source'))
import numpy as np
s=importlib.import_module('pileup.analysis.'+('alpha_pair_causal_controls' if is_single else 'alpha_pair_causal_shared_study'))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True);parser.add_argument('--ids',nargs='+',type=int)
    parser.add_argument('--workers',type=int,default=4);args=parser.parse_args()
    assert 1<=args.workers<=4
    output=Path(args.output);output.mkdir(parents=True,exist_ok=False)
    s.compatibility();s.verify()
    tasks=s.read(folder/('single-plan.json' if is_single else 'plan.json'))
    if args.ids is not None:
        tasks=[t for t in tasks if t['id'] in args.ids]
        assert {t['id'] for t in tasks}==set(args.ids)
    original=s.ART/('single' if is_single else 'cases')
    inputs=[{k:v for k,v in t.items() if k!='energy_dep_MeV'} for t in tasks] if is_single else tasks
    checks=[]
    with ProcessPoolExecutor(max_workers=args.workers,initializer=s.initialize) as pool:
        for result,arrays in pool.map(s.single if is_single else s.one,inputs):
            name=f"{result['id']:04d}"
            with np.load(original/(name+'.npz')) as expected:
                exact=all(np.array_equal(arrays[k],expected[k]) for k in arrays)
            np.savez_compressed(output/(name+'.npz'),**arrays)
            (output/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n')
            checks.append(dict(id=result['id'],physical_arrays_exact=exact))
    report=dict(cases=len(tasks),checks=checks,passed=all(c['physical_arrays_exact'] for c in checks))
    (output/'replay.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
    assert report['passed']

if __name__=='__main__':main()
