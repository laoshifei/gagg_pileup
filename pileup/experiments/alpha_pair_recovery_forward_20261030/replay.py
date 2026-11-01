"""Replay the frozen native-voltage fits into a new output directory."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import os
import sys
import time

folder=Path(__file__).resolve().parent
os.environ['GAGG_FORWARD_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'source'))
import numpy as np
from pileup.analysis import alpha_pair_recovery_forward_study as s


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True)
    parser.add_argument('--workers',type=int,default=6);parser.add_argument('--ids',nargs='+',type=int)
    args=parser.parse_args();output=Path(args.output).resolve()
    if output.exists():raise FileExistsError('Choose a new output directory')
    for freeze in ('source-freeze.json','model-freeze.json','input-freeze.json'):
        for name,expected in s.read(folder/freeze).items():
            p=s.ROOT/name
            if (folder/'source') in p.parents or p.parent==folder or p==s.ART/'inputs/observations.npz':
                assert s.sha(p)==expected,name
    s.build()
    with np.load(s.ART/'inputs/observations.npz') as f:
        items=[(i,v,t[:k]) for i,v,t,k in zip(f['ids'],f['voltage_mV'],f['known_times_ns'],f['known_counts'])
               if args.ids is None or int(i) in args.ids]
    if not items or (args.ids is not None and {int(x[0]) for x in items}!=set(args.ids)):
        raise ValueError('Requested IDs must be in this frozen diagnostic')
    # Comparison happens in the parent, never passed to worker inference.
    reference={int(x[0]):s.read(s.ART/'cases'/f'{int(x[0]):04d}.json')['fits'] for x in items}
    output.mkdir(parents=True);start=time.monotonic();equal=True
    with ProcessPoolExecutor(max_workers=args.workers,initializer=s.initialize,initargs=(True,)) as pool:
        for result in pool.map(s.job,items):
            s.write(output/f"{result['id']:04d}.json",result)
            equal &= result['fits']==reference[result['id']]
    result=dict(frames=len(items),fits_exactly_reproduced=bool(equal),seconds=time.monotonic()-start,workers=args.workers)
    s.write(output/'replay.json',result);print(result,flush=True)
    assert equal


if __name__=='__main__':main()
