"""Replay frozen inference into a new directory, preserving original results."""
from pathlib import Path
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import os
import sys
import time

folder = Path(__file__).resolve().parent
os.environ['GAGG_DECONV_PROJECT'] = str(folder.parents[2])
sys.path.insert(0, str(folder/'source'))

import numpy as np
from pileup.analysis.alpha_pair_deconvolution_study import (
    ARMS, ART, ROOT, initialize, inference_job, read_json, sha, write_json)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--ids', nargs='+', type=int)
    args = parser.parse_args()
    output = Path(args.output).resolve()
    if output.exists():
        raise FileExistsError('Choose a new output directory')
    for freeze in ('source-freeze.json', 'model-freeze.json', 'input-freeze.json'):
        for name, expected in read_json(folder/freeze).items():
            # Replaying does not depend on mutable analysis sources or the
            # old input producers; check frozen copies/current inference inputs.
            p = ROOT/name
            if (folder/'source') in p.parents or p.parent == folder or p == ART/'inputs/observations.npz':
                assert sha(p) == expected, name
    with np.load(ART/'inputs/observations.npz') as f:
        items = [(i, v, t[:k]) for i,v,t,k in zip(f['ids'],f['voltage_mV'],f['known_times_ns'],f['known_counts'])
                 if args.ids is None or int(i) in args.ids]
    if not items or (args.ids is not None and {int(x[0]) for x in items} != set(args.ids)):
        raise ValueError('Requested IDs must be in the frozen pilot')
    # References are checked by the parent after inference; never passed into
    # inference_job. Model initializer installs the data-read guard.
    reference = {int(x[0]): read_json(ART/'cases'/f'{int(x[0]):04d}.json')['fits'] for x in items}
    output.mkdir(parents=True)
    start = time.monotonic()
    def save_results(results):
        arrays = {arm: [] for arm in ARMS}; equal = True
        for result, transformed in results:
            write_json(output/f"{result['id']:04d}.json", result)
            equal &= result['fits'] == reference[result['id']]
            for arm in ARMS: arrays[arm].append(transformed[arm])
        np.savez_compressed(output/'transformed.npz', **{a:np.stack(x) for a,x in arrays.items()})
        record = dict(frames=len(items), ids=[int(x[0]) for x in items], fits_exactly_reproduced=bool(equal),
                      seconds=time.monotonic()-start, workers=args.workers, source='Frozen source/ and models; no Git checkout')
        write_json(output/'replay.json', record)
        print(json.dumps(record, indent=2), flush=True)
        assert equal
    if args.workers == 1:
        # Separate child keeps even parent-held references out of worker memory.
        with ProcessPoolExecutor(max_workers=1, initializer=initialize, initargs=(True,)) as pool:
            save_results(pool.map(inference_job, items))
    else:
        with ProcessPoolExecutor(max_workers=args.workers, initializer=initialize, initargs=(True,)) as pool:
            save_results(pool.map(inference_job, items))


if __name__ == '__main__':
    main()
