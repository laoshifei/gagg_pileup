"""Equivalent read-only ROOT extraction, overlapped with native production.

Only the first 200 donors are prefetched, leaving a native-production margin
before the frozen physical stage starts. No reconstruction outcome is read.
"""
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import json
import multiprocessing as mp
import os
import shutil
import sys
import tempfile
import time

FOLDER=Path(__file__).resolve().parent
PROJECT=FOLDER.parents[2]
OUTPUT=PROJECT/'jupyter/artifacts/alpha_pair_residual_20261023'
RECORD=PROJECT/'jupyter/run_records/alpha_pair_residual_20261023'
sys.path.insert(0,str(FOLDER/'inference-source'))


def extract(task):
    import numpy as np
    from pileup.analysis.sipm_response import SensorConfig,read_photon_hits,map_cells
    from pileup.analysis.alpha_pair_resolution_data import digest,load
    i,row=task;started=time.monotonic()
    path=OUTPUT/'photons'/f'{i:04}.npz'
    assert not path.exists(),'Do not overwrite a prepared photon cache'
    with tempfile.TemporaryDirectory(prefix='gagg-residual-prefill-') as temp:
        local=Path(temp)/Path(row['root']).name
        shutil.copyfile(PROJECT/row['root'],local)
        assert digest(local)==row['root_sha256']
        hits=read_photon_hits(local)
    config=load(OUTPUT/'native-inputs/metadata'/(Path(row['root']).stem+'.json'))['config']
    cells=map_cells(hits,SensorConfig(**config['sensor']))
    first=float(hits.time_ns[cells>=0].min())
    arrays=dict(time_ns=hits.time_ns-first,x_mm=hits.x_mm,y_mm=hits.y_mm,energy_eV=hits.energy_eV)
    temporary=path.with_suffix('.prefill-partial')
    with temporary.open('wb') as stream:
        np.savez_compressed(stream,**arrays);stream.flush();os.fsync(stream.fileno())
    with np.load(temporary) as f:
        assert set(f.files)==set(arrays)
        for k,v in arrays.items():np.testing.assert_array_equal(f[k],v)
    temporary.replace(path)
    return dict(index=i,root=row['root'],root_sha256=row['root_sha256'],
        cache_sha256=digest(path),seconds=time.monotonic()-started,passed=True)


def main():
    (OUTPUT/'photons').mkdir(exist_ok=True)
    log=RECORD/'photon-prefill.jsonl';done={}
    if log.exists():
        for line in log.read_text().splitlines():
            r=json.loads(line);done[r['index']]=r
    with ProcessPoolExecutor(max_workers=2,mp_context=mp.get_context('spawn')) as pool,log.open('a') as stream:
        while len(done)<200:
            try:manifest=json.loads((FOLDER/'manifest.json').read_text())
            except json.JSONDecodeError:
                time.sleep(1.);continue
            plans=sorted(manifest['planned_runs'],key=lambda r:r['seed'])
            index={p['seed']:i for i,p in enumerate(plans)}
            ready=[(index[r['seed']],r) for r in manifest['runs'] if r['status']=='completed'
                   and index[r['seed']]<200 and index[r['seed']] not in done]
            if not ready:
                time.sleep(5.);continue
            for future in as_completed([pool.submit(extract,r) for r in ready]):
                row=future.result();done[row['index']]=row
                stream.write(json.dumps(row)+'\n');stream.flush()
                if len(done)%20==0:print('PREFILL',len(done),'/200',flush=True)
    (RECORD/'photon-prefill-validation.json').write_text(json.dumps(dict(passed=True,
        files=len(done),workers=2,total_concurrent_native_and_prefill_processes_at_most=6,
        physics_changed=False,algorithm='Same read_photon_hits/map_cells/alignment as frozen donor_cache',
        root_hashes_verified=True,all_cached_arrays_readback_exact=True),indent=2)+'\n')
    print('PREFILL_COMPLETE',flush=True)


if __name__=='__main__':main()
