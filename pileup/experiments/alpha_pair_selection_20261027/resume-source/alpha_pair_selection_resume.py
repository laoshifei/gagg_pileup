"""Resume unchanged frozen fits while preserving every completed checkpoint row."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import multiprocessing as mp
import time

from .alpha_pair_selection_study import (EXPERIMENT, OUTPUT, RECORD, check_freeze,
    initialize, infer_one, now)
from .alpha_pair_information_study import load, save, digest, rows, serializable
from .alpha_pair_stability_study import copy_chunks


def resume(split, workers):
    check_freeze();assert 1<=workers<=8
    for n,h in load(EXPERIMENT/'resume-freeze.json')['files'].items():assert digest(EXPERIMENT/n)==h
    check=load(RECORD/(split+'-generation.json'));assert check['passed']
    for n,h in check['files'].items():assert digest(OUTPUT/split/n)==h
    target=OUTPUT/split/'predictions.jsonl';previous=digest(target)
    prior=rows(target);done={r['index'] for r in prior};assert len(prior)==len(done)
    expected=check['frames'];assert done<=set(range(expected))
    checkpoint=OUTPUT/split/'predictions.interrupted-01.jsonl'
    assert not checkpoint.exists();copy_chunks(target,checkpoint);assert digest(checkpoint)==previous
    save(RECORD/(split+'-interruption.json'),dict(created=now(),preserved_rows=len(prior),checkpoint_sha256=previous,
        checkpoint_path=str(checkpoint),reason='Calibration stopped advancing after200 rows; last write14:13:56UTC,confirmed unchanged at14:33:43UTC. Cause undetermined. Interrupted owned PTY withCtrl-C; exit1. Preserve all rows and finish unchanged tasks with2 workers.',
        changes='Only scheduling/concurrency/resumption; frozen fitter,model,observations,seeds and scientific protocol unchanged. No completed row recalculated or removed.'))
    pending=[i for i in range(expected) if i not in done];result=list(prior)
    resumed=OUTPUT/split/'predictions.resumed-01.jsonl';started=time.monotonic()
    with resumed.open('x') as stream,ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn'),initializer=initialize,initargs=(split,)) as pool:
        for r in prior:stream.write(json.dumps(r,default=serializable,allow_nan=False)+'\n')
        stream.flush()
        futures=[pool.submit(infer_one,i) for i in pending]
        for f in as_completed(futures):
            r=f.result();result.append(r);stream.write(json.dumps(r,default=serializable,allow_nan=False)+'\n');stream.flush()
            print('SELECTION_RESUME',split,len(result),expected,r['index'],round(r['seconds'],2),flush=True)
    assert len(result)==expected and {r['index'] for r in result}==set(range(expected))
    assert rows(resumed)[:len(prior)]==prior and digest(target)==previous
    resumed.replace(target)
    save(RECORD/(split+'-inference.json'),dict(passed=True,created=now(),frames=expected,sha256=digest(target),
        seconds=time.monotonic()-started,seconds_scope='Resumption only; interrupted run not timed completely',workers=workers,
        initial_workers=6,inherited_pilot=4 if split=='calibration' else 0,resumed_rows=len(prior),
        checkpoint_sha256=previous,truth_inputs=[],label_access_blocked=True,
        max_cost_replay_error=max(r['cost_replay_error'] for r in result)))
    print('SELECTION_RESUME_COMPLETE',split,expected,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('split',choices=('calibration','test'));p.add_argument('--workers',type=int,default=2)
    args=p.parse_args();resume(args.split,args.workers)
