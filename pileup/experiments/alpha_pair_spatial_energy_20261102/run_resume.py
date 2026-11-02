"""Resume missing frozen pair-fit IDs after an interrupted process.

No existing result is rewritten. This changes job scheduling only; all fitting
is delegated to the original frozen source and protocol.
"""
from pathlib import Path
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl

FOLDER=Path(__file__).resolve().parent
ROOT=FOLDER.parents[2]
os.environ['GAGG_SPATIAL_ENERGY_ROOT']=str(ROOT)
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
sys.path.insert(0,str(FOLDER/'source'))
from pileup.analysis import alpha_pair_spatial_energy_study as study


def main():
    for path,digest in study.read(FOLDER/'resume-freeze.json').items():
        assert study.sha(ROOT/path)==digest,path
    study.verify(True)
    assert study.read(study.REC/'numerical.json')['passed']
    assert study.read(study.EXP/'single-gate.json')['passed']
    assert study.read(study.REC/'pilot.json')['passed']
    lock=(study.REC/'resume.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    tasks=[t for t in study.read(study.EXP/'fit-plan.json') if t['kind']=='pair']
    folder=study.ART/'pairs'
    preserved={}
    pending=[]
    for task in tasks:
        path=folder/f"{task['id']:04d}.json"
        if not path.exists():
            pending.append(task)
            continue
        row=study.read(path)
        assert all(row[k]==v for k,v in task.items()),path
        assert row['read_inputs']==[task['waveform']] and not row['true_energy_input']
        assert len(row['energy_MeV'])==2 and all(0<e<100 for e in row['energy_MeV'])
        preserved[str(path.relative_to(ROOT))]=study.sha(path)
    original=study.read(study.REC/'interrupted-first-run.json')['completed_files']
    for path,digest in original.items():assert study.sha(ROOT/path)==digest,path
    stamp=datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    record=study.REC/f'resume-attempt-{stamp}.json'
    assert not record.exists()
    study.save(record,dict(started_utc=datetime.now(timezone.utc).isoformat(),
        pending_ids=[t['id'] for t in pending],preserved_files=preserved,
        runner_sha256=study.sha(Path(__file__)),workers=study.WORKERS))
    started=time.monotonic();finished=len(preserved)
    print('Resume',finished,'existing;',len(pending),'missing frozen tasks',flush=True)
    with ProcessPoolExecutor(max_workers=study.WORKERS,initializer=study.inference_initialize) as pool:
        futures=[pool.submit(study.fit_one,t) for t in pending]
        for future in as_completed(futures):
            row=future.result()
            path=folder/f"{row['id']:04d}.json"
            with path.open('x') as f:
                json.dump(row,f,indent=2,ensure_ascii=False);f.write('\n')
            finished+=1
            progress=dict(completed=finished,total=len(tasks),resume_seconds=time.monotonic()-started,
                          latest_id=row['id'],record=str(record.relative_to(ROOT)))
            study.save(study.REC/'resume-progress.json',progress)
            print('pairs',finished,'/',len(tasks),'resume elapsed',round(progress['resume_seconds'],1),flush=True)
    for path,digest in preserved.items():assert study.sha(ROOT/path)==digest,path
    assert len(list(folder.glob('*.json')))==len(tasks)
    study.save(study.REC/'pairs-run.json',dict(tasks=len(tasks),seconds=time.monotonic()-started,
        seconds_scope='Successful resume wall time only; interrupted attempt and discarded in-flight calculations excluded.',
        workers=study.WORKERS,preserved=len(preserved),computed=len(pending),
        original_completed_files_unchanged=True,resume_record=str(record.relative_to(ROOT))))
    print('All pair fits complete; previous results unchanged',flush=True)


if __name__=='__main__':main()
