"""Resume only missing frozen pair jobs; retain every completed fit unchanged."""
import fcntl
import os
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
import time
from . import alpha_current_pair_energy_study as study


def one(task):
    destination=study.ART/'fits'/f"{task['id']:04d}.json"
    assert not destination.exists()
    row=study.fit_one((task,'resume-completed'))
    staged=study.ART/'resume-completed'/destination.name
    assert not destination.exists()
    os.replace(staged,destination)
    return row


def main():
    with (study.REC/'resume.lock').open('w') as lock:
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        study.verify();assert study.read(study.EXP/'single-control.json')['passed']
        assert not (study.EXP/'pair-fits.json').exists()
        tasks=[t for t in study.read(study.EXP/'fit-plan.json') if t['kind']=='pair']
        pending=[];retained={}
        for task in tasks:
            path=study.ART/'fits'/f"{task['id']:04d}.json"
            if not path.exists():pending.append(task);continue
            row=study.read(path)
            assert all(row[k]==v for k,v in task.items()),task
            assert row['truth_guard'] and not row['true_E_input']
            retained[str(path.relative_to(study.ROOT))]=study.sha(path)
        started=time.monotonic()
        status=dict(state='running',pid=os.getpid(),started_utc=datetime.now(timezone.utc).isoformat(),
                    retained=len(retained),remaining=len(pending),completed=len(retained),workers=8)
        study.save(study.REC/'resume-status.json',status)
        study.save(study.EXP/'resume-retained-freeze.json',retained)
        print('Resume:',len(retained),'retained;',len(pending),'missing frozen jobs',flush=True)
        with ProcessPoolExecutor(max_workers=8,initializer=study.initialize,initargs=(True,)) as pool:
            futures=[pool.submit(one,t) for t in pending]
            for future in as_completed(futures):
                future.result();status['completed']+=1
                study.save(study.REC/'resume-status.json',status)
                if status['completed']%8==0:print('pair fits',status['completed'],'/112',flush=True)
        for path,digest in retained.items():assert study.sha(study.ROOT/path)==digest,path
        rows=[]
        for task in tasks:
            fit=study.read(study.ART/'fits'/f"{task['id']:04d}.json")
            rows.append(dict(id=task['id'],kind='pair',seconds=fit['seconds'],fits=1))
        study.save(study.EXP/'pair-fits.json',dict(rows=rows,fits=len(rows),seconds=time.monotonic()-started,
            timing_scope='Resume wall time only; original interrupted wall time unmeasured',
            sum_saved_fit_wall_seconds=sum(r['seconds'] for r in rows),retained_outputs=len(retained),resumed_fits=len(pending)))
        study.save(study.EXP/'pair-fit-freeze.json',{str((study.ART/'fits'/f"{t['id']:04d}.json").relative_to(study.ROOT)):
             study.sha(study.ART/'fits'/f"{t['id']:04d}.json") for t in tasks})
        status.update(state='completed',seconds=time.monotonic()-started)
        study.save(study.REC/'resume-status.json',status)
        print('All112pair fits completed; retained outputs unchanged',flush=True)


if __name__=='__main__':main()
