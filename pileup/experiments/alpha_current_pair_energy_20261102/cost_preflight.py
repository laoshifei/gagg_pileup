from pathlib import Path
import os,sys,json,time,hashlib
from concurrent.futures import ProcessPoolExecutor
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
ROOT=Path(__file__).resolve().parents[3]
PARENT=ROOT/'pileup/experiments/alpha_current_shared_20261102'
ART=ROOT/'jupyter/artifacts/alpha_current_shared_20261102'
REC=Path(__file__).resolve().parent
sys.path.insert(0,str(PARENT/'source'))
import numpy as np
from pileup.analysis import alpha_current_shared_study as shared
from pileup.analysis.alpha_pair_spatial_energy import common_inputs,fit_energy

def initialize():
    shared.ROOT,shared.EXP,shared.ART=ROOT,PARENT,ART
    shared.initialize()
    def audit(event,args):
        if event!='open' or not isinstance(args[0],(str,bytes,os.PathLike)):return
        p=Path(os.fsdecode(args[0])).resolve()
        if p.suffix=='.root' or p.name=='mothers.json' or p.is_relative_to(ART/'truth') or p.is_relative_to(ART/'photons') or p.is_relative_to(ART/'predictions'):
            raise PermissionError('Cost pilot cannot read event truth')
    sys.addaudithook(audit)

def one(job):
    i,method,gap=job
    path=ART/'pair-pilot/voltage/pair'/f'{i:04d}.npz'
    with np.load(path) as data:y=data['voltage_mV'][0]
    times=[128.,128.+gap]
    common=common_inputs(shared.frozen.OLD,y,times,shared.frozen.NOISE)
    model=shared.frozen.OLD if method=='old' else shared.frozen.NEW
    result=fit_energy(model,y,times,common)
    (REC/f'cost-fit-{i}-{method}.json').write_text(json.dumps(result,indent=2)+'\n')
    return dict(pilot_id=i,method=method,seconds=result['seconds'],evaluations=result['evaluations'],converged=result['converged'])

if __name__=='__main__':
    shared.ROOT,shared.EXP=ROOT,PARENT
    shared.verify()
    jobs=[(i,m,gap) for i,gap in [(0,8.),(2,40.),(4,1200.)] for m in ('old','corrected')]
    started=time.monotonic()
    with ProcessPoolExecutor(max_workers=6,initializer=initialize) as pool:rows=list(pool.map(one,jobs))
    result=dict(rows=rows,wall_seconds=time.monotonic()-started,pilot_fits=6,new_ROOTs=0,new_digitizations=0,
                inputs='Excluded L97 pilot mothers0/6; first saved noisy repetition at phase0',efficacy_used_for_design=False)
    (REC/'cost-preflight.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
