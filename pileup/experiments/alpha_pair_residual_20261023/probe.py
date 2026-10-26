"""Small training-parent diagnostic before a full validation scan."""
from pathlib import Path
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import time

PROJECT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(PROJECT))
for key in ('GAGG_CONDITIONAL_PROJECT','GAGG_RESOLUTION_PROJECT','GAGG_PAIR_PROJECT','GAGG_COMPONENT_PROJECT'):
    os.environ[key]=str(PROJECT)


def worker(task):
    import numpy as np
    from pileup.analysis.alpha_pair_conditional import ConditionalResponse
    from pileup.analysis.alpha_pair_conditional_study import forbid_labels
    from pileup.analysis.alpha_pair_residual import PrecisionBlend
    sys.addaudithook(forbid_labels)
    response=ConditionalResponse.load(PROJECT/'pileup/experiments/alpha_pair_conditional_20261022/response.npz')
    i,y,starts=task;answer=[]
    for fraction,rho in [(x,None) for x in (0.,.001,.003,.01,.03,.1,.3,1.)]+[(0.,r) for r in (.9,.95,.97)]:
        metric=PrecisionBlend(response,y,fraction,rho=rho)
        fits=[response.refine(y,s,metric=metric) for s in starts]
        # Choose among observed same-count seeds, never by event truth.
        count=max(len(s['times_ns']) for s in starts)
        fits=[f for f in fits if len(f['times_ns'])==count]
        fit=min(fits,key=lambda f:f['cost'])
        answer.append(dict(fraction=fraction,rho=rho,fit=fit))
    return dict(index=i,rows=answer)


def main():
    import numpy as np
    from pileup.analysis.alpha_pair_conditional_study import evaluate_fit
    previous=PROJECT/'jupyter/run_records/alpha_pair_conditional_20261022/local-cache.json'
    cache=Path(json.loads(previous.read_text())['directory'])
    cases=json.loads((cache/'development.cases.json').read_text())
    seeds=json.loads((cache/'development.seeds.json').read_text())
    with np.load(cache/'development.waveforms.npz') as f:y=f['voltage_mV']
    with np.load(cache/'development.truth.npz') as f:q=f['charge_pe'];t=f['times_ns'];train=f['train']
    used={};indices=[]
    for i,c in enumerate(cases):
        if not train[i] or c['kind']=='empty':continue
        if c['kind']=='pair' and c['gap_ns'] not in (8.,20.,80.,200.):continue
        key=(c['condition'],c['gap_ns']);limit=4 if c['kind']=='pair' else 8
        if used.get(key,0)>=limit:continue
        used[key]=used.get(key,0)+1;indices.append(i)
    print('TRAINING_PROBE',len(indices),'frames',flush=True)
    tasks=[(i,y[i],[seeds['ml_1'][i],seeds['ml_2'][i]]) for i in indices]
    start=time.monotonic();rows=[]
    with ProcessPoolExecutor(max_workers=4,mp_context=mp.get_context('spawn')) as pool:
        for result in pool.map(worker,tasks):
            i=result['index'];case=cases[i]
            for r in result['rows']:
                r.update(index=i,condition=case['condition'],gap_ns=case['gap_ns'],kind=case['kind'],
                         score=evaluate_fit(r['fit'],t[i],q[i],case))
                rows.append(r)
            if len(rows)%110==0:print('PROBE',len(rows)//11,'/',len(indices),flush=True)
    output=PROJECT/'jupyter/artifacts/alpha_pair_residual_20261023';output.mkdir(exist_ok=True,parents=True)
    (output/'training-precision-probe.json').write_text(json.dumps(dict(
        scope='Training parents only; no held-out selection or test',seconds=time.monotonic()-start,rows=rows),indent=2)+'\n')
    for fraction,rho in [(x,None) for x in (0.,.001,.003,.01,.03,.1,.3,1.)]+[(0.,r) for r in (.9,.95,.97)]:
        rr=[r for r in rows if r['fraction']==fraction and r['rho']==rho]
        print(fraction,rho,[(c,sum(r['score']['success'] for r in rr if r['condition']==c),
            int(sum(r['condition']==c for r in rr))) for c in ['low_low','mid_mid','high_high','low_high','high_low','low','mid','high']],flush=True)
    print('PROBE_COMPLETE',time.monotonic()-start,flush=True)


if __name__=='__main__':main()
