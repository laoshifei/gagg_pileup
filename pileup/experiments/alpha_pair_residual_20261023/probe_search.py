"""Training-only cost and accuracy check of residual candidate generation."""
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
RESPONSE=None


def initialize():
    global RESPONSE
    from pileup.analysis.alpha_pair_conditional import ConditionalResponse
    from pileup.analysis.alpha_pair_conditional_study import forbid_labels
    sys.addaudithook(forbid_labels)
    RESPONSE=ConditionalResponse.load(PROJECT/'pileup/experiments/alpha_pair_conditional_20261022/response.npz')


def worker(task):
    from pileup.analysis.alpha_pair_residual import fit_candidates
    i,y,seeds=task;start=time.monotonic()
    return dict(index=i,rows=[fit_candidates(RESPONSE,y,seeds,fraction=f,rho=rho) for f,rho in
        ((0.,None),(.1,None),(.3,None),(0.,.9))],seconds=time.monotonic()-start)


def main():
    import numpy as np
    from pileup.analysis.alpha_pair_conditional_study import evaluate_fit
    from pileup.analysis.alpha_pair_residual import choose
    cache=Path(json.loads((PROJECT/'jupyter/run_records/alpha_pair_conditional_20261022/local-cache.json').read_text())['directory'])
    cases=json.loads((cache/'development.cases.json').read_text())
    seeds=json.loads((cache/'development.seeds.json').read_text())
    with np.load(cache/'development.waveforms.npz') as f:y=f['voltage_mV']
    with np.load(cache/'development.truth.npz') as f:q=f['charge_pe'];t=f['times_ns'];train=f['train']
    used={};indices=[]
    for i,c in enumerate(cases):
        if not train[i]:continue
        if c['kind']=='pair' and c['gap_ns'] not in (8.,20.,80.,200.):continue
        key=(c['condition'],c['gap_ns']);limit=4 if c['kind']=='pair' else 8
        if used.get(key,0)>=limit:continue
        used[key]=used.get(key,0)+1;indices.append(i)
    print('TRAINING_SEARCH_PROBE',len(indices),flush=True)
    tasks=[(i,y[i],{k:v[i] for k,v in seeds.items()}) for i in indices]
    rows=[];start=time.monotonic()
    with ProcessPoolExecutor(max_workers=4,initializer=initialize,mp_context=mp.get_context('spawn')) as pool:
        for result in pool.map(worker,tasks):
            rows.append(result)
            if len(rows)%10==0:print('SEARCH_PROBE',len(rows),'/',len(indices),flush=True)
    output=PROJECT/'jupyter/artifacts/alpha_pair_residual_20261023'
    path=output/'training-search-probe.json'
    assert not path.exists()
    path.write_text(json.dumps(dict(scope='Training parents only',seconds=time.monotonic()-start,rows=rows),indent=2)+'\n')
    for j,(fraction,rho) in enumerate(((0.,None),(.1,None),(.3,None),(0.,.9))):
        for penalty in (.0001,.0003,.001,.003,.01):
            scores=[evaluate_fit(choose(r['rows'][j]['expanded'],penalty),t[r['index']],q[r['index']],cases[r['index']]) for r in rows]
            print(fraction,rho,penalty,[(c,sum(s['success'] for r,s in zip(rows,scores) if cases[r['index']]['condition']==c),
                sum(s['extra']>0 for r,s in zip(rows,scores) if cases[r['index']]['condition']==c)) for c in
                ['low_low','mid_mid','high_high','low_high','high_low','low','mid','high','all']],flush=True)
    print('SEARCH_PROBE_COMPLETE',time.monotonic()-start,flush=True)


if __name__=='__main__':main()
