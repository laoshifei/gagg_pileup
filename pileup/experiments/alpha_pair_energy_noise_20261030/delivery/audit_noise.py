"""Post-result diagnostic of the unexpectedly quiet weak-first Gaussian control.

No fit is changed. The alternative covariance explicitly uses true times and
model charges at true deposited energies, so it is not a blind prescription.
"""
from pathlib import Path
import csv
import json
import sys
import hashlib

import numpy as np

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'pileup').is_dir() and (p/'g4').is_dir())
EXP=ROOT/'pileup/experiments/alpha_pair_energy_noise_20261030'
OUT=ROOT/'jupyter/artifacts/alpha_pair_energy_noise_20261030'
sys.path.insert(0,str(EXP/'source-v2'))
from pileup.analysis.alpha_pair_recovery_energy import RecoveryEnergyResponse
from pileup.analysis.alpha_pair_energy_noise import diagonal_factor,local_energy_projection

response=RecoveryEnergyResponse.load(EXP)
results=[json.loads(p.read_text()) for p in sorted((OUT/'cases').glob('*.json'))]
audit=[]; qualified=[]; initial_failures=[]
for result in results:
    task=result['task']
    if not result['initial_from_reference_voltage']['converged']:initial_failures.append(task['id'])
    with np.load(OUT/'cases'/f"{task['id']:02}.npz") as f:
        physical=f['physical'];q=f['charges'];one=f['factor']
    nr=result['n_reference'];mask=np.all(q[nr:]>=500,axis=1)
    fits=[r for r in result['fits'] if r['domain']=='physical' and r['weight']=='gls']
    hits=[bool(np.all(abs(np.asarray(r['fit']['energy_MeV'])/task['energy_MeV']-1)<=.1)) for r in fits]
    qualified.append(dict(case=task['id'],condition=task['condition'],gap_ns=task['gap_ns'],
        all_frames=len(mask),qualified_frames=int(mask.sum()),qualified_successes=int(np.sum(np.asarray(hits)[mask])),
        unqualified_frames=int((~mask).sum())))
    if task['condition']!='weak_high':continue
    _,model_q=response.components(task['times_ns'],task['energy_MeV'])
    two=response.fixed_metric(dict(times_ns=task['times_ns'],charge_pe=model_q))
    transformed=response.transform(physical[nr:]);time=response.sample_grid[len(response.ar)-1:]
    region=(time>=task['times_ns'][0]+4)&(time<task['times_ns'][1]-4)
    observed=transformed.var(0,ddof=1)
    v1=diagonal_factor(one)[0]**2;v2=diagonal_factor(two)[0]**2
    projection=local_energy_projection(response,task['times_ns'],task['energy_MeV'],two,physical[:nr],physical[nr:])
    audit.append(dict(case=task['id'],block=task['block'],
        original_weight_time_ns=result['initial_from_reference_voltage']['times_ns'][0],
        true_times_ns=task['times_ns'],region_ns=[132.,164.],region_samples=int(region.sum()),
        original_early_variance=float(v1[region].mean()),
        two_component_early_variance=float(v2[region].mean()),
        observed_early_variance=float(observed[region].mean()),
        measured_over_original_variance=float(observed[region].mean()/v1[region].mean()),
        measured_over_two_component_variance=float(observed[region].mean()/v2[region].mean()),
        original_model_local_weak_logE_sd=result['projection']['gaussian_local_logE_sd'][0],
        two_component_model_local_weak_logE_sd=projection['gaussian_local_logE_sd'][0],
        two_component_projected_empirical_weak_logE_sd=projection['evaluation_projected_logE_sd'][0]))

record=dict(posthoc=True,trigger='Matched Gaussian weak-first fit was much narrower than physical repeats',
    scope='Covariance inspection only; alternative factor explicitly uses true K/t/E, no fit changed',
    uncertainty='8 repeated digitizations per optical pair; transformed sample variances are correlated; these ratios have no population CI',
    units='FIR-transformed mV squared; local scale is log-energy, not measured resolution or a bound',
    result=audit,physical_Q_ge_500_frames=qualified,
    initial_weight_optimizer_not_converged_cases=initial_failures,
    main_energy_fits_unchanged=True,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(EXP/'noise-audit.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
for name,rows in [('noise-audit',audit),('qualified-physical-frames',qualified)]:
    with (EXP/'tables'/f'{name}.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
print(json.dumps(record,ensure_ascii=False,indent=2))
