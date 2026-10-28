"""Frozen re-selection of previously computed, observation-only fit candidates."""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import importlib.metadata
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

from .alpha_pair_information_study import load,save,digest,rows
from .alpha_pair_physical_gate import assemble,replay_cost

PROJECT=Path(os.environ.get('GAGG_GATE_PROJECT',Path(__file__).resolve().parents[2])).resolve()
NAME='alpha_pair_physical_gate_20261027'
EXPERIMENT=PROJECT/'pileup/experiments'/NAME
OUTPUT=PROJECT/'jupyter/artifacts'/NAME
RECORD=PROJECT/'jupyter/run_records'/NAME
PREVIOUS=PROJECT/'pileup/experiments/alpha_pair_risk_training_v2_20261027'
INPUT=PROJECT/'jupyter/artifacts/alpha_pair_risk_training_v2_20261027'
OLD_RECORD=PROJECT/'jupyter/run_records/alpha_pair_risk_training_v2_20261027'


def models():
    from .alpha_pair_joint_likelihood import ConservingResponse,JointNoise
    return ConservingResponse.load(EXPERIMENT/'joint-response.npz'),JointNoise.load(EXPERIMENT/'joint_conserved_nll-noise.npz')


def check_freeze():
    for name,expected in load(EXPERIMENT/'freeze.json')['files'].items():assert digest(EXPERIMENT/name)==expected,name


def input_paths(domain):
    return dict(observation=INPUT/('native-observations.npz' if domain=='native' else 'evaluation/observations.npz'),
        baseline=INPUT/('native-predictions.baseline.jsonl' if domain=='native' else 'evaluation/predictions.baseline.jsonl'),
        predictions=INPUT/(domain+'-predictions.jsonl'),
        truth=INPUT/('native-labels.truth.npz' if domain=='native' else 'evaluation/labels.truth.npz'),
        cases=INPUT/('native-labels.cases.json' if domain=='native' else 'evaluation/labels.cases.json'))


def prepare():
    prior=load(OLD_RECORD/'final-validation.json');assert prior['passed']
    notebook=PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    assert digest(notebook)==prior['notebook_sha256']
    for name,expected in prior['source_sha256'].items():assert digest(PROJECT/name)==expected,name
    for folder in (EXPERIMENT,OUTPUT,RECORD):folder.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(notebook,RECORD/'notebook.before.ipynb')
    shutil.copyfile(PROJECT/'README.md',RECORD/'README.before.md')
    shutil.copyfile(PROJECT/'docs/literature-use.md',RECORD/'literature.before.md')
    for name in ('joint-response.npz','joint_conserved_nll-noise.npz','baseline-gate.json','settings.json'):
        shutil.copyfile(PREVIOUS/name,EXPERIMENT/name)
    shutil.copyfile(PREVIOUS/'gate-calibrations.json',EXPERIMENT/'previous-gates.json')
    save(EXPERIMENT/'previous-model-selection.json',load(PREVIOUS/'training-freeze.json'))
    protocol=dict(created=datetime.now(timezone.utc).isoformat(),scope='Pure alpha, fixed95/351ns, full2ns x1024 analog waveform, ample dynamic range, local unknown K<=2',
        question='Keep physical count decision; test whether the already frozen network contributes useful candidate starts beyond one extra local refinement',
        candidates=dict(baseline='Original K0/K1/K2',local='Preserve K0/K1; min(original K2, cached trial0 from baseline K2 with70 iterations)',
            mlp='Same already-computed unregularized MLP+NLL minimum from original/network initials and their two70-iteration trials; no learned-Q penalty',
            selection='Do not silently replace an ineligible global minimum with an eligible higher-cost candidate; preserve previous Q>=500 rule'),
        methods=dict(baseline='Original candidate + original physical gate',learned_gate='Same hybrid candidate + previous learned-logit gate, for exact paired replay',
            local_physical='Local-only candidate + recalibrated physical NLL gate',mlp_physical='Hybrid candidate + recalibrated physical NLL gate',
            mlp_fixed_gate='Hybrid candidate + original numerical physical gate; declared sensitivity analysis'),
        gate='max(NLL1-NLL2,0), candidate eachQ>=500; same280 native singles, observed-singleQ bins1500/5000/12500, plus-one tie-preserving rank<=.025; K0/K1 presence unchanged',
        budgets='Cached local control has one additional70-iteration start, hybrid has two; not equal-compute architecture comparison. No new optimization or training.',
        data='280 native calibration +760 exposed native evaluation +1000 previously exposed continuous Gaussian evaluation; no new independent data in this turn',
        counts=dict(all_frames=2040,calibration=280,evaluation=1760,primary_methods=4,sensitivity_methods=1,primary_scores=7040,all_scores=8800,forced_candidate_rows=2880),
        independence='All evaluation is now development reuse. No new optical parents, ROOTs, Gaussian draws, model fitting or network tuning. Prior seed/model provenance retained.',
        metrics='Each trueQ>=500; eachQwithin20%; eachtwithinmin(10ns,gap/4); no extra output. Report missing/extra/split/merge, matched Q bias/sigma68/time errors, forced K2 and exact rejection reasons.',
        diagnosis='For every correct K2 candidate classify accepted, subthreshold, no presence, nonpositive evidence or rank rejection. Counterfactual same-candidate gates isolate selection losses.',
        decision='Gate repair passes if native200ns high_weak >=20/40, native total>=241/360 and continuous total>=500/600, with all per-group controls extra<=5%,single success>=95%,abs medianQbias<=5%. Network adds practical value only if mlp_physical beats local_physical by>=5pp on both domains and passes controls. No automatic default replacement.',
        uncertainty='2000 paired bootstrap replicates; within each condition native optical-parent blocks, continuous draws independently. Descriptive exposed-data intervals, no physical uncertainty claim.',
        pilot=dict(native=[0,70,140,210,280,320,360,400,440,480,520,560,600,640,720,800,880,960],continuous=[0,40,80,120,160,200,240,280,320,360,600,680,717,999]),
        examples='Prespecified native high_weak replicate0 at40 and200ns; show full0..2048ns final selected outputs, absent second output=0',
        figures=['0004-93-1','0004-93-2','0004-94-1','0004-95-1','0004-96-1'],
        literature='Reuse Liddick2012 detection/parameter-estimation separation and false double controls; Morad2024 limited human-confirmed distinction of detection versus conditional classification. Current rank rule and ablation are project choices, not paper reproductions.')
    save(EXPERIMENT/'protocol.json',protocol)
    inputs={str(p.relative_to(PROJECT)):digest(p) for d in ('native','continuous') for p in input_paths(d).values()}
    for name in ('scores.jsonl','forced-K2.jsonl','report.json'):
        inputs[str((INPUT/name).relative_to(PROJECT))]=digest(INPUT/name)
    provenance=dict(prior_validation=prior,previous_source_freezes={name:digest(PREVIOUS/name) for name in ('freeze.json','evaluation-f64-freeze.json','training-freeze.json','report-v3-freeze.json')},
        inputs=inputs,mother_groups=load(OLD_RECORD/'input-provenance.json')['mother_groups'],
        notebook_before_sha256=digest(notebook),training_seeds=load(PREVIOUS/'protocol.json')['training']['seeds'],
        inherited_model_selection=load(PREVIOUS/'training-freeze.json'),new_ROOTs=0,new_random_draws=0)
    save(RECORD/'input-provenance.json',provenance)
    save(RECORD/'environment.json',dict(python=sys.version,executable=sys.executable,blas_threads=1,max_workers=1,
        dependencies={n:importlib.metadata.version(n) for n in ('numpy','scipy','matplotlib','pandas','nbformat')},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip(),
        uncommitted='Earlier information/allocation/risk studies retained; physical-gate sources added and frozen'))
    shutil.copytree(PREVIOUS/'evaluation-f64-source',EXPERIMENT/'inference-source')
    for name in ('alpha_pair_physical_gate.py','alpha_pair_physical_gate_study.py','test_alpha_pair_physical_gate.py'):
        shutil.copyfile(PROJECT/'pileup/analysis'/name,EXPERIMENT/'inference-source/pileup/analysis'/name)
    (EXPERIMENT/'run.py').write_text('''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
for name in ('GATE','RISK','INFORMATION'):os.environ['GAGG_'+name+'_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'inference-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_physical_gate_study import main
    main()
''')
    files={str(p.relative_to(EXPERIMENT)):digest(p) for p in (EXPERIMENT/'inference-source').rglob('*.py')}
    for name in ('run.py','protocol.json','previous-gates.json','previous-model-selection.json','baseline-gate.json','settings.json','joint-response.npz','joint_conserved_nll-noise.npz'):files[name]=digest(EXPERIMENT/name)
    save(EXPERIMENT/'freeze.json',dict(files=files,created=protocol['created']))
    print('PHYSICAL_GATE_PREPARED',len(files),flush=True)


def infer(*,pilot=False):
    import json
    from .alpha_pair_charge_allocation import forbid_labels
    check_freeze();provenance=load(RECORD/'input-provenance.json');protocol=load(EXPERIMENT/'protocol.json')
    if not pilot:assert load(RECORD/'pilot-validation.json')['passed'] and load(RECORD/'numerical-tests.json')['passed']
    response,noise=models();sys.addaudithook(forbid_labels)
    started=time.monotonic();maximum_error=0.;numbers={};failed=0
    for domain in ('native','continuous'):
        paths=input_paths(domain)
        for key in ('observation','baseline','predictions'):
            assert digest(paths[key])==provenance['inputs'][str(paths[key].relative_to(PROJECT))]
        with np.load(paths['observation']) as f:y=f['voltage_mV']
        base={r['index']:r['candidates'] for r in rows(paths['baseline'])}
        prediction={r['index']:r for r in rows(paths['predictions'])}
        assert set(base)==set(prediction)==set(range(len(y)))
        ids=protocol['pilot'][domain] if pilot else range(len(y));numbers[domain]=len(ids)
        path=OUTPUT/(domain+('-pilot' if pilot else '')+'-candidates.jsonl')
        with path.open('x') as stream:
            for count,i in enumerate(ids,1):
                sets,learned=assemble(base[i],prediction[i]['variants']['mlp_nll'])
                for bank in sets.values():
                    for candidate in bank:
                        actual=replay_cost(response,noise,y[i],candidate)
                        error=abs(actual-candidate['cost'])/max(1.,abs(actual));maximum_error=max(maximum_error,error)
                        assert error<1e-10,(domain,i,error)
                failed+=sum(not group[-1].get('converged',True) for group in sets.values())
                stream.write(json.dumps(dict(index=i,domain=domain,sets=sets,learned=learned),allow_nan=False)+'\n')
                if count%400==0:print('PHYSICAL_GATE_REPLAY',domain,count,len(ids),flush=True)
    outputs={p.name:digest(p) for p in OUTPUT.glob('*-pilot-candidates.jsonl' if pilot else '*-candidates.jsonl') if pilot or '-pilot-' not in p.name}
    save(RECORD/('pilot-validation.json' if pilot else 'inference.json'),dict(passed=True,frames=numbers,seconds=time.monotonic()-started,
        maximum_cost_replay_relative_error=maximum_error,retained_nonconverged_last_candidates=failed,truth_inputs=[],label_access_blocked=True,
        no_new_optimization=True,no_new_network_inference=True,outputs=outputs,source_sha256=digest(__file__)))
    print('PHYSICAL_GATE_REPLAY_COMPLETE',numbers,'pilot',pilot,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=('prepare','pilot','infer'));args=parser.parse_args()
    if args.command=='prepare':prepare()
    else:infer(pilot=args.command=='pilot')


if __name__=='__main__':main()
