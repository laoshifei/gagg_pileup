"""Frozen post-hoc diagnosis of L98 energy failures, without new inference."""
from __future__ import annotations
import argparse
import ast
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import time
import numpy as np
from . import alpha_mean_shape_study as base
from .alpha_pair_spatial_energy import common_inputs
from .alpha_current_energy_diagnosis import TERM_NAMES, residual_terms, local_operator, repeat_covariance

NAME = 'alpha_current_energy_diagnosis_20261102'
ROOT = Path(os.environ.get('GAGG_CURRENT_ENERGY_DIAGNOSIS_ROOT', Path(__file__).resolve().parents[2])).resolve()
EXP = ROOT/'pileup/experiments'/NAME
ART = ROOT/'jupyter/artifacts'/NAME
REC = ROOT/'jupyter/run_records'/NAME
PARENT = ROOT/'pileup/experiments/alpha_current_pair_energy_20261102'
SHARED = ROOT/'pileup/experiments/alpha_current_shared_20261102'
SHARED_ART = ROOT/'jupyter/artifacts/alpha_current_shared_20261102'
read, save, sha = base.read, base.save, base.sha


def prepare():
    assert not read(PARENT/'gate.json')['passed']
    for name in ('design-freeze.json', 'model-freeze.json', 'pair-fit-freeze.json', 'result-freeze.json'):
        for path, digest in read(PARENT/name).items():assert sha(ROOT/path) == digest, path
    for path in (EXP, ART, REC):path.mkdir(parents=True, exist_ok=False)
    for part in ('input', 'diagnostics', 'replay', 'figures'):(ART/part).mkdir()
    for part in ('bin', 'macros', 'original-runinfo', 'cases'):(EXP/part).mkdir()
    target = EXP/'source/pileup/analysis';target.mkdir(parents=True)
    (target.parent/'__init__.py').write_text('');(target/'__init__.py').write_text('')
    pending = ['alpha_current_energy_diagnosis_study', 'test_alpha_current_energy_diagnosis']
    copied = set();origins = {}
    while pending:
        name = pending.pop()
        if name in copied:continue
        new = name.startswith(('alpha_current_energy_diagnosis', 'test_alpha_current_energy_diagnosis'))
        source = (ROOT/'pileup/analysis' if new else PARENT/'source/pileup/analysis')/(name+'.py')
        assert source.exists(), source
        shutil.copyfile(source, target/source.name);copied.add(name)
        origins[str((target/source.name).relative_to(ROOT))] = dict(source=str(source.relative_to(ROOT)), sha256=sha(source))
        for node in ast.walk(ast.parse(source.read_text())):
            if isinstance(node, ast.ImportFrom) and node.level == 1:
                pending.extend([node.module] if node.module else [a.name for a in node.names])
    for name in ('sipm_age_mean.cpp', 'sipm_age_mean_jacobian.cpp'):
        shutil.copyfile(PARENT/'source/pileup/analysis'/name, target/name)
    shutil.copytree(PARENT/'source/g4', EXP/'source/g4')
    save(EXP/'source-origins.json', origins)
    models = [Path(p).name for p in read(PARENT/'model-freeze.json')]
    for name in models:shutil.copyfile(PARENT/name, EXP/name)
    save(EXP/'model-freeze.json', {str((EXP/n).relative_to(ROOT)):sha(EXP/n) for n in models})
    shutil.copyfile(PARENT/'bin/age.so', EXP/'bin/age.so')
    for name in ('build.json', 'analysis-environment.json', 'transport-environment.json', 'literature-evidence.json'):
        shutil.copyfile(PARENT/name, EXP/name)
    cases = [c for c in read(PARENT/'cases.json') if c['kind'] == 'pair']
    assert len(cases) == 40
    copies = [];originals = {}
    for c in cases:
        for label, original in [('voltage', ROOT/c['voltage']),
              ('truth', SHARED_ART/'truth/pair'/f"{c['source_id']:04d}.npz"),
              ('prediction', SHARED_ART/'predictions'/f"{c['source_id']:04d}.npz")]:
            destination = ART/'input'/f"{c['id']:04d}_{label}.npz"
            shutil.copyfile(original, destination)
            copies.append(dict(case=c['id'],kind=label,source=str(original.relative_to(ROOT)),
                               copy=str(destination.relative_to(ROOT)),sha256=sha(destination)))
            originals[str(original.relative_to(ROOT))] = sha(original)
    save(EXP/'input-origins.json', copies);save(EXP/'original-input-freeze.json', originals)
    save(EXP/'input-freeze.json', {r['copy']:r['sha256'] for r in copies})
    save(EXP/'case-plan.json', cases)
    fits = [read(ROOT/'jupyter/artifacts/alpha_current_pair_energy_20261102/fits'/f"{t['id']:04d}.json")
            for t in read(PARENT/'fit-plan.json') if t['kind']=='pair']
    assert len(fits)==112
    save(EXP/'reference-fits.json', fits)
    save(EXP/'preceding-gate.json', read(PARENT/'gate.json'))
    save(EXP/'preceding-result.json', read(PARENT/'result-summary.json'))
    save(EXP/'upstream-protocol.json', read(SHARED/'protocol.json'))
    ids = {p for c in cases for p in c['parents']};mothers = [r for r in read(PARENT/'input-mothers.json') if r['id'] in ids]
    assert len(mothers)==16
    for m in mothers:
        original = ROOT/m['snapshot_macro'];dest = EXP/'macros'/original.name
        shutil.copyfile(original, dest);m['diagnostic_macro'] = str(dest.relative_to(ROOT))
        shutil.copyfile(PARENT/'original-runinfo'/f"{m['id']:04d}.json", EXP/'original-runinfo'/f"{m['id']:04d}.json")
        assert sha(ROOT/m['root'])==m['root_sha256'] and sha(dest)==m['macro_sha256']
    save(EXP/'input-mothers.json', mothers)
    save(EXP/'provenance-freeze.json', {r['root']:r['root_sha256'] for r in mothers})
    save(EXP/'protocol.json', dict(created_utc=datetime.now(timezone.utc).isoformat(),
        base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        git_status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),
        purpose='Post-hoc truth-assisted diagnosis of exposed L98 failures, not inference or independent confirmation',
        sample='All40L98 pair cases:8disjoint mother pairs,5conditions,phase0;all8saved repetitions;no reserved mothers opened',
        methods='Frozen old and corrected L93 means;unchanged physical parameters,weights construction,energy bounds and original112fits',
        decomposition='y-background-mu(trueE)=strong isolated residual+weak isolated residual+shared-minus-isolated discrepancy+background/readout residual',
        stochastic_scope='Isolated and shared arrays come from separate digitizations with different random consumption;difference is not a noiseless causal interaction effect. Keep all cross covariances.',
        metrics='Individualrepeat0 and8repeatmean each have voltage-derived common metric,identical for2models;mean metric not separately calibrated. Reuse and hash-check L98 metric for112available referencefits.',
        local='TrueE Jacobian columns dmu/dlnE;SVD projection attrueE in inherited whitening;unit100delta_lnE,not percent energy error. Large shifts and boundaries invalidate local extrapolation.',
        repeated='For each observation metric,project all8saved repetitions with that one fixed metric;mean MC SE uses full4term covariance/8. Not total population or energy uncertainty.',
        scope_counts=dict(unique_cases=40,models=2,observation_metrics=2,projection_rows=160,reference_fits=112,
                          new_fits=0,new_ROOTs=0,new_digitizations=0,new_training=0),
        numerical='SVD/orthogonal weak projection identity;residual and covariance closure;first2fixedcases finite-difference gradient relative<3e-4;model identity<1e-8;reference cost relative<1e-8',
        cost='First fixed case used for numerical/cost preflight and retained;project40cases at8workers with1.5factor,limit600s;no efficacy selection',
        uncertainty='10000whole-mother bootstrap draws,seed742100000,common across conditions/methods/terms;only8clusters,descriptive conditional intervals;no parameter fitting or selection',
        replay='Repeat first2fixed cases into separate outputs;all numerical rows equal except timing',
        decision='Describe failure sources and local conditioning;no advancement gate,no parameter tuning,default replacement or new mother generation',
        literature='Reuse L98 saved Fixsen2004/Zhang2024 evidence and review boundaries;residual partition/SVD are project derivations,not claimed paper results;no new library/PDF/web review',
        workers=8,BLAS_threads=1))
    for name in ('README.md', 'docs/literature-use.md', 'jupyter/codex_0007_SiPM_Energy_Calibration.ipynb',
                 'jupyter/Research_Plan_20260921.ipynb', 'jupyter/Run_Catalog.ipynb'):
        p=ROOT/name;shutil.copyfile(p, REC/(p.name+'.before'))
    (EXP/'run.py').write_text("from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\n"
        "os.environ['GAGG_CURRENT_ENERGY_DIAGNOSIS_ROOT']=str(f.parents[2])\nos.environ['OPENBLAS_NUM_THREADS']='1'\nos.environ['OMP_NUM_THREADS']='1'\n"
        "sys.path.insert(0,str(f/'source'))\nif __name__=='__main__':\n    from pileup.analysis.alpha_current_energy_diagnosis_study import main\n    main()\n")
    files=[p for p in EXP.rglob('*') if p.is_file() and 'bin' not in p.relative_to(EXP).parts]
    save(EXP/'design-freeze.json', {str(p.relative_to(ROOT)):sha(p) for p in files})
    print('Frozen40cases,160truth-assisted projections,112referencefits;no new fits',flush=True)


def verify():
    counts={}
    for name in ('design-freeze.json','model-freeze.json','input-freeze.json'):
        data=read(EXP/name)
        for path,digest in data.items():assert sha(ROOT/path)==digest,path
        counts[name]=len(data)
    assert sha(EXP/'bin/age.so')==read(EXP/'build.json')['age_binary_sha256']
    return counts


def initialize():
    base.ROOT,base.EXP,base.ART,base.REC=ROOT,EXP,ART,REC;base.initialize(True)


def whiten(model,metric,values):
    x=model.anchor.transform(values)
    return metric.apply(x.reshape(-1,x.shape[-1]).T).T.reshape(x.shape)


def diagnose(job):
    case,folder=job;started=time.monotonic();cid=case['id']
    with np.load(ART/'input'/f'{cid:04d}_voltage.npz') as f:y=f['voltage_mV']
    with np.load(ART/'input'/f'{cid:04d}_truth.npz') as f:s=f['components_mV'];isolated=f['isolated_mV']
    with np.load(ART/'input'/f'{cid:04d}_prediction.npz') as f:pred={k:f[k] for k in f.files}
    refs={(f['repeat'],f['method']):f for f in read(EXP/'reference-fits.json') if f['case']==cid}
    t=case['times_ns'];energy=np.asarray(case['energy_dep_MeV']);weak=int(np.argmin(energy));order=[1-weak,weak]
    means={};gradients={};checks=[];arrays=dict(time_ns=base.OLD.readout.adc_centers,background_mV=base.OLD.background)
    for method,model in [('old',base.OLD),('corrected',base.NEW)]:
        means[method],gradients[method]=model.mean_and_gradient(t,np.log(energy),model.prepare(t),strength=1.)
        error=float(np.max(abs(means[method]-pred[method+'_joint'].sum(0))))
        assert error<1e-8
        check=dict(method=method,mean_identity_mV=error)
        if cid in (0,1):
            direction=np.array([.7,-.3]);h=2e-5
            finite=(model.components(t,energy*np.exp(h*direction)).sum(0)-model.components(t,energy*np.exp(-h*direction)).sum(0))/(2*h)
            relative=float(np.linalg.norm(gradients[method].T@direction-finite)/np.linalg.norm(finite))
            assert relative<3e-4;check['gradient_relative_error']=relative
        checks.append(check);arrays[method+'_gradient_mV']=gradients[method]
    rows=[]
    for kind,repeat in [('individual',0),('mean',-1)]:
        observation=y[0] if repeat==0 else y.mean(0)
        common=common_inputs(base.OLD,observation,t,base.NOISE);metric=common['metric']
        signature=hashlib.sha256(metric.factor.tobytes()+metric.modes.tobytes()).hexdigest()
        for method,model in [('old',base.OLD),('corrected',base.NEW)]:
            pieces,total=residual_terms(y,s,isolated,pred[method+'_joint'],pred[method+'_solo'],model.background,weak)
            d=whiten(model,metric,gradients[method][order]).T
            operator,geometry=local_operator(d)
            wp=whiten(model,metric,pieces);wt=whiten(model,metric,total)
            projections=100*(wp@operator.T);total_projection=100*(wt@operator.T)
            closure=float(np.max(abs(projections.sum(1)-total_projection)));assert closure<1e-8
            pp=projections[0] if repeat==0 else projections.mean(0)
            actual_r=wt[0] if repeat==0 else wt.mean(0)
            ptotal=total_projection[0] if repeat==0 else total_projection.mean(0)
            selected_parts=pieces[0] if repeat==0 else pieces.mean(0)
            covariance=repeat_covariance(projections[:,:,1])
            row=dict(case=cid,block=case['block'],condition=case['condition'],gap_ns=case['gap_ns'],method=method,kind=kind,
                true_E_MeV=energy[order].tolist(),truth_assisted=True,metric_signature=signature,geometry=geometry,
                projection_terms=pp.tolist(),projection_total=ptotal.tolist(),repeat_covariance=covariance,
                residual_rms_mV=float(np.sqrt(np.mean((observation-model.background-means[method])**2))),
                term_rms_mV=np.sqrt(np.mean(selected_parts**2,axis=1)).tolist(),true_cost=float(actual_r@actual_r),
                closure_mV=float(np.max(abs(pieces.sum(1)-total))),projection_closure=closure,reference_fit=None)
            if (repeat,method) in refs:
                fit=refs[repeat,method]
                assert fit['metric_signature']==signature
                assert fit['initial_MeV']==common['initial'].tolist() and fit['start_energy_MeV']==common['starts'].tolist()
                curve=model.components(t,fit['energy_MeV']).sum(0)
                fitted_r=whiten(model,metric,observation-model.background-curve)
                cost=float(fitted_r@fitted_r);costerror=abs(cost-fit['cost'])/max(1.,abs(fit['cost']))
                assert costerror<1e-8,(cid,method,kind,costerror)
                actual=100*np.log(np.asarray(fit['energy_MeV'])[order]/energy[order])
                row['reference_fit']=dict(id=fit['id'],energy_log_error=actual.tolist(),
                    linear_minus_fit=(ptotal-actual).tolist(),energy_relative_error=(100*(np.asarray(fit['energy_MeV'])[order]/energy[order]-1)).tolist(),
                    cost=cost,cost_replay_relative_error=costerror,cost_improvement=row['true_cost']-cost,
                    converged=fit['converged'],boundary=fit['boundary'],selected_start=fit['selected_start'])
                arrays[f'{method}_{kind}_fit_mean_mV']=curve
            arrays[f'{method}_{kind}_projected_repeats']=projections
            arrays[f'{method}_{kind}_jacobian']=d
            arrays[f'{method}_{kind}_weak_influence']=operator[1]
            arrays[method+'_terms_mV']=pieces
            rows.append(row)
    result=dict(case=cid,rows=rows,checks=checks,seconds=time.monotonic()-started)
    out=ART/folder;out.mkdir(exist_ok=True)
    np.savez_compressed(out/f'{cid:04d}.npz',**arrays);save(out/f'{cid:04d}.json',result)
    return dict(case=cid,seconds=result['seconds'])


def run(stage):
    verify();cases=read(EXP/'case-plan.json')
    if stage=='pilot':
        assert not (EXP/'cost-check.json').exists();initialize()
        assert not (ART/'diagnostics/0000.json').exists()
        row=diagnose((cases[0],'diagnostics'));estimate=row['seconds']*40/8*1.5
        save(EXP/'cost-check.json',dict(passed=estimate<600,first_case=row,projected_seconds=estimate))
        assert estimate<600;print('Pilot passed; projected seconds',estimate,flush=True);return
    assert read(EXP/'cost-check.json')['passed']
    assert not (EXP/'run-completion.json').exists()
    for c in cases[1:]:assert not (ART/'diagnostics'/f"{c['id']:04d}.json").exists()
    started=time.monotonic()
    with ProcessPoolExecutor(max_workers=8,initializer=initialize) as pool:
        pending=[pool.submit(diagnose,(c,'diagnostics')) for c in cases[1:]]
        for count,future in enumerate(as_completed(pending),2):
            future.result()
            if count%8==0:print('Diagnosed',count,'/40cases',flush=True)
    save(EXP/'run-completion.json',dict(cases=40,seconds_after_pilot=time.monotonic()-started,new_fits=0,new_ROOTs=0,new_digitizations=0))
    paths=sorted((ART/'diagnostics').glob('*'))
    save(EXP/'diagnostic-output-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in paths})


def score():
    verify();assert not (EXP/'result-summary.json').exists()
    for path,digest in read(EXP/'diagnostic-output-freeze.json').items():assert sha(ROOT/path)==digest
    results=[read(ART/'diagnostics'/f"{c['id']:04d}.json") for c in read(EXP/'case-plan.json')]
    rows=[r for result in results for r in result['rows']];assert len(rows)==160
    blocks=list(dict.fromkeys(r['block'] for r in rows))
    index=np.random.default_rng(742100000).integers(8,size=(10000,8));groups=[]
    for condition,gap in [('high_weak',8.),('high_weak',40.),('high_weak',200.),('high_weak',1200.),('weak_high',40.)]:
        for kind in ('individual','mean'):
            for method in ('old','corrected'):
                chosen=[r for r in rows if (r['condition'],r['gap_ns'],r['kind'],r['method'])==(condition,gap,kind,method)]
                chosen.sort(key=lambda r:blocks.index(r['block']));assert len(chosen)==8
                terms=np.array([r['projection_terms'] for r in chosen])[:,:,1];total=terms.sum(1)
                actual=[r['reference_fit'] for r in chosen];available=all(r is not None for r in actual)
                group=dict(condition=condition,gap_ns=gap,kind=kind,method=method,mother_pairs=8,
                    weak_terms_mean=terms.mean(0).tolist(),weak_terms_ci95=np.percentile(terms[index].mean(1),[2.5,97.5],axis=0).T.tolist(),
                    weak_projection_mean=float(total.mean()),weak_projection_ci95=np.percentile(total[index].mean(1),[2.5,97.5]).tolist(),
                    median_cosine=float(np.median([r['geometry']['cosine'] for r in chosen])),
                    median_normalized_condition=float(np.median([r['geometry']['normalized_condition'] for r in chosen])),
                    median_weak_retained_fraction=float(np.median([r['geometry']['weak_retained_fraction'] for r in chosen])),
                    residual_rms_mV=float(np.sqrt(np.mean([r['residual_rms_mV']**2 for r in chosen]))),
                    term_rms_mV=np.sqrt(np.mean(np.array([r['term_rms_mV'] for r in chosen])**2,axis=0)).tolist(),
                    median_mean_mc_se=float(np.median([r['repeat_covariance']['mean_mc_standard_error'] for r in chosen])),
                    mean_repeat_total_variance=float(np.mean([r['repeat_covariance']['total_variance'] for r in chosen])),
                    mean_repeat_diagonal_variance=float(np.mean([r['repeat_covariance']['diagonal_variance'] for r in chosen])),
                    mean_repeat_cross_variance=float(np.mean([r['repeat_covariance']['cross_variance'] for r in chosen])),
                    reference_available=available)
                if available:
                    fitted=np.array([r['energy_log_error'][1] for r in actual])
                    group.update(weak_fitted_log_mean=float(fitted.mean()),linear_minus_fit_rmse=float(np.sqrt(np.mean((total-fitted)**2))),
                        fit_better_than_truth=sum(r['cost_improvement']>=-1e-8 for r in actual),
                        boundary=sum(r['boundary'] for r in actual),nonconverged=sum(not r['converged'] for r in actual))
                groups.append(group)
    refs=[r['reference_fit'] for r in rows if r['reference_fit'] is not None];assert len(refs)==112
    checks=dict(max_mean_identity_mV=max(c['mean_identity_mV'] for r in results for c in r['checks']),
        max_gradient_relative_error=max(c.get('gradient_relative_error',0.) for r in results for c in r['checks']),
        max_cost_replay_relative_error=max(r['cost_replay_relative_error'] for r in refs),
        max_residual_closure_mV=max(r['closure_mV'] for r in rows),max_projection_closure=max(r['projection_closure'] for r in rows))
    summary=dict(groups=groups,checks=checks,blocks=blocks,cases=40,projection_rows=160,reference_fits=112,
        fit_better_than_truth=sum(r['cost_improvement']>=-1e-8 for r in refs),new_fits=0,new_ROOTs=0,new_digitizations=0,
        default_changed=False,diagnostic_only=True)
    save(EXP/'details.json',rows);save(EXP/'result-summary.json',summary)
    save(EXP/'result-freeze.json',{str((EXP/n).relative_to(ROOT)):sha(EXP/n) for n in ('details.json','result-summary.json')})
    print('Scored160diagnostic projections;',summary['fit_better_than_truth'],'/112saved fits cost <= truth cost',flush=True)


def replay():
    verify();assert not (EXP/'replay.json').exists()
    cases=read(EXP/'case-plan.json')[:2]
    with ProcessPoolExecutor(max_workers=2,initializer=initialize) as pool:list(pool.map(diagnose,[(c,'replay') for c in cases]))
    for c in cases:
        a=read(ART/'diagnostics'/f"{c['id']:04d}.json");b=read(ART/'replay'/f"{c['id']:04d}.json")
        assert {k:v for k,v in a.items() if k!='seconds'}=={k:v for k,v in b.items() if k!='seconds'}
        with np.load(ART/'diagnostics'/f"{c['id']:04d}.npz") as x,np.load(ART/'replay'/f"{c['id']:04d}.npz") as y:
            assert x.files==y.files and all(np.array_equal(x[k],y[k]) for k in x.files)
    save(EXP/'replay.json',dict(passed=True,cases=2,projection_rows=8,new_fits=0,exact_except_timing=True))
    print('Two diagnostic cases replay exactly;no new energy fits',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['prepare','verify','pilot','run','score','replay'])
    stage=parser.parse_args().stage
    if stage=='prepare':prepare()
    elif stage=='verify':print(verify())
    elif stage in ('pilot','run'):run(stage)
    elif stage=='score':score()
    else:replay()


if __name__=='__main__':main()
