"""Independent-mother mean/noise calibration and paired known-K/t ablation."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import hashlib
from importlib.metadata import version
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

from .alpha_pair_energy_calibrated import (METHODS, CalibratedEnergyResponse,
    calibrate_mean, calibrate_noise, compare_known_times, learned_metric)
from .alpha_pair_recovery_energy import RecoveryEnergyResponse

PROJECT = Path(os.environ.get("GAGG_ENERGY_CALIBRATED_PROJECT", Path(__file__).resolve().parents[2])).resolve()
NAME = "alpha_pair_energy_calibrated_20261030"
EXPERIMENT = PROJECT / "pileup/experiments" / NAME
OUTPUT = PROJECT / "jupyter/artifacts" / NAME
RECORD = PROJECT / "jupyter/run_records" / NAME
ART = PROJECT / "jupyter/artifacts"
CONDITIONS = [("weak", None), ("low", None), ("mid", None), ("high", None),
              ("high_weak", 8.), ("high_weak", 40.), ("high_weak", 200.), ("weak_high", 40.)]


def load(path):
    return json.loads(Path(path).read_text())


def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        while block:=f.read(1024*1024):h.update(block)
    return h.hexdigest()


def save(path, value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name(path.name+".partial")
    temporary.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    temporary.replace(path)


def snapshot():
    source=EXPERIMENT/"source/pileup/analysis"
    source.mkdir(parents=True)
    (source.parent/"__init__.py").write_text("")
    (source/"__init__.py").write_text("")
    for name in ("alpha_pair_energy_calibrated.py", "alpha_pair_energy_calibrated_study.py",
                 "test_alpha_pair_energy_calibrated.py", "alpha_pair_energy_covariance.py",
                 "alpha_pair_energy_noise.py", "alpha_pair_energy_oracle.py",
                 "alpha_pair_recovery_energy.py", "alpha_pair_model_sampler.py",
                 "test_alpha_pair_recovery_energy.py"):
        shutil.copyfile(PROJECT/"pileup/analysis"/name,source/name)
    (EXPERIMENT/"run.py").write_text(
        "from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\n"
        "os.environ['GAGG_ENERGY_CALIBRATED_PROJECT']=str(f.parents[2])\n"
        "sys.path.insert(0,str(f/'source'))\nif __name__=='__main__':\n"
        "    from pileup.analysis.alpha_pair_energy_calibrated_study import main\n    main()\n")


def prepare():
    assert not EXPERIMENT.exists(), "Keep existing experiment immutable"
    for p in (EXPERIMENT,OUTPUT/"inputs",OUTPUT/"calibration",RECORD):p.mkdir(parents=True,exist_ok=True)
    previous=PROJECT/"pileup/experiments/alpha_pair_energy_covariance_20261030"
    for name in ("energy-model.json","template.npz","response-parameters.json","sensor-config.json"):
        shutil.copyfile(previous/name,EXPERIMENT/name)
    snapshot()
    protocol=dict(created=datetime.now().astimezone().isoformat(),
        commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=PROJECT,text=True).strip(),
        uncommitted="Pre-existing README/literature changes and earlier untracked energy studies retained; actual source and configuration copied here",
        purpose="Separately calibrate energy-conditioned mean and effective residual covariance, then factorial2x2 energy-only ablation",
        scope="Pure alpha; TRUE COUNT/TIMES explicitly supplied; no blind event count, timing or particle classification",
        data="Mean: development fold0; noise: development fold1; evaluation: exposed20261026 confirmation mothers. Pairwise disjoint optical roots",
        exposure="Evaluation was inspected in previous research; retrospective holdout for the newly fitted parameters, NOT fresh untouched confirmation. Frozen inherited models were historically developed/selected",
        calibration="Mean1180frames/280mothers; noise1180signal+80empty/280different mothers; evaluation480frames/320different mothers. Repeated use of mothers is not independent sampling",
        mean="Frozen q0(E_dep) coordinate;9 causal shape features, ridge0.03(intercept0.0003),4ns smoothing; additive finite-support voltage correction and consistent model charge-area change; no realized Q regressor",
        noise="Native FIR grid, no interpolation; positive background scale from empty frames;3 bounded log-scale parameters, lag1/cross clipped±0.45; four positive centered fractional amplitude modes; residual moments can include mean error",
        ablation="Both old/noise metrics frozen at the same voltage-only seed E; new noise uses new mean. Same metric for both mean variants and same starts; fixed-C GLS, no parameter-dependent NLL search",
        methods=list(METHODS),primary="new_mean_new_noise versus old_mean_old_noise",
        seed="One-pulse voltage refinement then old fixed-K/t energy fit. Different from previous16-reference-repeat study; only same-wave current comparisons are paired",
        fitting="Common seed+linear+5charge-share starts for pairs, seed+linear for singles;120least-squares evaluations per start; fixed inherited energy bounds",
        conditions=CONDITIONS,frames=480,method_rows=1920,
        endpoint="All injected E_dep within±5/10/20%; no Q filter for primary. Bias,median,sigma68,RMSE,boundary,nonconvergence,paired gained/lost. Q>=500 subset separate",
        unknown_count_metrics="Misses/extras/split/merge/PSD efficiency are not measured",
        model="Frozen95/351ns dual decay,29ns recovery,14400cells,strength2.515148098; project assumptions, not literature alpha constants",
        observations="1024x2ns; native centers0.25+2n ns; existing full analog voltage/dynamic range; no clipping rejection",
        pilot="First weak/high single and first high_weak8/40 and weak_high40; cost and invariants only, no accuracy tuning",
        workers=6,BLAS_threads=1,cost_gate_seconds=1800,cost_safety_factor=1.5,
        no_new_ROOT=True,no_new_G4=True,default_changed=False)
    save(EXPERIMENT/"protocol.json",protocol)
    save(EXPERIMENT/"method-freeze.json",dict(before_calibration_and_evaluation=True,
        files={str(p.relative_to(EXPERIMENT)):digest(p) for p in EXPERIMENT.rglob("*") if p.is_file()}))
    dev=ART/"alpha_pair_joint_likelihood_20261026"
    old=ART/"alpha_pair_recovery_energy_20261029/inputs"
    confirm=ART/"alpha_pair_joint_confirmation_20261026"
    source_paths=[dev/"calibration.cases.json",dev/"calibration.observations.npz",dev/"calibration-ancestry.truth.npz",
        old/"development.energy.json",old/"evaluation.cases.json",old/"evaluation.labels.npz",old/"observations.npz",
        confirm/"physical-parents.json",confirm/"labels.truth.npz"]
    save(EXPERIMENT/"source-input-hashes.json",{str(p.relative_to(PROJECT)):digest(p) for p in source_paths})
    cases=load(dev/"calibration.cases.json")
    mothers=load(old/"development.energy.json");energy={r["root"]:r["edep_MeV"] for r in mothers}
    roots=[]
    with np.load(dev/"calibration-ancestry.truth.npz") as truth, np.load(dev/"calibration.observations.npz") as obs:
        components=truth["components_mV"];times=truth["times_ns"];voltage=obs["voltage_mV"]
        for fold,name in ((0,"mean"),(1,"noise")):
            selected=[i for i,c in enumerate(cases) if c["fold"]==fold and
                (c["kind"]=="single" or c["kind"]=="pair" and c["gap_ns"] in (8.,40.,200.))]
            rows=[];members=set()
            for i in selected:
                c=cases[i];k=len(c["mother_roots"]);members.update(c["mother_roots"])
                rows.append(dict(c,source_row=i,energy_MeV=[energy[r] for r in c["mother_roots"]],times_ns=times[i,:k].tolist()))
            assert len(rows)==1180 and len(members)==280
            roots.append(members);save(OUTPUT/"calibration"/f"{name}.json",rows)
            np.savez_compressed(OUTPUT/"calibration"/f"{name}.npz",components_mV=components[selected])
        empty=[i for i,c in enumerate(cases) if c["fold"]==1 and c["kind"]=="empty"]
        assert len(empty)==80
        np.savez_compressed(OUTPUT/"calibration/empty.npz",voltage_mV=voltage[empty])
        save(OUTPUT/"calibration/empty.json",[dict(cases[i],source_row=i) for i in empty])
    cases=load(old/"evaluation.cases.json");parents=load(confirm/"physical-parents.json")["records"]
    selected=[i for i,c in enumerate(cases) if (c["condition"],c["gap_ns"]) in CONDITIONS]
    assert len(selected)==480
    evaluation=set();labels=[];plan=[]
    with np.load(old/"evaluation.labels.npz") as truth,np.load(old/"observations.npz") as obs,np.load(confirm/"labels.truth.npz") as native:
        energies=truth["energy_MeV"];times=truth["times_ns"];charges=truth["charge_pe"]
        np.testing.assert_allclose(times[selected],native["times_ns"][selected],equal_nan=True)
        np.testing.assert_allclose(charges[selected],native["charge_pe"][selected],equal_nan=True)
        np.savez_compressed(OUTPUT/"inputs/observations.npz",voltage_mV=obs["voltage_mV"][selected],
            known_times_ns=times[selected],known_counts=np.array([len(cases[i]["parents"]) for i in selected]))
        np.savez_compressed(OUTPUT/"labels.components.npz",components_mV=native["components_mV"][selected])
        for index,i in enumerate(selected):
            c=cases[i];k=len(c["parents"]);rr=[parents[j]["root"] for j in c["parents"]];evaluation.update(rr)
            labels.append(dict(c,id=index,source_row=i,mother_roots=rr,energy_MeV=energies[i,:k].tolist(),
                times_ns=times[i,:k].tolist(),charge_pe=charges[i,:k].tolist(),block=int(truth["block"][i])))
            plan.append(dict(id=index))
    assert len(evaluation)==320 and not roots[0]&roots[1] and not (roots[0]|roots[1])&evaluation
    save(OUTPUT/"labels.json",labels);save(EXPERIMENT/"plan.json",plan)
    save(EXPERIMENT/"mother-split.json",dict(mean=sorted(roots[0]),noise=sorted(roots[1]),evaluation=sorted(evaluation),
        pairwise_disjoint=True,evaluation_historically_exposed=True))
    save(EXPERIMENT/"mother-provenance.json",dict(development=mothers,evaluation=parents))
    pilot=[]
    for condition in (CONDITIONS[0],CONDITIONS[3],CONDITIONS[4],CONDITIONS[5],CONDITIONS[7]):
        pilot.append(next(c["id"] for c in labels if (c["condition"],c["gap_ns"])==condition))
    save(EXPERIMENT/"pilot-ids.json",pilot)
    save(EXPERIMENT/"input-hashes.json",{str(p.relative_to(OUTPUT)):digest(p) for p in OUTPUT.rglob("*") if p.is_file()})
    save(EXPERIMENT/"data-freeze.json",dict(before_calibration=True,files={name:digest(EXPERIMENT/name) for name in
        ("method-freeze.json","input-hashes.json","mother-split.json","mother-provenance.json","source-input-hashes.json","plan.json","pilot-ids.json")}))
    save(RECORD/"environment.json",dict(python=sys.version,executable=sys.executable,workers=6,BLAS_threads=1,
        packages={n:version(n) for n in ("numpy","scipy","matplotlib","nbformat")}))
    print("PREPARED:280/280/320 disjoint mothers;1180 mean,1180+80noise,480evaluation frames",flush=True)


def verify(models=False):
    for name in ("method-freeze.json","data-freeze.json")+( ("model-freeze.json",) if models else () ):
        for rel,value in load(EXPERIMENT/name)["files"].items():assert digest(EXPERIMENT/rel)==value,rel
    for rel,value in load(EXPERIMENT/"input-hashes.json").items():assert digest(OUTPUT/rel)==value,rel


def records(name):
    rows=load(OUTPUT/"calibration"/f"{name}.json")
    with np.load(OUTPUT/"calibration"/f"{name}.npz") as f:components=f["components_mV"]
    return [dict(r,components_mV=components[i,:len(r["energy_MeV"])]) for i,r in enumerate(rows)]


def calibrate():
    verify();assert not (EXPERIMENT/"model-freeze.json").exists()
    start=time.monotonic();base=RecoveryEnergyResponse.load(EXPERIMENT)
    mean,mean_report=calibrate_mean(base,records("mean"))
    mean.save(EXPERIMENT/"mean-correction.npz")
    with np.load(OUTPUT/"calibration/empty.npz") as f:empty=f["voltage_mV"]
    noise,noise_report=calibrate_noise(mean,records("noise"),empty)
    assert noise_report["optimizer_converged"] and np.isfinite(noise["fractional_mode_variance"]).all()
    save(EXPERIMENT/"noise-calibration.json",noise)
    grid=np.geomspace(base.energy_bounds[0],base.energy_bounds[1],80)
    charge=[]
    for e in grid:
        curve,q=mean.components([128.],[e]);assert np.isfinite(curve).all() and (curve>=0).all()
        assert np.isfinite(learned_metric(mean,[128.],[e],noise).logdet)
        charge.append(float(q[0]))
    check=dict(energy_grid_MeV=grid.tolist(),new_isolated_charge_pe=charge,
        monotonic_isolated_charge=bool(np.all(np.diff(charge)>0)),finite_positive_covariance=True)
    assert check["monotonic_isolated_charge"]
    save(EXPERIMENT/"calibration-report.json",dict(mean=mean_report,noise=noise_report,checks=check,seconds=time.monotonic()-start))
    save(EXPERIMENT/"model-freeze.json",dict(before_evaluation_fits=True,files={name:digest(EXPERIMENT/name) for name in
        ("mean-correction.npz","noise-calibration.json","calibration-report.json")}))
    print(json.dumps(dict(mean=mean_report,noise=noise_report,parameters=noise,seconds=time.monotonic()-start),indent=2),flush=True)


def initialize():
    global BASE,MEAN,NOISE,VOLTAGE,TIMES,COUNTS
    BASE=RecoveryEnergyResponse.load(EXPERIMENT)
    MEAN=CalibratedEnergyResponse.load(BASE,EXPERIMENT/"mean-correction.npz")
    NOISE=load(EXPERIMENT/"noise-calibration.json")
    # Fitting may only open this study's waveform inputs, output files and frozen models.
    def audit(event,args):
        if event!="open" or not isinstance(args[0],(str,bytes,os.PathLike)):return
        path=Path(os.fsdecode(args[0])).resolve()
        forbidden=path.suffix.lower()==".root" or path.name.startswith("labels") or path.is_relative_to(OUTPUT/"calibration")
        forbidden |= path.is_relative_to(ART) and not path.is_relative_to(OUTPUT)
        forbidden |= path in (EXPERIMENT/"mother-provenance.json",EXPERIMENT/"mother-split.json")
        if forbidden:raise PermissionError(f"Inference cannot read per-event truth/calibration/old results: {path}")
    sys.addaudithook(audit)
    for p in (OUTPUT/"labels.json",OUTPUT/"calibration/mean.json",ART/"alpha_pair_recovery_energy_20261029/inputs/evaluation.labels.npz"):
        try:p.open("rb")
        except PermissionError:pass
        else:raise AssertionError("Truth-read guard failed")
    with np.load(OUTPUT/"inputs/observations.npz") as f:
        assert set(f.files)=={"voltage_mV","known_times_ns","known_counts"}
        VOLTAGE=f["voltage_mV"];TIMES=f["known_times_ns"];COUNTS=f["known_counts"]


def one(task,pilot=False):
    start=time.monotonic();i=task["id"];t=TIMES[i,:COUNTS[i]]
    fits=compare_known_times(BASE,MEAN,NOISE,VOLTAGE[i],t)
    assert tuple(fits["fits"])==METHODS
    for row in fits["fits"].values():
        assert row["times_ns"]==t.tolist() and np.isfinite(row["energy_MeV"]).all()
    result=dict(id=i,pilot=pilot,seconds=time.monotonic()-start,**fits,
        worker_truth_read_guard=True,known_K_t=True,true_E_Q_inference=False)
    save(OUTPUT/("pilot" if pilot else "cases")/f"{i:04}.json",result)
    return dict(id=i,seconds=result["seconds"])


def pilot():
    verify(True);initialize()
    results=[one(dict(id=i),True) for i in load(EXPERIMENT/"pilot-ids.json")]
    estimate=max(r["seconds"] for r in results)*480/6*1.5
    record=dict(passed=estimate<=1800,results=results,conservative_wall_seconds=estimate,workers=6,
        cost_only=True,accuracy_not_used=True)
    save(RECORD/"pilot-validation.json",record);print(json.dumps(record,indent=2),flush=True)
    assert record["passed"],"Cost gate exceeded; preserve original before revising execution"


def run(workers):
    verify(True);assert 1<=workers<=6 and load(RECORD/"pilot-validation.json")["passed"]
    started=time.monotonic();remaining=[]
    for task in load(EXPERIMENT/"plan.json"):
        p=OUTPUT/"cases"/f"{task['id']:04}.json"
        if p.exists():
            row=load(p);assert row["id"]==task["id"] and not row["pilot"] and tuple(row["fits"])==METHODS
        else:remaining.append(task)
    with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context("spawn"),initializer=initialize) as pool:
        futures=[pool.submit(one,task) for task in remaining]
        for i,future in enumerate(as_completed(futures),1):
            row=future.result()
            if i%20==0 or i==len(remaining):print(json.dumps(dict(done=i,total=len(remaining),**row)),flush=True)
    files={str(p.relative_to(OUTPUT)):digest(p) for p in sorted((OUTPUT/"cases").glob("*.json"))}
    assert len(files)==480
    save(RECORD/"result-freeze.json",dict(files=files,model_freeze_sha256=digest(EXPERIMENT/"model-freeze.json")))
    save(RECORD/"run-validation.json",dict(passed=True,frames=480,method_rows=1920,newly_run=len(remaining),
        workers=workers,seconds=time.monotonic()-started,known_K_t=True,true_E_Q_inference=False,default_changed=False))
    print("COMPLETE:480 frames,1920 method rows",flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument("action",choices=("prepare","calibrate","pilot","run"))
    parser.add_argument("--workers",type=int,default=6);args=parser.parse_args()
    if args.action=="run":run(args.workers)
    else:globals()[args.action]()


if __name__=="__main__":main()
