"""Known-input shared-cell mean-response test following isolated controls.

Truth E/t and photon ancestry are used for forward-model diagnostics only.
No energy fit, event detection or new transport validation is implied here.
"""
from __future__ import annotations

import argparse
import ast
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np

from .alpha_pair_causal_controls import compatibility, read, save, sha
from .alpha_pair_causal_recovery import CausalRecovery, DifferentiableAge, make_light_charge_map
from .alpha_pair_digitizer_optimized import simulate_avalanches
from .alpha_pair_energy_calibrated import CalibratedEnergyResponse
from .alpha_pair_recovery_energy import RecoveryEnergyResponse
from .sipm_response import SensorConfig, PhotonHits, AvalancheResult, map_cells
from .sipm_readout import ReadoutConfig, synthesize_readout

NAME = "alpha_pair_causal_shared_20261101"
ROOT = Path(os.environ.get("GAGG_CAUSAL_SHARED_ROOT", Path(__file__).resolve().parents[2])).resolve()
EXP = ROOT / "pileup/experiments" / NAME
ART = ROOT / "jupyter/artifacts" / NAME
REC = ROOT / "jupyter/run_records" / NAME
PREV = ROOT / "pileup/experiments/alpha_pair_causal_controls_20261101"
PREVART = ROOT / "jupyter/artifacts/alpha_pair_causal_controls_20261101"


def prepare():
    assert read(PREV / "single-gate.json")["single_gate_passed"]
    for folder in (EXP, ART, REC): folder.mkdir(parents=True, exist_ok=False)
    for folder in (ART / "photons", ART / "cases", EXP / "source/pileup/analysis", EXP / "bin"):
        folder.mkdir(parents=True)
    for name in ("template.npz", "energy-model.json", "response-parameters.json", "mean-correction.npz",
                 "noise-calibration.json", "sensor-config.json"):
        shutil.copyfile(PREV / name, EXP / name)
    donors = [p for p in read(PREV / "mother-provenance.json") if p["condition"] in ("weak", "high")]
    for donor in donors:
        name = f"{donor['parent']:04d}.npz"
        shutil.copyfile(PREVART / "photons" / name, ART / "photons" / name)
    save(EXP / "mother-provenance.json", donors)
    tasks = []
    for block in (0, 10, 20, 30):
        for condition, gap in (("high_weak",8.),("high_weak",40.),("high_weak",200.),("high_weak",1200.),("weak_high",40.)):
            parents = [block+240, block] if condition=="high_weak" else [block,block+240]
            tasks.append(dict(id=len(tasks), block=block, condition=condition, gap_ns=gap,
                parents=parents, times_ns=[128.,128.+gap], repeats=8,
                seed_base=712000000+block*1000, window_ns=4096.))
    save(EXP / "plan.json", tasks)
    save(EXP / "protocol.json", dict(created=datetime.now(timezone.utc).isoformat(),
        prerequisite="Isolated causal controls passed; original artifacts/hashes retained",
        commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        scope="Known-input forward diagnostic on eight historically exposed mothers, four mother pairs; not fitting or new transport confirmation",
        conditions="high then weak:8/40/200/1200ns; weak then high:40ns; four pairs, eight seeds",
        runs="Each case/seed:joint and two isolated digitizations, with recovery on and off; 960 digitizations, plus a separate two-repeat cost pilot",
        rng="Common seed by mother pair/repeat across controls; changed photon lists change draw consumption; no claim of identical avalanche draws",
        means="Ancestry event signals only, dark ancestry excluded but dark occupancy remains; observed first2048ns aperture averages and full4096ns ancestry charge both saved",
        model="Uniform14400-cell age expectation, fixed95/351ns light, fraction0.3; frozen isolated-voltage anchor plus delta(shared-isolated); strength1, no tuning",
        approximation="No per-event spatial truth is supplied to the uniform model. Cell-overlap calculation is post hoc mechanistic evidence only",
        grids_ns=[1.,.5],primary_grid_ns=1.,
        primary="Weak-component integrated-voltage suppression, relative to its isolated response in the SAME2048ns window",
        gates=dict(weak_suppression_error_percentage_points=5.,grid_change_percentage_points=1.,
                   recovery_off_interaction_percentage_points=3.,pilot_seconds=900.),
        gate_scope="Prespecified engineering point-estimate gates, not literature limits or independent statistical confirmation",
        expansion="Only after all response gates pass, design an energy fit comparison; otherwise stop and diagnose model assumptions",
        workers=4,BLAS_threads=1,no_new_ROOT=True,default_changed=False))
    source = EXP / "source/pileup/analysis"
    pending=["alpha_pair_causal_shared_study"];copied=set()
    while pending:
        name=pending.pop()
        if name in copied:continue
        p=ROOT/"pileup/analysis"/(name+".py");tree=ast.parse(p.read_text())
        shutil.copyfile(p,source/p.name);copied.add(name)
        pending.extend(n.module for n in ast.walk(tree) if isinstance(n,ast.ImportFrom) and n.level==1 and n.module)
    (source.parent/"__init__.py").write_text("");(source/"__init__.py").write_text("")
    for name in ("sipm_age_mean.cpp","sipm_age_mean_jacobian.cpp"):
        shutil.copyfile(ROOT/"pileup/analysis"/name,source/name)
    (EXP/"run.py").write_text("from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\n"
        "os.environ['GAGG_CAUSAL_SHARED_ROOT']=str(f.parents[2])\nsys.path.insert(0,str(f/'source'))\n"
        "if __name__=='__main__':\n    from pileup.analysis.alpha_pair_causal_shared_study import main\n    main()\n")
    command=["g++","-O3","-std=c++17","-fPIC","-shared",str(source/"sipm_age_mean.cpp"),str(source/"sipm_age_mean_jacobian.cpp"),"-o",str(EXP/"bin/age.so")]
    subprocess.run(command,check=True)
    save(EXP/"build.json",dict(command=command,compiler=subprocess.check_output(["g++","--version"],text=True).splitlines()[0],binary_sha256=sha(EXP/"bin/age.so")))
    initialize(False)
    for dt in (1.,.5):
        age=DifferentiableAge(CONFIG["sensor"],EXP/"bin/age.so",dt)
        np.savez_compressed(EXP/f"charge-map-{dt}.npz",**make_light_charge_map(age,ANCHOR))
    files=[p for p in EXP.rglob("*") if p.is_file() and "bin" not in p.relative_to(EXP).parts]
    save(EXP/"freeze.json",{str(p.relative_to(ROOT)):sha(p) for p in files})
    save(EXP/"input-freeze.json",{str(p.relative_to(ROOT)):sha(p) for p in (ART/"photons").glob("*.npz")})
    print("Prepared20 shared controls;960 digitizations; no fitting",flush=True)


def verify():
    for name in ("freeze.json","input-freeze.json"):
        for p,h in read(EXP/name).items():assert sha(ROOT/p)==h,p


def initialize(models=True):
    compatibility()
    global CONFIG,SENSOR,READOUT,DONORS,ANCHOR,MODELS
    CONFIG=read(EXP/"sensor-config.json")
    SENSOR=SensorConfig(**CONFIG["sensor"]);READOUT=ReadoutConfig(**CONFIG["readout"])
    ANCHOR=CalibratedEnergyResponse.load(RecoveryEnergyResponse.load(EXP),EXP/"mean-correction.npz")
    DONORS={}
    for p in (ART/"photons").glob("*.npz"):
        with np.load(p) as f:DONORS[int(p.stem)]={k:f[k].copy() for k in f.files}
    MODELS={}
    if models:
        for dt in (1.,.5):
            with np.load(EXP/f"charge-map-{dt}.npz") as f:mapping={k:f[k].copy() for k in f.files}
            MODELS[dt]=CausalRecovery(ANCHOR,CONFIG,EXP/"bin/age.so",mapping)


def digitize(parents,times,seed,recovery):
    fields={k:[] for k in ("time_ns","x_mm","y_mm","energy_eV")};owners=[]
    for j,(parent,t0) in enumerate(zip(parents,times)):
        donor=DONORS[parent];t=donor["time_ns"]+t0;keep=t<4096.
        for k in fields:fields[k].append((t if k=="time_ns" else donor[k])[keep])
        owners.append(np.full(int(keep.sum()),j))
    arrays={k:np.concatenate(v) for k,v in fields.items()};owner=np.concatenate(owners)
    order=np.argsort(arrays["time_ns"],kind="stable");owner=owner[order]
    hits=PhotonHits(**{k:v[order] for k,v in arrays.items()},window_ns=4096.)
    sensor=replace(SENSOR,recovery_enabled=recovery)
    av=simulate_avalanches(hits,sensor,seed);rec=av.records
    assignment=np.full(len(rec),-1,int);direct=rec["kind"]==0
    assignment[direct]=owner[rec["photon_index"][direct]]
    for j in np.flatnonzero(rec["parent_id"]>=0):assignment[j]=assignment[rec["parent_id"][j]]
    signals=[];charges=[]
    for j in range(len(parents)):
        subset=rec[assignment==j]
        event=AvalancheResult(subset,4096.,av.full_charge_C,np.array([]),{})
        signals.append(synthesize_readout(event,READOUT,0).arrays["signal_mV"][:1024])
        charges.append(float(subset["charge_C"].sum()/av.full_charge_C))
    return np.array(signals),np.array(charges)


def one(task):
    start=time.monotonic();result={}
    for recovery in (True,False):
        shared=[];isolated=[];qs=[];qi=[]
        for repeat in range(task["repeats"]):
            seed=task["seed_base"]+repeat
            curves,q=digitize(task["parents"],task["times_ns"],seed,recovery)
            singles=[digitize([p],[t],seed,recovery) for p,t in zip(task["parents"],task["times_ns"])]
            shared.append(curves);qs.append(q)
            isolated.append(np.concatenate([x[0] for x in singles]));qi.append(np.concatenate([x[1] for x in singles]))
        prefix="on" if recovery else "off"
        result[prefix+"_shared"]=np.asarray(shared);result[prefix+"_isolated"]=np.asarray(isolated)
        result[prefix+"_shared_charge"]=np.asarray(qs);result[prefix+"_isolated_charge"]=np.asarray(qi)
    return dict(id=task["id"],seconds=time.monotonic()-start,digitizations=task["repeats"]*6),result


def run(pilot=False):
    verify();tasks=read(EXP/"plan.json")
    if pilot:tasks=[dict(tasks[1],repeats=2)]
    else:assert read(REC/"pilot.json")["passed"]
    folder=ART/("pilot" if pilot else "cases");folder.mkdir(exist_ok=True)
    assert not list(folder.glob("*.npz")),"Keep previous output"
    start=time.monotonic()
    with ProcessPoolExecutor(max_workers=4,initializer=initialize) as pool:
        for i,(row,arrays) in enumerate(pool.map(one,tasks),1):
            np.savez_compressed(folder/f"{row['id']:04d}.npz",**arrays);save(folder/f"{row['id']:04d}.json",row)
            if pilot or i%4==0:print("shared controls",i,"/",len(tasks),"elapsed",round(time.monotonic()-start,1),flush=True)
    if pilot:
        estimate=row["seconds"]*4*20/4*1.5
        save(REC/"pilot.json",dict(passed=estimate<900.,projected_seconds=estimate,limit_seconds=900.,seconds=row["seconds"]))
    else:save(REC/"run.json",dict(cases=20,digitizations=960,seconds=time.monotonic()-start,workers=4))


def report():
    verify();initialize()
    provenance={r["parent"]:r for r in read(EXP/"mother-provenance.json")};rows=[];plots={}
    for task in read(EXP/"plan.json"):
        energies=[provenance[p]["energy_dep_MeV"] for p in task["parents"]]
        weak=task["parents"].index(task["block"])
        with np.load(ART/"cases"/f"{task['id']:04d}.npz") as f:
            shared=f["on_shared"].mean(0);isolated=f["on_isolated"].mean(0)
            offshared=f["off_shared"].mean(0);offisolated=f["off_isolated"].mean(0)
        actual=100*(shared[weak].sum()/isolated[weak].sum()-1)
        off=100*(offshared[weak].sum()/offisolated[weak].sum()-1)
        predictions={};relative={}
        for dt,model in MODELS.items():
            solo=model.components(task["times_ns"],energies,interaction=False)
            joint=model.components(task["times_ns"],energies)
            relative[dt]=float(100*(joint[weak].sum()/solo[weak].sum()-1))
            predictions[dt]=(solo,joint)
        # Independent-mother spatial overlap. This uses photon truth only for
        # diagnostics; the tested uniform model receives E/t, not these maps.
        distributions=[]
        for p in task["parents"]:
            h=PhotonHits(**DONORS[p],window_ns=20000.)
            cells=map_cells(h,SENSOR);counts=np.bincount(cells[cells>=0],minlength=SENSOR.n_cells).astype(float)
            distributions.append(counts/counts.sum())
        overlap=float(SENSOR.n_cells*np.dot(*distributions))
        row=dict(**task,energy_MeV=energies,actual_weak_area_change_percent=float(actual),
            predicted_weak_area_change_percent=relative[1.],prediction_error_percentage_points=relative[1.]-actual,
            fine_grid_prediction_percent=relative[.5],grid_change_percentage_points=relative[1.]-relative[.5],
            recovery_off_area_change_percent=float(off),spatial_overlap_relative_to_uniform=overlap,
            model_min_component_mV=float(predictions[1.][1].min()))
        rows.append(row)
        if task["block"]==0:
            plots[f"case{task['id']}_shared"]=shared;plots[f"case{task['id']}_isolated"]=isolated
            plots[f"case{task['id']}_model_shared"]=predictions[1.][1];plots[f"case{task['id']}_model_isolated"]=predictions[1.][0]
    groups=[]
    for condition,gap in (("high_weak",8.),("high_weak",40.),("high_weak",200.),("high_weak",1200.),("weak_high",40.)):
        group=[r for r in rows if (r["condition"],r["gap_ns"])==(condition,gap)]
        keys=["actual_weak_area_change_percent","predicted_weak_area_change_percent","prediction_error_percentage_points",
              "grid_change_percentage_points","recovery_off_area_change_percent","spatial_overlap_relative_to_uniform"]
        groups.append(dict(condition=condition,gap_ns=gap,mother_pairs=4,**{k:float(np.mean([r[k] for r in group])) for k in keys}))
    gates=read(EXP/"protocol.json")["gates"];failures=[]
    for group in groups:
        name=group["condition"]+"_"+str(int(group["gap_ns"]))
        if abs(group["prediction_error_percentage_points"])>gates["weak_suppression_error_percentage_points"]:failures.append(name+": shared response mismatch")
        if abs(group["grid_change_percentage_points"])>gates["grid_change_percentage_points"]:failures.append(name+": grid change")
        if abs(group["recovery_off_area_change_percent"])>gates["recovery_off_interaction_percentage_points"]:failures.append(name+": recovery-off control")
    if any(r["model_min_component_mV"] < -1e-10 for r in rows):failures.append("negative predicted component")
    summary=dict(groups=groups,response_gate_passed=not failures,failures=failures,
        inference_performed=False,default_changed=False,no_new_ROOT=True,
        limitations="Four optical pairs reused across gaps, orders and eight digitizer seeds; known-input diagnostics, engineering point-estimate gates; spatial overlap is explanatory evidence, not a fitted correction or proof of unique cause")
    save(EXP/"response-gate.json",summary);save(ART/"report.json",summary);save(ART/"rows.json",rows)
    np.savez_compressed(ART/"display-curves.npz",**plots)
    import csv
    (EXP/"tables").mkdir(exist_ok=True)
    with (EXP/"tables/shared-controls.csv").open("w") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)


def main():
    compatibility();parser=argparse.ArgumentParser();parser.add_argument("mode",choices=("prepare","pilot","run","report"))
    mode=parser.parse_args().mode
    if mode=="prepare":prepare()
    elif mode=="pilot":run(True)
    elif mode=="run":run(False)
    else:report()


if __name__=="__main__":main()
