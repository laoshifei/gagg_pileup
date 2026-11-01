"""Staged causal-response controls on existing optical mothers, not blind fits.

Run from the project root. The prepared experiment contains the actual sources,
configuration and input hashes. A failed isolated-response gate stops expansion
to pairs. Repeated digitizations do not add independent optical mothers.
"""
from __future__ import annotations

import argparse
import ast
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
from scipy.optimize import minimize_scalar

from .alpha_pair_causal_recovery import CausalRecovery, CausalReadout, DifferentiableAge, make_light_charge_map
from .alpha_pair_deconvolution import AgeRecovery
from .alpha_pair_digitizer_optimized import simulate_avalanches
from .alpha_pair_energy_calibrated import CalibratedEnergyResponse, learned_metric
from .alpha_pair_recovery_energy import RecoveryEnergyResponse
from .sipm_readout import ReadoutConfig, synthesize_readout
from .sipm_response import AvalancheResult, PhotonHits, SensorConfig

NAME = "alpha_pair_causal_controls_20261101"
ROOT = Path(os.environ.get("GAGG_CAUSAL_CONTROLS_ROOT", Path(__file__).resolve().parents[2])).resolve()
EXP = ROOT / "pileup/experiments" / NAME
ART = ROOT / "jupyter/artifacts" / NAME
REC = ROOT / "jupyter/run_records" / NAME
OLD = ROOT / "pileup/experiments/alpha_pair_recovery_forward_20261030"
PARENTS = ROOT / "jupyter/artifacts/alpha_pair_joint_confirmation_20261026"
HISTORY = ROOT / "jupyter/artifacts/sipm_energy_history_20261029/inputs"
BLOCKS = (0, 10, 20, 30)
PHASES = (0., .5, 1., 1.5)
NUMPY_TRAPEZOID_ADAPTER = not hasattr(np, "trapezoid")


def read(path):
    return json.loads(Path(path).read_text())


def save(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def compatibility():
    # Explicit runner-local adapter for the NumPy 2 rename. Frozen study files
    # stay unchanged. SciPy/NumPy optimizer differences are checked by replay.
    needed = not hasattr(np, "trapezoid")
    if needed:
        np.trapezoid = np.trapz
    assert np.trapezoid(np.array([0., 1., 2.]), np.array([0., 1., 2.])) == 2.
    return NUMPY_TRAPEZOID_ADAPTER


def snapshot():
    source = EXP / "source/pileup/analysis"
    source.mkdir(parents=True)
    (source.parent / "__init__.py").write_text("")
    (source / "__init__.py").write_text("")
    pending = ["alpha_pair_causal_controls"]
    copied = set()
    while pending:
        name = pending.pop()
        if name in copied:
            continue
        path = ROOT / "pileup/analysis" / (name + ".py")
        tree = ast.parse(path.read_text())
        shutil.copyfile(path, source / path.name)
        copied.add(name)
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.level == 1 and node.module:
                pending.append(node.module)
    for name in ("sipm_age_mean.cpp", "sipm_age_mean_jacobian.cpp"):
        shutil.copyfile(ROOT / "pileup/analysis" / name, source / name)
    (EXP / "run.py").write_text(
        "from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\n"
        "os.environ['GAGG_CAUSAL_CONTROLS_ROOT']=str(f.parents[2])\n"
        "sys.path.insert(0,str(f/'source'))\n"
        "if __name__=='__main__':\n"
        "    from pileup.analysis.alpha_pair_causal_controls import main\n    main()\n")


def prepare():
    for folder in (EXP, ART, REC):
        folder.mkdir(parents=True, exist_ok=False)
    (ART / "photons").mkdir()
    (ART / "single").mkdir()
    for name in ("template.npz", "energy-model.json", "response-parameters.json",
                 "mean-correction.npz", "noise-calibration.json", "sensor-config.json"):
        shutil.copyfile(OLD / name, EXP / name)
    parents = read(HISTORY / "parents.json")["records"]
    cases = read(HISTORY / "cases.json")
    with np.load(HISTORY / "labels.npz") as data:
        energies = data["energy_MeV"]
    truth = {c["parents"][0]: float(energies[i, 0]) for i, c in enumerate(cases) if c["kind"] == "single"}
    provenance, tasks = [], []
    for group, condition in enumerate(("weak", "low", "mid", "high")):
        for block in BLOCKS:
            parent = group * 80 + block
            src = PARENTS / "photons" / f"{parent:04d}.npz"
            dest = ART / "photons" / src.name
            shutil.copyfile(src, dest)
            provenance.append(dict(parent=parent, condition=condition, block=block,
                energy_dep_MeV=truth[parent], source=parents[parent], photon_sha256=sha(dest)))
            for phase in PHASES:
                tasks.append(dict(id=len(tasks), parent=parent, condition=condition,
                    block=block, phase_ns=phase, time_ns=128. + phase,
                    energy_dep_MeV=truth[parent], repeats=8,
                    seed_base=711000000 + parent * 100, window_ns=4096.))
    save(EXP / "single-plan.json", tasks)
    save(EXP / "mother-provenance.json", provenance)
    save(EXP / "protocol.json", dict(created=datetime.now(timezone.utc).isoformat(),
        commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        scope="Development causal/phase/mean control; 16 historically exposed optical mothers; no new G4, blind detection or independent transport confirmation",
        stages=["old three-frame replay", "C++ and readout numerical controls", "isolated physical mean controls", "pairs only after the isolated gate"],
        sensor="Unchanged S13360-6050PE proxy configuration", recovery_grid_ns=1.,
        optical_mothers=16, phases_ns=list(PHASES), repetitions_per_phase=8,
        digitizations=512, workers=4, BLAS_threads=1,
        clock="Earliest accepted-area optical arrival aligned to 128+phase ns, as in the historical donor cache",
        waveform="4096ns digitization; first 2048ns/1024 native aperture averages are the explicitly declared observed control window; tail charge recorded separately",
        phase_rng="Reuse each mother's eight digitizer seeds across phases; ordering and dark competition can change. Mothers, not repeats/phases, are independent sampling units",
        ancestry="Truth ancestry only separates event signal from dark signal for response validation. Full noisy analog voltage is retained",
        model="Frozen calibrated isolated voltage with causal support; uniform-age correction is not exercised for K=1 and is not established by a successful single test",
        fitting="Single-event E fit to each eight-repeat mean; known onset; frozen legacy covariance from voltage-only seed, no true E initialization",
        gates=dict(replay_energy_relative=2e-4, replay_cost_relative=1e-6,
            mean_energy_bias_percent=5., phase_bias_span_percentage_points=1.,
            pre_event_signal_mV=1e-10, pilot_projected_seconds=900.),
        gate_status="Prespecified engineering controls, not literature thresholds or significance tests",
        selection="Every parent/phase retained; failure stops pair expansion. No calibration strength tuning on these mothers",
        next_pairs="If singles pass, freeze a separate shared-versus-isolated response experiment before energy comparison",
        default_changed=False))
    snapshot()
    (EXP / "bin").mkdir()
    sources = [EXP / "source/pileup/analysis" / name for name in ("sipm_age_mean.cpp", "sipm_age_mean_jacobian.cpp")]
    command = ["g++", "-O3", "-std=c++17", "-fPIC", "-shared", *map(str, sources), "-o", str(EXP / "bin/age.so")]
    subprocess.run(command, check=True)
    save(EXP / "build.json", dict(command=command, compiler=subprocess.check_output(["g++", "--version"], text=True).splitlines()[0],
        source_hashes={str(p.relative_to(ROOT)): sha(p) for p in sources}, binary_sha256=sha(EXP / "bin/age.so"), fast_math=False))
    files = [p for p in EXP.rglob("*") if p.is_file() and "bin" not in p.relative_to(EXP).parts]
    save(EXP / "source-config-freeze.json", {str(p.relative_to(ROOT)): sha(p) for p in files})
    save(EXP / "input-freeze.json", {str(p.relative_to(ROOT)): sha(p) for p in (ART / "photons").glob("*.npz")})
    import scipy
    save(REC / "environment.json", dict(python=sys.version, executable=sys.executable,
        numpy=np.__version__, scipy=scipy.__version__, numpy_trapezoid_rename_adapter=compatibility()))
    print("Prepared 16 mothers, 64 parent/phase cases, 512 digitizations", flush=True)


def verify():
    for name in ("source-config-freeze.json", "input-freeze.json"):
        for path, expected in read(EXP / name).items():
            assert sha(ROOT / path) == expected, path


def baseline():
    verify()
    # Separate interpreter imports the *previous* frozen sources, without
    # modifying them or rebuilding their executable in place.
    program = '''import os,sys,json,numpy as np
from pathlib import Path
root=Path(sys.argv[1]); destination=Path(sys.argv[2])
if not hasattr(np,'trapezoid'):np.trapezoid=np.trapz
os.environ['GAGG_FORWARD_PROJECT']=str(root)
sys.path.insert(0,str(root/'pileup/experiments/alpha_pair_recovery_forward_20261030/source'))
from pileup.analysis import alpha_pair_recovery_forward_study as s
with np.load(s.ART/'inputs/observations.npz') as f:
    items=[(i,v,t[:k]) for i,v,t,k in zip(f['ids'],f['voltage_mV'],f['known_times_ns'],f['known_counts']) if int(i) in (40,160,400)]
references={int(i):s.read(s.ART/'cases'/f'{int(i):04d}.json')['fits'] for i,_,_ in items}
s.initialize(True)
rows=[]
for item in items:
    result=s.job(item); reference=references[int(item[0])]
    for method,fit in result['fits'].items():
        old=reference[method]
        rows.append(dict(id=int(item[0]),method=method,
            energy_relative=float(np.max(abs(np.array(fit['energy_MeV'])/old['energy_MeV']-1))),
            cost_relative=float(abs(fit['cost']/old['cost']-1)),
            convergence_equal=fit['converged']==old['converged']))
destination.write_text(json.dumps(dict(rows=rows,frames=len(items),truth_read_guard=True),indent=2)+'\\n')
'''
    assert not (REC / "baseline.json").exists()
    subprocess.run([sys.executable, "-B", "-c", program, str(ROOT), str(REC / "baseline.json")], check=True)
    report = read(REC / "baseline.json")
    gates = read(EXP / "protocol.json")["gates"]
    report["passed"] = all(r["energy_relative"] <= gates["replay_energy_relative"] and
        r["cost_relative"] <= gates["replay_cost_relative"] and r["convergence_equal"] for r in report["rows"])
    save(REC / "baseline.json", report)
    print("Baseline replay:", report["passed"], "max relative E", max(r["energy_relative"] for r in report["rows"]), flush=True)


def initialize():
    compatibility()
    global SENSOR, READOUT, ANCHOR, MODEL, NOISE, DONORS
    config = read(EXP / "sensor-config.json")
    SENSOR, READOUT = SensorConfig(**config["sensor"]), ReadoutConfig(**config["readout"])
    ANCHOR = CalibratedEnergyResponse.load(RecoveryEnergyResponse.load(EXP), EXP / "mean-correction.npz")
    NOISE = read(EXP / "noise-calibration.json")
    # K=1 does not use the charge map or interaction; load the frozen draft map
    # only for constructing the common interface, not to calibrate any labels.
    mapping = dict(light_pe=np.array([0., 50000.]), charge_pe=np.array([0., 50000.]), dt_ns=1., fast_fraction=.3)
    MODEL = CausalRecovery(ANCHOR, config, EXP / "bin/age.so", mapping)
    DONORS = {}
    for path in (ART / "photons").glob("*.npz"):
        with np.load(path) as f:
            DONORS[int(path.stem)] = {k: f[k].copy() for k in f.files}


def numerical():
    verify()
    initialize()
    config = read(EXP / "sensor-config.json")
    rng = np.random.default_rng(71101001)
    comparisons, gradients, readout_errors = [], [], []
    for dt in (.5, 1., 2.):
        age = DifferentiableAge(config["sensor"], EXP / "bin/age.so", dt)
        slow = AgeRecovery(config["sensor"], dt)
        x = rng.uniform(.1, 200., 160)
        directions = np.vstack([x, x * np.linspace(-.3, .3, len(x))])
        y, grad = age.with_gradient(x, directions)
        comparisons.append(float(abs(y - slow.forward(x)).max()))
        for d, g in zip(directions, grad):
            fd = (age.forward(x + 1e-5*d) - age.forward(x - 1e-5*d)) / 2e-5
            error = float(np.linalg.norm(g-fd) / max(np.linalg.norm(fd), 1e-12))
            gradients.append(error)
        readout = CausalReadout(config, dt)
        released = np.zeros(readout.n)
        released[[31, 62, 170]] = [1., 3., 2.]
        records = np.array([((i+1)*dt, q*SENSOR.full_charge_C) for i, q in enumerate(released) if q],
                           dtype=[("time_ns", "f8"), ("charge_C", "f8")])
        exact = synthesize_readout(AvalancheResult(records, 2048., SENSOR.full_charge_C, np.array([]), {}), READOUT, 0)
        readout_errors.append(float(abs(readout(released) - exact.arrays["signal_mV"]).max()))
    pre = []
    for e in (.25, 1.5, 4., 7.5):
        for phase in PHASES:
            t = 128. + phase
            curve = MODEL.isolated([t], [e])[0][0]
            pre.append(float(abs(curve[MODEL.readout.adc_centers + .25 <= t]).max()))
    report = dict(cpp_python_max_abs=max(comparisons), gradient_max_relative=max(gradients),
        native_readout_max_abs_mV=max(readout_errors), pre_event_max_mV=max(pre),
        controls="Discrete age sensitivities; direct original readout at exact right-edge charge times; four energies and four phases",
        limitations="Readout equivalence does not establish coarse-grid timing accuracy for distributed optical input; K=1 does not exercise shared-cell correction")
    report["passed"] = bool(max(comparisons)<1e-9 and max(gradients)<1e-5 and max(readout_errors)<1e-10 and max(pre)<1e-10)
    save(REC / "numerical.json", report)
    print(json.dumps(report), flush=True)


def single(task):
    assert "energy_dep_MeV" not in task, "Evaluation energy labels must stay in the parent/report"
    start = time.monotonic()
    donor = DONORS[task["parent"]]
    shifted = donor["time_ns"] + task["time_ns"]
    keep = shifted < task["window_ns"]
    hits = PhotonHits(shifted[keep], donor["x_mm"][keep], donor["y_mm"][keep], donor["energy_eV"][keep], task["window_ns"])
    physical, signals, charges = [], [], []
    for repeat in range(task["repeats"]):
        seed = task["seed_base"] + repeat
        av = simulate_avalanches(hits, SENSOR, seed)
        readout = synthesize_readout(av, READOUT, seed + 10000000)
        records = av.records
        is_event = records["kind"] == 0
        for j in np.flatnonzero(records["parent_id"] >= 0):
            is_event[j] = is_event[records["parent_id"][j]]
        event = AvalancheResult(records[is_event], task["window_ns"], av.full_charge_C, np.array([]), {})
        component = synthesize_readout(event, READOUT, 0).arrays["signal_mV"]
        physical.append(readout.arrays["analog_mV"][:1024])
        signals.append(component[:1024])
        charges.append(float(records[is_event]["charge_C"].sum() / av.full_charge_C))
    arrays = dict(analog_mV=np.asarray(physical), event_signal_mV=np.asarray(signals), event_charge_pe=np.asarray(charges))
    # Only voltage and declared onset enter this estimate; E truth is scored by
    # the parent/report later. Ancestry signal above is for physical diagnostics.
    voltage = arrays["analog_mV"].mean(0)
    initial = MODEL.initial_energy(voltage, [task["time_ns"]])
    metric = learned_metric(ANCHOR, [task["time_ns"]], initial, NOISE)
    target = metric.apply(ANCHOR.transform(voltage - MODEL.background))
    def objective(loge):
        curve = MODEL.isolated([task["time_ns"]], [np.exp(loge)])[0][0]
        residual = metric.apply(ANCHOR.transform(curve)) - target
        return float(residual @ residual)
    lo, hi = np.log(ANCHOR.energy_bounds)
    fit = minimize_scalar(objective, bounds=(lo, hi), method="bounded", options={"xatol":1e-10, "maxiter":150})
    result = dict(id=task["id"], parent=task["parent"], fitted_energy_MeV=float(np.exp(fit.x)),
        cost=float(fit.fun), converged=bool(fit.success), seconds=time.monotonic()-start,
        covariance_seed_energy_MeV=initial.tolist(), true_energy_input_to_fit=False,
        dropped_photons=int((~keep).sum()))
    return result, arrays


def simulate(pilot=False):
    verify()
    assert read(REC / "baseline.json")["passed"], "Baseline replay gate failed"
    assert read(REC / "numerical.json")["passed"], "Numerical gate failed"
    tasks = read(EXP / "single-plan.json")
    if pilot:
        selected = [dict(tasks[i], repeats=2) for i in (0, 48)]
        folder = ART / "pilot"
    else:
        assert read(REC / "pilot.json")["passed"]
        selected = tasks
        folder = ART / "single"
    folder.mkdir(exist_ok=True)
    assert not list(folder.glob("*.npz")), "Keep previous outputs"
    start = time.monotonic()
    with ProcessPoolExecutor(max_workers=4, initializer=initialize) as pool:
        inference_tasks = [{k: v for k, v in task.items() if k != "energy_dep_MeV"} for task in selected]
        for i, (result, arrays) in enumerate(pool.map(single, inference_tasks), 1):
            np.savez_compressed(folder / f"{result['id']:04d}.npz", **arrays)
            save(folder / f"{result['id']:04d}.json", result)
            if pilot or i % 8 == 0:
                print("isolated controls", i, "/", len(selected), "elapsed", round(time.monotonic()-start, 1), flush=True)
    if pilot:
        times = [read(folder / f"{task['id']:04d}.json")["seconds"] for task in selected]
        estimate = max(times)*4*64/4*1.5
        save(REC / "pilot.json", dict(passed=estimate<900., projected_seconds=estimate, limit_seconds=900.,
            purpose="Cost only; pilot uses two repeats, final planned controls use eight", seconds=times))
    else:
        save(REC / "single-run.json", dict(cases=64, digitizations=512, optical_mothers=16, seconds=time.monotonic()-start, workers=4))


def report():
    verify()
    initialize()
    rows = []
    tasks = read(EXP / "single-plan.json")
    for task in tasks:
        result = read(ART / "single" / f"{task['id']:04d}.json")
        with np.load(ART / "single" / f"{task['id']:04d}.npz") as data:
            event_mean = data["event_signal_mV"].mean(0)
            charge = float(data["event_charge_pe"].mean())
        curve = MODEL.isolated([task["time_ns"]], [task["energy_dep_MeV"]])[0][0]
        mask = MODEL.readout.adc_centers + .25 <= task["time_ns"]
        rows.append(dict(**task, fitted_energy_MeV=result["fitted_energy_MeV"], converged=result["converged"],
            energy_error_percent=100*(result["fitted_energy_MeV"]/task["energy_dep_MeV"]-1),
            model_area_error_percent=100*(curve.sum()/event_mean.sum()-1),
            normalized_mean_error=float(np.linalg.norm(curve-event_mean)/np.linalg.norm(event_mean)),
            pre_event_model_mV=float(abs(curve[mask]).max()),
            pre_event_physical_mV=float(abs(event_mean[mask]).max()), mean_charge_pe=charge))
    groups = []
    for condition in ("weak", "low", "mid", "high"):
        selected = [r for r in rows if r["condition"] == condition]
        phases = [float(np.mean([r["energy_error_percent"] for r in selected if r["phase_ns"]==phase])) for phase in PHASES]
        groups.append(dict(condition=condition, mothers=4, parent_phase_cases=16,
            energy_bias_percent=float(np.mean([r["energy_error_percent"] for r in selected])),
            phase_bias_percent=phases, phase_bias_span_percentage_points=float(np.ptp(phases)),
            model_area_error_percent=float(np.mean([r["model_area_error_percent"] for r in selected])),
            mean_waveform_relative_error=float(np.mean([r["normalized_mean_error"] for r in selected]))))
    gates = read(EXP / "protocol.json")["gates"]
    failures = []
    for group in groups:
        if abs(group["energy_bias_percent"])>gates["mean_energy_bias_percent"]:
            failures.append(group["condition"]+": isolated energy bias")
        if group["phase_bias_span_percentage_points"]>gates["phase_bias_span_percentage_points"]:
            failures.append(group["condition"]+": phase dependence")
    if any(not r["converged"] for r in rows): failures.append("single fit nonconvergence")
    if max(r["pre_event_model_mV"] for r in rows)>gates["pre_event_signal_mV"]: failures.append("pre-event model signal")
    if max(r["pre_event_physical_mV"] for r in rows)>gates["pre_event_signal_mV"]: failures.append("pre-event physical signal")
    summary = dict(groups=groups, single_gate_passed=not failures, failures=failures,
        no_pair_expansion=bool(failures), no_new_ROOT=True, default_changed=False,
        scope="16 selected historical mothers, 64 parent/phase cases, eight digital repeats each; response/mean control, not event-by-event energy resolution or blind efficiency")
    save(ART / "single-report.json", summary)
    save(ART / "single-rows.json", rows)
    save(EXP / "single-gate.json", summary)
    import csv
    (EXP / "tables").mkdir(exist_ok=True)
    with (EXP / "tables/single-controls.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


def main():
    compatibility()
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "baseline", "numerical", "pilot", "single", "report"))
    args = parser.parse_args()
    if args.mode == "prepare": prepare()
    elif args.mode == "baseline": baseline()
    elif args.mode == "numerical": numerical()
    elif args.mode == "pilot": simulate(True)
    elif args.mode == "single": simulate(False)
    else: report()


if __name__ == "__main__":
    main()
