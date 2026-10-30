"""Frozen paired covariance ablation on existing digitizations and new controls."""
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
import sys
import time

import numpy as np

from .alpha_pair_energy_covariance import METHODS, candidate_factor, compare_known_times
from .alpha_pair_energy_noise import gaussian_waveforms
from .alpha_pair_recovery_energy import RecoveryEnergyResponse

PROJECT = Path(os.environ.get("GAGG_ENERGY_COVARIANCE_PROJECT", Path(__file__).resolve().parents[2])).resolve()
NAME = "alpha_pair_energy_covariance_20261030"
EXPERIMENT = PROJECT / "pileup/experiments" / NAME
OUTPUT = PROJECT / "jupyter/artifacts" / NAME
RECORD = PROJECT / "jupyter/run_records" / NAME
OLD_NAME = "alpha_pair_energy_noise_20261030"
OLD_EXP = PROJECT / "pileup/experiments" / OLD_NAME
OLD_OUT = PROJECT / "jupyter/artifacts" / OLD_NAME
DOMAINS = ("physical", "energy_gaussian")


def load(path):
    return json.loads(Path(path).read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name+".partial")
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+"\n")
    tmp.replace(path)


def prepare():
    assert not (EXPERIMENT / "method-freeze.json").exists(), "Already frozen"
    for folder in (EXPERIMENT, OUTPUT / "inputs", RECORD):
        folder.mkdir(parents=True, exist_ok=True)
    for name in ("energy-model.json", "template.npz", "response-parameters.json", "sensor-config.json"):
        shutil.copyfile(OLD_EXP / name, EXPERIMENT / name)
    old_freeze = load(PROJECT / "jupyter/run_records" / OLD_NAME / "result-freeze.json")
    for rel, expected in old_freeze["files"].items():
        assert digest(OLD_OUT / rel) == expected, rel
    save(EXPERIMENT / "parent-result-hashes.json", old_freeze)
    source = EXPERIMENT / "source/pileup/analysis"
    source.mkdir(parents=True, exist_ok=True)
    (source.parent / "__init__.py").write_text("")
    (source / "__init__.py").write_text("")
    for name in ("alpha_pair_energy_covariance.py", "alpha_pair_energy_covariance_study.py",
                 "test_alpha_pair_energy_covariance.py", "alpha_pair_energy_noise.py",
                 "alpha_pair_energy_oracle.py", "alpha_pair_recovery_energy.py",
                 "alpha_pair_model_sampler.py", "test_alpha_pair_recovery_energy.py"):
        shutil.copyfile(PROJECT / "pileup/analysis" / name, source / name)
    (EXPERIMENT / "run.py").write_text(
        "from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\n"
        "os.environ['GAGG_ENERGY_COVARIANCE_PROJECT']=str(f.parents[2])\n"
        "sys.path.insert(0,str(f/'source'))\nif __name__=='__main__':\n"
        "    from pileup.analysis.alpha_pair_energy_covariance_study import main\n    main()\n")
    save(EXPERIMENT / "plan.json", [dict(id=i, domains=list(DOMAINS), repeats=8,
                                         gaussian_seed=713000000+1000*i) for i in range(24)])
    save(EXPERIMENT / "protocol.json", dict(
        created=datetime.now().astimezone().isoformat(),
        git_commit="67796208c81000563470044f2482e9622c582f77",
        git_status="Full/shared-filesystem scan timed out; scoped check preserved README/docs changes and prior untracked noise study; exact source is snapshotted",
        scope="Development paired diagnostic, known K/t and pure alpha; four reused exposed optical mother pairs",
        parent_experiment=OLD_NAME, physical_frames=192, new_model_gaussian_frames=192,
        total_scored_method_rows=1920, methods=list(METHODS), primary="energy_once versus reference_gls",
        response="Same frozen dual95/351ns, recovery strength2.515148098, scalar energy calibration; no mean/noise coefficient relearning",
        covariance="Same empirical banded formula with kappa=0. all_legacy_once evaluates old charge-conditioned mean curves; energy methods use the proposed recovery-energy curves and shared charges",
        reference="Old one-pulse covariance from16 independent reference repeats, common to every method as baseline/initialization; not a single-shot deployable calibration",
        energy_truth="Only generation and independent scoring. Worker inputs contain voltage, supplied times and old observed reference metric/initialization",
        gaussian="New Gaussian draws with corrected energy covariance at true model E/t; exact FIR-domain sampling, zero voltage-noise prefix; internal control, not new detector data",
        domains=list(DOMAINS), conditions="weak/high singles; high->weak8/40/200ns; weak->high40ns; blocks0,10,20,30",
        observation="1024x2ns, sample centers0.25+2n ns; first active photon128ns; source generation4.096us inherited; no ADC range rejection",
        fitting="Baseline unchanged6 starts; both once-weight variants use baseline energy plus linear+5share starts. FGLS up to4weight updates, subsequent steps warm-start only; stop max|delta ln E|<=0.001. NLL uses baseline/linear/5shares/once solution with120iterations each",
        likelihood="r^T C(E)^-1 r + logdet C(E), common nlog2pi omitted. Fixed GLS cost values from different C are not compared",
        end_points="All injected components E_dep within5/10/20%; signed bias,median,sigma68,RMSE,boundaries,solver and outer stability, paired gained/lost; physical Q>=500 subset separate",
        no_count_metrics="K/t fixed; misses/extras/splits/merges/type efficiency not measured",
        statistics="32 repeats per condition=4mother pairs x8; report every mother block; no population CI or independent optical-confirmation claim",
        pilot="Cases2 and5, both domains, first repeat, cost only; these waves remain in full development set",
        workers=4, blas_threads=1, cost_gate_seconds=1800., cost_safety_factor=1.5,
        next_gate="No default replacement on reused four mothers; independent mother data and blind K/t validation required",
        no_new_ROOT=True, no_new_G4=True, no_new_training=True, default_changed=False))
    save(EXPERIMENT / "method-freeze.json", dict(files={str(p.relative_to(EXPERIMENT)): digest(p)
        for p in EXPERIMENT.rglob("*") if p.is_file()}, before_new_gaussian_and_fits=True))
    response = RecoveryEnergyResponse.load(EXPERIMENT)
    labels = []
    for task in load(EXPERIMENT / "plan.json"):
        case = load(OLD_OUT / "cases" / f"{task['id']:02}.json")
        t = case["task"]
        with np.load(OLD_OUT / "cases" / f"{task['id']:02}.npz") as old:
            voltage = old["physical"][16:24].copy()
            reference_factor = old["factor"].copy()
            physical_charges = old["charges"][16:24].copy()
        curves, _ = response.components(t["times_ns"], t["energy_MeV"])
        mean = curves.sum(0)+response.baseline+response.background
        factor = candidate_factor(response, t["times_ns"], t["energy_MeV"])
        gaussian = gaussian_waveforms(response, mean, factor, task["gaussian_seed"], 8)
        np.savez_compressed(OUTPUT / "inputs" / f"{task['id']:02}.npz",
            physical=voltage, energy_gaussian=gaussian, known_times=np.asarray(t["times_ns"]),
            reference_factor=reference_factor)
        save(OUTPUT / "inputs" / f"{task['id']:02}.json", dict(
            reference_initial=case["initial_from_reference_voltage"],
            reference_source="Independent16-repeat mean voltage; same as parent study"))
        labels.append(dict(t, physical_component_Q_pe=physical_charges.tolist(),
                           gaussian_seed=task["gaussian_seed"], parent_npz_sha256=case["output_sha256"]))
    save(OUTPUT / "labels.json", labels)
    save(EXPERIMENT / "input-hashes.json", {str(p.relative_to(OUTPUT)): digest(p)
        for p in (OUTPUT / "inputs").glob("*") if p.is_file()})
    save(EXPERIMENT / "label-hash.json", {"labels.json": digest(OUTPUT / "labels.json")})
    save(EXPERIMENT / "data-freeze.json", dict(files={name: digest(EXPERIMENT / name)
        for name in ("method-freeze.json", "input-hashes.json", "label-hash.json")},
        before_first_fit=True))
    save(RECORD / "environment.json", dict(python=sys.version, executable=sys.executable,
        packages={n: version(n) for n in ("numpy", "scipy", "matplotlib", "nbformat")},
        workers=4, BLAS_threads=1))
    print("PREPARED:192 reused physical +192 new model Gaussian frames;1920 method rows", flush=True)


def verify():
    for rel, expected in load(EXPERIMENT / "data-freeze.json")["files"].items():
        assert digest(EXPERIMENT / rel) == expected, rel
    for rel, expected in load(EXPERIMENT / "method-freeze.json")["files"].items():
        assert digest(EXPERIMENT / rel) == expected, rel
    for rel, expected in load(EXPERIMENT / "input-hashes.json").items():
        assert digest(OUTPUT / rel) == expected, rel
    if (EXPERIMENT / "execution-freeze.json").exists():
        for rel, expected in load(EXPERIMENT / "execution-freeze.json")["files"].items():
            assert digest(EXPERIMENT / rel) == expected, rel


def initialize():
    global RESPONSE
    RESPONSE = RecoveryEnergyResponse.load(EXPERIMENT)
    # Inference explicitly cannot read per-event labels, raw ROOT or old results.
    def audit(event, args):
        if event != "open" or not isinstance(args[0], (str, bytes, os.PathLike)):
            return
        path = Path(os.fsdecode(args[0])).resolve()
        if path == OUTPUT / "labels.json" or path.suffix.lower() == ".root" or path.is_relative_to(OLD_OUT):
            raise PermissionError(f"Inference cannot read per-event truth/old results: {path}")
    sys.addaudithook(audit)


def one(task, pilot=False):
    started = time.monotonic()
    context = load(OUTPUT / "inputs" / f"{task['id']:02}.json")
    rows = []
    with np.load(OUTPUT / "inputs" / f"{task['id']:02}.npz") as data:
        assert set(data.files) == {"physical", "energy_gaussian", "known_times", "reference_factor"}
        for domain in DOMAINS:
            for repeat in range(1 if pilot else task["repeats"]):
                wave_started = time.monotonic()
                fits = compare_known_times(RESPONSE, data[domain][repeat], data["known_times"],
                    data["reference_factor"], context["reference_initial"])
                rows.append(dict(domain=domain, repeat=repeat, fits=fits,
                                 seconds=time.monotonic()-wave_started))
    result = dict(id=task["id"], pilot=pilot, rows=rows, seconds=time.monotonic()-started,
                  worker_truth_read_guard=True, supplied_K_t=True, true_energy_or_Q_input=False)
    save(OUTPUT / ("pilot" if pilot else "cases") / f"{task['id']:02}.json", result)
    return dict(id=task["id"], seconds=result["seconds"], waves=len(rows))


def pilot():
    verify(); initialize()
    results = [one(task, True) for task in load(EXPERIMENT / "plan.json") if task["id"] in (2, 5)]
    estimate = max(r["seconds"]/r["waves"] for r in results)*384/4*1.5
    record = dict(passed=estimate <= 1800., results=results,
                  conservative_wall_seconds=estimate, maximum_workers=4,
                  cost_only=True, accuracy_not_used_for_gate=True)
    save(RECORD / "pilot-validation.json", record)
    print(json.dumps(record, indent=2), flush=True)
    if not record["passed"]:
        raise RuntimeError("Cost gate failed; preserve frozen plan before revising computation")


def execution_plan():
    verify()
    assert not (EXPERIMENT / "execution-freeze.json").exists()
    assert not (OUTPUT / "cases").exists()
    prior = load(RECORD / "pilot-validation.json")
    assert not prior["passed"]
    estimate = prior["conservative_wall_seconds"]*4/6
    assert estimate <= 1800.
    plan = dict(workers=6, blas_threads=1, conservative_wall_seconds=estimate,
        reason="Cost-only concurrency change from4to6; no sample, method, seed or accuracy-based selection",
        physical_frames=192, gaussian_frames=192, method_rows=1920,
        pilot_accuracy_not_used=True, all_original_inputs_preserved=True)
    save(EXPERIMENT / "execution-plan.json", plan)
    shutil.copytree(EXPERIMENT / "source", EXPERIMENT / "source-v2")
    shutil.copyfile(PROJECT / "pileup/analysis/alpha_pair_energy_covariance_study.py",
                    EXPERIMENT / "source-v2/pileup/analysis/alpha_pair_energy_covariance_study.py")
    script = (EXPERIMENT / "run.py").read_text().replace("f/'source'", "f/'source-v2'")
    (EXPERIMENT / "run_v2.py").write_text(script)
    files = [p for p in (EXPERIMENT / "source-v2").rglob("*") if p.is_file()]
    files += [EXPERIMENT / "execution-plan.json", EXPERIMENT / "run_v2.py"]
    save(EXPERIMENT / "execution-freeze.json", dict(files={str(p.relative_to(EXPERIMENT)): digest(p) for p in files},
        full_evaluation_not_started=True, methods_unchanged=True))
    save(RECORD / "execution-plan-validation.json", dict(passed=True, **plan))
    print(json.dumps(plan, indent=2), flush=True)


def run(workers):
    verify()
    assert 1 <= workers <= 6
    assert load(RECORD / "execution-plan-validation.json")["passed"]
    start = time.monotonic()
    tasks = load(EXPERIMENT / "plan.json")
    remaining = []
    for task in tasks:
        path = OUTPUT / "cases" / f"{task['id']:02}.json"
        if path.exists():
            old = load(path)
            assert old["id"] == task["id"] and len(old["rows"]) == 16 and not old["pilot"]
        else:
            remaining.append(task)
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn"), initializer=initialize) as pool:
        pending = [pool.submit(one, task) for task in remaining]
        for i, future in enumerate(as_completed(pending), 1):
            print(json.dumps(dict(done=i, total=len(remaining), **future.result())), flush=True)
    files = {str(p.relative_to(OUTPUT)): digest(p) for p in sorted((OUTPUT / "cases").glob("*.json"))}
    assert len(files) == 24
    save(RECORD / "result-freeze.json", dict(files=files, method_freeze_sha256=digest(EXPERIMENT / "method-freeze.json")))
    save(RECORD / "run-validation.json", dict(passed=True, cases=24, newly_run=len(remaining),
        frames=384, method_rows=1920, workers=workers, seconds=time.monotonic()-start,
        known_K_t=True, true_E_Q_inference=False, no_new_ROOT=True, default_changed=False))
    print("COMPLETE:384 frames,1920 method rows", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("prepare", "pilot", "execution-plan", "run"))
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    if args.action == "prepare":
        prepare()
    elif args.action == "pilot":
        pilot()
    elif args.action == "execution-plan":
        execution_plan()
    else:
        run(args.workers)


if __name__ == "__main__":
    main()
