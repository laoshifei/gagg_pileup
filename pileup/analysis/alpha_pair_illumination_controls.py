"""Paired position interventions on frozen optical mothers (truth diagnostic).

Original, shuffled XY and uniform XY arms preserve every arrival time, photon
energy and active-area membership. No reconstruction or model fitting occurs.
Run prepare once, then use the frozen experiment's run.py from the project root.
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
import sys
import time

import numpy as np

from .alpha_pair_causal_controls import compatibility, read, save, sha
from .alpha_pair_causal_recovery import CausalRecovery, CausalReadout, DifferentiableAge
from .alpha_pair_digitizer_optimized import simulate_avalanches
from .alpha_pair_energy_calibrated import CalibratedEnergyResponse
from .alpha_pair_recovery_energy import RecoveryEnergyResponse
from .sipm_readout import ReadoutConfig, synthesize_readout
from .sipm_response import AvalancheResult, HC_EV_NM, PhotonHits, SensorConfig, map_cells

NAME = "alpha_pair_illumination_20261102"
ROOT = Path(os.environ.get("GAGG_ILLUMINATION_ROOT", Path(__file__).resolve().parents[2])).resolve()
EXP = ROOT / "pileup/experiments" / NAME
ART = ROOT / "jupyter/artifacts" / NAME
REC = ROOT / "jupyter/run_records" / NAME
ARMS = ("original", "shuffled", "uniform")
LAYOUTS = 4
REPEATS = 8
WORKERS = 4


def cells_for(donor, sensor):
    # Historical origin is first ACTIVE-area photon; outside photons can precede it.
    t = donor["time_ns"] - min(0., float(donor["time_ns"].min()))
    hits = PhotonHits(t, donor["x_mm"], donor["y_mm"], donor["energy_eV"],
                      window_ns=max(20000., float(t.max()) + 1.))
    return map_cells(hits, sensor)


def transform_positions(donor, sensor, parent, layout, arm):
    if arm not in ARMS:
        raise ValueError(arm)
    result = {k: v.copy() for k, v in donor.items()}
    inside = cells_for(donor, sensor) >= 0
    rng = np.random.default_rng(np.random.SeedSequence([713010001, parent, layout, ARMS.index(arm)]))
    if arm == "shuffled":
        permutation = rng.permutation(np.flatnonzero(inside))
        for k in ("x_mm", "y_mm"):
            result[k][inside] = donor[k][permutation]
    elif arm == "uniform":
        result["x_mm"][inside] = rng.uniform(sensor.x_min_mm, sensor.x_max_mm, int(inside.sum()))
        result["y_mm"][inside] = rng.uniform(sensor.y_min_mm, sensor.y_max_mm, int(inside.sum()))
    return result


def initialize():
    compatibility()
    global CONFIG, SENSOR, READOUT, DONORS
    CONFIG = read(EXP / "sensor-config.json")
    SENSOR = SensorConfig(**CONFIG["sensor"])
    READOUT = ReadoutConfig(**CONFIG["readout"])
    DONORS = {}
    for arm in ARMS:
        for layout in range(LAYOUTS):
            for p in (ART / "photons" / arm / str(layout)).glob("*.npz"):
                with np.load(p) as f:
                    DONORS[arm, layout, int(p.stem)] = {k: f[k].copy() for k in f.files}


def digitize(parents, times, seed, recovery, arm, layout):
    fields = {k: [] for k in ("time_ns", "x_mm", "y_mm", "energy_eV")}
    owners = []
    for j, (parent, onset) in enumerate(zip(parents, times)):
        donor = DONORS[arm, layout, parent]
        t = donor["time_ns"] + onset
        keep = t < 4096.
        for k in fields:
            fields[k].append((t if k == "time_ns" else donor[k])[keep])
        owners.append(np.full(int(keep.sum()), j))
    arrays = {k: np.concatenate(v) for k, v in fields.items()}
    order = np.argsort(arrays["time_ns"], kind="stable")
    owner = np.concatenate(owners)[order]
    hits = PhotonHits(**{k: v[order] for k, v in arrays.items()}, window_ns=4096.)
    av = simulate_avalanches(hits, replace(SENSOR, recovery_enabled=recovery), seed)
    rec = av.records
    assignment = np.full(len(rec), -1, int)
    direct = rec["kind"] == 0
    assignment[direct] = owner[rec["photon_index"][direct]]
    for j in np.flatnonzero(rec["parent_id"] >= 0):
        assignment[j] = assignment[rec["parent_id"][j]]
    signals, charges = [], []
    for j in range(len(parents)):
        subset = rec[assignment == j]
        event = AvalancheResult(subset, 4096., av.full_charge_C, np.array([]), {})
        signals.append(synthesize_readout(event, READOUT, 0).arrays["signal_mV"][:1024])
        charges.append(float(subset["charge_C"].sum() / av.full_charge_C))
    dark_charge = rec[assignment < 0]["charge_C"].sum() / av.full_charge_C
    assert np.isclose(sum(charges) + dark_charge, rec["charge_C"].sum() / av.full_charge_C, rtol=1e-12)
    analog = synthesize_readout(av, READOUT, seed + 900000000).arrays["analog_mV"][:1024]
    return np.asarray(signals), np.asarray(charges), analog


def one(task):
    start = time.monotonic()
    result = {}
    arms = task.get("arms", ARMS)
    for arm in arms:
        for recovery in (True, False):
            shared, isolated, qs, qi, analog = [], [], [], [], []
            for repeat in range(task["repeats"]):
                seed = task["seed_base"] + repeat
                curves, q, full = digitize(task["parents"], task["times_ns"], seed, recovery, arm, task["layout"])
                singles = [digitize([p], [t], seed, recovery, arm, task["layout"])
                           for p, t in zip(task["parents"], task["times_ns"])]
                shared.append(curves)
                qs.append(q)
                analog.append(full)
                isolated.append(np.concatenate([v[0] for v in singles]))
                qi.append(np.concatenate([v[1] for v in singles]))
            prefix = arm + ("_on" if recovery else "_off")
            for key, values in (("shared", shared), ("isolated", isolated),
                                ("shared_charge", qs), ("isolated_charge", qi), ("analog", analog)):
                result[prefix + "_" + key] = np.asarray(values)
    assert all(np.isfinite(v).all() for v in result.values())
    return dict(id=task["id"], case=task["case"], layout=task["layout"],
                seconds=time.monotonic() - start, digitizations=len(arms) * task["repeats"] * 6), result


def snapshot():
    source = EXP / "source/pileup/analysis"
    source.mkdir(parents=True)
    (source.parent / "__init__.py").write_text("")
    (source / "__init__.py").write_text("")
    pending, copied = ["alpha_pair_illumination_controls"], set()
    while pending:
        name = pending.pop()
        if name in copied:
            continue
        p = ROOT / "pileup/analysis" / (name + ".py")
        tree = ast.parse(p.read_text())
        shutil.copyfile(p, source / p.name)
        copied.add(name)
        pending.extend(n.module for n in ast.walk(tree)
                       if isinstance(n, ast.ImportFrom) and n.level == 1 and n.module)
    for name in ("sipm_age_mean.cpp", "sipm_age_mean_jacobian.cpp"):
        shutil.copyfile(ROOT / "pileup/analysis" / name, source / name)
    (EXP / "run.py").write_text(
        "from pathlib import Path\nimport os,sys\nf=Path(__file__).resolve().parent\n"
        "os.environ['GAGG_ILLUMINATION_ROOT']=str(f.parents[2])\n"
        "os.environ['OPENBLAS_NUM_THREADS']='1'\nos.environ['OMP_NUM_THREADS']='1'\n"
        "sys.path.insert(0,str(f/'source'))\n"
        "if __name__=='__main__':\n"
        "    from pileup.analysis.alpha_pair_illumination_controls import main\n    main()\n")
    (EXP / "bin").mkdir()
    command = ["g++", "-O3", "-std=c++17", "-fPIC", "-shared",
               str((source / "sipm_age_mean.cpp").relative_to(ROOT)),
               str((source / "sipm_age_mean_jacobian.cpp").relative_to(ROOT)),
               "-o", str((EXP / "bin/age.so").relative_to(ROOT))]
    subprocess.run(command, cwd=ROOT, check=True)
    save(EXP / "build.json", dict(command=command,
         compiler=subprocess.check_output(["g++", "--version"], text=True).splitlines()[0],
         binary_sha256=sha(EXP / "bin/age.so")))


def predict():
    """Freeze the legacy predictor and a truth-photon diagnostic BEFORE runs."""
    initialize()
    anchor = CalibratedEnergyResponse.load(RecoveryEnergyResponse.load(EXP), EXP / "mean-correction.npz")
    models = {}
    for dt in (1., .5):
        with np.load(EXP / f"charge-map-{dt}.npz") as f:
            mapping = {k: f[k].copy() for k in f.files}
        models[dt] = CausalRecovery(anchor, CONFIG, EXP / "bin/age.so", mapping)
    provenance = {r["parent"]: r for r in read(EXP / "mother-provenance.json")}
    rows, curves = [], {}
    for task in read(EXP / "cases.json"):
        weak = task["parents"].index(task["block"])
        energies = [provenance[p]["energy_dep_MeV"] for p in task["parents"]]
        row = dict(case=task["id"], condition=task["condition"], gap_ns=task["gap_ns"],
                   block=task["block"], energy_dep_MeV=energies)
        for dt, model in models.items():
            solo = model.components(task["times_ns"], energies, interaction=False)
            joint = model.components(task["times_ns"], energies)
            row[f"legacy_{dt}_percent"] = float(100 * (joint[weak].sum() / solo[weak].sum() - 1))
            light = []
            for parent, onset in zip(task["parents"], task["times_ns"]):
                donor = DONORS["original", 0, parent]
                t = donor["time_ns"] + onset
                keep = (cells_for(donor, SENSOR) >= 0) & (t >= 0) & (t < 2048.)
                weights = np.interp(HC_EV_NM / donor["energy_eV"][keep],
                                    SENSOR.pde_wavelength_nm, SENSOR.pde_probabilities[-1])
                histogram = np.bincount(np.floor(t[keep] / dt).astype(int), weights=weights,
                                        minlength=int(2048 / dt))
                assert np.isclose(histogram.sum(), weights.sum(), rtol=1e-13)
                light.append(histogram)
            light = np.array(light)
            total = light.sum(0)
            shares = np.divide(light, total, out=np.zeros_like(light), where=total > 0)
            response = model.age.forward(total)
            empirical_solo = np.array([model.readout(model.age.forward(v)) for v in light])
            empirical_joint = np.array([model.readout(response * share) for share in shares])
            row[f"photon_{dt}_percent"] = float(100 * (empirical_joint[weak].sum() / empirical_solo[weak].sum() - 1))
            if dt == 1.:
                for key, value in (("legacy_solo", solo), ("legacy_joint", joint),
                                   ("photon_solo", empirical_solo), ("photon_joint", empirical_joint)):
                    curves[f"case{task['id']}_{key}"] = value
        rows.append(row)
    old_rows = read(EXP / "previous-rows.json")
    difference = max(abs(r["legacy_1.0_percent"] - old_rows[r["case"]]["predicted_weak_area_change_percent"]) for r in rows)
    assert difference < 1e-9, difference
    save(EXP / "predictions.json", dict(rows=rows, legacy_max_difference_percentage_points=difference,
         photon_model="Truth on-area photon counts/time/energy, full-recovery PDE weighting, uniform age expectation, no fitted scale or isolated anchor. Diagnostic only; dark history/neighbour XT closure remain approximate."))
    np.savez_compressed(EXP / "prediction-curves.npz", **curves)


def prepare():
    compatibility()
    previous = ROOT / "pileup/experiments/alpha_pair_causal_shared_20261101"
    previous_art = ROOT / "jupyter/artifacts/alpha_pair_causal_shared_20261101"
    for folder in (EXP, ART, REC):
        folder.mkdir(parents=True, exist_ok=False)
    for name in ("template.npz", "energy-model.json", "response-parameters.json", "mean-correction.npz",
                 "noise-calibration.json", "sensor-config.json", "charge-map-1.0.npz", "charge-map-0.5.npz",
                 "mother-provenance.json"):
        shutil.copyfile(previous / name, EXP / name)
    shutil.copyfile(previous_art / "rows.json", EXP / "previous-rows.json")
    old_cases = read(previous / "plan.json")
    save(EXP / "cases.json", old_cases)
    tasks = []
    for case in old_cases:
        for layout in range(LAYOUTS):
            tasks.append(dict(case, id=len(tasks), case=case["id"], layout=layout, repeats=REPEATS,
                              seed_base=713100000 + case["block"] * 10000 + layout * 100))
    save(EXP / "plan.json", tasks)
    save(EXP / "protocol.json", dict(created_utc=datetime.now(timezone.utc).isoformat(),
        commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        git_status=subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True),
        scope="Known-input mechanistic diagnostic; same eight historical optical mothers and four pairs. No inference, new optical transport or independent performance confirmation.",
        intervention=dict(original="Unchanged input", shuffled="Joint XY permutation within each mother's active-area photons; exact marginal cell counts retained; removes time/energy-position association",
                          uniform="Independent uniform XY inside active area; time/energy pair and active-area membership retained, spatial marginal and association changed",
                          outside="All fields of outside-area photons retained", layout_rng="SeedSequence([713010001,parent,layout,arm_index])"),
        design=dict(mother_pairs=4, conditions=5, arms=list(ARMS), layouts=LAYOUTS, digitizer_seeds_per_layout=REPEATS,
                    production_digitizations=len(tasks) * REPEATS * len(ARMS) * 6,
                    controls="joint plus both singles, recovery on/off for each arm; independent digital seeds by pair/layout/repeat; same seeds across arms/orders/gaps",
                    coupling="Photon uniforms remain paired while altered firings change secondary RNG consumption; no identical-avalanche assumption",
                    time_window_ns=4096, observed_window_ns=2048, workers=WORKERS, blas_threads=1),
        endpoint="Equal-mother mean of 100*(weak shared area / same-arm isolated area - 1), integrated first2048ns aperture signals. Own-arm denominator prevents isolated-response changes masquerading as shared suppression.",
        primary="Original vs uniform at8/40ns; shuffled,200/1200ns and reverse40ns are mechanistic controls, all retained.",
        models="Frozen previous anchored uniform prediction; additionally prespecified truth-photon-time/PDE weighted uniform age diagnostic at1/.5ns, never used as blind input or tuned to simulations.",
        gates=dict(pilot_projected_seconds=900., recovery_off_abs_percentage_points=3., uniform_effect_at_short_gaps_percentage_points=3.,
                   inherited_uniform_model_abs_error_percentage_points=5., empirical_grid_change_percentage_points=1.),
        uncertainty="10000 paired nested resamples of mother pair, layout block and digital seed; four exposed mother pairs only. Descriptive intervals, no independent transport/systematic coverage.",
        decision="Spatial contribution supported on selected mothers if uniform reduces short-gap suppression by>=3pp with paired interval above zero and recovery-off controls pass. Passing transformed inputs does not validate the model on original inputs. No automatic expansion to E fits.",
        no_new_ROOT=True, default_changed=False, library="Configured shared paper library unavailable; reuse project L86 restrictions, no new external literature claim."))
    sensor = SensorConfig(**read(EXP / "sensor-config.json")["sensor"])
    for row in read(EXP / "mother-provenance.json"):
        parent = row["parent"]
        with np.load(previous_art / "photons" / f"{parent:04d}.npz") as f:
            donor = {k: f[k].copy() for k in f.files}
        for arm in ARMS:
            for layout in range(LAYOUTS):
                dest = ART / "photons" / arm / str(layout)
                dest.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(dest / f"{parent:04d}.npz", **transform_positions(donor, sensor, parent, layout, arm))
    (ART / "baseline-reference").mkdir()
    for case in (1, 4):
        shutil.copyfile(previous_art / "cases" / f"{case:04d}.npz", ART / "baseline-reference" / f"{case:04d}.npz")
    snapshot()
    predict()
    import scipy
    save(EXP / "environment.json", dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                                        numpy_trapezoid_adapter=compatibility()))
    save(EXP / "source-config-freeze.json", {str(p.relative_to(ROOT)): sha(p)
         for p in EXP.rglob("*") if p.is_file() and "bin" not in p.relative_to(EXP).parts})
    save(EXP / "input-freeze.json", {str(p.relative_to(ROOT)): sha(p) for p in ART.rglob("*.npz")})
    print("Prepared80 tasks;11520 production digitizations; inputs and predictions frozen", flush=True)


def verify():
    count = 0
    for name in ("source-config-freeze.json", "input-freeze.json"):
        for path, expected in read(EXP / name).items():
            assert sha(ROOT / path) == expected, path
            count += 1
    return count


def check():
    hashes = verify()
    initialize()
    spatial, checked = [], 0
    for key, donor in DONORS.items():
        arm, layout, parent = key
        reference = DONORS["original", 0, parent]
        old_cells, new_cells = cells_for(reference, SENSOR), cells_for(donor, SENSOR)
        inside = old_cells >= 0
        assert np.array_equal(new_cells >= 0, inside)
        assert all(np.array_equal(donor[k], reference[k]) for k in ("time_ns", "energy_eV"))
        assert all(np.array_equal(donor[k][~inside], reference[k][~inside]) for k in donor)
        if arm in ("original", "shuffled"):
            assert np.array_equal(np.bincount(old_cells[inside], minlength=SENSOR.n_cells),
                                  np.bincount(new_cells[inside], minlength=SENSOR.n_cells))
        expected = transform_positions(reference, SENSOR, parent, layout, arm)
        assert all(np.array_equal(donor[k], expected[k]) for k in donor)
        checked += 1
    for block in (0, 10, 20, 30):
        for arm in ARMS:
            for layout in range(LAYOUTS):
                probabilities = []
                for parent in (block + 240, block):
                    cells = cells_for(DONORS[arm, layout, parent], SENSOR)
                    counts = np.bincount(cells[cells >= 0], minlength=SENSOR.n_cells).astype(float)
                    probabilities.append(counts / counts.sum())
                spatial.append(dict(block=block, arm=arm, layout=layout,
                    overlap=float(SENSOR.n_cells * np.dot(*probabilities))))
    checks = []
    folder = ART / "baseline-replay"
    folder.mkdir(exist_ok=False)
    cases = read(EXP / "cases.json")
    for case_id in (1, 4):
        task = dict(cases[case_id], case=case_id, layout=0, repeats=2, arms=["original"])
        metadata, arrays = one(task)
        with np.load(ART / "baseline-reference" / f"{case_id:04d}.npz") as old:
            exact = all(np.array_equal(arrays["original_" + k], old[k][:2]) for k in old.files)
        assert exact
        np.savez_compressed(folder / f"{case_id:04d}.npz", **arrays)
        checks.append(dict(case=case_id, first_two_repeats_exact=exact))
    save(EXP / "input-check.json", dict(passed=True, frozen_hashes=hashes, photon_arrays_checked=checked,
         preserved="Times, energies, all outside photons and per-photon active membership exact; shuffled cell counts exact; seeds reproduce transforms",
         baseline_replays=checks, baseline_digitizations=24, spatial=spatial))
    print("Input intervention checks and two old physical replays passed", flush=True)


def run(mode, output=None, ids=None):
    verify()
    assert read(EXP / "input-check.json")["passed"]
    tasks = read(EXP / "plan.json")
    if mode == "pilot":
        tasks = [dict(tasks[4], repeats=2)]  # first pair,40ns,layout0
        folder = ART / "pilot"
    elif mode == "run":
        assert read(REC / "pilot.json")["passed"]
        folder = ART / "cases"
    else:
        if not output or not ids:
            raise ValueError("Replay requires --output and --ids")
        tasks = [t for t in tasks if t["id"] in ids]
        assert {t["id"] for t in tasks} == set(ids)
        folder = Path(output)
    folder.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    rows = []
    with ProcessPoolExecutor(max_workers=WORKERS, initializer=initialize) as pool:
        for row, arrays in pool.map(one, tasks):
            if mode == "replay":
                with np.load(ART / "cases" / f"{row['id']:04d}.npz") as f:
                    row["physical_arrays_exact"] = all(np.array_equal(arrays[k], f[k]) for k in arrays)
                assert row["physical_arrays_exact"]
            np.savez_compressed(folder / f"{row['id']:04d}.npz", **arrays)
            save(folder / f"{row['id']:04d}.json", row)
            rows.append(row)
            if len(rows) % 4 == 0 or mode != "run":
                print(mode, len(rows), "/", len(tasks), "elapsed", round(time.monotonic() - start, 1), flush=True)
    summary = dict(tasks=len(rows), digitizations=sum(r["digitizations"] for r in rows),
                   seconds=time.monotonic() - start, workers=WORKERS)
    if mode == "pilot":
        estimate = rows[0]["seconds"] * 4 * 80 / WORKERS * 1.5
        summary.update(projected_seconds=estimate, limit_seconds=900., passed=estimate < 900.)
        save(REC / "pilot.json", summary)
        assert summary["passed"], summary
    elif mode == "run":
        save(REC / "run.json", summary)
    else:
        summary.update(checks=rows, passed=all(r["physical_arrays_exact"] for r in rows))
        save(folder / "replay.json", summary)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "check", "pilot", "run", "replay"))
    parser.add_argument("--output")
    parser.add_argument("--ids", type=int, nargs="+")
    args = parser.parse_args()
    if args.mode == "prepare":
        prepare()
    elif args.mode == "check":
        check()
    else:
        run(args.mode, args.output, args.ids)


if __name__ == "__main__":
    main()
