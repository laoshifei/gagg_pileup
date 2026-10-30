"""Independent scoring and fixed-number figures for covariance ablation."""
from __future__ import annotations

from collections import defaultdict
import csv
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .alpha_pair_energy_covariance import METHODS, candidate_factor, gaussian_terms
from .alpha_pair_energy_covariance_study import (PROJECT, EXPERIMENT, OUTPUT, RECORD,
    OLD_OUT, DOMAINS, load, save, digest, verify)
from .alpha_pair_energy_noise import diagonal_factor
from .alpha_pair_recovery_energy import RecoveryEnergyResponse

CONDITIONS = [("weak", None), ("high", None), ("high_weak", 8.),
              ("high_weak", 40.), ("high_weak", 200.), ("weak_high", 40.)]
LABELS = ["Weak single", "Strong single", "Strong→weak\n8 ns", "Strong→weak\n40 ns",
          "Strong→weak\n200 ns", "Weak→strong\n40 ns"]
NAMES = {"reference_gls":"Reference GLS", "all_legacy_once":"All / legacy", "energy_once":"Energy / once",
         "energy_fgls":"Energy / FGLS", "energy_nll":"Energy / NLL"}
COLORS = dict(zip(METHODS, ["#777777", "#d4a800", "#287cb5", "#2d9b72", "#ce553b"]))


def csv_table(name, rows):
    path = EXPERIMENT / "tables" / f"{name}.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def finish(fig, name):
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(OUTPUT / f"{name}.{ext}", dpi=170, bbox_inches="tight")
    plt.close(fig)


def main():
    verify()
    for rel, expected in load(RECORD / "result-freeze.json")["files"].items():
        assert digest(OUTPUT / rel) == expected
    for rel, expected in load(EXPERIMENT / "label-hash.json").items():
        assert digest(OUTPUT / rel) == expected
    model = RecoveryEnergyResponse.load(EXPERIMENT)
    truth = {r["id"]: r for r in load(OUTPUT / "labels.json")}
    results = {i: load(OUTPUT / "cases" / f"{i:02}.json") for i in truth}
    frames, components = [], []
    lookup = {}
    baseline_checks = 0
    timings = []
    for i, result in results.items():
        label = truth[i]
        energy = np.asarray(label["energy_MeV"])
        old = load(OLD_OUT / "cases" / f"{i:02}.json")
        old_fits = {r["repeat"]: r["fit"] for r in old["fits"] if r["domain"] == "physical" and r["weight"] == "gls"}
        for row in result["rows"]:
            assert tuple(row["fits"]) == METHODS
            repeat, domain = row["repeat"], row["domain"]
            timings.append(row["seconds"])
            qualified = (bool(np.all(np.asarray(label["physical_component_Q_pe"][repeat]) >= 500.))
                         if domain == "physical" else None)
            if domain == "physical":
                np.testing.assert_allclose(row["fits"]["reference_gls"]["energy_MeV"], old_fits[repeat]["energy_MeV"], rtol=1e-10, atol=1e-12)
                baseline_checks += 1
            for method, fit in row["fits"].items():
                assert fit["times_ns"] == label["times_ns"]
                error = np.asarray(fit["energy_MeV"])/energy-1.
                assert np.isfinite(error).all()
                key = (i, domain, repeat, method)
                assert key not in lookup
                common = dict(case=i, block=label["block"], condition=label["condition"], gap_ns=label["gap_ns"],
                              domain=domain, repeat=repeat, method=method)
                frame = dict(common, success5=int(np.all(abs(error) <= .05)),
                    success10=int(np.all(abs(error) <= .10)), success20=int(np.all(abs(error) <= .20)),
                    not_converged=int(not fit["converged"]), boundary=int(fit["boundary"]),
                    outer_not_stable=int(method == "energy_fgls" and not fit["outer_stable"]),
                    physical_all_Q_ge_500=qualified,
                    cost=fit["cost"], total_evaluations=fit["total_evaluations"])
                frames.append(frame); lookup[key] = frame
                for j, value in enumerate(error):
                    components.append(dict(common, component=j, label=label["labels"][j],
                        truth_energy_MeV=float(energy[j]), predicted_energy_MeV=fit["energy_MeV"][j],
                        relative_error=float(value), success10=int(abs(value) <= .10)))
    assert len(frames) == 1920 and baseline_checks == 192
    groups, comp_groups = defaultdict(list), defaultdict(list)
    for row in frames:
        groups[row["condition"], row["gap_ns"], row["domain"], row["method"]].append(row)
    for row in components:
        comp_groups[row["condition"], row["gap_ns"], row["domain"], row["method"]].append(row)
    summaries, parents, paired = [], [], []
    for condition, gap in CONDITIONS:
        for domain in DOMAINS:
            for method in METHODS:
                key = condition, gap, domain, method
                rows = groups[key]; cr = comp_groups[key]
                weak = np.array([r["relative_error"] for r in cr if r["label"] == "weak"])
                common = dict(condition=condition, gap_ns=gap, domain=domain, method=method)
                summary = dict(common, draws=len(rows), mother_blocks=len({r["block"] for r in rows}),
                    **{k:sum(r[k] for r in rows) for k in
                       ("success5", "success10", "success20", "not_converged", "boundary", "outer_not_stable")},
                    weak_bias=float(weak.mean()) if len(weak) else None,
                    weak_median=float(np.median(weak)) if len(weak) else None,
                    weak_sigma68=float(np.diff(np.percentile(weak, [16, 84]))[0]/2) if len(weak) else None,
                    weak_rmse=float(np.sqrt(np.mean(weak*weak))) if len(weak) else None,
                    qualified_frames=sum(bool(r["physical_all_Q_ge_500"]) for r in rows) if domain == "physical" else None,
                    qualified_success10=sum(r["success10"] for r in rows if r["physical_all_Q_ge_500"]) if domain == "physical" else None)
                assert summary["draws"] == 32
                summaries.append(summary)
                for block in (0, 10, 20, 30):
                    br = [r for r in rows if r["block"] == block]
                    be = [r["relative_error"] for r in cr if r["block"] == block and r["label"] == "weak"]
                    parents.append(dict(common, block=block, draws=len(br), success10=sum(r["success10"] for r in br),
                                        weak_bias=float(np.mean(be)) if be else None))
                if method != "reference_gls":
                    base = [lookup[r["case"], domain, r["repeat"], "reference_gls"] for r in rows]
                    paired.append(dict(common, draws=len(rows),
                        reference_success=sum(r["success10"] for r in base), new_success=summary["success10"],
                        gained=sum(not a["success10"] and b["success10"] for a,b in zip(base, rows)),
                        lost=sum(a["success10"] and not b["success10"] for a,b in zip(base, rows))))

    noise_rows = []
    audit_curves = {}
    time_grid = model.sample_grid[len(model.ar)-1:]
    early = (time_grid >= 132.) & (time_grid < 164.)
    for i, label in truth.items():
        if label["condition"] != "weak_high":
            continue
        with np.load(OUTPUT / "inputs" / f"{i:02}.npz") as data:
            physical = data["physical"].copy(); original = data["reference_factor"].copy()
        measured = model.transform(physical).var(0, ddof=1)
        mean_variances = {}
        for method in METHODS:
            variances = []
            for row in results[i]["rows"]:
                if row["domain"] != "physical":
                    continue
                fit = row["fits"][method]
                if method == "reference_gls":
                    factor = original
                else:
                    energy = fit["energy_MeV"] if method == "energy_nll" else fit["covariance_energy_MeV"]
                    factor = candidate_factor(model, label["times_ns"], energy, legacy=method == "all_legacy_once")
                variances.append(diagonal_factor(factor)[0]**2)
            predicted = np.mean(variances, axis=0)
            mean_variances[method] = predicted
            noise_rows.append(dict(case=i, block=label["block"], method=method, physical_repeats=8,
                interval_start_ns=132., interval_end_exclusive_ns=164.,
                measured_variance=float(measured[early].mean()), predicted_variance=float(predicted[early].mean()),
                measured_over_predicted=float(measured[early].mean()/predicted[early].mean()),
                candidate_energy_from_waveform=method != "reference_gls"))
        audit_curves[i] = (measured, mean_variances)

    for name, rows in (("frames", frames), ("components", components), ("summary", summaries),
                       ("parent-summary", parents), ("paired", paired), ("noise-audit", noise_rows)):
        csv_table(name, rows)
    by_group = {(r["condition"],r["gap_ns"],r["domain"],r["method"]):r for r in summaries}
    plt.rcParams.update({"font.size":10, "axes.spines.top":False, "axes.spines.right":False})
    fig, axes = plt.subplots(2, 1, figsize=(12.5, 8), sharex=True)
    x = np.arange(6); width = .15
    for panel, (ax, domain) in enumerate(zip(axes, DOMAINS)):
        for j, method in enumerate(METHODS):
            count = [by_group[c,g,domain,method]["success10"] for c,g in CONDITIONS]
            bars = ax.bar(x+(j-2)*width, np.asarray(count)/32*100, width, color=COLORS[method], label=NAMES[method])
            ax.bar_label(bars, labels=[str(n) for n in count], fontsize=8, padding=2)
        ax.set_ylim(0, 116); ax.set_ylabel("All energies within ±10% (%)")
        ax.set_title(f"({'ab'[panel]}) "+("Reused shared-cell digitizations" if panel == 0 else "New energy-covariance Gaussian controls"))
        ax.grid(axis="y", alpha=.2); ax.set_axisbelow(True)
    axes[0].legend(ncol=5, loc="upper center", bbox_to_anchor=(.5,1.24), fontsize=9)
    axes[1].set_xticks(x, LABELS)
    fig.suptitle("Figure 0007-46-1: Known K/t energy recovery; labels are successes out of 32", y=1.03)
    finish(fig, "0007-46-1_energy_recovery")

    chosen = ("reference_gls", "energy_once", "energy_nll")
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 7.8))
    for index, (ax, (condition, gap)) in enumerate(zip(axes.flat, CONDITIONS[2:])):
        distributions = [[r["relative_error"]*100 for r in comp_groups[condition,gap,"physical",m] if r["label"] == "weak"] for m in chosen]
        boxes = ax.boxplot(distributions, tick_labels=[NAMES[m] for m in chosen], patch_artist=True, showfliers=True)
        for box, m in zip(boxes["boxes"], chosen):
            box.set_facecolor(COLORS[m]); box.set_alpha(.55)
        for j, m in enumerate(chosen, 1):
            block_means = [r["weak_bias"]*100 for r in parents if r["condition"]==condition and r["gap_ns"]==gap and r["domain"]=="physical" and r["method"]==m]
            ax.scatter(j+np.linspace(-.10,.10,4), block_means, c="black", s=17, zorder=3)
        ax.axhspan(-10,10,color="green",alpha=.10); ax.axhline(0,c="black",lw=.6)
        ax.set_title(f"({'abcd'[index]}) {LABELS[index+2].replace(chr(10),' ')}")
        ax.set_ylabel("Weak deposited-energy error (%)")
    fig.suptitle("Figure 0007-47-1: All physical weak-energy errors; dots = four mother-block means; separate y scales")
    finish(fig, "0007-47-1_weak_energy")

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.2))
    for j, method in enumerate(METHODS):
        values = [r["measured_over_predicted"] for r in noise_rows if r["method"] == method]
        axes[0].plot([0,10,20,30], values, "o-", color=COLORS[method], label=NAMES[method])
    axes[0].axhline(1,c="black",ls="--",lw=.8); axes[0].set_xlabel("Optical mother block")
    axes[0].set_ylabel("Measured / predicted early FIR variance")
    axes[0].set_title("(a) Weak→strong 40 ns, 132 ≤ t < 164 ns")
    measured, means = audit_curves[5]
    zoom = (time_grid >= 108) & (time_grid <= 210)
    axes[1].plot(time_grid[zoom], measured[zoom], c="black", lw=1, label="8-repeat measured variance")
    for method in chosen:
        axes[1].plot(time_grid[zoom], means[method][zoom], color=COLORS[method], label=NAMES[method])
    axes[1].set_yscale("log"); axes[1].set_xlabel("Time from window start (ns)")
    axes[1].set_ylabel("FIR variance (mV²)"); axes[1].set_title("(b) Block 0: explicit front-region zoom")
    axes[0].legend(fontsize=8); axes[1].legend(fontsize=8)
    fig.suptitle("Figure 0007-48-1: Noise coverage using fitted energies; eight repeats per block, no population CI")
    finish(fig, "0007-48-1_noise_coverage")

    # Fixed before scoring: block0, weak-first40ns, first physical repeat.
    i, repeat = 5, 0
    label = truth[i]
    row = next(r for r in results[i]["rows"] if r["domain"] == "physical" and r["repeat"] == repeat)
    with np.load(OUTPUT / "inputs" / f"{i:02}.npz") as data:
        observation = data["physical"][repeat].copy()
    with np.load(OLD_OUT / "cases" / f"{i:02}.npz") as data:
        true_components = data["components"][16+repeat].copy()
    fig, axes = plt.subplots(4, 1, figsize=(12.5, 10), sharex=True)
    axes[0].plot(model.sample_grid, observation, c="black", lw=.9, label="Observed total")
    predictions = {}
    for method in chosen:
        curves, _ = model.components(label["times_ns"], row["fits"][method]["energy_MeV"])
        predictions[method] = curves
        total = curves.sum(0)+model.baseline+model.background
        axes[0].plot(model.sample_grid, total, c=COLORS[method], lw=1, label=NAMES[method])
        axes[3].plot(model.sample_grid, observation-total, c=COLORS[method], lw=.8, label=NAMES[method])
    for j, ax in enumerate(axes[1:3]):
        ax.plot(model.sample_grid, true_components[j], c="black", lw=.9, label="True contribution in this digitization")
        for method in chosen:
            ax.plot(model.sample_grid, predictions[method][j], c=COLORS[method], lw=1,
                    label=f"{NAMES[method]}: {row['fits'][method]['energy_MeV'][j]:.4f} MeV")
        ax.set_title(f"({'bc'[j]}) {label['labels'][j].capitalize()} component; true deposited E={label['energy_MeV'][j]:.4f} MeV")
    axes[0].set_title("(a) Total waveform"); axes[3].set_title("(d) Observed total minus fitted total")
    for ax in axes:
        ax.set_xlim(0,2048); ax.set_ylabel("Voltage (mV)"); ax.legend(fontsize=8)
    axes[-1].set_xlabel("Time from window start (ns)")
    fig.suptitle("Figure 0007-49-1: Full 0–2.048 µs window; preselected block0, weak→strong40ns, repeat0")
    finish(fig, "0007-49-1_full_waveform")

    summary = dict(scope=load(EXPERIMENT/"protocol.json")["scope"],
        scored_frames=384, scored_method_rows=len(frames), physical_baseline_reproduced=baseline_checks,
        optical_pairs=4, new_gaussian_frames=192, true_E_Q_inference=False, known_K_t=True,
        count_and_PSD_metrics_assessed=False, statistics="Four reused mother blocks; no generalization/CI claim",
        primary="energy_once versus reference_gls", summary=summaries, paired=paired, noise_audit=noise_rows,
        not_converged={m:sum(r["not_converged"] for r in frames if r["method"]==m) for m in METHODS},
        fgls_outer_not_stable=sum(r["outer_not_stable"] for r in frames),
        median_seconds_per_wave_all_methods=float(np.median(timings)),
        figures=["0007-46-1", "0007-47-1", "0007-48-1", "0007-49-1"])
    save(EXPERIMENT / "result-summary.json", summary)
    save(RECORD / "report-validation.json", dict(passed=True, scored_method_rows=len(frames),
        summary_groups=len(summaries), physical_baseline_reproduced=baseline_checks,
        all_known_times_preserved=True, no_fit_failure_dropped=True,
        figures={p.name:digest(p) for p in OUTPUT.glob("0007-*") if p.suffix in (".png", ".svg")}))
    print(json.dumps(dict(scored=len(frames), baseline_reproduced=baseline_checks,
                         not_converged=summary["not_converged"],
                         physical=[r for r in summaries if r["domain"]=="physical"]), indent=2))


if __name__ == "__main__":
    main()
