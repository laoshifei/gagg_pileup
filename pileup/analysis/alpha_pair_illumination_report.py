"""Response summaries and paired uncertainty for position interventions.

Photon truth is used only for this mechanistic diagnostic. Resampling keeps
arms paired and treats mother pairs, layouts and digital seeds as nested units.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from . import alpha_pair_illumination_controls as study
from .alpha_pair_causal_controls import read, save, sha


def percentile(values):
    return np.percentile(values, [2.5, 97.5]).tolist()


def bootstrap(data, predictions, seed):
    """data[mother,arm,recovery,shared_or_solo,layout,digital_repeat]."""
    rng = np.random.default_rng(seed)
    n = 10000
    mothers = rng.integers(0, 4, (n, 4))
    layouts = rng.integers(0, 4, (n, 4, 4))
    repeats = rng.integers(0, 8, (n, 4, 4, 8))
    result = np.empty((n, 3, 2))
    for arm in range(3):
        for recovery in range(2):
            sample = []
            for kind in range(2):
                a = data[:, arm, recovery, kind]
                sampled = a[mothers[:, :, None, None], layouts[:, :, :, None], repeats]
                sample.append(sampled.mean(axis=(2, 3)))
            result[:, arm, recovery] = (100 * (sample[0] / sample[1] - 1)).mean(axis=1)
    return result, {k: np.asarray(v)[mothers].mean(1) for k, v in predictions.items()}


def validate_resampling():
    data = np.full((4,3,2,2,4,8),100.)
    perturbation = np.random.default_rng(17).normal(size=(4,4,8))
    for arm, level in enumerate((85.,88.,92.)):
        data[:,arm,0,0] = level+perturbation
    boot, pred = bootstrap(data, {"constant":np.full(4,-8.)}, 29)
    assert boot.shape == (10000,3,2)
    assert np.max(np.abs(boot[:,1,0]-boot[:,0,0]-3.)) < 1e-12
    assert np.max(np.abs(boot[:,2,0]-boot[:,0,0]-7.)) < 1e-12
    assert np.array_equal(boot[:,:,1],np.zeros((10000,3)))
    assert np.array_equal(pred["constant"],np.full(10000,-8.))
    save(study.REC / "resampling-check.json",dict(passed=True, replicates=10000,
         paired_arm_offsets_preserved_pp=[3.,7.], neutral_control_exact_zero=True,
         constant_model_predictions_preserved=True))
    print("Nested paired resampling controls passed",flush=True)


def report():
    count = study.verify()
    cases = read(study.EXP / "cases.json")
    predicted = read(study.EXP / "predictions.json")["rows"]
    # State0:on,state1:off; type0:shared,type1:isolated. Weak event only.
    areas = np.empty((20, 3, 2, 2, 4, 8))
    rows, display = [], {}
    files = []
    for case in cases:
        cid = case["id"]
        weak = case["parents"].index(case["block"])
        strong = 1 - weak
        layouts = []
        for layout in range(4):
            p = study.ART / "cases" / f"{cid*4+layout:04d}.npz"
            files.append(p)
            with np.load(p) as f:
                arrays = {k: f[k].copy() for k in f.files}
            assert all(np.isfinite(v).all() for v in arrays.values())
            layouts.append(arrays)
        original_solo = np.stack([f["original_on_isolated"] for f in layouts]).mean((0, 1))
        for arm_index, arm in enumerate(study.ARMS):
            for state, prefix in enumerate(("on", "off")):
                for kind, suffix in enumerate(("shared", "isolated")):
                    signals = np.stack([f[f"{arm}_{prefix}_{suffix}"] for f in layouts])
                    assert signals.shape == (4, 8, 2, 1024)
                    areas[cid, arm_index, state, kind] = signals[:, :, weak].sum(-1)
                    if cid < 5:
                        display[f"case{cid}_{arm}_{prefix}_{suffix}"] = signals.mean((0, 1))
            on = areas[cid, arm_index, 0].mean((1, 2))
            off = areas[cid, arm_index, 1].mean((1, 2))
            assert on[1] > 0 and off[1] > 0
            solo = np.stack([f[f"{arm}_on_isolated"] for f in layouts]).mean((0, 1))
            shared_q = np.stack([f[f"{arm}_on_shared_charge"] for f in layouts]).mean((0, 1))[weak]
            isolated_q = np.stack([f[f"{arm}_on_isolated_charge"] for f in layouts]).mean((0, 1))[weak]
            point = float(100 * (on[0] / on[1] - 1))
            rows.append(dict(case=cid, block=case["block"], condition=case["condition"], gap_ns=case["gap_ns"],
                arm=arm, recovery_on_percent=point, recovery_off_percent=float(100*(off[0]/off[1]-1)),
                full_charge_change_percent=float(100*(shared_q/isolated_q-1)),
                isolated_weak_area_gain_percent=float(100*(solo[weak].sum()/original_solo[weak].sum()-1)),
                isolated_strong_area_gain_percent=float(100*(solo[strong].sum()/original_solo[strong].sum()-1)),
                legacy_error_pp=float(predicted[cid]["legacy_1.0_percent"]-point),
                photon_error_pp=float(predicted[cid]["photon_1.0_percent"]-point)))
            if cid < 5:
                display[f"case{cid}_{arm}_analog"] = np.stack([f[f"{arm}_on_analog"] for f in layouts]).mean((0, 1))
    np.savez_compressed(study.ART / "weak-areas.npz", areas=areas)
    np.savez_compressed(study.ART / "display-curves.npz", **display)
    groups = []
    for group_id, (condition, gap) in enumerate((("high_weak",8.), ("high_weak",40.),
                                               ("high_weak",200.), ("high_weak",1200.), ("weak_high",40.))):
        ids = [c["id"] for c in cases if c["condition"] == condition and c["gap_ns"] == gap]
        predictions = {key: [predicted[i][key] for i in ids]
                       for key in ("legacy_1.0_percent", "photon_1.0_percent", "photon_0.5_percent")}
        boot, pred_boot = bootstrap(areas[ids], predictions, 713900000 + group_id)
        means = areas[ids].mean(axis=(-1, -2))
        change = (100*(means[:, :, :, 0]/means[:, :, :, 1]-1)).mean(0)
        arms = {}
        for a, arm in enumerate(study.ARMS):
            arm_rows = [r for r in rows if r["case"] in ids and r["arm"] == arm]
            arms[arm] = dict(on_percent=float(change[a, 0]), on_ci95=percentile(boot[:, a, 0]),
                off_percent=float(change[a, 1]), off_ci95=percentile(boot[:, a, 1]),
                legacy_error_pp=float(np.mean(predictions["legacy_1.0_percent"])-change[a, 0]),
                legacy_error_ci95=percentile(pred_boot["legacy_1.0_percent"]-boot[:, a, 0]),
                photon_error_pp=float(np.mean(predictions["photon_1.0_percent"])-change[a, 0]),
                photon_error_ci95=percentile(pred_boot["photon_1.0_percent"]-boot[:, a, 0]),
                isolated_weak_gain_percent=float(np.mean([r["isolated_weak_area_gain_percent"] for r in arm_rows])),
                isolated_strong_gain_percent=float(np.mean([r["isolated_strong_area_gain_percent"] for r in arm_rows])))
        groups.append(dict(condition=condition, gap_ns=gap, mothers=4, layouts=4, repeats=8, arms=arms,
            uniform_minus_original_pp=float(change[2,0]-change[0,0]),
            uniform_minus_original_ci95=percentile(boot[:,2,0]-boot[:,0,0]),
            shuffled_minus_original_pp=float(change[1,0]-change[0,0]),
            shuffled_minus_original_ci95=percentile(boot[:,1,0]-boot[:,0,0]),
            uniform_minus_shuffled_pp=float(change[2,0]-change[1,0]),
            uniform_minus_shuffled_ci95=percentile(boot[:,2,0]-boot[:,1,0]),
            uniform_effect_minus_off_effect_pp=float(change[2,0]-change[0,0]-change[2,1]+change[0,1]),
            uniform_effect_minus_off_effect_ci95=percentile(boot[:,2,0]-boot[:,0,0]-boot[:,2,1]+boot[:,0,1]),
            legacy_percent=float(np.mean(predictions["legacy_1.0_percent"])),
            photon_percent=float(np.mean(predictions["photon_1.0_percent"])),
            photon_grid_change_pp=float(np.mean(predictions["photon_1.0_percent"])-np.mean(predictions["photon_0.5_percent"])) ))
    gates = read(study.EXP / "protocol.json")["gates"]
    off_pass = all(abs(a["off_percent"]) <= gates["recovery_off_abs_percentage_points"]
                   for g in groups for a in g["arms"].values())
    short_pass = all(g["uniform_minus_original_pp"] >= gates["uniform_effect_at_short_gaps_percentage_points"]
                     and g["uniform_minus_original_ci95"][0] > 0 for g in groups[:2])
    grid_pass = all(abs(g["photon_grid_change_pp"]) <= gates["empirical_grid_change_percentage_points"] for g in groups)
    summary = dict(groups=groups, recovery_off_control_passed=off_pass, empirical_grid_passed=grid_pass,
        spatial_contribution_supported_on_selected_mothers=off_pass and short_pass,
        inherited_model_passes_uniform_point_gate=all(abs(g["arms"]["uniform"]["legacy_error_pp"]) <= gates["inherited_uniform_model_abs_error_percentage_points"] for g in groups),
        scope="Four historically exposed mother pairs, paired nested descriptive intervals; transformed inputs do not validate original inputs or establish independent reconstruction performance.",
        default_changed=False, energy_inference=False, new_ROOT=0,
        production_digitizations=read(study.REC / "run.json")["digitizations"])
    save(study.EXP / "response-summary.json", summary)
    save(study.ART / "report.json", summary)
    save(study.ART / "rows.json", rows)
    tables = study.EXP / "tables"
    tables.mkdir(exist_ok=True)
    with (tables / "position-controls.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with np.load(study.ART / "pilot/0004.npz") as pilot, np.load(study.ART / "cases/0004.npz") as full:
        pilot_exact = all(np.array_equal(pilot[k], full[k][:2]) for k in pilot.files)
    assert pilot_exact
    save(study.REC / "result-audit.json", dict(frozen_hashes_checked=count, finite_arrays=True,
         production_task_count=len(files), output_hashes={str(p.relative_to(study.ROOT)):sha(p) for p in files},
         pilot_first_two_repeats_exact=pilot_exact, bootstrap_resamples=10000,
         point_estimator="Equal-mother mean of each mother's ratio of means; bootstrap re-computes the same estimator with paired arms"))
    for g in groups:
        print(g["condition"], g["gap_ns"], {a: round(v["on_percent"],4) for a,v in g["arms"].items()},
              "uniform-original pp", round(g["uniform_minus_original_pp"],4), g["uniform_minus_original_ci95"], flush=True)
    print("spatial contribution supported:", summary["spatial_contribution_supported_on_selected_mothers"], flush=True)


def plots():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    study.verify()
    study.initialize()
    summary = read(study.ART / "report.json")
    groups = summary["groups"]
    colors = dict(original="#355C9B", shuffled="#D17B27", uniform="#238B67")
    names = dict(original="Original XY", shuffled="Shuffled XY", uniform="Uniform XY")
    plt.rcParams.update({"font.size":10, "axes.spines.top":False, "axes.spines.right":False})
    saved = []

    def finish(fig, name):
        for extension in ("png", "svg"):
            path = study.ART / (name + "." + extension)
            fig.savefig(path, dpi=180, bbox_inches="tight")
            saved.append(dict(path=str(path.relative_to(study.ROOT)), sha256=sha(path)))
        plt.close(fig)

    # Equal-mother spatial maps, pooled across the four randomized layouts.
    maps = {}
    for group, offset in (("Weak",0), ("Strong",240)):
        for arm in study.ARMS:
            densities = []
            for block in (0,10,20,30):
                for layout in range(4):
                    cells = study.cells_for(study.DONORS[arm,layout,block+offset], study.SENSOR)
                    counts = np.bincount(cells[cells>=0], minlength=14400).astype(float)
                    coarse = (counts/counts.sum()).reshape(20,6,20,6).sum((1,3)) / .09
                    densities.append(coarse)
            maps[group,arm] = np.mean(densities,axis=0)
    fig, axes = plt.subplots(2,3,figsize=(11,7),layout="constrained")
    vmax = max(v.max() for v in maps.values())
    spatial = read(study.EXP / "input-check.json")["spatial"]
    for row, group in enumerate(("Weak","Strong")):
        for col, arm in enumerate(study.ARMS):
            ax = axes[row,col]
            im = ax.imshow(maps[group,arm],origin="lower",extent=(-3,3,-3,3),vmin=0,vmax=vmax,cmap="magma")
            overlap = np.mean([r["overlap"] for r in spatial if r["arm"]==arm])
            ax.set(title=f"({chr(97+row*3+col)}) {group}: {names[arm]}",xlabel="x (mm)",ylabel="y (mm)")
            if row==0:
                ax.text(.04,.04,f"pair overlap = {overlap:.3f}",transform=ax.transAxes,color="white",fontsize=9)
    fig.colorbar(im,ax=axes,label="Normalized photon density (mm$^{-2}$)",shrink=.85)
    fig.suptitle("Fig. 0007-76-1: Position interventions on the same optical mothers")
    finish(fig,"0007-76-1_spatial_maps")

    fig, axes = plt.subplots(1,2,figsize=(12,4.6),layout="constrained")
    x = np.arange(5)
    labels = ["H-W 8", "H-W 40", "H-W 200", "H-W 1200", "W-H 40"]
    for index, arm in enumerate(study.ARMS):
        point = np.array([g["arms"][arm]["on_percent"] for g in groups])
        intervals = np.array([g["arms"][arm]["on_ci95"] for g in groups])
        xx = x+(index-1)*.18
        axes[0].vlines(xx,intervals[:,0],intervals[:,1],color=colors[arm],lw=2)
        axes[0].plot(xx,point,"o",color=colors[arm],label=names[arm])
    axes[0].plot(x,[g["legacy_percent"] for g in groups],"k--",label="Frozen anchored model",lw=1.2)
    axes[0].plot(x,[g["photon_percent"] for g in groups],":",color="#8F4A98",label="Truth-time uniform model",lw=1.5)
    for arm,offset in (("shuffled",-.08),("uniform",.08)):
        field=arm+"_minus_original"
        point=np.array([g[field+"_pp"] for g in groups])
        ci=np.array([g[field+"_ci95"] for g in groups])
        axes[1].vlines(x+offset,ci[:,0],ci[:,1],color=colors[arm],lw=2)
        axes[1].plot(x+offset,point,"o",color=colors[arm],label=names[arm]+" minus original")
    for ax in axes:
        ax.axhline(0,color=".6",lw=.8)
        ax.set_xticks(x,labels,rotation=20)
        ax.set_xlabel("Order and gap (ns)")
        ax.grid(axis="y",alpha=.18)
        ax.legend(fontsize=8)
    axes[0].set(title="(a) Weak-component area change",ylabel="Shared / isolated - 1 (%)")
    axes[1].set(title="(b) Paired intervention effect",ylabel="Change from original (percentage points)")
    axes[0].legend(loc="lower center",bbox_to_anchor=(.58,.01),fontsize=8)
    fig.suptitle("Fig. 0007-77-1: Recovery-on response; descriptive nested 95% intervals")
    finish(fig,"0007-77-1_response_effect")

    fig, axes = plt.subplots(1,2,figsize=(11,4),layout="constrained")
    edges = np.arange(0.,1200.1,50.)
    for ax,offset,group in zip(axes,(0,240),("Weak","Strong")):
        for arm in study.ARMS:
            numerator=np.zeros(len(edges)-1)
            denominator=numerator.copy()
            for block in (0,10,20,30):
                for layout in range(4):
                    donor=study.DONORS[arm,layout,block+offset]
                    active=study.cells_for(donor,study.SENSOR)>=0
                    center=donor["x_mm"]**2+donor["y_mm"]**2<=1.5**2
                    denominator+=np.histogram(donor["time_ns"][active],edges)[0]
                    numerator+=np.histogram(donor["time_ns"][active & center],edges)[0]
            probability=np.divide(numerator,denominator,out=np.full_like(numerator,np.nan),where=denominator>0)
            ax.plot((edges[1:]+edges[:-1])/2,probability,color=colors[arm],label=names[arm])
        ax.axhline(np.pi*1.5**2/36,color=".4",ls=":",label="Uniform area fraction")
        ax.set(xlabel="Time from first active-area photon (ns)",ylabel="Fraction inside radius 1.5 mm",ylim=(0,1),
               title=f"({'a' if offset==0 else 'b'}) {group} mothers pooled")
        ax.legend(fontsize=8)
        ax.grid(alpha=.15)
    fig.suptitle("Fig. 0007-78-1: Space-time association; arrival times unchanged")
    finish(fig,"0007-78-1_time_position")

    with np.load(study.ART / "display-curves.npz") as f:
        display={k:f[k].copy() for k in f.files}
    with np.load(study.EXP / "prediction-curves.npz") as f:
        models={k:f[k].copy() for k in f.files}
    fig,axes=plt.subplots(2,2,figsize=(12,6.8),layout="constrained")
    t=np.arange(1024)*2+.25
    for column,(case,weak,label) in enumerate(((1,1,"H-W 40 ns"),(4,0,"W-H 40 ns"))):
        for arm in study.ARMS:
            joint=display[f"case{case}_{arm}_on_shared"]
            solo=display[f"case{case}_{arm}_on_isolated"]
            axes[0,column].plot(t,joint.sum(0),color=colors[arm],label=names[arm],lw=1.2)
            axes[1,column].plot(t,joint[weak]-solo[weak],color=colors[arm],label=names[arm],lw=1.2)
        for model,color,ls,label_model in (("legacy","black","--","Frozen anchored model"),
                                           ("photon","#8F4A98",":","Truth-time uniform model")):
            joint=models[f"case{case}_{model}_joint"]
            solo=models[f"case{case}_{model}_solo"]
            axes[0,column].plot(t,joint.sum(0),color=color,ls=ls,label=label_model,lw=1)
            axes[1,column].plot(t,joint[weak]-solo[weak],color=color,ls=ls,label=label_model,lw=1)
        for row in range(2):
            axes[row,column].set(xlim=(0,2048),xlabel="Time (ns)",ylabel="Voltage (mV)")
            axes[row,column].grid(alpha=.15)
        axes[0,column].set_title(f"({chr(97+column)}) {label}: total event signal")
        axes[1,column].set_title(f"({chr(99+column)}) {label}: weak shared minus isolated")
        axes[1,column].axhline(0,color=".6",lw=.7)
    axes[0,0].legend(fontsize=8)
    fig.suptitle("Fig. 0007-79-1: Complete observed window, first fixed mother pair, 32 repeats")
    finish(fig,"0007-79-1_waveforms")
    save(study.REC / "figures.json",saved)
    print("Saved four numbered figures as PNG and SVG",flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("report", "plots", "validate"))
    args = parser.parse_args()
    if args.mode == "validate":
        validate_resampling()
    elif args.mode == "report":
        report()
    else:
        plots()


if __name__ == "__main__":
    main()
