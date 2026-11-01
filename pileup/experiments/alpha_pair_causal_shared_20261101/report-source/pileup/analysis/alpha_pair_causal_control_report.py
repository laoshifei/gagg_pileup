"""Standalone figures for the isolated and shared causal-response controls."""
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .alpha_pair_causal_controls import compatibility
from .alpha_pair_causal_recovery import CausalRecovery
from .alpha_pair_energy_calibrated import CalibratedEnergyResponse
from .alpha_pair_recovery_energy import RecoveryEnergyResponse

ROOT = Path(os.environ.get("GAGG_CAUSAL_REPORT_ROOT", Path(__file__).resolve().parents[2])).resolve()
SINGLE = ROOT / "jupyter/artifacts/alpha_pair_causal_controls_20261101"
SHARED = ROOT / "jupyter/artifacts/alpha_pair_causal_shared_20261101"
EXPERIMENT = ROOT / "pileup/experiments/alpha_pair_causal_shared_20261101"


def read(path): return json.loads(Path(path).read_text())


def save(fig, name):
    fig.savefig(SHARED / (name + ".png"), dpi=170, bbox_inches="tight")
    fig.savefig(SHARED / (name + ".svg"), bbox_inches="tight")
    plt.close(fig)


def main():
    compatibility()
    plt.rcParams.update({"font.size":10, "axes.spines.top":False, "axes.spines.right":False})
    groups=read(SINGLE/"single-report.json")["groups"]
    rows=read(SINGLE/"single-rows.json")
    conditions=[g["condition"] for g in groups]
    phases=np.array([0.,.5,1.,1.5])
    fig, axes=plt.subplots(1,2,figsize=(10,3.8))
    ax=axes[0];ax.axhspan(-5,5,color="#d9eedb",alpha=.5)
    ax.bar(conditions,[g["energy_bias_percent"] for g in groups],color="#35699b",alpha=.8)
    for i,condition in enumerate(conditions):
        values=[np.mean([r["energy_error_percent"] for r in rows if r["condition"]==condition and r["block"]==b]) for b in (0,10,20,30)]
        ax.scatter(i+np.linspace(-.13,.13,4),values,color="#343434",s=20,zorder=4)
    ax.axhline(0,color="black",lw=.7);ax.set_ylabel("Mean-waveform energy bias (%)")
    ax.set_title("(a) Four optical mothers per group")
    for g in groups:axes[1].plot(phases,g["phase_bias_percent"],"o-",label=g["condition"])
    axes[1].set(xlabel="Arrival phase within 2 ns sampling period (ns)",ylabel="Mean energy bias (%)",title="(b) Same mothers and seeds across phases")
    axes[1].legend(ncol=2);fig.suptitle("Figure 0007-70-1: Isolated response and sampling-phase controls")
    fig.tight_layout();save(fig,"0007-70-1_isolated_energy")

    anchor=CalibratedEnergyResponse.load(RecoveryEnergyResponse.load(EXPERIMENT),EXPERIMENT/"mean-correction.npz")
    with np.load(EXPERIMENT/"charge-map-1.0.npz") as f:mapping={k:f[k] for k in f.files}
    model=CausalRecovery(anchor,read(EXPERIMENT/"sensor-config.json"),EXPERIMENT/"bin/age.so",mapping)
    time=model.readout.adc_centers
    fig,axes=plt.subplots(4,2,figsize=(11,10))
    for index,condition in enumerate(conditions):
        selected=[r for r in rows if r["condition"]==condition and r["phase_ns"]==0.]
        physical=[];predicted=[]
        for r in selected:
            with np.load(SINGLE/"single"/f"{r['id']:04d}.npz") as f:physical.append(f["event_signal_mV"].mean(0))
            predicted.append(model.isolated([r["time_ns"]],[r["energy_dep_MeV"]])[0][0])
        physical=np.array(physical);predicted=np.array(predicted)
        mean=physical.mean(0);sem=physical.std(0,ddof=1)/2
        for ax in axes[index]:
            ax.plot(time,mean,label="Digitized event mean",color="#246b88")
            ax.fill_between(time,mean-sem,mean+sem,alpha=.18,color="#246b88",label="Mother mean +/- SEM")
            ax.plot(time,predicted.mean(0),"--",label="Frozen causal anchor",color="#c65b28")
            ax.axvline(128,color="gray",ls=":",lw=.8);ax.set_ylabel(condition+" (mV)")
        axes[index,0].set_xlim(0,2048);axes[index,1].set_xlim(110,230)
    axes[0,0].set_title("(a) Full declared observation window")
    axes[0,1].set_title("(b) Leading-edge detail")
    for ax in axes[-1]:ax.set_xlabel("Time (ns)")
    axes[0,0].legend(fontsize=8)
    fig.suptitle("Figure 0007-71-1: Isolated mean voltage at true calibration coordinates")
    fig.tight_layout();save(fig,"0007-71-1_isolated_waveforms")

    groups=read(SHARED/"report.json")["groups"]
    rows=read(SHARED/"rows.json")
    labels=[("H-W" if g["condition"]=="high_weak" else "W-H")+f"\n{g['gap_ns']:g} ns" for g in groups]
    x=np.arange(len(groups));fig,axes=plt.subplots(1,2,figsize=(11,4))
    axes[0].bar(x-.17,[g["actual_weak_area_change_percent"] for g in groups],.32,label="Digitized mean",color="#246b88")
    axes[0].bar(x+.17,[g["predicted_weak_area_change_percent"] for g in groups],.32,label="Uniform model",color="#c65b28")
    for i,g in enumerate(groups):
        selected=[r for r in rows if (r["condition"],r["gap_ns"])==(g["condition"],g["gap_ns"])]
        axes[0].scatter(i-.17+np.linspace(-.07,.07,4),[r["actual_weak_area_change_percent"] for r in selected],s=15,color="black",zorder=5)
    axes[0].plot(x,[g["recovery_off_area_change_percent"] for g in groups],"o",label="Recovery off",color="#37965a")
    axes[0].set_xticks(x,labels);axes[0].set_ylabel("Weak-component area change (%)")
    axes[0].set_title("(a) Shared relative to isolated, same window");axes[0].legend(fontsize=8)
    axes[1].axhspan(-5,5,color="#d9eedb",alpha=.5)
    axes[1].plot(x,[g["prediction_error_percentage_points"] for g in groups],"o-",label="Model minus physical",color="#c65b28")
    axes[1].plot(x,[g["grid_change_percentage_points"] for g in groups],"s-",label="1 ns minus 0.5 ns grid",color="#71529c")
    axes[1].set_xticks(x,labels);axes[1].set_ylabel("Difference (percentage points)")
    axes[1].set_title("(b) Prespecified +/-5 point mean-response gate")
    axes[1].legend(fontsize=8)
    for ax in axes:ax.axhline(0,color="black",lw=.7)
    fig.suptitle("Figure 0007-72-1: Shared-cell suppression and numerical-grid check")
    fig.tight_layout();save(fig,"0007-72-1_shared_response")

    data=np.load(SHARED/"display-curves.npz")
    fig,axes=plt.subplots(2,2,figsize=(11,7))
    for row,(identifier,weak,label) in enumerate([(1,1,"High then weak, 40 ns"),(4,0,"Weak then high, 40 ns")]):
        actual=data[f"case{identifier}_shared"];solo=data[f"case{identifier}_isolated"]
        predicted=data[f"case{identifier}_model_shared"];isolated=data[f"case{identifier}_model_isolated"]
        axes[row,0].plot(time,actual.sum(0),label="Digitized shared mean",color="#246b88")
        axes[row,0].plot(time,predicted.sum(0),"--",label="Causal uniform model",color="#c65b28")
        axes[row,1].plot(time,actual[weak]-solo[weak],label="Digitized shared - isolated",color="#246b88")
        axes[row,1].plot(time,predicted[weak]-isolated[weak],"--",label="Model correction",color="#c65b28")
        for ax in axes[row]:ax.set_xlim(0,2048);ax.set_xlabel("Time (ns)");ax.legend(fontsize=8)
        axes[row,0].set_ylabel(label+"\ntotal event signal (mV)")
        axes[row,1].set_ylabel("Weak-component change (mV)")
    axes[0,0].set_title("(a) Full signal window")
    axes[0,1].set_title("(b) Full weak-component correction")
    fig.suptitle("Figure 0007-73-1: Forward diagnostics, known E/t and ancestry")
    fig.tight_layout();save(fig,"0007-73-1_shared_waveforms")
    print("Created four numbered PNG/SVG figures",flush=True)


if __name__=="__main__":main()
