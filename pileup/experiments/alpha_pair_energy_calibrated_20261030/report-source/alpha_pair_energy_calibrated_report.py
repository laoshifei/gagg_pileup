"""Independent scoring of frozen energy-mean/noise ablation; truth audit only."""
from __future__ import annotations

from collections import defaultdict
import csv
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .alpha_pair_energy_calibrated import (METHODS, BandedModes, CalibratedEnergyResponse, learned_metric)
from .alpha_pair_energy_covariance import candidate_factor
from .alpha_pair_energy_calibrated_study import (EXPERIMENT,OUTPUT,RECORD,CONDITIONS,load,save,digest,verify)
from .alpha_pair_recovery_energy import RecoveryEnergyResponse

NAMES=dict(zip(METHODS,["Old mean / old noise","New mean / old noise","Old mean / new noise","New mean / new noise"]))
COLORS=dict(zip(METHODS,["#707070","#e49427","#3984bb","#27845b"]))
LABELS=["Weak single","Low single","Mid single","High single","High→weak\n8 ns",
        "High→weak\n40 ns","High→weak\n200 ns","Weak→high\n40 ns"]


def table(name,rows):
    p=EXPERIMENT/"tables"/f"{name}.csv";p.parent.mkdir(parents=True,exist_ok=True)
    with p.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def finish(fig,name):
    fig.tight_layout()
    for ext in ("png","svg"):fig.savefig(OUTPUT/f"{name}.{ext}",dpi=170,bbox_inches="tight")
    plt.close(fig)


def oracle_diagnostic(model,metric,voltage,t,e):
    curves,_=model.components(t,e);total=curves.sum(0)+model.baseline+model.background
    residual=metric.apply(model.transform(voltage-total));columns=[]
    for j in range(len(e)):
        up=np.asarray(e).copy();down=up.copy();up[j]*=np.exp(1e-4);down[j]*=np.exp(-1e-4)
        derivative=(model.components(t,up)[0].sum(0)-model.components(t,down)[0].sum(0))/2e-4
        columns.append(metric.apply(model.transform(derivative)))
    jac=np.array(columns).T;inverse=np.linalg.pinv(jac,rcond=1e-10)
    projection=inverse@residual;covariance=inverse@inverse.T
    singular=np.linalg.svd(jac,compute_uv=False)
    return dict(white_second_moment=float(np.mean(residual**2)),
        local_condition=float(singular[0]/singular[-1]),
        projected_logE=projection.tolist(),model_local_logE_sd=np.sqrt(np.diag(covariance)).tolist())


def main():
    verify(True)
    for rel,value in load(RECORD/"result-freeze.json")["files"].items():assert digest(OUTPUT/rel)==value
    labels=load(OUTPUT/"labels.json");base=RecoveryEnergyResponse.load(EXPERIMENT)
    mean=CalibratedEnergyResponse.load(base,EXPERIMENT/"mean-correction.npz")
    noise=load(EXPERIMENT/"noise-calibration.json")
    with np.load(OUTPUT/"inputs/observations.npz") as f:voltage=f["voltage_mV"]
    with np.load(OUTPUT/"labels.components.npz") as f:truth_components=f["components_mV"]
    frames=[];components=[];diagnostics=[];lookup={};results={}
    for label in labels:
        i=label["id"];result=load(OUTPUT/"cases"/f"{i:04}.json");results[i]=result
        assert result["id"]==i and not result["pilot"] and result["worker_truth_read_guard"]
        e=np.array(label["energy_MeV"]);t=label["times_ns"]
        metric_old=BandedModes(candidate_factor(base,t,e));metric_new=learned_metric(mean,t,e,noise)
        for method,model,metric in ((METHODS[0],base,metric_old),(METHODS[1],mean,metric_old),
                                   (METHODS[2],base,metric_new),(METHODS[3],mean,metric_new)):
            fit=result["fits"][method];assert fit["times_ns"]==t
            error=np.asarray(fit["energy_MeV"])/e-1.;assert np.isfinite(error).all()
            common=dict(id=i,condition=label["condition"],gap_ns=label["gap_ns"],block=label["block"],method=method)
            frame=dict(common,success5=int(np.all(abs(error)<=.05)),success10=int(np.all(abs(error)<=.1)),
                success20=int(np.all(abs(error)<=.2)),boundary=int(fit["boundary"]),not_converged=int(not fit["converged"]),
                all_Q_ge_500=bool(np.all(np.asarray(label["charge_pe"])>=500.)),
                evaluations=fit["total_evaluations"],cost=fit["cost"])
            frames.append(frame);lookup[i,method]=frame
            names=label["condition"].split("_")
            for j,err in enumerate(error):
                components.append(dict(common,component=j,particle="alpha",strength_group=names[j],
                    true_E_dep_MeV=float(e[j]),predicted_E_dep_MeV=fit["energy_MeV"][j],relative_error=float(err)))
            # True E/t are only used by this separate post-fit diagnostic.
            d=oracle_diagnostic(model,metric,voltage[i],t,e)
            for j in range(len(e)):
                diagnostics.append(dict(common,component=j,strength_group=names[j],
                    white_second_moment=d["white_second_moment"],local_condition=d["local_condition"],
                    projected_logE=d["projected_logE"][j],model_local_logE_sd=d["model_local_logE_sd"][j],
                    oracle_true_E_t=True))
        if (i+1)%80==0:print(f"Scored {i+1}/480",flush=True)
    assert len(frames)==1920 and len(components)==2560
    groups=defaultdict(list);cg=defaultdict(list);dg=defaultdict(list)
    for r in frames:groups[r["condition"],r["gap_ns"],r["method"]].append(r)
    for r in components:cg[r["condition"],r["gap_ns"],r["method"],r["strength_group"]].append(r)
    for r in diagnostics:dg[r["condition"],r["gap_ns"],r["method"],r["strength_group"]].append(r)
    summaries=[];paired=[];errors=[];audit=[]
    for condition,gap in CONDITIONS:
        for method in METHODS:
            rows=groups[condition,gap,method]
            common=dict(condition=condition,gap_ns=gap,method=method)
            assert len(rows)==(80 if gap is None else 40)
            summary=dict(common,frames=len(rows),
                **{k:sum(r[k] for r in rows) for k in ("success5","success10","success20","boundary","not_converged")},
                qualified_frames=sum(r["all_Q_ge_500"] for r in rows),
                qualified_success10=sum(r["success10"] for r in rows if r["all_Q_ge_500"]))
            summaries.append(summary)
            if method!=METHODS[0]:
                before=[lookup[r["id"],METHODS[0]] for r in rows]
                paired.append(dict(common,frames=len(rows),old_success=sum(r["success10"] for r in before),new_success=summary["success10"],
                    gained=sum(not a["success10"] and b["success10"] for a,b in zip(before,rows)),
                    lost=sum(a["success10"] and not b["success10"] for a,b in zip(before,rows))))
            for strength in condition.split("_"):
                er=np.array([r["relative_error"] for r in cg[condition,gap,method,strength]])
                errors.append(dict(common,strength_group=strength,components=len(er),bias=float(er.mean()),
                    median=float(np.median(er)),sigma68=float(np.diff(np.percentile(er,[16,84]))[0]/2),
                    rmse=float(np.sqrt(np.mean(er**2)))))
                dr=dg[condition,gap,method,strength];projections=np.array([r["projected_logE"] for r in dr])
                sd=np.array([r["model_local_logE_sd"] for r in dr])
                audit.append(dict(common,strength_group=strength,frames=len(dr),
                    white_second_moment=float(np.mean([r["white_second_moment"] for r in dr])),
                    projection_bias=float(projections.mean()),projection_sd=float(projections.std(ddof=1)),
                    predicted_sd_rms=float(np.sqrt(np.mean(sd**2))),
                    measured_to_predicted_projection_variance=float(projections.var(ddof=1)/np.mean(sd**2)),
                    interpretation="True E/t local linear diagnostic across mothers; includes residual mean mismatch; not repeated fixed-mother noise or confidence coverage"))
    for name,rows in (("frames",frames),("components",components),("summary",summaries),("paired",paired),
                       ("energy-errors",errors),("oracle-diagnostics",diagnostics),("oracle-summary",audit)):table(name,rows)
    by={(r["condition"],r["gap_ns"],r["method"]):r for r in summaries}
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    fig,axes=plt.subplots(2,1,figsize=(12,7.5))
    for panel,ax in enumerate(axes):
        subset=CONDITIONS[panel*4:panel*4+4];n=80 if panel==0 else 40;x=np.arange(4);width=.19
        for j,method in enumerate(METHODS):
            counts=[by[c,g,method]["success10"] for c,g in subset]
            bars=ax.bar(x+(j-1.5)*width,np.array(counts)/n*100,width,color=COLORS[method],label=NAMES[method])
            ax.bar_label(bars,labels=[f"{v}/{n}" for v in counts],fontsize=8,padding=2)
        ax.set_xticks(x,LABELS[panel*4:panel*4+4]);ax.set_ylim(0,115);ax.set_ylabel("All energies within ±10% (%)")
        ax.set_title(f"({'ab'[panel]}) "+("Single-pulse controls" if panel==0 else "Known count/time pairs"));ax.grid(axis="y",alpha=.2)
    axes[0].legend(ncol=2,fontsize=9,loc="upper right")
    fig.suptitle("Figure 0007-52-1: Paired mean/noise ablation on retrospective disjoint mothers")
    finish(fig,"0007-52-1_energy_recovery")

    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for j,(ax,(condition,gap)) in enumerate(zip(axes.flat,CONDITIONS[4:])):
        data=[[r["relative_error"]*100 for r in cg[condition,gap,m,"weak"]] for m in METHODS]
        boxes=ax.boxplot(data,tick_labels=["Old/old","New/old","Old/new","New/new"],patch_artist=True,showfliers=True)
        for box,m in zip(boxes["boxes"],METHODS):box.set_facecolor(COLORS[m]);box.set_alpha(.6)
        ax.axhspan(-10,10,color="green",alpha=.12);ax.axhline(0,c="black",lw=.6)
        ax.set_ylabel("Weak E_dep error (%)");ax.set_title(f"({'abcd'[j]}) {LABELS[j+4].replace(chr(10),' ')}")
    fig.suptitle("Figure 0007-53-1: All weak-alpha errors, including boundaries; independent y scales")
    finish(fig,"0007-53-1_weak_errors")

    fig,axes=plt.subplots(1,2,figsize=(12,4.5));x=np.arange(4)
    for method in METHODS:
        rows=[next(r for r in audit if (r["condition"],r["gap_ns"],r["method"],r["strength_group"])==(c,g,method,"weak")) for c,g in CONDITIONS[4:]]
        axes[0].plot(x,[r["projection_bias"]*100 for r in rows],"o-",label=NAMES[method],c=COLORS[method])
        axes[1].plot(x,[r["measured_to_predicted_projection_variance"] for r in rows],"o-",c=COLORS[method])
    axes[0].axhline(0,c="black",ls="--",lw=.8);axes[1].axhline(1,c="black",ls="--",lw=.8)
    axes[0].set_title("(a) Local weak-energy mean projection");axes[0].set_ylabel("Mean projected delta log(E) ×100")
    axes[1].set_title("(b) Projection variance / model prediction");axes[1].set_ylabel("Ratio (log scale)");axes[1].set_yscale("log")
    for ax in axes:ax.set_xticks(x,LABELS[4:]);ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle("Figure 0007-54-1: Post-fit oracle E/t diagnostics; linear projections are not fitted errors")
    finish(fig,"0007-54-1_oracle_projection")

    fig,axes=plt.subplots(3,2,figsize=(13,10),sharex=True)
    examples=[]
    for col,condition in enumerate(("high_weak","weak_high")):
        label=next(r for r in labels if r["condition"]==condition and r["gap_ns"]==40.);i=label["id"]
        examples.append(dict(id=i,condition=condition,selection="First manifest row at40ns; preselected before accuracy scoring"))
        grid=base.sample_grid;k=len(label["energy_MeV"]);weak=condition.split("_").index("weak")
        axes[0,col].plot(grid,voltage[i],c="black",lw=.8,label="Observed voltage")
        axes[1,col].plot(grid,truth_components[i,weak],c="black",lw=1,label="True weak component")
        for method,model in ((METHODS[0],base),(METHODS[3],mean)):
            fit=results[i]["fits"][method];curves,_=model.components(fit["times_ns"],fit["energy_MeV"])
            full=curves.sum(0)+model.baseline+model.background
            axes[0,col].plot(grid,full,c=COLORS[method],ls="--",label=NAMES[method])
            axes[1,col].plot(grid,curves[weak],c=COLORS[method],ls="--",label=NAMES[method])
            axes[2,col].plot(grid,voltage[i]-full,c=COLORS[method],lw=.8,label=NAMES[method])
        for row in range(3):
            axes[row,col].set_ylabel("Voltage (mV)");axes[row,col].grid(alpha=.2);axes[row,col].legend(fontsize=8)
            axes[row,col].set_title(f"({'abcdef'[row*2+col]}) {condition.replace('_','→')},40 ns; "+["full waveform","weak component","residual"][row])
        axes[2,col].set_xlabel("Time (ns)");axes[2,col].set_xlim(0,2048)
    fig.suptitle("Figure 0007-55-1: Full observation windows; fixed manifest examples, known K/t")
    finish(fig,"0007-55-1_full_waveforms")
    save(OUTPUT/"report.json",dict(summary=summaries,paired=paired,energy_errors=errors,oracle_summary=audit,examples=examples))
    save(RECORD/"report-validation.json",dict(passed=True,frames=480,method_rows=len(frames),component_rows=len(components),
        independent_new_calibration_mothers=True,retrospective_evaluation=True,known_K_t=True,
        nonconverged=sum(r["not_converged"] for r in frames),boundary=sum(r["boundary"] for r in frames),
        figures=["0007-52-1","0007-53-1","0007-54-1","0007-55-1"],count_and_type_metrics="Not measured"))
    print(json.dumps(summaries,ensure_ascii=False,indent=2),flush=True)


if __name__=="__main__":main()
