# Detection tradeoffs and conditional later-charge profiles

网络起点＋物理判决在两个域所有空／单组额外输出帧≤5%的开发约束下，原生40 ns严格成功为0/40；即使允许所有合资格候选通过也只有3/40。放宽判决只能接纳已有正确候选，不能修正其余候选的电荷／时刻错误。ΔJ≤4的160条条件剖面中，13条触及50 pe下界、13条存在多个分离区段；0条真实Q2评分探针把参考最小代价再降低超过0.01。这些数值及有限搜索误差需与包含率一起解释，不将该层级称为95%置信区间。

| 方法 | 开发控制允许的最大α | 原生40 ns：默认→放宽 | 原生200 ns：默认→放宽 | 原生总成功 | 高斯总成功 |
| --- | ---: | ---: | ---: | ---: | ---: |
| baseline | 0.05000 | 0→0/40 | 20→20/40 | 241/360 | 501/600 |
| learned_gate | 0.04878 | 0→0/40 | 7→12/40 | 232/360 | 491/600 |
| local_physical | 0.05000 | 0→0/40 | 20→20/40 | 241/360 | 501/600 |
| mlp_physical | 0.05000 | 0→0/40 | 20→20/40 | 241/360 | 502/600 |

| 输入／间隔 | ΔJ≤4包含真实Q2 | 范围宽度／真实Q2中位 | 触及50 pe下界 | 分离区段 | 真Q探针降低参考最小值 |
| --- | ---: | ---: | ---: | ---: | ---: |
| native／40 | 36/40 | 2.693 | 7/40 | 12/40 | 0/40 |
| native／200 | 39/40 | 0.751 | 0/40 | 0/40 | 0/40 |
| continuous／middle | 37/40 | 1.246 | 6/40 | 1/40 | 0/40 |
| continuous／long | 35/40 | 0.563 | 0/40 | 0/40 | 0/40 |

## Experiment definition

Pure alpha, fixed95/351 ns, full1024×2 ns analog waveform, ample dynamic range. Original blind count fits allow K≤2; profiles explicitly assume K=2. Charge is equivalent pe, not calibrated deposited energy. All samples are previously exposed development reuse:280 native gate singles,760 native evaluation frames and1000 conditional-Gaussian frames. No new ROOT, training or random draws.

The threshold sweep reuses all four prior candidate/gate workflows.127 alpha points cover the attainable plus-one upper ranks with40/80/80/80 calibration samples. Presence,Q≥500 and positive physical evidence guards stay fixed. At alpha=.025 all7040 scored decisions exactly match the previous experiment. The observed ≤5% worst-control envelope is retrospective; no default is changed. Full errors and controls appear in threshold-summary.csv.

Profiles use160 prespecified strong-then-weak frames: native40/200 ns and Gaussian24–64/128–256 ns,40 each. For fixed chronological Q2, optimize Q1 and ordered t1/t2. Both times lie in the observed K1 trigger region[−96,+512]ns clipped to[0,2046]. Flag cached candidates outside this region.33 logarithmic50–50000 pe points plus observed-candidate points; two best observation-derived starts,50 iterations, forward/reverse continuation. The original candidate-dependent covariance and logdet are retained. Frozen-nuisance slices are included in the finite envelope. Disconnected ranges and boundary truncation are reported; global optimality is not claimed.

Levels1/4/9 are compatibility contours, not confidence levels. After primary profiles are frozen, a separate scoring probe fixes true later Q2 and refits nuisance using observed starts only. It never updates profiles, gates or default fits. Negative probe delta below−.01 flags a better point missed by the finite grid/search. Inclusion and interpolated width are empirical diagnostics, not calibrated coverage. All failed convergence flags are retained.

## Reproduction and provenance

Run from the project root with the project Python environment and OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1. Actual inference uses7 workers; the maximum is8. Input observations, candidates and truth come from the risk-training-v2/physical-gate studies; hashes, prior model/parent provenance and Git HEAD are in local input-provenance.json/environment.json. No new macro or stochastic seed is needed for this deterministic diagnostic. Earlier uncommitted work is retained.

run.py uses inference-source; report.py uses report-source; deliver.py uses delivery-source. Source, models, protocol and task indices are frozen with hashes, so no Git checkout is needed. Local arrays/models are ignored by Git and must accompany a reproducibility bundle.

```bash
bash jupyter/run_gagg_python.sh pileup/experiments/alpha_pair_evidence_profile_20261027/report.py sweep
bash jupyter/run_gagg_python.sh pileup/experiments/alpha_pair_evidence_profile_20261027/report.py analyze
bash jupyter/run_gagg_python.sh pileup/experiments/alpha_pair_evidence_profile_20261027/report.py figures
bash jupyter/run_gagg_python.sh pileup/experiments/alpha_pair_evidence_profile_20261027/deliver.py validate
```

prepare/pilot/infer/probe are one-time commands that reject existing outputs. To recompute inference, copy the frozen experiment and its inputs into a separate project tree. Notebook/document delivery refuses to replace subsequent user edits. Notebook writes use256 KiB chunks and complete readback before replacement for shared-folder compatibility.

Local artifacts: jupyter/artifacts/alpha_pair_evidence_profile_20261027/. Local records: jupyter/run_records/alpha_pair_evidence_profile_20261027/. Notebook0004 sections98–102, figures0004-98-1/98-2/99-1/100-1/101-1. All original475 cells retained; new code cells executed. Full0–2048 ns examples are prespecified; plotted profile components are conditional K2 explanations, not final blind reconstructions.

6 numerical tests,4 pilot profiles,160 primary profiles,160 true-Q2 scoring probes,7040 default-score replays and127-point sweep verified. final-validation.json records hashes, label-access probes, executed Notebook and image checks. Gate sweep involves894080 cached choices, not that many independent samples or fresh fits.

## Evidence and limitations

Liddick2012 (NIMA66970–78; DOI10.1016/j.nima.2011.12.034, siliconDSSD40MSPS) motivates separate detection/accuracy/false-double evaluation; its numerical limits do not transfer. Baker/Cousins1984 (NIM221437–442; DOI unverified in existing library note), PDF1/4/5, distinguishes estimation, coverage and goodness of fit under Poisson/multinomial histograms. No histogram likelihood is transplanted to voltage. Both scientific notes remain AI drafts pending review; source conditions and unresolved details are preserved in docs/literature-use.md L63. Only local existing notes were consulted; no PDF reread or library edits.

Before confidence claims, check grid/search stability and independent repeated-sampling calibration. No new optical parents, detector measurements, free-tau, alpha/gamma or full high-rate K≥3 validation; keep the current default reconstruction.
