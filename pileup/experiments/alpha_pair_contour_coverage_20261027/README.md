# Adaptive charge contours and fixed-model coverage pilot

5条定向数值诊断中，5/5通过新增的边界括区与更强局部重拟合检查。随后冻结方案，生成32条新条件高斯双脉冲和12条空／单控制；32/32条双脉冲通过数值筛查，全部32条都保留在主要包含率分母中。40ns组ΔJ≤4范围包含真实Q2为14/16，已知K2严格整对恢复为3/16。200ns组ΔJ≤4范围包含真实Q2为16/16，已知K2严格整对恢复为10/16。这检验的是冻结模型内、固定参数点的重复抽样表现；每组16条不足以认证95%覆盖，更不能推广为原生SiPM的盲恢复或普适最短间隔。

两条的变化来源并不一样：第160条仅同网格加强局部搜索，ΔJ≤9端点就移动6.72%，原节点最大代价下降3.021；后续加点与复核再改变约0.071%。第190条同网格端点变化可忽略，后续加点与复核改变约3.45%。因此前者明确存在原节点的局部搜索不足；后者主要在后续细化中变化，但该步骤也重新优化了未知参数，不能全归因于画图插值。

5条最终括区均小于1%，更强的末次局部搜索通过筛查；最小代价的总改善仍很小，5条的严格整对恢复结果均未改变。这是定向数值诊断的结论，不是总体恢复率或全局最优证明。

| 间隔/ns | ΔJ上限 | 范围包含真Q2 | 该比例的精确二项95%区间 | 跨度/真Q2中位 | 分离区段帧数 |
| --- | --- | --- | --- | --- | --- |
| 40 | 1 | 9/16 | 29.9%–80.2% | 0.780 | 2/16 |
| 40 | 4 | 14/16 | 61.7%–98.4% | 1.605 | 5/16 |
| 40 | 9 | 16/16 | 79.4%–100.0% | 4.049 | 6/16 |
| 200 | 1 | 10/16 | 35.4%–84.8% | 0.345 | 0/16 |
| 200 | 4 | 16/16 | 79.4%–100.0% | 0.690 | 0/16 |
| 200 | 9 | 16/16 | 79.4%–100.0% | 1.035 | 0/16 |

**范围宽度也要看**：ΔJ≤4时，40／200ns两组的总跨度／真Q2中位分别为1.605和0.690；40ns有5/16条出现分离区段，3/16条碰到50pe下界。宽范围可以包含真值，同时仍无法精确分配两个事件的电荷。

**事件数选择会改变分母**：40ns原盲基线只在6/16条输出两个事件，其中ΔJ≤4区段包含真值为4/6；全体已知K2样本的14/16不能直接当作最终盲输出的覆盖率。200ns原盲基线严格整对恢复12/16，条件剖面最优点为10/16；这批小样本未显示逐事件恢复的普遍改善，也未检验差异显著性。

**数值一致性**：32条×3层级的区段与真实Q2探针包含判定有0个不一致；0个探针低于主参考代价超过0.01，0个探针未收敛。插值ΔJ与实际探针的最大差仍为0.508，不能把成员判定一致扩大为整条连续剖面已全局精确。

The numerical protocol targets the two previously changing level-9 Gaussian boundaries plus three prespecified controls. Actual algorithms, bounds, tolerances, budgets and selection are frozen in protocol.json; source and model hashes are in freeze.json. Initial serialization failure, original source in inference-source-v1, original freeze-v1.json and the empty interrupted output are retained. The fix only changes result serialization, not the physical model or numerical settings.

After all five diagnostics pass, coverage-protocol.json and coverage-freeze.json freeze a new, complete algorithm and independent-noise design before generation. Sixteen repeated pairs at each40/200ns gap have fixed Q=[18000,800]pe and t1=200.3ns. Twelve empty/single controls use the original blind baseline/gate. All32 pairs, including numerical failures, stay in the primary denominator. The main inclusion metric is membership in the computed, interpolated interval set; downstream fixed-trueQ2 probes are separate diagnostics and never alter the profiles. Exact binomial95% intervals quantify Monte Carlo uncertainty of inclusion proportions, not nominal charge intervals. Small per-condition samples cannot certify95% coverage. No threshold tuning, new optical parents, ROOTs, training or default change.

Run from the project root with the project Python environment, OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1. Up to8 workers; refinement5, pilot2, main7. Commands use frozen source/model snapshots rather than live code, with hashes checked; no Git checkout needed. Arrays, weights, Notebook and run records are local ignored inputs and require a separate reproducibility bundle.

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_contour_coverage_20261027/run.py tests
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_contour_coverage_20261027/coverage.py tests
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_contour_coverage_20261027/report.py analyze
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_contour_coverage_20261027/report.py figures
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_contour_coverage_20261027/deliver.py validate
```

The run/generate/pilot/infer/probe commands refuse existing outputs. Recompute in a separate project tree with all frozen inputs. coverage.py freeze requires the completed five-frame numerical gate. Source and seeds are fixed before any new-noise result, and pilot frames are reused exactly in the full sample. Full-waveform examples0/1/16/17 are prespecified. Inference workers cannot open truth/cases/ROOT; scoring-only trueQ2 probes run after primary hashes freeze. Full cost replay, sampler covariance/FIR checks, sixteenth-sample denominator checks, old-file hashes, Notebook preservation and embedded figures are validated. Notebook writes use256KiB chunks and full readback on shared storage.

Scope: pure alpha, fixed95/351ns, full1024x2ns analog voltage, ample dynamic range, knownK2 ordered later charge, observed-trigger-local time domain and50-50000pe charge box. This is frozen-model repeated sampling, not independent nativeSiPM validation, deposited-energy resolution, a universal minimum spacing or high-rateK>=3 recovery. Numerical stability, interval width/inclusion and individual-event accuracy are reported separately. Existing Baker/Cousins1984 and Liddick2012 notes are reused with detector/statistical conditions and AI-pending-review status; see docs/literature-use.md L65. Expand parameter points and independent samples after resolving numeric diagnostics; use separate calibration and confirmation samples for future threshold selection.

Local outputs: jupyter/artifacts/alpha_pair_contour_coverage_20261027/; records: jupyter/run_records/alpha_pair_contour_coverage_20261027/. Notebook0004 sections108–113 retain all518 preceding cells and add executed diagnostics and numbered figures. No Git commit or push in this study.
