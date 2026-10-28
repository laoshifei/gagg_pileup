# Numerical stability of conditional later-charge profiles

32条预定波形中，32条在重点ΔJ≤4层级通过末次加密的数值筛查；30条在ΔJ≤1／4／9三层级同时通过。增加搜索和加密网格没有改善这批样本的已知K2严格恢复数，原生40ns组3/8条的分离区段仍存在。0条真实Q2评分探针把最终参考代价降低超过0.01；该探针只用于评分，不反馈主剖面。仍有数值变化或搜索遗漏，不能把当前范围当作已稳定的置信区间。本轮完成数值审核，独立覆盖率校准留待处理这些变化后开展。

| 条件 | 最后加密三层级均稳定 | ΔJ≤4宽度中位：旧→新／真实Q2 | 多区段：旧→新 | 已知K2严格正确：旧→新 |
| --- | ---: | ---: | ---: | ---: |
| native／40 | 8/8 | 12.809→14.237 | 3→3/8 | 1→1/8 |
| native／200 | 8/8 | 0.928→0.950 | 0→0/8 | 1→1/8 |
| continuous／middle | 6/8 | 1.423→1.394 | 0→0/8 | 1→1/8 |
| continuous／long | 8/8 | 0.736→0.764 | 0→0/8 | 6→6/8 |

末次加密在ΔJ≤1／4分别为32/32和32/32通过；未通过的是continuous第160条（ΔJ≤9，端点变化6.61%）、continuous第190条（ΔJ≤9，端点变化5.40%）。这两条的区段数、边界标志及参考代价不变，只是端点变化超过预定5%。

**必须按同一批样本比较**：原生40ns组本轮8条的范围跨度／真实Q2中位数为12.81→14.24；上一轮全部40条的旧中位数为2.69。前两者才是加密前后比较，不能说从2.69恶化到14.24。系统抽取8条不保证统计代表性；跨度也包含分离区段之间的空隙。四组已知K2的严格恢复数仍分别为1／1／1／6（各8条），增加搜索没有改善这批样本的逐事件准确度。

本轮额外27180次局部优化中，1248次保留未收敛标志；这不是失败波形条数。新真实Q2评分探针有1条未收敛，取旧／新较低代价后实际选用的探针有0条未收敛。有限搜索仍可能遗漏更低解。

另做了**评分后的一致性检查**：图中区段来自网格的lnQ插值，真Q包含计数来自实际拟合探针，两者不能定义混用。本轮32条×3层级共96个包含判定有0个不一致；插值ΔJ与探针ΔJ的最大绝对差为0.369。包含判定相同不等于整条连续曲线已收敛。这项检查不调整主曲线，也不是新的独立覆盖率检验。

## Fixed experiment

Pure alpha, fixed95/351ns,1024x2ns analog waveform, ample dynamic range. Known K=2, chronological later Q2, equivalent pe rather than calibrated deposited energy. All32 frames are previously exposed development data: systematic positions0,5,...,35 in each of the four40-frame native40/200ns and conditional-Gaussian24-64/128-256ns groups. The four first cases are pilot and plotting examples, retained in the final32. No new ROOT, training, seed or random draw. Inherited model/data/seed provenance is recorded locally.

Four stages: unchanged original finite profile; same grid with6 separated observed starts and150 iterations per direction; every geometric midpoint with4 starts/150 iterations; further midpoints at1/4/9 crossings and adjacent to the minimum,4 starts/150 iterations. Retain better old nodes and a fixed-nuisance slice from each stage minimum. Thus grid stages can also lower old-node costs through the slice. This is not an equal-compute algorithm comparison. The fixed response, parameter-dependent covariance and logdet remain unchanged, as do the charge50-50000pe and observed-trigger-local time domains. Details and numeric seed-selection rules are in protocol.json and the frozen source.

The last-step numerical screen requires reference improvement<=.05, unchanged interval count and boundary flags, and every matched interval endpoint within a factor1.05, separately at levels1/4/9. These are prespecified engineering tolerances, not certified error bounds or confidence levels. Finite multistart, incomplete convergence and interpolated contour endpoints remain limitations. No generalization to unexamined128 frames or physical independent tests is claimed.

After primary hashes are frozen, a separate true-Q2-only scoring probe refits nuisance from observed starts. Keep the lower old/new probe cost for scoring, never as a primary input. A cost below the primary reference by>.01 flags a missed point. Full known-K2 recovery, Q bias/sigma68, timing, missing/extra and split/merge are saved separately from interval inclusion. Original blind count candidates/gates and all160 earlier profiles remain untouched.

## Reproduction

Run from the project root using the project Python environment with OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1. The limit is8 workers; pilot uses4 and main7. Frozen source and model hashes are in freeze.json, report-freeze.json and delivery-freeze.json. run.py/report.py/deliver.py load these snapshots rather than current live code; no Git checkout is needed. Models/arrays/local records are ignored by Git and must be supplied alongside a reproducibility bundle.

```bash
bash jupyter/run_gagg_python.sh pileup/experiments/alpha_pair_profile_stability_20261027/run.py tests
bash jupyter/run_gagg_python.sh pileup/experiments/alpha_pair_profile_stability_20261027/report.py analyze
bash jupyter/run_gagg_python.sh pileup/experiments/alpha_pair_profile_stability_20261027/report.py figures
bash jupyter/run_gagg_python.sh pileup/experiments/alpha_pair_profile_stability_20261027/deliver.py validate
```

pilot/infer/probe are one-time commands refusing existing output. To recompute the complete run, use a separate project tree with the frozen inputs. Preparation interruption from an unscoped Git diff on shared storage was recorded, then resumed with identical protocol and verified backups. Scoped Git stat has a60s timeout; exact source and prior-document backups preserve provenance. Notebook writes use256KiB chunks, full readback and concurrent-edit checks.

Local outputs: jupyter/artifacts/alpha_pair_profile_stability_20261027/; records: jupyter/run_records/alpha_pair_profile_stability_20261027/. Notebook0004 sections103-107, five fixed figures including full0-2048ns comparisons. Original497 cells and outputs preserved; new cells executed. Validation checks6 analytic/numerical tests,32 four-stage profiles,32 true-Q2 probes,128 conditional-fit scores,384 contour scores,288 paired comparisons, actual cost replay, source/data hashes, truth-open guards, figures and Notebook preservation. No new default, Git commit or push.

## Evidence and next step

Only existing local scientific notes were reused: Baker/Cousins1984,NIM221437-442(PDF1/4/5; DOI unverified) distinguishes estimation, coverage and goodness of fit for Poisson/multinomial histograms, not this correlated voltage. Liddick2012,NIMA66970-78,DOI10.1016/j.nima.2011.12.034(PDF3-8;65umDSSD,40MSPS/12bit,25us,artificial doubles) motivates separate detection/accuracy/background reporting; its numerical limits are not transferred. Both notes remain AI drafts pending human scientific review. See docs/literature-use.md L64 for precise adoption and unresolved issues.

Address remaining grid/search/interpolation changes before independent repeated-sampling coverage calibration; distinguish model-internal calibration and native transfer. No free-tau,alpha/gamma,high-rateK>=3 or universal timing-limit result is claimed.
