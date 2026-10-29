# Independent blind joint-parameter compatibility study

用校准集中盲选为双事件的36条真实双脉冲，经验95%上分位得到联合代价阈值c=11.864；40/200ns分别贡献6/30条。这是当前选门和设计下的混合样本标定，不是统一95%置信保证。

40ns：原盲流程输出双事件14/64，严格整对恢复0/64；选出双事件后的联合区域包含真实Q1/Q2/t1/t2为13/14，按全部真实双脉冲计算则为13/64。时间漏检57个，额外输出7个。

200ns：原盲流程输出双事件63/64，严格整对恢复38/64；选出双事件后的联合区域包含真实Q1/Q2/t1/t2为58/63，按全部真实双脉冲计算则为58/64。时间漏检15个，额外输出14个。

找到达到或超出容差边界的替代解：79条盲选双事件，其中真实双脉冲77条、控制样本误报双事件2条；严格正确为38/79。分母为0时准确率未定义。

有限搜索未找到替代解：0条盲选双事件，其中真实双脉冲0条、控制样本误报双事件0条；严格正确为0/0。分母为0时准确率未定义。

质量标记没有改变原拟合或增加正确恢复事件；找到替代参数说明该代价区域内存在歧义，有限搜索未找到不证明全局唯一或物理可辨识。

| 间隔/ns | 盲选双事件 | 严格整对恢复 | 选双后联合包含 | 全体中选双且包含 | 时间漏检/128 | 额外输出 |
| --- | --- | --- | --- | --- | --- | --- |
| 40 | 14/64 | 0/64 | 13/14 | 13/64 | 57/128 | 7 |
| 200 | 63/64 | 38/64 | 58/63 | 58/64 | 15/128 | 14 |

This is a frozen-model diagnostic, not a new reconstruction algorithm or a theorem of nominal confidence coverage. The original blind candidate search and native count gate alpha0.025 are unchanged. Joint cost regions use the original K2 candidate as reference, not a claimed global optimum. Eight constrained exterior searches look for actual alternative Q/t sets at or beyond the20%/10ns accuracy boundary; absent finite-search witnesses do not prove containment or identifiability.

Calibration144 frames precedes independently generated confirmation288 frames. See protocol.json for the fixed stratified designs, SeedSequence words, complete numerical settings and stopping rules. The selected-pair empirical95% quantile is frozen before confirmation, with domains retained and no test tuning. Per-gap exact binomial sampling intervals condition on this frozen calibration; pooled fixed strata do not receive iid-binomial intervals. Truth is blocked in inference workers and used only for simulation and downstream calibration/scoring. Numerical failures, missed/extra/split/merge and quality-acceptance losses remain explicit.

Purealpha,fixed95/351ns,full1024x2ns float64 voltage,ample dynamic range,localK<=2. No new ROOT,optical mothers,training,free-tau,mixedalpha/gamma or rate-scan claims. Q is equivalentpe,not depositedenergy. No marginal profile interval widths or complete topology are computed. Source/model/protocol snapshots avoid a required Git checkout; local ignored arrays/models/Notebook must be supplied separately.

From the project root, using the project environment and BLAS1:

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_blind_uncertainty_v2_20261028/run.py tests
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_blind_uncertainty_v2_20261028/run.py generate-calibration
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_blind_uncertainty_v2_20261028/run.py pilot --workers 2
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_blind_uncertainty_v2_20261028/run.py infer-calibration --workers 4
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_blind_uncertainty_v2_20261028/run.py calibrate
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_blind_uncertainty_v2_20261028/run.py generate-test
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_blind_uncertainty_v2_20261028/run.py infer-test --workers 4
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_blind_uncertainty_v2_20261028/report.py analyze
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_blind_uncertainty_v2_20261028/report.py figures
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_blind_uncertainty_v2_20261028/deliver.py validate
```

Existing output paths are protected. For clean replay use an isolated project copy and the recorded input assets; the pilot's nine rows are reused unchanged. Primary and report/delivery sources have separate frozen manifests. The sibling preparation without_v2 failed to import an omitted absolute analysis dependency before tests/generation and is preserved as an unsuccessful preparation, not a result. Local Notebook0004 sections127–132 append to600 preserved cells and include eight explicitly numbered figures. No commit or push.
