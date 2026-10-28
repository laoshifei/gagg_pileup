# Event selection and independent matched-Gaussian gate confirmation

40ns：同为0.025门槛参数，改用新单脉冲校准后，输出两个事件由17/64变为33/64，严格整对恢复由0/64变为1/64。即使直接给这批真实双脉冲选K=2，当前候选中严格成功也只有10/64；这是同一有限候选的诊断上限，不是物理可辨识极限。

40ns被选中候选的有序Q2误差中位由+61.1%变为+39.5%，这是选择群体改变，不是同一条拟合参数变准；包含时间未匹配候选，不称能量分辨率。时间未匹配的额外输出由7增为17个，多检出没有转化成相应的准确恢复。额外输出表示10ns内找不到匹配真事件，不必然表示整窗原本无信号。

200ns：同为0.025门槛参数，改用新单脉冲校准后，输出两个事件由63/64变为63/64，严格整对恢复由46/64变为46/64。即使直接给这批真实双脉冲选K=2，当前候选中严格成功也只有46/64；这是同一有限候选的诊断上限，不是物理可辨识极限。

四档单脉冲共256条，原／新主门槛误报双事件分别为1/1条；含任意时间未匹配输出的帧分别为1/1条。各Q档结果另列，不能把这个等权混合比例当作任意输入谱的误报率。

| 间隔/ns | 门槛校准来源 | 输出双事件 | 严格整对恢复 | 时间漏检事件 | 额外输出事件 |
| --- | --- | --- | --- | --- | --- |
| 40 | 原门槛 | 17/64 | 0/64 | 54/128 | 7 |
| 40 | 新条件高斯单脉冲 | 33/64 | 1/64 | 48/128 | 17 |
| 200 | 原门槛 | 63/64 | 46/64 | 8/128 | 7 |
| 200 | 新条件高斯单脉冲 | 63/64 | 46/64 | 8/128 | 7 |

The physical waveform fitter is unchanged. All six decision rules are prespecified: original calibration or fresh matched-Gaussian calibration, each with rank parameter0.025/0.05/0.10. The primary contrast is0.025. Use256 independent single pulses for calibration (64 at each800/3000/8000/18800pe); freeze gates before generating448 new confirmation waves (64 each empty,the four singles,and18000/800pe pairs at40/200ns). First time200.3ns. Actual source/model/design snapshots are recorded in freeze.json,with report/delivery manifests separate. Inference workers cannot open truth/cases/ROOT. Test truth is used only after prediction hashes freeze; no test-based gate tuning or default update.

Previous32 charge profiles are diagnostic data only. Retain all interval components and distinguish P(C),P(S),P(C|S),P(C and S),P(C|not S). No new charge interval profiles or nominal coverage claims. Counterfactual forcedK2 uses the same observed candidates and true-pair evaluation only; its strict-success count bounds gate-only improvement for these candidates,not physical resolvability.

Purealpha,fixed95/351ns,full1024x2ns float64 voltage,ample dynamic range,localK<=2. Exact frozen JointNoise Gaussian in FIR space,conditioned zero first-three residuals. SeedSequence[202610272001,split,index,995],split1calibration/2confirmation; all704 draws independent. No optical mothers,new ROOTs or training. Fitted-Q-bin calibration at four fixed true charges is not a uniform nuisance-parameter false-positive guarantee. Calibration sizes also differ: at parameter0.025,n64 accepts only rank1/65 whereas n80 can accept2/81. Report sample sizes,attainable ranks and actual evidence cutoffs; differences are not solely due to noise-model matching. Exact binomial confirmation bars condition on the frozen calibration and do not include calibration uncertainty. Q is equivalentpe,not deposited energy; controlled gaps are not injection rates.

Run from the project root with the project environment and BLAS1. Frozen scripts use snapshots rather than requiring Git checkout. Pilot2workers,initial calibration6,resumed calibration2,confirmation2; maximum8. Calibration stopped writing after200 rows for about20minutes,cause undetermined. The owned run was interrupted; all200 rows preserved,then only56 missing tasks recomputed using the unchanged fitter with2workers. The separate resume.py/source/manifest and interruption record retain this history. Source/config/run scripts are shareable; ignored local arrays,models and Notebook require a separate reproducibility bundle. Numeric phases refuse existing outputs; replay generation/inference in a separate tree with all frozen inputs. Pilot rows are reused byte-for-byte in calibration.

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_selection_20261027/run.py tests
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_selection_20261027/run.py generate-calibration
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_selection_20261027/run.py pilot --workers 2
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_selection_20261027/run.py infer-calibration --workers 2
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_selection_20261027/run.py gates
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_selection_20261027/run.py generate-test
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_selection_20261027/run.py infer-test --workers 2
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_selection_20261027/report.py analyze
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_selection_20261027/report.py figures
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_selection_20261027/deliver.py validate
```

Local outputs jupyter/artifacts/alpha_pair_selection_20261027/; records jupyter/run_records/alpha_pair_selection_20261027/. Notebook0004 sections114–119 preserve543 prior cells and add executed tables and six numbered figures. Full-waveform examples320/321/384/385 are prespecified. Check relativeQ bias/sigma68,time errors,missed/extra/split/merge,strict pair recovery and per-control false outputs separately. Existing Liddick2012 notes support separating detection,parametric accuracy and false pairs; detector conditions and pending-review status are retained in docs/literature-use.md L66. Future work: improve accepted candidates,then independent selection-aware interval calibration/confirmation and native transfer. No commit or push.
