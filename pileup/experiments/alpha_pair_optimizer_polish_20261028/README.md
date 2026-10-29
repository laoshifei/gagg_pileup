# Count-preserving optimizer ablation

独立288条确认采用同一波形的配对比较，原事件数门和K判决保持不变；主方案为预先固定的多起点后Powell，另外两种仅作过程消融。

40ns：整对严格恢复0/64→0/64；改善0、退化0，配对精确检验p=1。时间漏检50→51，额外输出4→5。

200ns：整对严格恢复36/64→34/64；改善1、退化3，配对精确检验p=0.625。时间漏检18→20，额外输出15→17。

80条被盲选为K2的确认波形中，20条代价降低超过0.05。控制组原误报双事件1/160，K误报数保持不变是设计约束，不能作为新筛选器有效的证据。

较低代价与逐事件准确度须分别解释；有限搜索未找到更优解不证明物理不可辨识。本轮不替换默认算法，也不外推为原生SiPM、α/γ、free-τ或高计数率性能。

| 间隔/ns | 原流程 | 宽域重启 | 多起点 | 再加Powell（主方案） | 改善/退化 |
| --- | --- | --- | --- | --- | --- |
| 40 | 0/64 | 0/64 | 0/64 | 0/64 | 0/0 |
| 200 | 36/64 | 36/64 | 34/64 | 34/64 | 1/3 |

The original waveform-only candidate generator and native count gate remain unchanged. K=0/1 outputs are exact copies; only selected K=2 parameters are polished. The primary configuration was frozen before fresh test generation: observed starts, six L-BFGS-B starts and bounded Powell, retaining earlier lower-cost fits. Wide restart changes the parameterization/domain as well as iteration budget. No post-polish Q>=500 cut, true count, truth parameters or prior exterior witness enters the optimizer. Improved objective is not guaranteed improved event recovery, and finite searches do not prove a global optimum.

Fifteen previously exposed development frames include six cases selected for suboptimal original objective. They are not an unbiased performance sample. Independent confirmation uses 288 fresh model draws and same-waveform paired comparisons. See protocol.json and test-design.truth.json for fixed design and random seeds. Full per-event and denominator-aware scoring is saved locally. No default algorithm promotion, training or new ROOT files.

Pure alpha, fixed95/351ns, full1024x2ns float64 voltage, ample dynamic range; conditional Gaussian model with first three FIR residuals fixed at zero. Q is equivalent photoelectrons, not deposited energy. No native SiPM, mixed alpha/gamma, free-tau or continuous-rate claim. Models and arrays are ignored local assets, required separately for replay. Source snapshots contain actual dependencies; no Git checkout is required.

From project root using the project environment and one BLAS thread:

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_optimizer_polish_20261028/run.py tests
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_optimizer_polish_20261028/run.py development-inputs
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_optimizer_polish_20261028/run.py infer-development --workers 4
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_optimizer_polish_20261028/run.py score-development
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_optimizer_polish_20261028/run.py generate
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_optimizer_polish_20261028/run.py infer-test --workers 4
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_optimizer_polish_20261028/run.py score-test
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_optimizer_polish_20261028/report.py analyze
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_optimizer_polish_20261028/report.py figures
```

Existing outputs are protected. Replay in a separate project copy with original local development assets and models. Primary and reporting source freezes are separate manifests. Notebook delivery also requires the recorded before-version of the local Notebook, README and literature-use document; old 627 cells and outputs are preserved. Full commands/environment/hashes and validations are in jupyter/run_records/alpha_pair_optimizer_polish_20261028/. No commit or push was performed.
