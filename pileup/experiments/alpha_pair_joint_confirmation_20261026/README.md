# 联合条件似然的冻结后独立确认

开发实验 `alpha_pair_joint_likelihood_20261026` 按预定条件选中 `joint_conserved_nll`。本实验保持其均值、候选依赖协方差、面积约束、搜索和增峰门槛，另取新光学母事件确认；不在新数据上重新选择参数。

## 固定设计

320个新原生母ROOT：1.625、2.5、5.486、8 MeV源能量各80个。源种子从276000000起，联合数字化基种子276500000；源、雪崩和读出随机流在生产前核对旧manifest及案例记录。实际物理源码及G4二进制由旧实验复制并冻结，没有调整探测器或发光参数。

每ROOT一个20 μs单α原生窗。保留光子共同推进SiPM微胞，构造十种电荷组合×8/40/200 ns×40对，加80空窗与320原生单事件，共1600帧。两事件不是独立电压直接相加，首对优化数字化结果须与原标量实现逐数组一致。分析输入为充分动态范围的2 ns×1024点模拟电压，不按ADC削顶拒绝窗口。

固定95/351 ns双τ，未知局部K≤2，Q为实现雪崩等效pe，t为有效SiPM面首光子。整对成功需两真实Q≥500 pe、各Q误差≤20%、各t误差≤min(10 ns,Δ/4)、无额外输出。另报常规10 ns匹配、漏检、误拆／误并与条件Q/t误差；Q精度不是沉积能量分辨率。

比较冻结的新似然、原观测GLS、原模型GLS及同280门槛母事件的原模型GLS。前两类组合平均各15点，沿用开发前规定的恢复与控制要求。逐点Wilson95区间外，用40个母事件块的分层配对bootstrap（2000次，种子276990001）估计平均差区间；前20／后20块分开重抽，跨能量、条件、间隔和方法共享权重，不把复用母事件当独立帧。

## 运行与保留

`freeze.json` 保留实际物理源码、二进制、分析源码、模型、门槛和协议。各ROOT的实际宏、种子、哈希及状态写入manifest、ROOT RunInfo和本地Run_Catalog。

初次4进程生产收到SIGTERM：2个完整核验ROOT保留，4次未完成尝试的宏、日志及可能产生的ROOT全部保留；原因未确定。随后降为单进程，按原预定种子重跑中断项，未根据物理结果筛选。`interruption-1.json` 记录对应关系。

从项目根目录执行，已有结果支持断点继续：

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_joint_confirmation_20261026/run.py native --workers 1
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_joint_confirmation_20261026/run.py physical --workers 4
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_joint_confirmation_20261026/run.py infer --method joint --workers 4
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_joint_confirmation_20261026/run.py infer --method observed --workers 4
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_joint_confirmation_20261026/run.py infer --method model --workers 4
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_joint_confirmation_20261026/run.py report
```

已完成ROOT不覆盖。中断后若存在尚未登记的实际宏／日志，应先核验并登记中断或恢复状态，再安排相同种子的明确重试；不删除历史结果。

完整结果与验证分别保存于本地 `jupyter/artifacts/alpha_pair_joint_confirmation_20261026/` 和同名 `jupyter/run_records/`。Notebook0004第63–65节使用固定图0004-63-1、63-2、64-1。Notebook及jupyter内容只在本地维护。

评分复核与图表脚本另存于 `report-source/`，由 `report-freeze.json` 记录哈希及所依赖的推断快照，不修改已经冻结的物理／推断源码。结果重放命令：

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_joint_confirmation_20261026/report.py validate
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_joint_confirmation_20261026/report.py figures
```

同一入口的 `notebook` 和 `catalog` 分别追加分析笔记和运行目录；本轮均已执行，检测到相同交付标识时会停止，避免重复追加。

本实验为同一固定模型的独立光学样本，不是实物验证、连续能量／位置覆盖、free-τ最短窗或随机高率K≥3性能结论。母事件复用与有限控制数量限制均保留。

## 完成结果

320个新ROOT和1600帧均完成，1200对的两真实Q均≥500 pe。全部方法保持冻结；三套候选搜索共4800次盲帧分析，四种选择规则共6400条评分。

| 方法 | 弱组合15点平均 / % | 较强组合15点平均 / % | 空／单控制通过 |
| --- | ---: | ---: | --- |
| 原观测GLS | 60.50 | 92.83 | 否：weak Q偏差−9.04% |
| 原模型GLS，原560例门槛 | 58.83 | 92.67 | 是 |
| 原模型GLS，共同280例门槛 | 58.50 | 92.67 | 是 |
| 冻结的joint_conserved_nll | 64.50 | 93.33 | 是 |

相对原观测GLS，弱组合平均增加4.00个百分点，母事件分层配对bootstrap95区间为[1.00,7.17]个百分点；较强组合增加0.50个百分点，区间[−1.17,2.33]，不据此宣称较强组合有明确提高。完整方案空／weak／low／mid／high额外输出为0/0/0/1/0（每组80例），单事件Q偏差中位依次+0.076%／−0.007%／−0.021%／−0.043%。

收益主要来自先弱后强：weak_high、8 ns由原观测／模型GLS的0/40提高到27/40。但high_weak的8和40 ns仍0/40；200 ns新20/40，低于原观测23/40和原模型26/40。仅在评分端指定真实K后，新high_weak在8/40/200 ns为0/3/20（各40对）。通过整体预定规则不意味着所有条件改善，不按真实强弱顺序逐帧挑选方法。

320个ROOT／宏／种子、物理与推断快照、门槛、母事件独立性、6400条评分及140条汇总全部核验。新方法41/1600最终候选未获优化器成功状态，原观测7、原模型2；保留并计入结果，未选择性重跑。三个新图和7个代码单元已执行；Notebook0004旧327个单元保留，Run_Catalog登记320完成ROOT及4次中断尝试。
