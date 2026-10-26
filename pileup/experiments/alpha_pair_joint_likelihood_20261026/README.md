# 联合分量诊断与归一化拟合似然

本实验继续检验强α尾部中弱α事件的恢复，包含条件均值、实际分量互相关、候选依赖协方差和 `log det C`，并修正经验模板的完整面积不守恒问题。仍为固定95/351 ns双τ、充分动态范围、2 ns采样、完整2.048 μs局部模拟电压、未知K≤2。Q为实际雪崩等效pe，t为有效SiPM面首光子时刻；没有新ROOT或网络训练。

## 数据与标定

- 汇总既有三批6120帧：5400共同微胞双事件、560原生单事件、160空窗。按母事件 `replicate//2%2` 分折，训练／留出各280母事件，跨间隔和条件复用母事件。
- 留出折的280原生单事件只用于增峰门槛，另320个此前公开母事件的1600帧用于开发评价；三组母ROOT不交叉。旧测试已暴露，不能再称独立测试。
- 原生 `signal_mV` 包含暗计数。初次汇总 `calibration.truth.npz` 保留但未用于标定，正式输入为按雪崩归属重放的 `calibration-ancestry.truth.npz`。560例粒子＋暗分量闭合误差最大1.48×10⁻¹² mV，实际历史metadata配置逐例一致。
- 双事件真实分量来自联合光子／共同微胞数字化，不是把独立电压简单相加。真分量只用于标定、明确的诊断和事后评价。

## 两种新盲流程

电子学极点为2/20/2 ns，固定AR3变换后保留1021维边缘数据。该变换不改变闪烁双τ。噪声近似包含候选强度、Q及邻居依赖的方差、相邻点相关、分量交叉项、电子噪声带和暗噪声带。

| 流程 | 均值 | 总电荷条件化 |
| --- | --- | --- |
| `old_mean_nll` | 原条件均值 | 不使用 |
| `joint_conserved_nll` | 训练折重新标定＋每分量完整面积归一化 | 固定份额0.95的秩≤2修正 |

后者每分量完整面积为Q×13.618501389 mV·ns/pe；不把截断观测积分当作完整Q。高斯近似目标为 `rᵀ C⁻¹ r + log det C`，可为负。增峰用原始NLL差，按拟合Q的1500/5000/12500 pe边界和独立单事件经验秩p≤0.025选择。首个输出需Q≥500 pe且NLL低于空模型。

两流程都改变代价、候选排序及门槛，完整方案还改变均值和积分约束，因此不是单一因素的因果对照。保留两种旧基线，并将旧模型GLS按相同280门槛母事件重新标定。推断worker禁止读取ROOT、truth和case标签。

## 运行和快照

`development-source/` 与 `development-freeze.json` 为初次冻结。小样本在K=0的空分量轴上触发SciPy FIR错误，尚未写入任何预测；该版本完整保留。`development-source-v2/`、`development-freeze-v2.json` 和 `run-v2.py` 修正空轴滤波，物理、标定参数和搜索预算不变，全部正式预测使用v2。

从项目根目录执行；输入数据和本地Python环境需存在：

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_joint_likelihood_20261026/run-v2.py infer --variant joint_conserved_nll
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_joint_likelihood_20261026/run-v2.py infer --variant old_mean_nll
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_joint_likelihood_20261026/run-v2.py select
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_joint_likelihood_20261026/report.py figures
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_joint_likelihood_20261026/report.py validate
```

推断支持按逐帧JSONL断点继续，已有预测不会覆盖。初次数据与均值标定由 `alpha_pair_joint_data.py`、`alpha_pair_joint_calibration.py` 和 `alpha_pair_joint_study.py` 完成，输入哈希、母ROOT、种子、原始宏来源和依赖版本保存在本地同名run_records。重建数据应使用新的输出位置并保留原记录。

## 评价与边界

整对成功：两真实Q≥500 pe、各Q误差≤20%、各t误差≤min(10 ns,Δ/4)、没有额外输出。另报常规10 ns匹配、漏检、误拆／误并和条件Q/t误差；时间失配输出也可能计入额外，不能全部视为虚构事件。

预定选用要求：各空／单控制额外≤5%，各单严格成功≥95%、Q偏差中位绝对值≤5%；弱组合15点平均超过两种旧基线，较强组合15点平均相对原观测权重下降不超过2个百分点。跨点重用母事件，均值是描述性指标，不是独立曝光的合并估计。

文献复用Fixsen2004（DOI 10.1016/j.nima.2003.11.313）的条件均值／非平稳协方差思想，原条件为TES模拟、每能量5000次，共享库AI笔记待人工审核。原文加权最小二乘不等同本项目的含logdet高斯近似；滤波、守恒、低秩修正及GAGG性能是项目推导和计算。未修改共享库。

结果位置：本地 `jupyter/artifacts/alpha_pair_joint_likelihood_20261026/`；运行与验证位置：`jupyter/run_records/alpha_pair_joint_likelihood_20261026/`。Notebook0004新增第58–62节和固定图0004-58-1、59-1、60-1、60-2、61-1。Notebook和全部jupyter内容只在本地维护。

## 实际开发结果

| 方法 | 弱组合15点平均 / % | 较强组合15点平均 / % | 控制通过 |
| --- | ---: | ---: | --- |
| 原观测GLS | 61.00 | 93.17 | 否，保留旧结果 |
| 原模型GLS，原560例门槛 | 58.83 | 93.33 | 是 |
| 原模型GLS，本轮280例门槛 | 58.33 | 93.33 | 是 |
| old_mean_nll | 65.00 | 92.33 | 是 |
| joint_conserved_nll | 65.17 | 94.17 | 是 |

按预定规则选中 `joint_conserved_nll`。空／weak／low／mid／high控制额外输出分别0/0/0/2/0（每组80例），四档单事件Q偏差中位为−0.069%／+0.024%／−0.044%／−0.051%。high_weak的8/40/200 ns盲成功为0/0/22（每点40对）；已知K直接选现有双事件候选后为0/3/22，40 ns的失败仍主要不能由放宽增峰门槛消除。

4项数值测试通过，4800条开发评分和105条汇总重放一致；两版冻结源码和原输入哈希保持。完整方案1600帧中有40个最终候选未获优化器成功状态，保留并计入全部结果，没有选择性重跑；局部停止／有限预算不能证明全局最优。五幅图和10个新增代码单元已执行，原306个Notebook单元保留。独立确认单独位于 `alpha_pair_joint_confirmation_20261026`，不把开发指标混作新样本性能。
