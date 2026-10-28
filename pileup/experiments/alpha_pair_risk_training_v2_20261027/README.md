# 连续电荷监督学习与物理拟合

2026-10-27执行，承接同一似然下电荷分配搜索未明显改善的结果。比较三种电荷损失的小型MLP，并实际运行MLP＋原NLL和MLP＋学习型Q惩罚NLL。通过预定开发要求的方案：无；不自动替换默认方法。

## 数据、物理量与冻结顺序

纯α，固定95/351 ns双τ，2 ns×1024完整模拟电压，充分动态范围，局部未知K≤2。Q使用pe，原生数据为实现雪崩等效电荷；新高斯Q是均值／协方差的条件参数，不作新实现雪崩电荷或沉积能量。t为有效面首光子。本轮无新ROOT或新光学母事件，未改变探测器参数。

9600训练＋960验证，事件数概率0.1/0.3/0.6。单Q为500–25000 pe对数均匀，双事件混合独立对数均匀、强弱随机先后和近等强，个体Q同域、总Q≤48000。间隔4–320 ns的均匀／对数均匀各一半。它是声明的模拟设计分布，非实测源能谱或旧四档测试答案。使用冻结的ConservingResponse与JointNoise条件高斯抽样；FIR前三点的残差前缀设零，原物理NLL评价完整1021个创新。

训练／验证／评价参数种子基数依次270310001/270310002/270310003，每帧seed=base×100000+index；噪声使用SeedSequence([seed,993])独立流。网络初始化270310010，小批次270310011；每种损失只有一个训练初始化。每500步按公共验证指标选择：强制K2严格整对正确最高，其次Q绝对ln误差中位最低，再次count BCE最低。各模型及混合所选损失/宽度先冻结，随后生成新评价。

新连续评价1000帧：360对为近等强／先强后弱／先弱后强×6–12/24–64/128–256 ns×40对，弱600–1500、强12000–24000 pe；另240对两Q独立均匀500–23000、间隔均匀4–320 ns，检验训练分布改变；80空＋四档连续单Q各80。它是模型内新参数／噪声样本，不是未见物理确认。

原生760帧沿用已公开确认批（三组合×8/40/200 ns×40对＋400控制），另280个原生单事件校准门槛。母事件分组与旧模型标定／门槛／评价保持不交叉，条件间复用。旧评价本轮只作公开迁移诊断。所有方法同一评价帧，基线逐项重放核对。

## 实现与先导修订

观测单事件初拟合生成1408维特征：完整1021个白化FIR残差＋385个观测初拟合对齐残差＋lnQ与时刻；asinh(z/3)、训练集均值／标准差归一化后截到±8。这是网络数值预处理，不是ADC裁剪。MLP隐藏层128/96、SiLU，输出两个时间槽（−384至384 ns/2 ns网格、峰附近五点平均）、总Q和分配份额，以及K2分数。最终Q候选clip到50–50000 pe，盲输出仍要求每Q≥500。训练目标为模拟事件数/Q/t参数，没有将测试真值或实际分量波形送入推断。

三损失为固定尺度绝对平方((Qhat−Q)/5000)²、对数平方ln(Qhat/Q)²、相对平方(Qhat/Q−1)²。总损失count BCE＋3×双事件平均Q损失＋0.4×双事件时间CE，时间软标签σ=2 ns；3000步、batch64、Adam0.0007、相同初始化及小批次流。NumPy实现，未新增ML依赖。

初版`alpha_pair_risk_training_20261027`完整保留。24帧先导及训练输入发现单事件初拟合总Q偏差：5757个训练双事件有740个所需|ln(Qsum/Qpilot)|>0.1，范围−0.758758至0.425039。因此正式训练及评价前另建v2，把总量头从Qpilot×exp(0.1tanh(z))改成Qpilot×exp(tanh(z))。生成函数AST相同、数组逐字节哈希相同，v2复用初版训练/验证数据，并重做梯度、过拟合先导和物理混合检查；没有重抽样选容易样本。`preflight-head-amendment.json`记录依据与顺序。

验证选择：absolute第3000步，验证强制K2成功30.33%，总训练49.5s；logarithmic第2500步，验证强制K2成功38.00%，总训练48.4s；relative第3000步，验证强制K2成功41.00%，总训练53.6s。两个混合流程都用`relative`，lnQ惩罚宽度为[0.20801263325687447, 0.21970885078116778]。网络和原最佳K2两起点，各最多70次L-BFGS-B、时刻起点±24 ns，保留最佳初值；regularized加入Σ[(lnQ−lnQnetwork)/σ]²，σ=max(0.2,验证全部双事件的68%绝对lnQ误差)，包含错时刻。网络与拟合看同一波形，所以是数据依赖正则化，不是独立先验或新增观测证据。

网络分数不当似然、不再乘事件率先验。各输出使用同一280单事件、拟合单Q边界1500/5000/12500 pe、加一并保留并列的上尾秩≤0.025，各自校准增峰；保留原K0/K1及每Q≥500。整体比较同时改变候选和检测分数，额外报告强制K2诊断，不能把所有差异归因于拟合。推断禁止打开truth/case/ROOT；失败与未收敛输出均计入。

## 实际结果与控制

严格整对正确：两个真实Q≥500，每个Q误差≤20%，每个时刻误差≤min(10 ns,gap/4)，无额外输出。合计数的条件配比是人为设计，不能解释为某注入率的效率。

| 方法 | 新高斯600对 | 其中分布改变240对 | 已公开原生360对 |
| --- | ---: | ---: | ---: |
| baseline | 500/600 (83.3%) | 236/240 (98.3%) | 241/360 (66.9%) |
| nn_absolute | 232/600 (38.7%) | 133/240 (55.4%) | 80/360 (22.2%) |
| nn_logarithmic | 225/600 (37.5%) | 127/240 (52.9%) | 71/360 (19.7%) |
| nn_relative | 257/600 (42.8%) | 135/240 (56.2%) | 71/360 (19.7%) |
| mlp_nll | 482/600 (80.3%) | 232/240 (96.7%) | 226/360 (62.8%) |
| mlp_regularized | 465/600 (77.5%) | 232/240 (96.7%) | 198/360 (55.0%) |

| 方法 | 原生先强后弱8 ns | 40 ns | 200 ns |
| --- | ---: | ---: | ---: |
| baseline | 0/40 | 0/40 | 20/40 |
| nn_absolute | 0/40 | 1/40 | 2/40 |
| nn_logarithmic | 0/40 | 1/40 | 6/40 |
| nn_relative | 0/40 | 1/40 | 2/40 |
| mlp_nll | 0/40 | 0/40 | 7/40 |
| mlp_regularized | 0/40 | 0/40 | 5/40 |

各组空／单控制、Q偏差/σ68、时间误差、漏检、额外输出、误拆／误并、K2输出及强制K2候选均记录在CSV。误差宽度是时间匹配条件下的(p84−p16)/2，没有沉积能标，不能当能谱FWHM。额外输出可能为有信号区域的错时刻或误拆。

预定开发要求：新连续双事件总成功率至少+5个百分点；连续及原生各空／单组额外输出≤5%，单组成功≥95%，匹配Q偏差中位绝对值≤5%；原生双事件总数不下降。

| 方法 | 开发要求全部通过 | 未通过项 |
| --- | --- | --- |
| nn_absolute | False | native_extra、native_single_success、continuous_pair_gain、native_pair_no_drop |
| nn_logarithmic | False | continuous_pair_gain、native_pair_no_drop |
| nn_relative | False | native_extra、native_single_success、continuous_pair_gain、native_pair_no_drop |
| mlp_nll | False | continuous_pair_gain、native_pair_no_drop |
| mlp_regularized | False | native_extra、native_single_success、continuous_pair_gain、native_pair_no_drop |

通过预定开发要求的方案：无；不自动替换默认方法。 成对区间2000次bootstrap（种子270310099）：连续按独立帧，原生每条件按母事件块（强弱20块×2对，中／中40块）；没有把原生合并数按独立帧计算区间。区间限当前模型/公开资料，退化区间不等于物理确定性。无训练多种子、连续能量新母样本、位置失配或实物确认，不能推断普适间隔极限、最高率、α/γ、free-τ或高率K≥3。

## 验证、输出与重放

6项数值测试通过：三风险标量最优值、三损失输出梯度、网络梯度、时间解码/总Q、门槛并列与阈下保护、保存/Adam下降。24帧先导后完成三次正式训练；四核上BLAS单线程。1760评价帧×六方法＝10560主评分，960双事件×六方法＝5760强制K2诊断；原生基线逐项一致，K0/K1不变，NLL混合保留不差于原K2的值。网络重放最大差0、NLL最大相对差0。网络耗时表只包括batch forward/解码，排除了初拟合特征和基线；混合耗时也是增量，不声称实时吞吐。

本地Notebook0004第87–92节新增24单元、执行12代码单元；原431单元及输出保持。图0004-87-1、88-1、89-1、90-1、90-2、91-1固定编号。预定native先强后弱40/200 ns第0对展示全0–2048 ns、原生真实分量与条件均值；最终只选K1时第二输出画零。

源码在`inference-source/`及`evaluation-source/`，`freeze.json`、`evaluation-freeze.json`、`training-freeze.json`与`report-freeze.json`分别保护计算、推断、权重与报告。首次生成基线发现冻结包缺少既有`alpha_pair_residual`依赖，尚未写出任何基线预测；另加`dependency-source/`与`dependency-freeze.json`、`evaluate_complete.py`补齐，来自旧确认实验的相同源码，未改算法、权重或测试配置，原快照保留。实际HEAD为`bcb27007fadf5df7dbb1d2447638b1f371069a8c`加冻结代码；之前信息/分配研究未提交改动保留。

从项目根目录运行，保留本地输入、模型及运行记录，无需切换Git：

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_risk_training_v2_20261027/evaluate_f64.py baseline --workers 4
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_risk_training_v2_20261027/evaluate_f64.py predict --domain continuous --workers 4
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_risk_training_v2_20261027/evaluate_f64.py predict --domain native --workers 4
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_risk_training_v2_20261027/report_v3.py analyze
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_risk_training_v2_20261027/report_v3.py figures
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_risk_training_v2_20261027/report_v3.py notebook
```

推断按index续跑。数据生成/训练不覆盖现有目录；`run.py generate --role train|validation|evaluation`、`train --risk absolute|logarithmic|relative`与`freeze-training`供新独立输出路径使用，先做对应prepare/先导并保留冻结顺序。v2数据来源、生成记录、先导修订及全部文件哈希在`jupyter/run_records/alpha_pair_risk_training_v2_20261027/`；数组、候选、评分/CSV/图在`jupyter/artifacts/alpha_pair_risk_training_v2_20261027/`。Notebook、模型和数据按约定仅本地维护，单独克隆代码不包含输入或训练权重。

报告排版v2补齐各子图(a)–(f)标记。随后数值检查发现初始候选的偏置扣除留在float32，而原基线/精修采用float64，部分单脉冲的代价比较差约10⁻⁵，触发不劣化断言；未放宽检查，统一为float64后重算全部1040原生/1000高斯混合推断。旧推断记录保存在artifacts/run_records下的`inference-f32-before/`，4例回归包含实际失败帧，训练权重/数据/物理目标/优化预算均保持。最终推断为`evaluation-f64-source/`、`evaluation-f64-freeze.json`、`evaluate_f64.py`；最终报告为`report-source-v3/`、`report-v3-freeze.json`及`report_v3.py`，同步以float64重放代价。此前源码和报告快照均保留。

文献来源和采用边界见`docs/literature-use.md` L61：Morad2024的分阶段评价和Kim2023的监督/弱信号边界仅作为设计依据；本轮MLP/三风险/惩罚不是其原算法复现。
