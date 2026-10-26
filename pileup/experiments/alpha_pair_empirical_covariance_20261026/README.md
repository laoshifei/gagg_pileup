# Empirical pulse covariance for blind weak-alpha reconstruction

以孤立单脉冲残差估计非平稳相关噪声，比较四种预先登记的权重。均值条件响应、95/351 ns双衰减时间、物理配置、候选搜索和优化预算不变；无网络训练。输入完整模拟电压、充分动态范围，局部未知K为0/1/2。Q为实际雪崩等效pe，不直接当作能量。

## 数据划分和方法

560个已有原生单事件按replicate奇偶分为280个协方差训练母样本与280个增峰门槛标定母样本。训练的weak/low/mid/high数量40/80/80/80，标定同数；两组与开发验证的320个母ROOT均无种子交叉。验证复用上一轮已公开资料中的十组合×8/40/200 ns×40对及400空／单控制，共1600帧，**不是新的独立测试**。

训练读取原20 μs记录，真实Q/t对齐后的−64到2112 ns、2 ns网格完全被观测覆盖。每档残差中心化后做SVD，保存经验方差和前32个相关模式；均值残差没有加回信号模型。插值带来的相关性及每档仅40或80个母样本的统计限制保留。

协方差为`C=diag(max(v-rowSum(U²),0.01 mV²))+UUᵀ`。四种预定候选：仅对角方差、16模式×0.5份额、16模式×0.9份额、32模式×0.9份额，模式幅度乘份额平方根。正定白化已与稠密逆矩阵二次型作数值对照。

盲推断只用电压的单事件OLS初估选择Q档和平移C，每帧所有K/Q/t候选共享同一C。没有用真实事件数、源能量、分量、Q/t或微胞状态作推断输入。worker审计钩子禁止读取ROOT、truth和case文件。

每个候选使用独立的280个单事件重新标定`(C1-C2)/C1`上尾秩，固定`p<=0.025`，Q分档边界1500/5000/12500 pe。另用相同280例重标定旧AR模型权重，区分噪声权重与标定样本量的影响。原先560例标定的旧流程同时保留。

选择规则在运行前写入`development-plan.json`：各空／单控制额外输出≤5%、各单成功≥95%、各单Q误差中位绝对值≤5%；弱组合15点等权平均超过两种旧盲基线，较强组合15点相对原观测权重最多下降2个百分点。不同间隔及组合复用母样本，平均仅是描述性量，不累计独立曝光。逐对成功要求两个真实Q≥500 pe，各Q误差≤20%、各t误差≤min(10 ns,Δ/4)，且无额外输出；t指有效SiPM面首光子时刻。

文献仅借用Fixsen等2004的条件均值／协方差思想（DOI 10.1016/j.nima.2003.11.313；TES模型、每能量5000次，既有AI精读笔记待审）。该文没有这里的低秩截断、正则化或GAGG双事件性能，这些均为项目实现。

## 重放和记录

在项目根目录、已有本地数组和项目Python环境下，分别执行：

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_empirical_covariance_20261026/run.py infer --variant diagonal
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_empirical_covariance_20261026/run.py infer --variant rank16_half
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_empirical_covariance_20261026/run.py infer --variant rank16_full
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_empirical_covariance_20261026/run.py infer --variant rank32_full
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_empirical_covariance_20261026/run.py select
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_empirical_covariance_20261026/report.py all
```

每种1880帧，4个worker，逐帧检查点；已完成帧保留。`development-freeze.json`保留229个实际分析／物理／配置文件哈希，`development-source/`和`source/`可重放而不需切换Git提交。ROOT、数组、二进制和本地Notebook不随GitHub分发；实际文件哈希和输入来源在本地`jupyter/run_records/alpha_pair_empirical_covariance_20261026/`。

`report-source/`及`report-freeze.json`另外保存185个报告脚本／依赖文件，未改变推断快照。`report.py all`核验冻结输入、重算评分并生成六幅图；若本地旧Notebook还未添加这一阶段，可单独执行`report.py notebook`追加并实际执行新单元，旧内容及图号保留。已交付过的Notebook不会重复追加同一章节。

报告交付后执行`annotate_results.py --notebook`可重放补充的局部收敛核验，并在第54／56节插入收敛次数和局部收益／整体取舍的阅读说明；保留已有全部计算和图像，不改变冻结拟合结果。

本地Notebook0004第55–57节记录训练划分、全部候选及控制，固定图号0004-55-1、56-1、56-2。本轮没有新ROOT，Run_Catalog保留原有登记。

## 开发结果

四种新权重共完成7520次盲帧分析。下表只评价1600个开发验证帧，各类平均是五种组合×三个间隔的15点等权描述性平均；与旧七间隔35点平均覆盖不同。

| 方法 | 弱组合平均 / % | 较强组合平均 / % | 空／单控制达标 | high_weak 200 ns |
| --- | ---: | ---: | --- | ---: |
| 原观测权重＋固定惩罚 | 61.00 | 93.17 | 否：弱Q偏差−9.74%、low额外输出6/80 | 21/40 |
| 原模型权重＋经验筛选，560例门槛标定 | 58.83 | 93.33 | 是 | 19/40 |
| 原模型权重，改为共同280例门槛标定 | 58.83 | 93.33 | 是 | 19/40 |
| 经验对角方差 | 59.67 | 89.50 | 否：mid额外输出6/80，成功74/80 | 20/40 |
| 16模式，份额0.5 | 57.83 | 86.83 | 是 | 9/40 |
| 16模式，份额0.9 | 56.17 | 87.33 | 是 | 8/40 |
| 32模式，份额0.9 | 56.67 | 87.50 | 否：high额外输出6/80，成功74/80 | 8/40 |

四种新方法的high_weak、8及40 ns均0/40；全部没有达到预定选用条件，`development-selection.json`的chosen为null。对角方差相对旧模型权重的弱组合平均略增，但仍低于原观测权重，且较强组合与控制退化。16模式两方案通过单事件控制也未改善双事件恢复。

局部收益保留：weak_high、8 ns的旧AR共同门槛标定对照为0/40，四个新候选依次25/40、24/40、24/40、23/40；这一方向明显改善，同时存在其他方向的损失，不能只据此升级整体流程。

本轮不升级默认流程，未启动新的独立测试。完整175条开发汇总、8000条逐帧开发评分和120组同门槛样本量的配对变化保存在本地结果目录。逐帧评价同时保留通常10 ns匹配下的漏检、额外输出、误拆／误并及条件Q/t误差；额外输出包括未在容限内匹配的输出，不必然是凭空产生的物理事件。

结论限于当前四种孤立脉冲协方差近似及有限校准样本。下一步检验双事件的条件均值、联合方差和互相关；若C随候选Q/t变化，高斯近似目标须包含log det C。新增盲方法通过开发约束后再冻结并生成新母样本，不将本轮已公开数据重新称作独立测试。
