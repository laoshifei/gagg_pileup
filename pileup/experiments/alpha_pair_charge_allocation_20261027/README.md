# 总电荷与强弱分配的盲搜索对照

2026-10-27执行，承接局部信息诊断。本轮改变候选搜索和数值坐标，检验强弱Q反相关方向的有限优化；均值、协方差、logdet及物理模型保持。四种变体均未同时通过预定原生开发约束，不选择新默认流程。

## 数据与预定方法

280个独立原生单事件用于门槛，760帧已公开原生数据与760帧同Q/t高斯draw0用于开发比较，共1800帧；高斯draw1本轮未运行。每种输入包含三种组合×8/40/200 ns×40对，加空／四档单事件各80帧。训练280、门槛280和评价320个母事件互不重叠；跨条件复用母事件。没有新ROOT或网络训练。

保持纯α、固定95/351 ns双τ、2 ns×1024完整模拟电压、充分动态范围、未知局部K≤2。原生Q为实现雪崩等效pe，高斯Q为相同数值的条件参数，不作沉积能量；t为有效面首光子。

四种变体为`component_local`、`allocation_local`、`component_profile`、`allocation_profile`，均保留基线K0/K1和原最佳K2。local对原K2候选增加最多70次迭代；profile另外使用同一组最多8个起点，各筛选20次，最好的两个精修70次，并复用local结果。

新坐标为`t1/10,t2/10,100*ln(S),10*r`，S=Q1+Q2，Q1=L+(U−L)r，Q2=S−Q1；L=max(50,S−50000)，U=min(50000,S−50)，r∈[0,1]。完整保持每个Q在50–50000 pe的原域。时间范围仍为起点±24 ns，L-BFGS-B设置相同。数值缩放和有限差分的物理步长也改变；不声称增加信息，不给MLE目标加入Jacobian。

共同profile以观测单事件时刻为锚，前／后4/8/12/20/40/80/160/320 ns配第二候选；第一事件份额0.01/0.03/0.05/0.10/0.15/0.5/0.85/0.90/0.95/0.97/0.99。每点在可行域内数值优化lnS，最多24次、容差10⁻⁴，并保留原值。每方向选4个间隔的最佳份额作共同起点，不读取真实时刻、Q或K。有限剖面与局部优化均无全局最优保证。

各方法在相同280个原生单事件上重新标定K2增峰分布；原单Q分档1500/5000/12500 pe、加一／并列上尾秩p≤0.025和Q≥500保护保持。另存直接沿用原数值门槛的敏感性表。推断worker阻止打开truth、case、ROOT；所有失败和未收敛候选保留。

## 实际结果

严格整对成功要求两真实Q≥500、各Q误差≤20%、各时刻误差≤min(10 ns,间隔/4)、无额外输出。每点40对：

| 输入／方法 | 先强后弱8 ns | 40 ns | 200 ns | 先弱后强8 ns |
| --- | ---: | ---: | ---: | ---: |
| native / baseline | 0/40 | 0/40 | 20/40 | 27/40 |
| native / component_local | 0/40 | 0/40 | 20/40 | 27/40 |
| native / allocation_local | 0/40 | 0/40 | 21/40 | 27/40 |
| native / component_profile | 0/40 | 0/40 | 20/40 | 27/40 |
| native / allocation_profile | 0/40 | 0/40 | 21/40 | 27/40 |
| matched / baseline | 0/40 | 0/40 | 22/40 | 26/40 |
| matched / component_local | 0/40 | 0/40 | 22/40 | 26/40 |
| matched / allocation_local | 0/40 | 0/40 | 22/40 | 26/40 |
| matched / component_profile | 0/40 | 0/40 | 22/40 | 27/40 |
| matched / allocation_profile | 0/40 | 0/40 | 22/40 | 27/40 |

原生40 ns、allocation_profile的K2候选中，31/40代价下降，但强制K2严格成功3/40，盲结果0/40。不能以似然下降冒充逐事件恢复。

- `baseline`：原生空／weak／low／mid／high额外输出帧依次为0/0/0/1/0（每组80帧）。
- `component_local`：原生空／weak／low／mid／high额外输出帧依次为0/0/0/1/0（每组80帧）。
- `allocation_local`：原生空／weak／low／mid／high额外输出帧依次为0/0/0/1/0（每组80帧）。
- `component_profile`：原生空／weak／low／mid／high额外输出帧依次为0/0/0/1/0（每组80帧）。
- `allocation_profile`：原生空／weak／low／mid／high额外输出帧依次为0/0/0/1/0（每组80帧）。

全部条件的Q偏差、sigma68误差宽度、时间误差、漏检、额外输出、误拆／误并、K2输出和时间匹配，以及原门槛对照均保存在CSV。sigma68为匹配后Q相对误差的(p84−p16)/2，未标定沉积能量，不能当能谱FWHM。

预定原生约束包括空／单额外输出≤5%、单事件严格成功≥95%、单Q偏差中位绝对值≤5%；先强后弱40 ns至少提高10个百分点、200 ns不下降；先弱后强8 ns最多下降5个百分点，中／中各点≥95%，六个强弱组平均不下降。

| 变体 | 开发要求全部通过 | 未通过项 |
| --- | --- | --- |
| component_local | False | high_weak40_gain |
| allocation_local | False | high_weak40_gain |
| component_profile | False | high_weak40_gain |
| allocation_profile | False | high_weak40_gain |

四种变体均未同时通过预定原生开发约束，不选择新默认流程。

差值区间用40个母事件块、前20／后20分层、2000次共同权重bootstrap。中／中每组涉及40块，强弱各涉及前20块、每块2对，不能把全部帧独立重抽。它是已公开开发样本的条件性描述；退化[0,0]不表示总体差异确定，高斯draw0不是新光学母事件。

## 成本、验证和边界

原生每帧平均增量：component_local 0.094 s、73次函数调用；allocation_local 0.095 s、71次函数调用；component_profile 3.979 s、3039次函数调用；allocation_profile 4.049 s、3080次函数调用。profile单独使用时包含共同剖面成本；本次共同工作只计算一次。数字排除已缓存基线，逐帧并发wall time不是实时吞吐保证，迭代上限相同不等于函数调用数相同。

7项测试通过：完整电荷域往返、总Q守恒、目标不变、可辨认相关电荷玩具问题、剖面保留原值、空候选保留、真值访问阻断。32帧先导后4进程运行，正式新增1768帧用时3251.3s，先导32.3s。实际环境与源码哈希保留；源HEAD为`bcb27007fadf5df7dbb1d2447638b1f371069a8c`加冻结新增代码，上一轮信息分析改动完整保留。

1800帧／7200组变体候选，7600条主评分／13680条含原门槛评分，252行汇总，3600条强制K2诊断核验通过；基线逐条与原缓存一致，K0/K1不变，原K2代价不被劣化，重放最大相对差0。标量剖面266981点，其中950点优化器未报成功，仍保留最佳已算值；完整局部拟合未收敛统计另列。

不读取真实Q/t的四种方法都属于盲推断；强制K2图和实际分量仅用于诊断／显示。图85-1为预定第0对、40和200 ns、全长0–2048 ns，分别画总波形、强分量和弱分量；若最终只选一个，弱输出画零，绝不拿被拒绝K2候选冒充成功解堆积。

既有Fixsen2004支持条件均值／协方差原则，Liddick2012支持检测／估计分别评价；据AI精读稿待人工审核，本轮未重查PDF。坐标和剖面为项目实现，非论文已验证的GAGG改进；详见`docs/literature-use.md` L60。

局部K≤2、固定双τ、已公开三种组合和三档间隔限制结论。没有验证连续Q／位置失配、混合α/γ、free-τ、随机高率K≥3或实物性能。若同一似然搜索收益不足，后续应检验估计目标、响应／噪声模型和来自训练集的先验；禁止用四档真实能量约束测试答案。

## 文件与重放

本地Notebook为`jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb`第82–86节；旧410个单元及输出保留，新增21单元，11代码单元已执行。固定图0004-82-1、83-1、83-2、84-1、85-1。

`jupyter/artifacts/alpha_pair_charge_allocation_20261027/`保存全部候选、剖面、评分、CSV和图；`jupyter/run_records/alpha_pair_charge_allocation_20261027/`保存输入来源、环境、验证与备份。`freeze.json`保存28个计算文件；报告源码与运行入口另由`report-freeze.json`保护。方法门槛是实际280校准帧导出的`gate-calibrations.json`。

从项目根目录运行，先保留全部本地输入／模型和运行记录：

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_charge_allocation_20261027/run.py infer --workers 4
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_charge_allocation_20261027/report.py analyze
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_charge_allocation_20261027/report.py figures
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_charge_allocation_20261027/report.py notebook
```

计算按index续跑，正式运行要求先导验证；不要在已冻结目录再次prepare。报告用冻结源码，无需切换Git。输入在`input-provenance.json`中逐一记录，图中真实分量另见`plot-truth-provenance.json`；仅克隆仓库不能替代本地模型NPZ、数据数组和运行记录。Jupyter／Notebook／模型数组按项目约定不上传。Notebook重放保留后续用户修改；文档生成另检查已交付文件哈希，避免覆盖新内容。
