# 当前 SiPM 能标的保护能点与独立确认

上一轮硬反演边界使高能端正常电荷涨落被拒绝，选模惩罚反而偏向有系统偏差的简单比例能标。本实验扩大物理标定支持范围，并将内部目标能区与标定保护点分开；保留原四候选和上一轮失败结果。

每个新 ROOT 保存一个 20 μs 原生单 α 窗口，中央入射，源在 5 μs 发射。G4 内部串行，最多四个独立任务；启动至少间隔一秒，输出为项目根目录 `g4/data/YYYYMMDD_HHhMMmSSs.root`。这不是注入率扫描。物理源码、SiPM 参数、固定双 τ 95/351 ns、2 ns 采样保持原实验配置，使用充分动态范围的完整模拟电压，不以 ADC 触顶拒绝。

## 样本与评价协议

| 阶段 | 新 ROOT | 源能量 / MeV | 用途 |
| --- | ---: | --- | --- |
| 先导 | 16 | 1.5、2、3.5、9，各4 | 配置、实际沉积能量和计算成本；不参加拟合/选择 |
| 拟合 | 96 | 1.5、2、3.5、4.5、6.5、9，各16 | 加旧280个fold0单事件，共376 |
| 选择 | 64 | 2、3.5、4.5、6.5，各16 | 加旧280个fold1单事件，共344 |
| 冻结后新确认 | 320 | 1.625、1.8、2、2.5、3.5、4.5、5.486、6.5、7、8，各32 | 1.8/7仅在确认中出现；检查插值和新母事件迁移 |

先导必须通过预定沉积能量覆盖要求才能继续。源1.5/9 MeV为拟合保护点，不进入内部目标选模损失。目标按预先规定的源能组定义，保留组内全部事件，不根据观测Q或真实E删掉边缘样本。

四候选为Q=aE、Q=aL(E)、归一化累计光产额的Q=ax/(1+bx)、单调PCHIP。累计α光产额已含淬灭，不重复Birks。选择只用344个选择母事件：八个内部目标组相对能量MSE等权平均，域外罚1；反演仍限于拟合集实际E范围，不外推或裁到端点。选择集不再混入拟合。

新单事件主要控制是：每个源能组至少95%有限能量估计，能量相对偏差中位数绝对值≤5%。这是32例/组的点估计判据，不提供总体95%有效率保证。另报误差sigma68、能量±5/10/20%比例及所有域外事件。sigma68不是能谱FWHM。

新单事件测试输入为已知初级粒子的真实分量电荷，排除暗计数祖先，故是Q→E_dep能标验证，不是盲事件检出性能。E_source、E_dep及Q/等效pe分别记录。波形算法比较使用已公开历史1600帧的冻结模型GLS/联合NLL输出，只更新能量换算，Q/t/事件数及所有门保持不变。历史重放不能冒充新堆积确认。

## 运行顺序

使用项目 Python 环境，在项目根目录执行：

```bash
python pileup/experiments/sipm_energy_guardband_20261029/run.py pilot --workers 4
python pileup/experiments/sipm_energy_guardband_20261029/run.py pilot-check
python pileup/experiments/sipm_energy_guardband_20261029/run.py fit --workers 4
python pileup/experiments/sipm_energy_guardband_20261029/run.py select --workers 4
python pileup/experiments/sipm_energy_guardband_20261029/energy.py calibrate
python pileup/experiments/sipm_energy_guardband_20261029/run.py confirmation --workers 4
python pileup/experiments/sipm_energy_guardband_20261029/energy.py export
python pileup/experiments/sipm_energy_guardband_20261029/energy.py predict
python pileup/experiments/sipm_energy_guardband_20261029/energy.py evaluate
```

`energy.py`及分析源码在先导通过后、补充标定生产前冻结；能标冻结是确认生产的前置条件。执行下一步前必须完成前一步。已完成ROOT会跳过；失败尝试保留并用新秒级输出名重试同一种子。不得同时运行两个生产阶段，因为它们共同更新运行清单。运行中的清单可能正在写入，结果分析在阶段结束后读取。

实际C++源码与配置在`source/`，二进制身份在`freeze.json`及运行清单；Python运行器用`analysis-source/`中的真实源码副本。重建不依赖切换Git提交，仍需可用的G4/ROOT依赖。旧标定数据与历史拟合输出属于本地输入，不随GitHub分发。

本地输出：`jupyter/artifacts/sipm_energy_guardband_20261029/`。运行记录：`jupyter/run_records/sipm_energy_guardband_20261029/`。研究过程追加到codex_0007，所有新ROOT/实际宏/种子追加到Run_Catalog；已有笔记输出及历史数据保留。

## 适用范围

中央纯α模拟、当前代理器件参数、固定双τ和充分动态范围。单事件静态能标通过不保证共享微胞堆积中仍无偏，亦不验证α/γ混合、free-τ、未知K≥3、高率能谱或实物装置。原理解释复用Rosado2018既有待审笔记，但本项目有理饱和候选不是该论文完整公式的复现。


## 实际完成结果

496个计划母事件全部完成并核验；实际创建499份ROOT，其中3份处理未完成尝试保留但不纳入统计，按原种子重试关系见retry-audit.json。选定`empirical_pchip`。新确认十个源能组通过预设控制，有限逆变换320/320，能量±10%内318/320，最大组偏差中位数绝对值3.48%。各组条件误差宽度、点估计和独立性限制见结果摘要。

历史同一波形/拟合的弱组合能量±10%整对恢复为模型GLS 31.67%、联合NLL 38.17%，配对差+6.50个百分点，母事件块95%区间[+3.00,+10.17]。原电荷、定时与事件数未改变；这不是新堆积确认。默认拟合与选门保持，下一步检查共享微胞条件下的能量响应及高率完整窗口。
