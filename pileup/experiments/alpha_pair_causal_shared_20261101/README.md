# 因果共享微胞响应：均匀近似未通过短间隔控制

2026-11-01完成因果响应的分阶段控制：16个历史母事件、4相位、512次孤立数字化通过预定的组均值能标与相位检查；随后4对历史母事件的960次共享/孤立、恢复开关对照未通过共享响应门。强→弱8/40 ns的弱分量面积变化分别为−14.48%/−14.47%，均匀微胞模型仅预测−8.41%/−8.27%，相差6.06/6.20个百分点。按预先协议停在响应诊断，没有双脉冲能量拟合、新ROOT或默认算法替换。

| 顺序/间隔 | 数字化弱分量面积变化/% | 模型预测/% | 差值/百分点 | 恢复关闭/% |
| --- | ---: | ---: | ---: | ---: |
| high_weak，8 ns | -14.475 | -8.411 | +6.064 | +0.762 |
| high_weak，40 ns | -14.472 | -8.269 | +6.203 | +0.938 |
| high_weak，200 ns | -7.946 | -4.295 | +3.652 | +0.453 |
| high_weak，1200 ns | +1.050 | -0.236 | -1.286 | +1.555 |
| weak_high，40 ns | -10.862 | -6.493 | +4.369 | +0.956 |

此处端点是同一2048 ns观测窗内、按祖先标签分解的弱分量电压面积，相对其孤立响应的变化；不是逐事件能量误差或盲恢复率。4对母事件跨顺序/间隔复用，每条件八个数字化种子。正式960次，另12次成本先导。

参数在评价前冻结：14400微胞、29 ns恢复、95/351 ns形状、快分量0.3、均匀照明与恢复修正强度1。模型以已标定的孤立电压为基础，加上年龄模型的共享减孤立修正；这是项目近似。1→0.5 ns网格对组均值抑制预测的最大变化约0.00042个百分点，无法解释6个百分点的偏差。恢复关闭下0.45–1.56%的有限样本差值保留，不能写成严格零。

四对母事件的无权重空间重叠指标Ncell×Σp强p弱平均2.158（均匀为1）。它提示照明非均匀性值得检验，未证明是唯一原因，也未用于调整本轮参数。下一步可保持光子时刻/能量不变，仅随机化有效面内位置作对照；随后用独立开发母事件标定空间分布，冻结后再检验共享响应。

嵌套母事件/数字化重采样的模型误差95%百分位范围：8 ns为[4.42,7.63]、40 ns为[4.81,7.56]个百分点。只有四对历史母事件，该范围仅描述选中样本；本轮门按预设点估计判定，不作普适分辨率或全体算法结论。

报告第一次读取空间分布时，负的“有效面首光子”相对时间触发输入检查。模拟输入已加128 ns且未受影响；仅在空间映射诊断中平移时间原点。原source、freeze及所有模拟输出保留；修正报告在report-source，入口run_report.py，原因与哈希另存。

重放与报告（项目根目录，保留本地输入/模型）：

```bash
python3 -B pileup/experiments/alpha_pair_causal_shared_20261101/replay.py --ids 1 --output jupyter/artifacts/causal_shared_replay_new --workers 1
python3 -B pileup/experiments/alpha_pair_causal_shared_20261101/run_report.py
python3 -B pileup/experiments/alpha_pair_causal_shared_20261101/run_plots.py
python3 -B pileup/experiments/alpha_pair_causal_shared_20261101/run_audit.py
```

四幅图0007-70-1至73-1与完整过程在本地codex_0007第69–74节。原电压、逐分量、种子/母ROOT映射与源哈希保留。此分支停在已知输入响应检查，未开展未知K/t、α/γ或高率谱验证。

## 首次重放的编译准备

Git不包含编译库。保留本地输入与模型数组后，从项目根目录用本机C++17编译器生成同一源码的库：

```bash
mkdir -p pileup/experiments/alpha_pair_causal_shared_20261101/bin
g++ -O3 -std=c++17 -fPIC -shared pileup/experiments/alpha_pair_causal_shared_20261101/source/pileup/analysis/sipm_age_mean.cpp pileup/experiments/alpha_pair_causal_shared_20261101/source/pileup/analysis/sipm_age_mean_jacobian.cpp -o pileup/experiments/alpha_pair_causal_shared_20261101/bin/age.so
```

源文件及模型哈希保持不变；编译器和平台差异需要先通过数值检查。旧基线重放还需前一实验自身的编译库与本地输入。
