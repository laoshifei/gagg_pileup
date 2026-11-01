# 因果孤立响应与采样相位控制

16个历史光学母事件、四档信号，每母四个到达相位、每相位八次数字化，共512次正式数字化，另4次成本先导。组均值能量偏差与相位跨度通过预定工程门；未测盲检出、逐事件能量分辨率或共享恢复效果。

| 信号档 | 平均能量偏差/% | 四相位跨度/百分点 | 平均面积偏差/% |
| --- | ---: | ---: | ---: |
| weak | -2.5260 | 0.0106 | +2.178 |
| low | -0.1496 | 0.0066 | +0.171 |
| mid | -0.3675 | 0.0051 | +0.498 |
| high | +0.0033 | 0.0163 | -0.037 |

每档只有四个独立光学母事件。弱档逐母平均偏差范围为-11.32%至+5.05%，组均值通过不表示所有母事件都通过±5%。同母相位和数字化重复不能计为独立输运。

原生读出为0.25 ns内部步长、0.5 ns采样孔径、2 ns周期；4096 ns数字化，事先声明的观测窗为前2048 ns。事件分量通过仿真祖先标签仅用于响应检查，完整含噪电压另存。单事件能量拟合使用八次含噪电压平均、已知时刻和固定权重，不读取评价能量真值。

本机Python3.10.12/NumPy1.26.4/SciPy1.15.3；运行器将缺少的np.trapezoid名称显式接到np.trapz，旧实验源码不改。旧三帧九拟合重放通过，最大相对E差1.30e-7；C++/Python最大差8.53e-14、梯度相对差1.55e-10、原读出差8.89e-16 mV。事件前均值为零是必要控制，不能证明整体均值或双事件响应准确。

实际源码、配置、编译命令、协议、母事件/种子与哈希随实验保留，数组和Notebook在本地。重放到新目录（项目根目录运行，需本地输入与模型npz）：

```bash
python3 -B pileup/experiments/alpha_pair_causal_controls_20261101/replay.py --ids 0 48 --output jupyter/artifacts/causal_single_replay_new --workers 2
```

完整结果在本地codex_0007第69–74节。预定组均值偏差门±5%、相位跨度1个百分点是本项目工程判据，不是文献阈值。后续共享控制另建实验目录，保留本实验。

## 首次重放的编译准备

Git不包含编译库。保留本地输入与模型数组后，从项目根目录用本机C++17编译器生成同一源码的库：

```bash
mkdir -p pileup/experiments/alpha_pair_causal_controls_20261101/bin
g++ -O3 -std=c++17 -fPIC -shared pileup/experiments/alpha_pair_causal_controls_20261101/source/pileup/analysis/sipm_age_mean.cpp pileup/experiments/alpha_pair_causal_controls_20261101/source/pileup/analysis/sipm_age_mean_jacobian.cpp -o pileup/experiments/alpha_pair_causal_controls_20261101/bin/age.so
```

源文件及模型哈希保持不变；编译器和平台差异需要先通过数值检查。旧基线重放还需前一实验自身的编译库与本地输入。
