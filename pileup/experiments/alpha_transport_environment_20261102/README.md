# 当前Geant4环境重建与同种子检查

当前 Ubuntu/WSL 的 Geant4 11.3.2、ROOT 6.34.10 已完成主工程和实验快照构建。16条跨环境同种子对照全部正常输出；主工程在当前环境重跑四档各一例，28组光子/电压数组逐项一致。文档中的单α命令也已通过，共新增21个有效原生ROOT。

旧 Fedora 环境为 Geant4 11.4.2、ROOT 6.40.04。16对初始随机状态全部相同，但沉积能量、产生光子数、排序后的完整到达光子记录和原生信号均未逐事件精确复现。最大沉积能量差为0.10252 keV，到达光子数相对差为[-1.449%, +2.297%]。这些是同种子工程对照，不能据此断定总体响应有偏差，或把差异单独归因于Geant4版本；数据集、编译器和平台也发生了变化。

| 检查 | 结果 |
| --- | --- |
| 冻结G4源码与主工程 | 39文件内容完全相同；未改物理参数 |
| 主工程／快照CTest | 各1/1通过 |
| 宏路径及光子记录针对性测试 | 2项通过 |
| 跨环境初始随机状态 | 16/16相同 |
| 跨环境E_dep／产生光子数／完整到达记录 | 分别0/16逐项相同 |
| 当前环境主工程重跑 | 4/4通过，28组数组完全一致 |
| 文档单α命令 | 1/1通过 |
| 新ROOT／新增独立光学母样本 | 21／0 |

从项目根目录运行主工程：

```bash
source /opt/root/root63410/bin/thisroot.sh
source /opt/geant4/geant4-v1132/bin/geant4.sh
cmake -S g4 -B g4/build -DCMAKE_BUILD_TYPE=Release \
  -DGeant4_DIR=/opt/geant4/geant4-v1132/lib/cmake/Geant4 \
  -DROOT_DIR=/opt/root/root63410/cmake
cmake --build g4/build -j 4
ctest --test-dir g4/build --output-on-failure
cd g4/build
./gagg ../macros/run_single_alpha_sipm.mac
```

旧build完整保留在.local_archive/build_migration_20261102/g4_build_from_previous_machine。
当前二进制在g4/build/gagg，冻结实验二进制在本目录bin/gagg，均需在相同外部依赖下使用；构建产物、ROOT及Notebook不随Git分发。

本实验source/g4保留39份实际源码/配置，analysis-source保留检查与报告源码；build.json记录独立快照构建命令。外部依赖本身不包含在快照中。原16条控制的协议、种子及哈希在protocol.json、plan.json、design-freeze.json，main-replay/main-smoke各有单独计划、结果和源码冻结。tables/root_catalog.csv给出21个ROOT与实际宏/日志/哈希，runinfo保存全部RunInfo；同环境比较为数组逐项相同，不要求ROOT文件哈希一致。

```bash
python3 -B pileup/experiments/alpha_transport_environment_20261102/run.py verify
python3 -B pileup/experiments/alpha_transport_environment_20261102/validate.py
```

validate.py只读核对已有数据，不重复生产。run.py的preflight/complete用于原16条任务的续跑，已完成时拒绝重复；replay_main.py与smoke_main.py也保护已完成结果。新生产应另建实验与输出，不能重写既有快照。

用户明确要求使用当前环境，因此后续以本机11.3.2/6.34.10继续。336个正式强母＋8个成本先导的草案保留，尚未运行；后续须另建实验快照，将问题明确为旧冻结模型在当前环境的新母事件迁移检验。先检查四档孤立响应与读出一致性，再根据冻结方案验证逐母投影MSE、均值不确定度和孤立E控制。不把跨环境差异当作独立同环境的均值修正收益；不使用已暴露对照调修正。

本阶段为运行环境检查，没有新增独立母样本、离线数字化或能量拟合。旧的跨环境精确复现门仍记为未通过；用户随后指定当前环境，这一选择和21条最终状态分别记录在current-environment-decision.json与environment-status.json。未检验总体分布等效或实际探测器性能。本地codex_0007第119–122节保留两幅固定编号图；L95记录研究影响。
