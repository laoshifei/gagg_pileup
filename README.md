# GAGG 探测器响应与脉冲堆积模拟

本项目使用 Geant4 模拟 GAGG(Ce) 中的能量沉积、闪烁发光及光子输运，记录到达虚拟 SiPM 的光子时间序列，再用 Python 分析波形。研究对象是 p–¹¹B 相关的 α 粒子、热轫致辐射背景和 ¹²C 俘获/退激 γ，目标是研究计数率升高时的随机堆积，以及 α 数量和能谱的恢复能力。

Geant4 保留到达光子，并可通过 `G4VDigitizerModule` / `G4DigiManager` 在窗末生成 SiPM 雪崩、电流、电压和 ADC，保存到同一 ROOT 的 `SiPMEvent`。当前单粒子默认器件为滨松 **S13360-6050PE**（25°C、3 V 过压），采用数据表典型参数和近似读图 PDE；恢复时间使用同系列代理值，光耦合与读出仍需标定。[SiPM 响应说明](docs/sipm-response.md) 包含统一的原生与离线入口、逐项器件依据和验证结果。原有光子计数波形与数字化电压分别保留；已有分粒子波形来自仿真真值。当前已完成纯 α 和 α/γ 的盲拟合与机器学习对照、弱 α 诊断，以及纯 α 单事件能标的独立确认；盲重建仅输入观测波形，真值辅助诊断单独标明。共享微胞下的能量迁移、高率完整能谱、模型失配与实物验证仍待完成。

## 项目目录

| 目录 | 内容 |
| --- | --- |
| [g4](g4/) | Geant4 工程，包含源码、宏、分析脚本和结果图。 |
| [pileup](pileup/) | 纯 α 与 α/γ 混合计数率实验：配置、可独立重跑的源码快照、运行器和数据检查。 |
| [analysis](analysis/) | 单粒子 SiPM 分析、孤立 α/γ 波形拟合和拟合窗口扫描。 |
| [tools](tools/) | 各系统独立 Python 环境、Notebook 启动及本机内核注册工具。 |
| [11BH_alpha_spectra](11BH_alpha_spectra/) | 简化 p–¹¹B 单 α 能谱、生成脚本、CSV 和检查图。 |
| [bremsstrahlung](bremsstrahlung/) | 参数化热轫致辐射光子谱、生成脚本、CSV 和检查图。 |
| [stopping_power](stopping_power/) | 阻止本领表、Birks 光产额曲线和生成程序。 |

ROOT 原始数据、训练权重、分析数组、构建产物和本地研究笔记不随仓库分发，排除规则见 [.gitignore](.gitignore)。实验源码、配置、宏和运行清单随仓库保留；重放既有结果仍需另行提供对应的本地输入与权重，并在本机编译可执行程序。

## 当前模型

模拟流程：入射粒子 → GAGG 内能量沉积 → 闪烁光子产生与传播 → 完整时间窗的光子排序 → 可选 SiPM 数字化 → ROOT。后续可分析电信号，也可用保留的光子输入重新离线扫描器件参数。

- GAGG 尺寸为 10 × 10 × 1 mm，入口位于 z = 250 mm；前 Mylar 膜厚 5 μm，侧膜厚 20 μm。
- 背面虚拟光子记录层为 10 × 10 × 0.01 mm，记录进入的光子并终止其轨迹；S13360-6050PE 数字化只接收中央 6 × 6 mm 有效面内的光子，区域外光子仍保留在 ROOT。当前记录层采用理想折射率匹配，尚未加入真实封装及光学脂。
- 堆积模式下，一个 Geant4 event 是一个时间窗。各源独立抽样 N ~ Poisson(R × T)，再在窗内均匀抽样发射时间。
- 默认 T = 10 μs，可用 `/gagg/run/windowLength` 配置。纯 α 首轮实验使用完整 20 μs 窗和中央 10 μs 评价区。部分历史 ROOT 树标题写着 1 ms；解释历史数据时使用其 RunInfo.t_length_ps。
- 常用宏启用圆盘源并将粒子瞄准晶体前表面；速率表示有效注入粒子/光子数率，不直接等于等离子体总产额或探测事件率。
- α 使用有效单粒子能谱；¹²C γ 使用单光子谱线混合，当前没有生成完整三 α 关联事件或 γ 级联符合对。

几何、时间窗、光学与固定能谱参数见 [Constants.hh](g4/include/Constants.hh)；源抽样见 [PrimaryGeneratorAction.cc](g4/src/PrimaryGeneratorAction.cc)；光产额与边界设置见 [DetectorConstruction.cc](g4/src/DetectorConstruction.cc)。

## 编译与运行

需要 CMake、C++17 编译器、带 UI/可视化组件的 Geant4、ROOT 的 Core/RIO/Tree，以及 nlohmann_json ≥ 3.2（读取 SiPM 配置；优先使用独立 CMake 包，也支持 ROOT 自带头文件）。既有生产数据使用 Geant4 11.3.2、ROOT 6.34.10；2026-10-09 本机使用 Geant4 11.4.2、ROOT 6.40.04，主工程、两个冻结实验和阻止本领工具已重新构建。环境检查记录保存在本地 `jupyter/Run_Catalog.ipynb`，历史结果仍按各自的运行记录解释。

从项目根目录运行，先按本机安装方式加载 Geant4/ROOT 环境：

```bash
cmake -S g4 -B g4/build
cmake --build g4/build -j 2
cd g4/build
./gagg ../macros/run_all_pileup.mac
```

未指定输出文件时，程序使用相对目录 `../data` 和秒级文件名 `YYYYMMDD_HHhMMmSSs.root`，因此上述运行方式将新结果写入 `g4/data/`；同名时等待下一秒。宏可用 `/gagg/run/outputFile` 指定相对于启动目录的 ROOT 路径；显式指定的已有同名文件会被拒绝覆盖。无参数运行 `./gagg` 会进入交互可视化模式。

当前 [run_all_pileup.mac](g4/macros/run_all_pileup.mac) 配置为轫致辐射 10⁶ Hz、¹²C γ 10 Hz、α 10⁶ Hz，生成 10 个时间窗。每窗期望初级粒子数约 20。另有单源宏与 [低速率三源宏](g4/macros/run_all_p11b_sources_pileup.mac)。大规模运行前应先用独立输出目录和少量事件检查配置。

单 α / γ 电信号的入口如下。两个宏分别发射一个 5.486 MeV α 和一个 661.657 keV γ，20 μs 窗内在 5 μs 发射，中央正入射；每次只产生一个 ROOT。单γ宏选用固定32例样本中有信号的一例种子，便于查看电信号；该宏不是探测效率样本。固定样本运行器仍从预定种子生成并保留所有γ，包括仅含暗计数和电子学噪声的无沉积事件。

```bash
ctest --test-dir g4/build --output-on-failure
cd g4/build
./gagg ../macros/run_single_alpha_sipm.mac
./gagg ../macros/run_single_gamma_sipm.mac
```

从项目根目录运行 `bash tools/run_gagg_python.sh -m analysis.single_particle_sipm run --alpha 4 --gamma 32` 可保存固定样本及三组随机数对照、实际源码与宏、ROOT 索引、读出交叉检查和图。默认数字化关闭，原生 SiPM 运行启用后要求 `/run/beamOn 1`。具体参数、削顶标记与模型限制见 [原生数字化](docs/sipm-response.md#native-sipm)。 本次完整改动与实际结果另见本地独立笔记 [codex_0001_SiPM_Digitization.ipynb](jupyter/codex_0001_SiPM_Digitization.ipynb)。

## 纯 α 堆积实验

首轮沿用当前 p–¹¹B 有效单 α 能谱，扫描注入率 10⁵、10⁶、10⁷、10⁸ Hz，关闭光子背景。准备与运行方法见 [pileup/README.md](pileup/README.md)。每个实验目录保存自己的源码、宏、运行器和配置；以后修改主工程时，旧实验仍能自行编译重跑，无需切换 Git 提交。Geant4、ROOT 等外部依赖须保持可用。

运行接口支持以下命令，均应在 `/run/beamOn` 前设置；`preWindow` 和 `postWindow` 定义中央评价区，不会删除保护区内的粒子或光子：

```text
/gagg/run/seed 260923101
/gagg/run/windowLength 20 us
/gagg/run/preWindow 5 us
/gagg/run/postWindow 5 us
/gagg/run/outputFile waveform.root
```

当前生产采用 G4 串行、每个 ROOT 一个波形窗、最多 8 个独立任务。ROOT 统一存入 `g4/data/`，使用原有秒级文件名，任务启动至少间隔 1 秒。具体入口与运行记录见 [pileup 说明](pileup/README.md)。

截至 2026-09-30，首轮四档各 100 个独立时间窗已全部完成并通过逐窗数据检查，生产 ROOT 合计约 54.27 GB。2026-10-13 另建 [S13360 纯 α 计数率先导实验](pileup/experiments/alpha_s13360_rate_fit_20261013/README.md)，保留源码、配置、宏、种子和运行清单；四档比较样本分别为 16、4、3、1 个独立窗，未使用首轮保留测试集。2026-10-14 补充 20、40、60、80 MHz 各 3 窗，并按用户指定的“采集动态范围足够”前提，对全部 36 窗进行 [完整电压统一重分析](pileup/experiments/alpha_full_signal_20261014/README.md)。当前 [盲重建基线](pileup/analysis/alpha_pileup_fit.py) 从完整含噪电压估计脉冲数、到达时间和电荷；本轮不考虑 ADC 位数、量化和削顶，旧 ADC 整窗拒绝结果只作为历史记录。同日完成 [时序 CNN 与拟合的比较](pileup/experiments/alpha_ml_comparison_20261014/README.md)：26 窗训练、10 窗验证、17 个新种子测试窗，20–100 MHz 的事件恢复率改善，误检增加；另保存近邻/幅度外推失败及高率偶然时间匹配对照。完整过程记录在本地 [独立 ML Notebook](jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb)。这些小样本结果不构成最终能谱或实物性能验证。下一步仍需条件能量标定、模型失配与扩大独立样本评价。方法依据与适用条件见 [文献采用记录](docs/literature-use.md)；运行明细与研究计划在本地 Jupyter 笔记中维护。

<a id="ml-multimethod-20261014"></a>
同日扩展为[五种流程的共同测试](pileup/experiments/alpha_ml_multimethod_20261014/README.md)：新增U-Net、局部Transformer及CNN＋物理拟合，保留原CNN与盲拟合；两网络各训练2400步，22个新独立窗口覆盖八档源率，原17窗作回归。算法参数仅由验证集选择，完整电压输入和足够动态范围前提保持。不同率的改善并不一致，受控近邻、幅度外推失败及高密度偶然匹配均原样记录；详见[独立Notebook第二轮](jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb#ml-multimethod-20261014)和[采用记录L21](docs/literature-use.md#ml-multimethod-20261014)。

2026-10-22 完成[条件响应与相关重拟合的独立检验](pileup/experiments/alpha_pair_conditional_20261022/README.md)：利用旧母样本标定电荷／邻居相关响应，冻结参数后另生成240个单α光学母ROOT，构造2320帧共同微胞双事件及空／单控制。中／中电荷、8 ns间隔的整对成功由旧固定模板混合的0/40提高至新GLS混合的40/40，两个冻结网络种子均如此；低Q偏差和先强后弱漏检仍明显，不能推广成普适8 ns分辨率或最高计数率。完整比较、误报控制及限制见[采用记录L50](docs/literature-use.md#alpha-pair-conditional-20261022)，本地codex_0004第41–44节保留带编号图和完整波形。

2026-10-23 完成[相关权重正则化与残差新增候选](pileup/experiments/alpha_pair_residual_20261023/README.md)的独立检验：240个新光学母ROOT、2320帧、14条处理流程，网络权重保持不变。同批低Q单事件的匹配后电荷偏差中位数由原GLS的-13.59%变为-3.82%；相对原来同样允许未知K选择的条件拟合，先强后弱40 ns整对成功由0/40变为39/40、200 ns由37/40变为40/40。低Q单事件额外输出8/80；主方案未通过全部预设控制要求，完整新流程暂不替换默认方案。结果限于纯α、固定双τ、局部K≤2，详见[采用记录L51](docs/literature-use.md#alpha-pair-residual-20261023)；本地codex_0004第45–48节保存编号图与完整波形。

2026-10-23 完成[真实弱α与经验增峰筛选](pileup/experiments/alpha_pair_null_gate_20261023/README.md)的独立检验：另生成320个新母ROOT、3200帧；约700 pe原生单事件Q偏差中位由-9.74%变为-0.89%，低Q单事件额外输出6/80→0/80。本批空／单控制全部通过预定要求，真实弱双事件仍有明显失败区间；原方法、仅权重调整、经验筛选及已知K诊断全部保留。范围仍限纯α、固定双τ和局部K≤2，详见[采用记录L52](docs/literature-use.md#alpha-pair-null-gate-20261023)；本地codex_0004第49–52节保留编号图与分量分轴的完整波形。

<a id="alpha-gamma-mixed-20261015"></a>
2026-10-15 新增 [α/轫致辐射 γ](pileup/experiments/alpha_gamma_brems_fit_20261015/README.md) 与 [α/661.657 keV γ](pileup/experiments/alpha_gamma_mono662_fit_20261015/README.md) 的盲拟合：α 注入 1/10 MHz、γ 均为 10 MHz，各条件 8 个独立原生共享微胞窗口。拟合未知脉冲数、时刻、电荷和类型，保持完整含噪电压输入；另有147次空/单/双脉冲控制。以下恢复要求真值Q≥500 p.e.、时间误差≤10 ns、类型正确及电荷误差≤20%，不是能量分辨率。

| γ 来源 | α 注入率 / MHz | α 完整恢复 | γ 完整恢复 |
| --- | ---: | ---: | ---: |
| 轫致辐射谱 | 1 | 44/64（68.8%） | 2/30（6.7%） |
| 轫致辐射谱 | 10 | 270/571（47.3%） | 3/42（7.1%） |
| 661.657 keV | 1 | 49/58（84.5%） | 29/32（90.6%） |
| 661.657 keV | 10 | 282/605（46.6%） | 12/31（38.7%） |

详细的检出/分型、α样本纯度、漏检及误拆分/误合并、条件电荷误差和偶然匹配对照见本地 [新 Notebook](jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb)。固定模板仅来自两条旧孤立示例，结果受阈下信号和模型失配影响，不能推出通用最大计数率。采用依据与待验证事项见 [L22](docs/literature-use.md#alpha-gamma-mixed-20261015)。

2026-10-15–16 实施 [独立标定与稀疏加权拟合](pileup/experiments/alpha_gamma_regularized_fit_20261015/README.md)：另生成56个标定、18个单粒子开发检查、4个混合验证及8个冻结后测试窗，加入电荷相关模板、固定相关噪声权重、稀疏选择与无L1幅度重拟合，并保留弱分量。相同新测试窗中，α源10 MHz时，α完整恢复在轫致辐射背景下由81/153增至97/153，在单能γ背景下由69/155增至91/155；α源1 MHz时两种背景均退步，γ未稳定改善。额外预测53→19、漏检60→61、误合并标志24→31，不能称为全面提高。实际结果、失败控制、纯度和条件电荷误差见 [Notebook第13节](jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-regularized-20261015)及 [L24](docs/literature-use.md#alpha-gamma-regularized-20261015)；每条件仅两窗，尚无占用依赖标定、沉积能谱或实物验证。

<a id="alpha-gamma-model-selection-20261016"></a>
2026-10-16 继续[联合时刻/类型搜索与残差模型比较](pileup/experiments/alpha_gamma_model_selection_20261016/README.md)，按α、γ同等重要选择方案并另生成8个测试窗。同模板控制24/40→40/40；新原生低α率轫致辐射条件相对上一版为alpha 12/16→14/16、gamma 7/12→11/12。高率为brems：α 81/142→78/142，γ 0/6→1/6；mono662：α 85/143→86/143，γ 0/6→1/6，仍不足以可靠恢复两类。真值辅助诊断和标定覆盖检查说明候选搜索之外的问题仍需处理；完整同窗三方法比较见[Notebook第14节](jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-model-selection-20261016)和[L25](docs/literature-use.md#alpha-gamma-model-selection-20261016)。

<a id="alpha-gamma-ml-20261016"></a>
2026-10-16—17 完成[α＋γ多任务CNN、U-Net与网络＋物理拟合](pileup/experiments/alpha_gamma_ml_20261016/README.md)，另生成16个冻结后测试窗。混合流程相对扩展标定拟合的α完整恢复374/641→416/641，γ35/75→28/75，额外预测41→64；未达到两类同时改善，保留现有拟合默认流程。弱γ错分仍突出，匹配模板控制也显示网络尚未掌握稳定分离。完整结果、失败分层及独立性核验见[Notebook第8节](jupyter/codex_0006_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb#alpha-gamma-ml-implementation-20261016)和[L29](docs/literature-use.md#alpha-gamma-ml-20261016)。

<a id="alpha-gamma-prior-ablation-20261019"></a>
2026-10-19 完成[电荷/类别先验消融](pileup/experiments/alpha_gamma_prior_ablation_20261019/README.md)：在20个已公开开发窗固定候选、时刻、模板及权重，拆开两个显式网络惩罚。无/仅电荷/仅类别/两者的α完整恢复为483/791、517/791、503/791、520/791，γ为42/86、37/86、42/86、38/86；仅电荷项使γ错分增加。该结果支持优先检查跨类别共用网络Q的约束，仍需完整流程及新独立测试验证。过程追加到[ML Notebook第9节](jupyter/codex_0006_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb#alpha-gamma-prior-ablation-20261019)与[L30](docs/literature-use.md#alpha-gamma-prior-ablation-20261019)。

<a id="alpha-gamma-full-prior-20261019"></a>
2026-10-19 将[仅类别惩罚接回完整混合流程](pileup/experiments/alpha_gamma_full_prior_20261019/README.md)：20个已公开开发窗通过预定逐类恢复/precision/误报筛选后，固定λQ=0、λC=10及旧网络，另生成16个新独立原生窗口。新测试α完整恢复345/643→382/643，γ30/79→29/79，额外预测40→54。本批独立数据未通过与开发相同的筛选，不替换默认拟合流程。未通过项：gamma_pooled_recovery, cell_brems_1000000_alpha, cell_brems_10000000_gamma。详见[Notebook第10节](jupyter/codex_0006_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb#alpha-gamma-full-prior-20261019)与[L31](docs/literature-use.md#alpha-gamma-full-prior-20261019)。

## 分析已有数据

Python 分析及 Notebook 依赖见 [requirements-analysis.txt](requirements-analysis.txt)，包括 numpy、matplotlib、uproot、awkward、pandas、scipy 和 Notebook 执行组件。[requirements-analysis.lock.txt](requirements-analysis.lock.txt) 固定已验证的依赖版本；各系统分别安装对应的二进制包，不共享虚拟环境。环境工具保存在随仓库分发的 [tools](tools/) 中；`jupyter/` 下的研究笔记、图像和运行记录仍只在本地维护，GitHub 克隆不包含这些文件。

项目位于共享文件夹时，**Debian 和 Fedora 各用本机用户目录中的环境**。目录按项目路径、发行版、架构和 Python 次版本区分，默认位于 `~/.local/share/gagg-pileup/envs/`，也支持本机 `XDG_DATA_HOME`。七本当前研究笔记统一选择 **GAGG pileup (本机环境)** 内核（`gagg-pileup`）。切换系统时，Notebook 与 ROOT 路径保持相对引用；Python、依赖、内核注册和临时连接文件留在本机。

已接入项目内核的本机，日常继续使用原有 Notebook 启动入口，无需切换到项目目录：

```bash
bash ~/start_jupyter.sh
```

本机入口先检查 GAGG 环境，自动准备依赖，并把 `gagg-pileup` 内核注册到原有 Notebook 环境，再按原设置打开 `/media/sf_works`。直接打开 `GAGG_pileup/jupyter/` 下的笔记，即可使用已保存的项目内核设置；原有其他内核仍可使用。环境可用时直接复用，首次使用、依赖清单更新或系统 Python 次版本变化时自动准备，无需提前执行 setup。首次准备需要下载依赖；系统需提供 Python 3 和 venv/ensurepip 组件。安装记录保存在 `jupyter/run_records/`，中断后可重试同一个本机启动入口。

新机器只需在本机启动脚本执行 Jupyter 之前接入一次以下调用。`notebook_dir` 指向共享目录，`venv_dir` 指向该启动脚本原来使用的本机 Jupyter 环境；按该机器实际目录调整这两个变量。此后环境检查与注册由原脚本自动完成：

```bash
gagg_environment="$notebook_dir/GAGG_pileup/tools/gagg_environment.py"
if [[ -f "$gagg_environment" ]]; then
    "${GAGG_BASE_PYTHON:-/usr/bin/python3}" "$gagg_environment" register --prefix "$venv_dir"
fi
```

2026-10-14 已接入当前 Fedora 上的 `~/start_jupyter.sh` 和 `~/software/jupyter/start_jupyter.sh`；原文件备份和实际启动验证见本地 [入口记录](jupyter/run_records/local_jupyter_entry_20261014_093600/local_server_checks.json)。旧的 `jupyter/gagg_environment.py` 入口继续转交给 `tools/`，保留现有启动方式。Debian 的本机脚本尚需切换系统后接入并实测。

从项目根目录运行 `bash tools/start_gagg_jupyter.sh` 可独立启动经典 NbClassic 界面，也适用于没有原启动脚本的新机器。存在本地 SiPM 笔记时默认打开该笔记，否则显示项目文件列表。两种入口均可传入 `--no-browser`、`--ServerApp.port=8890` 等 Jupyter 参数。查看本机环境位置或执行分析脚本（分析入口同样会自动准备环境）：

```bash
python3 tools/gagg_environment.py info
bash tools/run_gagg_python.sh -m unittest pileup.analysis.test_alpha_response -v
```

项目内备用启动器按脚本位置定位项目，允许从其他目录启动；分析脚本仍按各自说明选择工作目录。自定义解释器可设置本机 `GAGG_BASE_PYTHON`，自定义环境可设置本机 `GAGG_ENV_DIR`；自定义环境应放在本机文件系统中。本机入口保留原来的 Notebook 配置，备用 NbClassic 入口使用 [classic_notebook_config.json](tools/classic_notebook_config.json)，如存在本地 `jupyter/classic_notebook_config.json` 则优先使用它。

共享 `.venv` 曾分别基于 Python 3.14.4 和 3.13.5 创建，切换系统后出现缺少标准库的内核启动错误；它不再是 Jupyter 默认入口。2026-10-09/10 与 2026-10-13 的环境记录仍按当时条件保留，用于追溯历史结果。2026-10-13 的旧环境和笔记备份位于 `jupyter/run_records/notebook_repair_20261013_112957/`；本次跨系统环境调整、原笔记和输出备份位于 `jupyter/run_records/notebook_repair_20261014_084256/`。本次只调整运行环境与笔记执行入口，不据此推断新旧 Geant4 的物理响应一致性。 本次已在 Fedora 44 / Python 3.14.3 实际执行七本 Notebook，35项分析测试及经典页面/内核检查通过；Debian目录选择已检查，实际执行待切换系统确认。[验证结果](jupyter/run_records/notebook_repair_20261014_084256/final_validation.json) 保留本次环境和逐本执行状态。

取得对应 ROOT 数据后，从项目根目录运行：

```bash
cd g4/analysis
bash ../../tools/run_gagg_python.sh plot_pileup_waveform.py ../data/20260624_16h28m34s.root --bin-width-ns 1.0
```

脚本只画第一个 event，PNG 保存到启动命令所在目录。它使用 `uproot.source.file.MemmapSource`，只加载首窗绘图需要的光子时间和粒子标签；不加载无关的位置、能量或后续窗光子数组。高计数率下绘制全部真值分量仍需预留内存。

单粒子比较示例与分析工具说明见 [分析目录说明](g4/analysis/README.md)。

纯 α 响应核验与条件能标模块见 [pileup/analysis/alpha_response.py](pileup/analysis/alpha_response.py)，使用方法见 [pileup 说明](pileup/README.md#α-响应核验与能量标定)。首批 100 个校准/验证窗的检查、长尾模板及能标比较已完成，实际输出记录在本地 `jupyter/Alpha_Response_Calibration.ipynb`；当前电荷重建基线尚未完成与此条件能标的连接。

孤立 α/γ 波形拟合与短窗口扫描分别见 [alpha_gamma_waveform_fit.py](analysis/alpha_gamma_waveform_fit.py) 和 [alpha_gamma_fit_window_scan.py](analysis/alpha_gamma_fit_window_scan.py)。默认示例需要另行提供对应 ROOT 和冻结运行记录；相关本地路径在模块中列明。两条孤立模拟波形只能展示模型行为，不能据此推断分类效率。

2026-10-19 新增[因果短前缀的独立检验](pileup/experiments/alpha_gamma_short_window_20261019/README.md)：149个标定窗、45个开发窗，另生成30α＋144γ测试窗，比较固定95/351 ns两幅度及完整类别模板的OLS/GLS。开发选择的“触发前20 ns＋触发后300 ns”在新测试中使19个α和12个γ的Q≥500 pe目标全部达到电荷误差≤5%且类型正确；100 ns时分别仅14/19和7/12。该结果约束完整响应形状，只说明这批孤立模拟的电荷恢复；小样本不能认证普适最短窗，也不回答自由双τ或高率堆积性能。[本地Notebook第10节](jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#short-window-implementation-20261019)保留逐类效率、误差、零信号及停止失败结果；[采用记录L35](docs/literature-use.md#short-window-implementation-20261019)记录文献边界。

2026-10-20 完成[自洽短前缀拟合的七流程比较](pileup/experiments/alpha_gamma_advanced_fit_20261020/README.md)：另生成168标定、84开发和168冻结后原生测试。开发选中`expanded_gls`；同批新测试在触发前20＋后150 ns时，“Q误差≤5%且类型正确”由原标定GLS的alpha16/16、gamma11/15，经扩充标定GLS的alpha16/16、gamma9/15，到所选方法的alpha16/16、gamma9/15。全部100/150/200/300/500 ns、局部区间覆盖及双脉冲假设失效保存在[本地Notebook第11节](jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#advanced-fit-implementation-20261020)和[采用记录L37](docs/literature-use.md#advanced-fit-implementation-20261020)。结论限于同模型孤立固定完整形状的实现电荷，不是free双τ或高率堆积性能。

2026-10-20 实现[跨事件共享与标定约束双τ](pileup/experiments/alpha_gamma_shared_tau_20261020/README.md)：113条独立标定联合得到有效响应τ为91.74/364.82 ns，与固定95/351、逐事件约束/free双τ及完整模板GLS共同扫描50–1000 ns。开发冻结约束宽度后新增168个独立原生窗口；触发前20＋后150 ns时，带标定约束的Q±5%且类型正确为alpha 9/17、gamma 3/12，完整模板GLS为alpha 15/17、gamma 8/12。匹配均值＋白噪声对照另分开评价电荷和双τ精度；不能移作原生SiPM或材料寿命的普适最短窗。详见[本地Notebook第12节](jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#shared-tau-implementation-20261020)及[L38](docs/literature-use.md#shared-tau-implementation-20261020)。

离线 SiPM 数字化入口为 [sipm_response.py](pileup/analysis/sipm_response.py)，使用同一 [S13360-6050PE 配置](g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json)，当前示例、验证入口和笔记均采用该器件。用完整时间窗推进所有粒子共享的 pixel 状态，另存雪崩表、波形、种子、配置及源码记录；不修改输入 ROOT。命令、概率定义、实测单雪崩模板接口和限制见 [使用说明](docs/sipm-response.md)。历史无能量字段的 ROOT 仅在明确启用常数有效 PDE 近似时接受。

## ROOT 数据与运行记录

| 树 | 每行含义 | 主要字段 |
| --- | --- | --- |
| RunInfo | 一次运行的信息 | 时间窗和中央评价区 ps、事件数、实际种子与随机引擎状态、步长、宏路径与内容、各源启用状态及注入率。 |
| WaveformEvent | 一个完整时间窗 | 初级粒子数、总沉积能量 MeV、全部到达光子的时间 ps、位置 mm 与光子能量 eV。 |
| PrimaryPhoton | 一个初级粒子及其后代贡献，包含 α | 粒子/源标签、源能量和发射时间、沉积能量、闪烁光子数、到达光子的时间、位置与光子能量 eV。 |
| SiPMEvent（启用时） | 整个时间窗对应的一个 Geant4 digi | 雪崩时间 ns、微胞、类型、父编号、原光子索引、电荷 C、恢复系数；电流 A、电压 mV、ADC码与削顶标记。 |

`RunInfo.sipm_config_json` 保存实际完整参数；另记录模块版本、随机数引擎说明和独立基础种子，`SiPMEvent` 记录该事件派生后的两个实际种子。新增电信号的粒子标签不参与数字化计算。

ps 除以 1000 得到 ns。`edep_total_MeV` 是沉积能量真值，不是淬灭及光收集后的可见能量。模板形状、绝对光输出和能量标定须分别处理。

2026-10-10 起，主工程在两棵树中新增 `photon_energy_eV`（`vector<float>`）分支，记录进入虚拟 SiPM 时 pre-step 的光子能量。`WaveformEvent` 中它与 `all_photon_arrival_t_ps/x_mm/y_mm` 按索引对应，`PrimaryPhoton` 中与 `photon_arrival_t_ps/x_mm/y_mm` 对应；能量随同一光子一起按到达时间排序。真空波长可用 `wavelength_nm = 1239.842 / photon_energy_eV` 换算。历史 ROOT 和冻结实验没有该字段，分析时先检查分支是否存在；`primary_energy_MeV` 仍表示初级粒子能量。

新运行将发射与光子到达时间向下量化到 1 ps，记录范围为 `[0,T)`；旧版使用最近整数舍入。相对连续时间的量化误差小于 1 ps，与旧版整数记录最多相差 1 ps；物理输运和发光模型不变。`primary_time_ps` 是源发射时间，并非进入晶体或首次沉积的时刻。旧文件可能没有新增的 RunInfo 字段。

本地现有多源样本均只有一个 10 μs 窗；下表来自各自 ROOT 的 RunInfo 和 PrimaryPhoton，不由当前宏推断。完整文件名为 `<时间戳>.root`，这些 ROOT 文件不随 Git 仓库分发。

| 时间戳 | brems / ¹²C γ / α 注入率（Hz） | 初级粒子数 | 实际 brems / ¹²C γ / α 数 |
| --- | --- | ---: | --- |
| 20260624_16h28m34s | 10⁶ / 10 / 10⁶ | 19 | 10 / 0 / 9 |
| 20260624_16h28m48s | 10⁷ / 10² / 10⁷ | 206 | 102 / 0 / 104 |
| 20260624_16h29m17s | 10⁸ / 10³ / 10⁸ | 1993 | 977 / 0 / 1016 |

单窗示例不能给出可靠的重建效率或误差统计。后续主实验应使用多段独立随机流，处理窗前历史和窗末截断，评价 α 效率、漏检、误拆分/误合并、假 α 率、能量偏差/分辨率和整体能谱偏差。

每次运行记录代码提交号及未提交改动、完整宏、随机种子、事件数、源谱和注入率、时间窗、关键参数、依赖环境和输出路径。实际宏和 manifest 保存在对应提交目录，ROOT 与宏的索引记录在本地 `jupyter/Run_Catalog.ipynb`。原 `data.log` 已归档到 `jupyter/run_records/`，新代码不再生成这个重复日志；RunInfo 继续保存运行条件和实际宏内容。现有历史数据的几何、光学参数及环境仍有追溯缺口。


2026-10-20 实现[纯α主动拆峰和局部模型搜索](pileup/experiments/alpha_local_model_selection_20261020/README.md)：固定95/351 ns模板，未知K/t/Q，增加单→双拆分、增删/合并及多初值联合重拟合。独立标定/开发后冻结，在23个新原生窗口上与原基线同窗比较，20/40 MHz完整恢复分别121/291→150/291、90/526→140/526，10 MHz为90/136→88/136；40 MHz额外预测50→66。八档完整表、误差、失败和时间错位对照见[本地Notebook0003第17节](jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#pure-alpha-local-model-selection-20261020)及[采用记录L40](docs/literature-use.md#pure-alpha-local-model-selection-20261020)。低残差不能替代逐事件恢复；高率每档仅2窗，电荷尚未连接沉积能标。

2026-10-21 完成[相关GLS、局部种群搜索与轻量GRU的独立对照](pileup/experiments/alpha_correlated_global_20261020/README.md)：固定双τ和充分动态范围，三种新方法与旧SSE/局部K/CNN在23个冻结后新原生窗同窗比较。同一批23个新原生窗中，旧局部K、对角WLS、相关GLS、GLS＋局部DE、旧CNN和轻量GRU的完整恢复数依次为726、776、636、635、1038、982（共同目标4333个）。这是非均衡曝光总计，逐率结果及额外输出须同时阅读。 协方差白化仍不充分，已知K种群搜索诊断与未知K盲效率分开报告；不构成通用最高率或能量分辨率结论。[L42采用记录与验证边界](docs/literature-use.md#alpha-correlated-global-gru-20261021)对应本地Notebook0003第19节及固定图19-1至19-5，原输出保留。

<a id="alpha-ml-refresh-20261021"></a>
2026-10-21 按当前方法[重做纯α机器学习比较](pileup/experiments/alpha_ml_refresh_20261021/README.md)：四个冻结网络与当前WLS、旧混合及新CNN＋残差WLS统一评价，另生成16个冻结后新种子窗口，已有23窗作为公开共同复评分别报告。新确认当前WLS恢复419/2098，新混合为446/2098，额外输出112→133。完整恢复增加伴随额外输出增加，属于取舍，未达到两项同时改善；不据此自动替换默认方法。未重训网络、未改物理参数，保持固定双τ／充分动态范围；逐率、误差、负对照及八档完整波形见本地Notebook0004第15–21节和[L43采用记录](docs/literature-use.md#alpha-ml-refresh-20261021)。

<a id="alpha-component-ml-20261021"></a>
2026-10-21完成[分量监督与U-Net–GRU的新训练比较](pileup/experiments/alpha_component_ml_20261021/README.md)：四方案各两个训练种子，49训练/26验证ROOT，冻结后60个新窗口覆盖0.1–100 MHz。U-Net–GRU原目标完整恢复1124/1164/4253、额外输出462/455，相对同批新GRU小幅改善，尚未明显超越旧U-Net；两个结构、两个种子加入分量损失后完整恢复数均下降，当前不采用该损失。WLS为812/4253、额外300；逐率错误、62软件控制和oracle诊断一并记录，不能只按合并数排名。实现保留未知K密集解码，固定95/351 ns、完整模拟电压及充分动态范围。实际源码/配置/模型/宏、数值检查和8幅全窗图保留，本地Notebook0004第23–28节与Run_Catalog记录过程；高率仍每档仅两短窗，Q未作沉积能标。

2026-10-22完成[单脉冲校准与局部双分量ML](pileup/experiments/alpha_pair_curriculum_20261021/README.md)：三种单脉冲目标×两种子、无序/有序×直接/物理四种分量方案×两种子，共14次训练；24个冻结后新原生母ROOT。共同微胞控制盲拟合恢复179/267、额外35；有序物理头＋refit两种子分别为184/267和191/267、额外14和13。原生孤立、均值、线性叠加与联合光子重数字化分域报告，保留严格近邻时间和额外输出。原生随机触发裁帧与真时刻居中诊断共享24母ROOT。局部两槽模型尚不能直接替代完整20µs高率未知K研究；过程见本地Notebook0004第30–34节及[L47](docs/literature-use.md#alpha-pair-curriculum-20261021)。

2026-10-26完成[强→弱事件失败诊断](pileup/experiments/alpha_pair_weak_diagnostics_20261026/README.md)与[四种经验协方差权重对照](pileup/experiments/alpha_pair_empirical_covariance_20261026/README.md)。旧测试明确转为开发数据：high_weak、40 ns的原盲恢复0/40，固定真实时刻后13/40，真值启动联合GLS为8/40，精确减去真实强分量后40/40；后三项为真值辅助诊断。四种新权重均未同时满足恢复与误拆约束，未升级默认方法或生成新独立测试。完整结果记录在本地Notebook0004第53–57节、六幅固定编号图及[L53–L54](docs/literature-use.md#alpha-pair-weak-diagnosis-20261026)；保持固定双τ、充分动态范围及局部未知K≤2，不能据此给出普适最小间隔或最高计数率。

2026-10-26进一步实现[联合分量与归一化似然](pileup/experiments/alpha_pair_joint_likelihood_20261026/README.md)：用独立母事件分折核验条件均值和互相关，加入完整模板面积守恒、候选依赖协方差及logdet。两种新流程各分析1880帧；其中1600帧开发评价中，完整方案弱组合15点平均65.17%、较强组合94.17%，原观测GLS分别61.00%、93.17%，空／单控制通过，四档单事件Q偏差中位绝对值均小于0.07%。开发选中该方案，但high_weak的8/40 ns仍0/40、200 ns为22/40。旧数据收益不能当独立确认；过程与五幅固定编号图见本地Notebook0004第58–62节及[L55](docs/literature-use.md#alpha-pair-joint-likelihood-20261026)。

随后完成[320个新母事件的独立确认](pileup/experiments/alpha_pair_joint_confirmation_20261026/README.md)：冻结后1600帧中，弱组合平均由原观测GLS的60.50%提高至64.50%，母事件分块配对95区间为+1.00至+7.17个百分点；较强组合92.83%→93.33%，差值区间跨零。新方法400个空／单控制仅mid有1帧额外输出，单事件Q偏差中位绝对值均小于0.08%。收益主要来自weak_high、8 ns的0/40→27/40；high_weak、40 ns仍0/40，200 ns反由23/40降至20/40，不能视为全面改善。Notebook0004第63–65节和[L56](docs/literature-use.md#alpha-pair-joint-confirmation-20261026)保留完整指标、三个固定图、有限搜索及条件取舍；320个完成ROOT和4次中断尝试已登记本地Run_Catalog。

<a id="alpha-pair-weak-search-20261026"></a>
2026-10-26完成[弱分量起点与局部剖面对照](pileup/experiments/alpha_pair_weak_search_20261026/README.md)：两种新搜索各3480帧，原开发与已公开确认批均作为本轮开发资料。最新批弱组合平均由原NLL的64.50%变为增加原份额起点的64.33%、补充弱份额起点的64.17%；先强后弱8/40/200 ns均仍为0/0/20，各40对。两种方法的空／单控制通过，但均未通过恢复率约束，不升级方法或新增母ROOT。真值辅助剖面、候选排序与门槛诊断说明降低似然代价不等于恢复正确参数，尚不能推断物理分辨极限。9600条评分重放通过，六幅编号图与完整波形记录在本地Notebook0004第66–70节及[L57](docs/literature-use.md#alpha-pair-weak-search-20261026)。

2026-10-26完成[原生与匹配高斯模型对照](pileup/experiments/alpha_pair_model_closure_20261026/README.md)：760帧已公开原生数据与相同Q/t的1520帧高斯抽样，另280帧标定门槛。冻结算法、共用原生门槛时，先强后弱8/40/200 ns分别为原生0/0/20（每点40对）、高斯0/0/42（每点80对）；匹配模型仍重现短间隔失败。改用高斯门槛未增加双事件成功，却使原生强单事件额外输出0/80→7/80。40 ns高斯即使固定真实时刻也仅29/80严格恢复，仍属真值辅助诊断，不能推断所有算法的分辨极限。5项采样测试、4560条盲评分及2880条辅助评分核验通过，不升级默认方案。五幅固定编号图、全窗波形与原生非高斯残差记录在本地Notebook0004第71–76节及[L58](docs/literature-use.md#alpha-pair-model-closure-20261026)。

2026-10-27完成[恢复阶段与局部信息诊断](pileup/experiments/alpha_pair_information_20261027/README.md)：拆解1080帧已有双事件结果，计算444个真实／名义参数点，并对2160条已有辅助拟合误差作比较。先强后弱40 ns时，弱lnQ的局部尺度中位为0.478；固定真实时刻为0.474，固定强事件Q/t为0.049，强弱Q局部相关约−0.995，提示电荷分配歧义。上述值基于已知K和真参数，不是实际相对电荷标准差或盲恢复保证；漏检与定位损失也仍存在。6项数值测试和全部参数点的差分检查通过，原评分与输入保持一致。本地Notebook0004第77–81节保存五幅固定编号图及误差分布，详见[L59](docs/literature-use.md#alpha-pair-information-20261027)。

<a id="alpha-pair-charge-allocation-20261027"></a>
2026-10-27完成[总电荷／分配坐标与共同剖面搜索对照](pileup/experiments/alpha_pair_charge_allocation_20261027/README.md)：四种变体在1800帧上比较，保持同一物理似然与电荷域，以相同280个原生单事件分别校准增峰门槛。原生先强后弱40 ns严格恢复由基线0/40变为allocation_profile的0/40，200 ns由20/40变为21/40。四种变体均未同时通过预定原生开发约束，不选择新默认流程。七项数值测试、13680条评分和3600条已知K诊断核验通过；五幅固定编号图、全长分量波形及完整控制记录在本地Notebook0004第82–86节与[L60](docs/literature-use.md#alpha-pair-charge-allocation-20261027)。

<a id="alpha-pair-risk-training-20261027"></a>
2026-10-27完成[连续Q监督学习与物理拟合](pileup/experiments/alpha_pair_risk_training_v2_20261027/README.md)：9600训练、960验证，比较绝对／对数／相对Q三损失；按验证选择`relative`用于两种混合拟合，冻结后另生成1000个连续Q高斯评价样本。带学习型Q惩罚的NLL相对原物理基线，严格双事件成功为新高斯500/600→465/600、已公开原生迁移241/360→198/360；完整六流程及分条件结果保留。通过预定开发要求的方案：无；不自动替换默认方法。没有新光学母事件，不能当物理独立确认。6项数值检查、10560主评分与5760已知K诊断核验通过；Notebook0004第87–92节、六幅固定编号图及[L61](docs/literature-use.md#alpha-pair-risk-training-20261027)记录训练、控制和全窗分量。

<a id="alpha-pair-physical-gate-20261027"></a>
2026-10-27完成[物理事件数判决与网络起点对照](pileup/experiments/alpha_pair_physical_gate_20261027/README.md)：在相同缓存拟合候选上分离事件数判决：原生先强后弱200 ns的学习型门槛拒绝13个严格正确K2，改用物理似然判决后严格恢复7/40→20/40；40 ns仍0/40。追加局部精修／网络起点＋物理判决，原生均241/360，高斯分别500/600、501/600。门槛修复通过开发控制，但网络相对局部精修的增益未达两域各5个百分点要求；保持原默认流程。 全部数据本轮均为已公开开发复用，没有新训练／拟合／ROOT。6项数值检查、2040帧代价重放、7040主评分／8800含敏感性评分和2880已知K2诊断核验；Notebook0004第93–97节保存五幅固定图、完整波形与候选拒绝去向，依据及边界见[L62](docs/literature-use.md#alpha-pair-physical-gate-20261027)。

<a id="alpha-pair-evidence-profile-20261027"></a>
2026-10-27完成[事件数门槛与条件电荷剖面](pileup/experiments/alpha_pair_evidence_profile_20261027/README.md)：固定四流程候选，扫描127个可取秩门槛；两域十个空／单控制组均≤5%额外输出时，四流程原生40 ns先强后弱仍0/40。另计算160条已知K2的第二槽Q2剖面及单独真实Q2评分探针，报告范围宽度、分离区段、边界与数值搜索遗漏，不把ΔJ层级称为置信水平。6项数值测试、7040条默认评分精确重放及全量剖面代价检查通过。全部数据为开发复用，不改默认算法；Notebook0004第98–102节、五幅固定图及[L63](docs/literature-use.md#alpha-pair-evidence-profile-20261027)保存完整结果。

<a id="alpha-pair-profile-stability-20261027"></a>
2026-10-27完成[条件电荷剖面的数值稳定性审核](pileup/experiments/alpha_pair_profile_stability_20261027/README.md)：四组各系统抽取8条，共32条已公开开发波形，对比原曲线、增加观测起点、全域中点加密和相容边缘再加密。最后加密在1／4／9三层级同时通过预定端点／区段／参考代价筛查的有30/32；0条真实Q2评分探针将最终参考再降低>0.01。结果为已知K2的有限数值诊断，不是盲效率改善或覆盖率校准；原160条剖面及默认流程保留。6项数值检查、32条四阶段曲线和独立标记的评分探针核验完成；本地Notebook0004第103–107节、五幅固定图及[L64](docs/literature-use.md#alpha-pair-profile-stability-20261027)保留逐条变化和完整波形。

<a id="alpha-pair-contour-coverage-20261027"></a>
2026-10-27完成[边界细化与新条件高斯重复抽样](pileup/experiments/alpha_pair_contour_coverage_20261027/README.md)：5条定向数值诊断中，5/5通过新增的边界括区与更强局部重拟合检查。随后冻结方案，生成32条新条件高斯双脉冲和12条空／单控制；32/32条双脉冲通过数值筛查，全部32条都保留在主要包含率分母中。40ns组ΔJ≤4范围包含真实Q2为14/16，已知K2严格整对恢复为3/16。200ns组ΔJ≤4范围包含真实Q2为16/16，已知K2严格整对恢复为10/16。这检验的是冻结模型内、固定参数点的重复抽样表现；每组16条不足以认证95%覆盖，更不能推广为原生SiPM的盲恢复或普适最短间隔。 本地Notebook0004第108–113节保存固定图号、全部区段与完整波形；方法条件和文献采用见[L65](docs/literature-use.md#alpha-pair-contour-coverage-20261027)。原默认盲算法保持。

<a id="alpha-pair-selection-20261027"></a>
2026-10-27完成[事件选择与条件高斯门槛确认](pileup/experiments/alpha_pair_selection_20261027/README.md)：固定拟合器，256条新单脉冲校准后另生成448条独立空／单／双确认波形；同时分解旧32条区间的选择分母。40ns：同为0.025门槛参数，改用新单脉冲校准后，输出两个事件由17/64变为33/64，严格整对恢复由0/64变为1/64。即使直接给这批真实双脉冲选K=2，当前候选中严格成功也只有10/64；这是同一有限候选的诊断上限，不是物理可辨识极限。

40ns被选中候选的有序Q2误差中位由+61.1%变为+39.5%，这是选择群体改变，不是同一条拟合参数变准；包含时间未匹配候选，不称能量分辨率。时间未匹配的额外输出由7增为17个，多检出没有转化成相应的准确恢复。额外输出表示10ns内找不到匹配真事件，不必然表示整窗原本无信号。

200ns：同为0.025门槛参数，改用新单脉冲校准后，输出两个事件由63/64变为63/64，严格整对恢复由46/64变为46/64。即使直接给这批真实双脉冲选K=2，当前候选中严格成功也只有46/64；这是同一有限候选的诊断上限，不是物理可辨识极限。

四档单脉冲共256条，原／新主门槛误报双事件分别为1/1条；含任意时间未匹配输出的帧分别为1/1条。各Q档结果另列，不能把这个等权混合比例当作任意输入谱的误报率。 本轮未计算新的区间覆盖，默认算法保持；本地Notebook0004第114–119节及[L66](docs/literature-use.md#alpha-pair-selection-20261027)记录条件、全部指标和固定编号完整波形。


### 2026-10-28：弱脉冲Q/t与均值估计方程诊断

40ns的严格整对恢复，原NLL／扩大搜索NLL／固定一次权重GLS／迭代GLS分别为2/16、2/16、2/16、2/16。所有失败与未收敛样本均保留在分母中。

同一组在固定真实时刻／固定真实电荷／从全部真实参数启动后自由拟合时，严格成功分别为4/16、7/16、6/16；这三项使用额外真值信息，只是误差来源诊断。

200ns的严格整对恢复，原NLL／扩大搜索NLL／固定一次权重GLS／迭代GLS分别为10/16、10/16、11/16、11/16。所有失败与未收敛样本均保留在分母中。

同一组在固定真实时刻／固定真实电荷／从全部真实参数启动后自由拟合时，严格成功分别为13/16、14/16、13/16；这三项使用额外真值信息，只是误差来源诊断。

扩大搜索使4/32条的J降低超过0.05；从真值出发仍有4/32条找到比扩大搜索更低超过0.05的J。有限起点和有限预算不能证明全局最优；J更低也不等于Q和t更接近真值。0.05只用于描述数值差异，不是统计显著性阈值。

共有7/32条出现：真值启动后的自由拟合严格成功，扩大搜索候选失败，但成功解的J反而高出0.05以上。其中40ns为4条、200ns为3条。这些反例说明同一代价的排名与本次逐事件误差并不一致；增加搜索起点无法保证选中更准确的那个解，也不能据此判定物理上一定不可分辨。

没有方法同时通过预先冻结的四项推进条件，因此本轮不增加112条确认样本、不改变默认重建流程。这是44条已暴露模型样本上的开发筛查结论，不是物理不可分辨极限，也没有证明其他优化器无效。

结果为已知K、暴露过的44条条件高斯开发样本；未作为盲效率或物理分辨极限。笔记0004第120–125节含六幅固定编号图，旧570个单元及输出保留。源码及模型快照：`pileup/experiments/alpha_pair_mean_score_20261028/`；本地结果与运行记录：`jupyter/artifacts/alpha_pair_mean_score_20261028/`、`jupyter/run_records/alpha_pair_mean_score_20261028/`。本轮不改变默认流程，无新ROOT、训练、提交或推送。


### 2026-10-28：盲重建的联合参数可靠性独立验证

用校准集中盲选为双事件的36条真实双脉冲，经验95%上分位得到联合代价阈值c=11.864；40/200ns分别贡献6/30条。这是当前选门和设计下的混合样本标定，不是统一95%置信保证。

40ns：原盲流程输出双事件14/64，严格整对恢复0/64；选出双事件后的联合区域包含真实Q1/Q2/t1/t2为13/14，按全部真实双脉冲计算则为13/64。时间漏检57个，额外输出7个。

200ns：原盲流程输出双事件63/64，严格整对恢复38/64；选出双事件后的联合区域包含真实Q1/Q2/t1/t2为58/63，按全部真实双脉冲计算则为58/64。时间漏检15个，额外输出14个。

找到达到或超出容差边界的替代解：79条盲选双事件，其中真实双脉冲77条、控制样本误报双事件2条；严格正确为38/79。分母为0时准确率未定义。

有限搜索未找到替代解：0条盲选双事件，其中真实双脉冲0条、控制样本误报双事件0条；严格正确为0/0。分母为0时准确率未定义。

质量标记没有改变原拟合或增加正确恢复事件；找到替代参数说明该代价区域内存在歧义，有限搜索未找到不证明全局唯一或物理可辨识。

校准144条、独立确认288条；原NLL和事件数门不变，使用完整四参数代价区域及有限容差外搜索。模型条件、固定分层及数值限制详见 [实验记录](pileup/experiments/alpha_pair_blind_uncertainty_v2_20261028/README.md)。本地0004笔记第127–132节保留原600单元，并给出八幅固定图号的新结果。没有新ROOT、训练或默认算法替换。


### 2026-10-28：保持事件数判决的优化器独立消融

独立288条确认采用同一波形的配对比较，原事件数门和K判决保持不变；主方案为预先固定的多起点后Powell，另外两种仅作过程消融。

40ns：整对严格恢复0/64→0/64；改善0、退化0，配对精确检验p=1。时间漏检50→51，额外输出4→5。

200ns：整对严格恢复36/64→34/64；改善1、退化3，配对精确检验p=0.625。时间漏检18→20，额外输出15→17。

80条被盲选为K2的确认波形中，20条代价降低超过0.05。控制组原误报双事件1/160，K误报数保持不变是设计约束，不能作为新筛选器有效的证据。

较低代价与逐事件准确度须分别解释；有限搜索未找到更优解不证明物理不可辨识。本轮不替换默认算法，也不外推为原生SiPM、α/γ、free-τ或高计数率性能。

| 间隔/ns | 原流程 | 宽域重启 | 多起点 | 再加Powell（主方案） | 改善/退化 |
| --- | --- | --- | --- | --- | --- |
| 40 | 0/64 | 0/64 | 0/64 | 0/64 | 0/0 |
| 200 | 36/64 | 36/64 | 34/64 | 34/64 | 1/3 |

先检查15条已暴露开发波形，再冻结方案生成288条新确认波形；比较原流程、宽域重启、多起点和Powell。完整记录见[实验说明](pileup/experiments/alpha_pair_optimizer_polish_20261028/README.md)；本地0004笔记第133–137节新增5幅固定编号图，原627个单元与输出保留。全部结果仍限固定双τ的纯α条件高斯模型。

### 2026-10-28：当前 SiPM 电荷到沉积能量的独立验证

完成[独立能标与三间隔验证](pileup/experiments/sipm_energy_calibration_v2_20261028/README.md)：280/280个母事件分别拟合和选择四种能标，320个不重叠的既有确认母事件提供1600帧。盲推断只将冻结的拟合电荷换算成沉积能量，未调整原拟合或选门；该批波形已用于电荷/定时评价，不称作新未见输运样本。

本轮能标未通过预设单事件控制。所选光产额比例能标在源能量1.625、2.5、5.486、8 MeV四组上，真实分量电荷换算的沉积能量偏差中位数为+13.52%、+5.91%、−1.57%、−5.95%。选择集显示：考虑饱和的能标偏差较小，但高能端正常电荷涨落超出硬反演域；拒绝惩罚因此选中了有系统偏差的比例能标。需补目标能区外的标定保护点和中间能点，再独立确认。

在当前未通过控制的能标下，含弱α的15个组合/间隔中，能量±10%严格整对恢复均值为模型GLS 21.33%、联合NLL 24.83%；配对差+3.50个百分点，母事件块95%区间[+0.67,+6.50]。这是能标与拟合误差共同作用的诊断，不能作为最终能谱性能，也不与原电荷±20%端点直接比较。固定图号的条件图、绝对计数谱、漏检/额外输出表保存在本地`codex_0007_SiPM_Energy_Calibration.ipynb`；协议及限制见[文献采用L72](docs/literature-use.md#sipm-energy-transfer-20261028)。8项数值测试、母事件划分和3200条旧盲输出一致性核验通过；没有新ROOT或默认算法替换。


### 2026-10-29：能标保护点与冻结后新单事件确认

[扩充能标实验](pileup/experiments/sipm_energy_guardband_20261029/README.md)完成496个有效原生单α ROOT，包括16先导、96拟合、64选模和320冻结后新确认；另保留3份未纳入分析的尝试ROOT，按原种子重试未增加独立样本数。保护源能点1.5/9MeV拓宽实际沉积能量支持，确认含两档未参与标定的1.8/7MeV。选择`empirical_pchip`，十组单事件控制通过，最大组能量偏差中位数绝对值3.48%，有限换算320/320。这里用已知分量Q验证能标，不能作为盲事件检出率；32例/组的通过是点估计判据。

历史盲拟合仅更新能量换算，含弱α的能量±10%严格整对恢复为模型GLS 31.67%、联合NLL 38.17%；原Q/t、检出数和门保持不变。受控条件图与绝对计数谱已更新，仍需检查共享微胞前史、连续高率和模型失配。过程见本地codex_0007第8–14节及[采用记录L73](docs/literature-use.md#sipm-energy-guardband-20261029)，Run_Catalog新增实际496个ROOT/宏/种子映射。
