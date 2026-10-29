# 本项目文献采用记录

本文件记录 GAGG(Ce) α 能谱与堆积重建项目实际采用的文献结论、证据位置、适用条件和待验证事项。推导、计算与图像继续保存在 [研究计划](../jupyter/Research_Plan_20260921.ipynb) 和 [仿真与分析笔记](../jupyter/Note.ipynb)，在具体采用条目中链接相关章节。

## 知识库入口与使用顺序

共享知识库根目录为 `/media/sf_works/papers`，从本项目只读访问。以下链接按本项目与 `papers/` 为同级目录组织。

1. 从 [topic-index.md](../../papers/topic-index.md) 定位主题，并从 [review-index.md](../../papers/review-index.md) 查询已有审核意见。
2. 优先读取与当前问题相关的主题页：
   - [数字脉冲处理与堆积恢复](../../papers/methods/数字脉冲处理与堆积恢复.md)
   - [脉冲形状甄别与粒子鉴别](../../papers/methods/脉冲形状甄别与粒子鉴别.md)
   - [闪烁体的光产额与能量分辨率](../../papers/concepts/闪烁体的光产额与能量分辨率.md)
3. 沿主题页链接读取相关 `sources/` 精读笔记；已有记录足以回答的问题直接复用。
4. 关键公式、参数或引用证据有缺口时，按具体问题定点回查原文，保留 DOI、页码、公式、图表等定位信息。

## 采用记录

### 2026-09-30：纯 α 随机堆积第一版方法与评价口径

**研究问题**：沿用当前 p–¹¹B 有效单 α 谱，如何从总光子到达波形中估计未知的脉冲数、时间、光输出及沉积能谱，并比较四档注入率？已完成的 400 窗只证明数据准备和一致性检查完成，尚无重建性能结果。

**本次阅读范围**：从主题与审核索引进入三个优先主题页，再复用下列来源笔记的相关方法、结果、限制和用户意见；核对了 Morad 第 553 项实际人工意见。本轮没有重新阅读全文 PDF、运行论文代码或验证论文性能。旧 PLAN 第 9 节的“只读公开摘录”属于 2026-09-22 的历史阅读范围，现在补入的是共享库已有的全文精读记录，不冒充本轮原文复核。

**本项目决策**：先核验已有 α 响应和光输出标定，再做纯 α 的候选发现、非负幅度/时间联合拟合与残差增删候选。先在 10⁶ Hz 独立测试窗闭环，10⁵ Hz 作低率对照，再比较 10⁷、10⁸ Hz。算法结构是本项目待验证方案；不直接移植其他探测器的性能数字。执行顺序见 [PLAN 第 12 节](../jupyter/Research_Plan_20260921.ipynb#alpha-next-20260930)。

#### L01：Zhang 2024——模板误差与隐藏脉冲

- **来源位置**：[数字脉冲主题 B23-052](../../papers/methods/数字脉冲处理与堆积恢复.md)；[来源笔记第 3、4、6 节](../../papers/sources/2024_Zhang_NIMA_Development-of-a-pile-up-pulse-recovery-algorithm-for-the-LaBr3-detector.md)；原文 PDF2–5，式 2–5、图 2–7；[DOI](https://doi.org/10.1016/j.nima.2024.169273)。
- **文献条件**：2 英寸 × 2 英寸 LaBr₃(Ce)、R6233 PMT、500 MHz/14 bit、32 μs 离线记录；模拟采用与拟合同一模板、417/1042/1667 keV 三离散能量、高斯噪声 σ=7 ADC 道（约 1.5 keV）、指数间隔、10⁶–10⁸ s⁻¹。实验快中子区域不超过 4×10⁶ s⁻¹。
- **结论与采用**：固定模板失配可诱发假脉冲；作者按残差尝试增峰并重拟合。本项目先查 α 形状随能量、对齐及尾部截断的变化，再决定是否需要分能量模板；增峰阈值和停止规则用独立验证集校准。幅度依赖、非负约束及光子计数的噪声模型均需按本项目定义，不称逐式复现。
- **验证限制**：论文约 10 ns 间隔及建议不超过 6×10⁷ s⁻¹来自特定模拟；实验“恢复超过 98%”借助模拟效率推算。均不能作为 GAGG 的分辨间隔、最大计数率或已测效率。式 5 修正系数坐标及增峰标准仍有缺口。
- **审核/项目状态**：来源为“AI 精读草稿，待人工审核”；本项目仅采用方法方向，尚未实现或验证。

#### L02：Adrich 2021——峰数不等于脉冲数

- **来源位置**：[数字脉冲主题 B23-044](../../papers/methods/数字脉冲处理与堆积恢复.md)；[来源笔记第 3、4、6 节](../../papers/sources/2021_Adrich_NIMA_Least-squares-fitting-algorithm-for-peak-pile-up-correction-in-gamma-ray-spectroscopy.md)；原文 PDF2–5，式 1–2、表 1、图 3–6；[DOI](https://doi.org/10.1016/j.nima.2020.164962)。
- **文献条件**：1 英寸 × 1 英寸 LaBr₃:Ce、R6231-100 PMT、有源底座、CAEN DT5730，¹³⁷Cs γ；约 10⁵ 个低率单脉冲建模板，高率实测约 80 kcps、400 万事例。当前复用笔记未列采样率；平滑宽度、噪声/协方差和 p 值阈值未完整给出。另将实测单脉冲两两叠加，保存能量与偏移作真值。
- **结论与采用**：一个观测峰也可能包含多个脉冲；单模板拟合失败后允许增加分量。项目拟合的 K 为未知数，不固定为 2，也不等于 R×T；用独立单脉冲检验假拆分，用少量双脉冲辅助定位失败。
- **验证限制**：表 1 的 99.68% 是已识别双候选的接受结果；92.3(10)% 是求和峰变化估计；8.5 ns 约 50%、17 ns 超过 90% 是两能量均通过 ±3σ 的合成测试。分母、成功定义和样本不同，不能互换。高拟合优度不等于正确恢复真实脉冲数。
- **审核/项目状态**：“AI 精读草稿，待人工审核”；采用候选数与真值数分离的原则，项目阈值、匹配和能量检验尚待实现。

#### L03：Lopatin 2012——发现候选与估计能量分开

- **来源位置**：[数字脉冲主题 B23-024](../../papers/methods/数字脉冲处理与堆积恢复.md)；[来源笔记第 3 节](../../papers/sources/2012_Lopatin_IEEE_Pileup-attenuation-for-spectroscopic-signals-using-a-sparse-reconstruction.md)；原文 PDF2–4，式 1–6、算法 1。
- **文献条件**：HPGe、¹³⁷Cs/²⁴¹Am 混合源；模型使用单位积分字典、Poisson 到达和独立高斯噪声，实际采样率、源强、计数率及噪声幅度未给出完整信息。证据主要为实测能谱的定性改善，仍有未消除的求和结构。
- **结论与采用**：论文使用较强正则化寻找起点、较弱正则化估计能量。若项目后续使用稀疏候选，拟合能量时另做无 L1 收缩的非负幅度重拟合，并检验偏差；这是项目调整，不是论文原算法。
- **验证限制**：原式 4 未约束非负，式 5 对同一项列出相互冲突的阈值条件。保留疑点，不自行修成“作者算法”。目前仅列为候选发现的后备方法，不与第一版模板法同时铺开。
- **审核/项目状态**：“AI 精读草稿，待人工审核”；未复现，未选择为当前必需模块。

#### L04：GAGG 响应——形状模板与能量标定分别核验

- **来源位置**：[PSD 主题 B23-021/B23-031](../../papers/methods/脉冲形状甄别与粒子鉴别.md)、[光产额主题 B23-031](../../papers/concepts/闪烁体的光产额与能量分辨率.md)；[Kobayashi 2012 第 3–4 节](../../papers/sources/2012_Kobayashi_NIMA_Significantly-different-pulse-shapes-for-gamma-and-alpha-rays-in-Gd3Al2Ga3O12-Ce3%2B-scintillating-crystals.md)，PDF2–3、表 2，[DOI](https://doi.org/10.1016/j.nima.2012.07.055)；[Tamagawa 2015 第 3 节](../../papers/sources/2015_Tamagawa_NIMA_Alpha-gamma-pulse-shape-discrimination-in-Gd3Al2Ga3O12-GAGG-Ce3%2B-crystal-scintillator-using-shape-indicator.md)，PDF2–3、图 2–3，[DOI](https://doi.org/10.1016/j.nima.2015.05.052)。
- **文献条件**：Kobayashi 的 5 mm GAGG 立方体/H2256 PMT、5 GHz 采样、512 波形平均、0–1000 ns 衰减拟合，比较 0.662 MeV γ 与约 5.47 MeV α；另用 1 cm 立方体/双 PMT、净 1 μs 积分测相对光量。Tamagawa 为 1×1×3 cm³ GAGG/R6233K PMT、8 bit/500 MHz、2048 ns，外部约 5.48 MeV α、γ 标能范围 0.66–1.84 MeV。未提供可移植到本项目的噪声模型或高率 pileup 性能。
- **结论与采用**：平均形状差异不是逐事件堆积效率；有限积分、表面收光及 α 淬灭影响能标，γ 积分线性也不证明 α 动能线性。旧 α 模板的首光子对齐、1 μs 截断及逐事件归一化须核验；Qα(E_dep) 需保留绝对光输出与响应展宽，不能仅由单位面积模板得到。
- **项目依据**：[Note 第 4.2 节](../jupyter/Note.ipynb#single-alpha-gamma-waveforms) 及 [模板复算记录](../jupyter/artifacts/note_single_particle_4c4ltqti/run_record.json)。现模型 α 为共享 95/351 ns、光量权重 0.30/0.70；Kobayashi 表 2 α 为 145/764 ns、0.83/0.17。这是已知模型差异，本轮不替换冻结实验参数；先验证现模型下的算法，再另建物理失配对照。
- **审核/项目状态**：两篇来源均为“AI 精读草稿，待人工审核”，状态未改变。2026-09-30 已完成第一批长尾、能量分组与条件能标检查，见下面的项目验证；真实 GAGG 标定、未知脉冲数重建与物理失配仍待验证。

#### L05：Morad 2024——把漏检和条件性能分开报告

- **来源位置**：[数字脉冲主题中的 Morad 2024 条目](../../papers/methods/数字脉冲处理与堆积恢复.md)；[来源笔记第 3–5 节](../../papers/sources/2024_Morad_TNS_Model-Based-Deep-Learning-Algorithm-for-Detection-and-Classification-at-High-Event-Rates.md)；原文 PDF2–3、7、10 表 II；[DOI](https://doi.org/10.1109/TNS.2024.3371573)。
- **文献条件**：25 mm × 25 mm CLYC、PMT 800 V、1 GS/s、17 mV 触发；实测单脉冲库和模拟库，Poisson 到达及线性叠加构造高率 n/γ 信号，模拟加 AWGN。文中“0 dB”采用其本底噪声幅度定义，不是通用信号功率定义。6 Mcps 是构造信号的测试条件，不是持续在线 DAQ 证明。
- **结论与采用**：第 553 项确认的 6 Mcps 结果是条件分类 89.1%、含漏检的总准确率 72.9%、另列虚警 14.1%。项目只借鉴分母分离：匹配成功子集的能量误差必须与全体目标 α 的恢复效率一起报告，假脉冲数/有效时间另列；不把这些数值作为 α 重建性能，也不据此启动深度学习。
- **审核状态**：**“第553项的限定采用意见已获人工确认，整篇其余内容仍待审核”**；已读取 [实际人工意见](../../papers/global-review.md#human-review-553-morad2024-confirmed)，确认范围没有扩大。论文虚警分母、训练/测试独立性、公式及其他问题仍保留。
- **项目状态**：评价原则纳入 PLAN；α 的匹配规则、信号定义、误拆分/误合并统计与置信区间尚待实现。

<a id="sipm-use-20261010"></a>
### 2026-10-10：SiPM 光子响应与数字化接口（通用方案）

**问题与决定**：用户尚未选型。采用“光子到达 → 概率触发/微胞恢复 → 单雪崩电荷与波形 → 电子学/ADC”的分层。完整公式、字段和边界探测注意事项见 [PLAN 第 13 节](../jupyter/Research_Plan_20260921.ipynb#sipm-response-20261010)。2026-10-10 用户授权实现后已建立离线模型与假设参数验证，实际器件尚未标定，冻结实验未改动。

**实现与验证（2026-10-10）**：建议的演变见 [PLAN 第13.9节](../jupyter/Research_Plan_20260921.ipynb#sipm-implementation-plan-20261010)，完成记录见 [第13.10节](../jupyter/Research_Plan_20260921.ipynb#sipm-implementation-results-20261010)。[离线实现](../pileup/analysis/sipm_response.py) 使用共享微胞队列，保留线性基线，支持二维 PDE 表、恢复及相关噪声；读出支持参数化电流或实测电压模板。软件/合成极限检查和已有单窗后处理不等同于实物验证。单候选串扰/AP分支、其 q/q0 缩放、Gamma 触发延迟和 ADC 孔径是明确标注的项目子模型，未声称直接复现某篇论文的完整代码。配置中所有数值均标为假设，概率定义和适用限制见 [说明](sipm-response.md)。

**证据范围**：只读检索共享库主题与审核索引、主题页和三篇来源笔记；定点回查下列 PDF 页。三篇均保留“AI 精读草稿，待人工审核”，未找到针对这三篇的后续人工确认。没有重新全文精读或递归阅读其引用；另核对官方 Geant4 11.4 文档与本机 11.4.2 源码。

#### L06：Acerbi 与 Gundacker 2019——PDE、单雪崩电荷与响应叠加

- **来源**：*Understanding and simulating SiPMs*，NIMA 926, 16–35；[DOI](https://doi.org/10.1016/j.nima.2018.11.118)；[笔记](../../papers/sources/2019_Acerbi_NIMA_Understanding-and-simulating-SiPMs.md)；[原文](../../papers/raw/2019_Acerbi_NIMA_Understanding-and-simulating-SiPMs.pdf)。本轮核 PDF4 式3、PDF5 式5、PDF7 恢复/电路、PDF13 式18、PDF15 式19–24；PDF15另看图像。
- **采用**：PDE 含触发概率和微胞填充因子，雪崩电荷约为过压乘电容；低占用可叠加单次雪崩响应。本项目拟用面积归一双指数核记账电荷；原式19按峰高归一，该转换是项目建模选择。
- **条件与缺口**：综述的 TOF-PET 前沿计时示例在特定低暗率下省略部分长尾等效应，不能据此忽略高率 GAGG 的恢复、基线及相关噪声。核和噪声参数须由选定器件/读出链标定；当前95/351 ns是已包含在到达时间中的闪烁参数，不作为 SiPM 常数。

- **流程补核（2026-10-10）**：本次实际查看 PDF16 图26（刊页31），核实其包含 PDE/LTE 筛选、SPTR、串扰/暗事件、单胞信号叠加、低通与阈值定时；适用范围仍是指定 TOF-PET 示例。此前的 SiPMDigitizer 是项目拟建模块名，尚无已验证的同名实现；后续新增文献已定位G4SiPM专门方法论文，见L09；流程依据与项目接口分界见 [PLAN 第13.8节](../jupyter/Research_Plan_20260921.ipynb#sipm-literature-flow-20261010)。

#### L07：Jaime Rosado 2019——恢复同时影响触发概率和电荷

- **来源**：*Modeling the Nonlinear Response of Silicon Photomultipliers*，IEEE Sensors Journal 19(24), 12031–12039；[DOI](https://doi.org/10.1109/JSEN.2019.2938018)；[笔记](../../papers/sources/2019_Rosado_IEEESensorsJ_Modeling-the-Nonlinear-Response-of-Silicon-Photomultipliers.md)；[原文](../../papers/raw/2019_Rosado_IEEESensorsJ_Modeling-the-Nonlinear-Response-of-Silicon-Photomultipliers.pdf)。本轮核 PDF1 范围/式1–3、PDF2 式4–6与种子定义；PDF2另看图像。
- **采用**：逐微胞按上次触发后的实际过压计算 PDE 和电荷；所有粒子的光子共享同一器件历史。短脉冲占用式只作极限校验，不能对长 GAGG 脉冲逐 bin 独立施加硬饱和。
- **条件与缺口**：据既有笔记，实验为 S13360-1325CS/1350CS、3×3×20 mm³ LYSO、400 ns积分、约300–2100 keV及连续 LED；采样率/全部测试温度未完整列出。平均电荷式含两种子等近似，不替代强非线性逐事件波形/方差；恢复时间、PDE参数和相关噪声不移植成项目默认值。

#### L08：Ogawara 与 Ishikawa 2016——光学输运与实测读出响应分层

- **来源**：*Signal pulse emulation for scintillation detectors using Geant4 Monte Carlo with light tracking simulation*，RSI 87, 075114；[DOI](https://doi.org/10.1063/1.4959186)；[笔记](../../papers/sources/2016_Ogawara_RSI_Signal-pulse-emulation-for-scintillation-detectors-using-Geant4-Monte-Carlo-with-light-tracking-simulation.md)；[原文](../../papers/raw/2016_Ogawara_RSI_Signal-pulse-emulation-for-scintillation-detectors-using-Geant4-Monte-Carlo-with-light-tracking-simulation.pdf)。本轮核 PDF4 图1/方法/仪器、PDF6 §III.B。
- **采用及条件**：借鉴 Geant4 到达时间与读出响应卷积分层。原装置为 R6231-100 PMT、5 GS/s、600 MHz；BC408 提供包含有限闪烁和读出链的有效响应，并非 SiPM SPE。指定晶体平均662 keV脉冲验证不能作为当前高率 pileup 性能。

**工程核查**：调研时的历史 ROOT 只有光子时间和位置，支持假设常数有效 PDE 的首版后处理。2026-10-10 随后已在主工程两棵数据树增加 photon_energy_eV，记录光子进入 SiPM 时的能量，与时间/位置同步排序；新输出可据此查光谱 PDE，历史数据仍缺此字段。构建及同种子小样本回归通过，见 [PLAN 第13.6节](../jupyter/Research_Plan_20260921.ipynb#photon-energy-20261010)。真实封装反射/返回路径也不能靠离线缩放补出。Geant4 EFFICIENCY 是边界吸收后的探测概率；Digitize 需要显式调用。依据：[官方边界文档](https://geant4.web.cern.ch/documentation/dev/bfad_html/ForApplicationDevelopers/TrackingAndPhysics/physicsProcess.html#boundary-process)、[数字化文档](https://geant4.web.cern.ch/documentation/dev/bfad_html/ForApplicationDevelopers/Detector/digitization.html)及本机 G4OpBoundaryProcess 的 DoAbsorption/InvokeSD、G4DigiManager 接口。加入器件响应后重新核验模板和光输出能标；旧 η=0.78348 不作 PDE。

<a id="sipm-p10-use-20261010"></a>
### 2026-10-10：新增 P10 SiPM 文献与已实现模型的对照

**范围与状态**：新增18篇题录/导读已筛查，重点回读以下4篇完整来源笔记并定点核PDF；其余仅作导航。SiPM审核检索及作者/题名补查未定位这4篇后续人工意见，仍保留AI草稿待人工审核。完整公式、证据表、代码对照与后续优先级见 [PLAN第13.11节](../jupyter/Research_Plan_20260921.ipynb#sipm-p10-20261010)。库内文件只读；本次更新文献记录，不修改器件模型、假设参数或既有仿真结果。

#### L09：Niggemann等2015——G4SiPM的逐胞时间队列

- **来源**：*G4SiPM: A novel silicon photomultiplier simulation package for Geant4*；[DOI 10.1016/j.nima.2015.01.067](https://doi.org/10.1016/j.nima.2015.01.067)；[笔记](../../papers/sources/2015_Niggemann_NIMA_G4SiPM-A-novel-silicon-photomultiplier-simulation-package-for-Geant4.md)；[原文](../../papers/raw/2015_Niggemann_NIMA_G4SiPM-A-novel-silicon-photomultiplier-simulation-package-for-Geant4.pdf)
- **核验/采用**：回查PDF1–3，另看PDF2图1/式3–9原页。直接支持光子输运后建立时间—微胞队列、加入暗事件、递归相关触发和恢复。四邻居Binomial串扰、双分量Poisson后脉冲是候选扩展；当前单子候选模型并非该论文算法的完整复现。
- **条件/限制**：作者的5%符合限于S10362-33-050、3600胞、485 nm、12 ns光脉冲和250 ns积分；过压/温度与全部参数未列全。式4的PDE参考面、式5–6的概率语义及触发恢复假设不能无条件移植。未获取/运行作者包，当前兼容性未核。

#### L10：Garutti、Klanner、Rolph、Schwandt 2021——低占用Geiger Array与单雪崩电流

- **来源**：*Simulation of the response of SiPMs; Part I: Without saturation effects*；[DOI 10.1016/j.nima.2021.165853](https://doi.org/10.1016/j.nima.2021.165853)；[笔记](../../papers/sources/2021_Garutti_NIMA_Simulation-of-the-response-of-SiPMs%3B-Part-I-Without-saturation-effects.md)；[原文](../../papers/raw/2021_Garutti_NIMA_Simulation-of-the-response-of-SiPMs%3B-Part-I-Without-saturation-effects.pdf)
- **核验/采用**：PDF2–4，PDF4原页核式1–3/AP段。采用“放电数组→单位面积单胞核→积分/电路卷积”的分层及门前历史要求。其快慢衰减核是指数之和；现有双指数差核是另一现象学选择，未替换代码。
- **条件/限制**：低占用，不自行计算PDE、绝对增益或电路推导响应；输入主放电不能与到面光子混用。AP的0/1候选也是作者明确的简化。KETEK谱对照条件与证据层级见PLAN，未把其定性比较当作当前GAGG高率验证。

#### L11：Klanner 2024——饱和、概率/电荷恢复与公共负载

- **来源**：*Simulation of the response of SiPMs Part II: With saturation effects*；[DOI 10.1016/j.nima.2023.169018](https://doi.org/10.1016/j.nima.2023.169018)；[笔记](../../papers/sources/2024_Klanner_NIMA_Simulation-of-the-response-of-SiPMs-Part-II-With-saturation-effects.md)；[原文](../../papers/raw/2024_Klanner_NIMA_Simulation-of-the-response-of-SiPMs-Part-II-With-saturation-effects.pdf)
- **核验/采用**：PDF3–6，另看PDF4式3–6与接受步骤。采用高占用必须保留恢复历史和检查负载压降的物理依据；继续保留当前x/y逐胞模拟，未改成论文的均匀单胞平均场近似。
- **条件/限制**：平均响应验证不能自动支持方差、时间分辨或高DCR性能；据笔记PDF7–11仍有亚ns/12 ns区失配及高DCR验证缺口，本轮未重查这些结果页。式5恢复基线和AP条件化存在需源码澄清之处，不照抄有歧义的式子。当前跨阻输出与ADC削顶不包含公共偏置反馈。

#### L12：Corsi等2007——单胞电荷与前端电路

- **来源**：*Modelling a silicon photomultiplier (SiPM) as a signal source for optimum front-end design*；[DOI 10.1016/j.nima.2006.10.219](https://doi.org/10.1016/j.nima.2006.10.219)；[笔记](../../papers/sources/2007_Corsi_NIMA_Modelling-a-silicon-photomultiplier-%28SiPM%29-as-a-signal-source-for-optimum-front-end-design.md)；[原文](../../papers/raw/2007_Corsi_NIMA_Modelling-a-silicon-photomultiplier-%28SiPM%29-as-a-signal-source-for-optimum-front-end-design.pdf)
- **核验/采用**：PDF2–3文本核电荷关系与前放描述，图形细节沿用完整笔记。单胞电荷约为过压乘有效电容，端口脉形取决于器件电容、淬灭和前端阻抗/带宽。优先以同读出链SPE或可验证电路约束波形。
- **条件/限制**：单胞暗脉冲电路拟合，不包含完整探测概率或高占用噪声；未移植旧器件参数、未使用有疑点的反演式。

**本次检索时的项目判断**：保留Geant4→ROOT→离线逐胞雪崩→读出的接口。优先器件/SPE/双脉冲/暗态标定，再比较独立prompt/delayed串扰、多个后脉冲时间分量与增益涨落，检查前史长度及高率负载。各项仍是待验证的模型选择，未因新增文献改写已有35项软件检查的意义。记录核查结果见 [清单](../jupyter/run_records/sipm_literature_p10_20261010.json)。

<a id="sipm-native-use-20261010"></a>
### 2026-10-10：原生 Geant4 数字化的采用记录

用户随后选择从官方 `G4VDigitizerModule` / `G4DigiManager` 实现，先检查单α和单γ电信号。现已增加 [SiPMDigitizer](../g4/src/SiPMDigitizer.cc) 和C++数值核心，在完整窗光子排序后显式调用，保存当前G4Event的digi集合并写入 `SiPMEvent`；Python离线入口保留用于参数扫描与独立读出核对。此实现采用L06/L07的PDE—恢复—电荷关系、L09的统一时序队列方向和L10的雪崩数组/电路分层。L09四邻居多候选、双分量后脉冲以及L11公共偏置反馈尚未实现；参数仍为本项目假设，文献AI待审核状态不变。

官方机制依据：[Geant4 digitization](https://geant4.web.cern.ch/documentation/dev/bfad_html/ForApplicationDevelopers/Detector/digitization.html)及本机11.4.2头文件。未引入第三方G4SiPM，亦不将项目 `SiPMDigitizer` 称为Geant4提供的完整SiPM物理模型。单粒子条件、配置和验证入口见 [响应说明](sipm-response.md#native-sipm)，完整改动与实际结果见独立笔记 [codex_0001_SiPM_Digitization.ipynb](../jupyter/codex_0001_SiPM_Digitization.ipynb)，研究路线见 [PLAN第13.12节](../jupyter/Research_Plan_20260921.ipynb#sipm-native-20261010)。这些软件/合成核验不替代真实器件或高率实验验证。

<a id="hamamatsu-s13360-20261013"></a>
## L13：滨松 S13360-6050PE 的数据表配置（2026-10-13）

用户要求将假设 SiPM 换成一个具体滨松型号；本次选择单通道 **S13360-6050PE，25°C、3 V 过压**。最新入口为 [v2 配置](../g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json)，完整依据、对照结果和限制写在 [codex_0001_SiPM_Digitization.ipynb 器件参数节](../jupyter/codex_0001_SiPM_Digitization.ipynb#hamamatsu-s13360-20261013)。`datasheet_based` 表示典型参数与注明的近似，不能称为实物标定。

- **厂家原始资料**：Hamamatsu Photonics, *MPPC S13360 series*，文件名 `s13360_series_kapd1052e.pdf`，Cat.No.KAPD1052E07，Jul.2025；无DOI。[官方产品页](https://www.hamamatsu.com.cn/cn/zh-cn/product/optical-sensors/mppc/mppc_mppc-array/S13360-6050PE.html)与[官方PDF](https://www.hamamatsu.com.cn/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s13360_series_kapd1052e.pdf)。项目下载副本、SHA256和页面图保存在本地 `jupyter/run_records/codex_hamamatsu_s13360_sources/`，未写入共享文献库。核对PDF文本及2/11页型号/几何；实际查看3–5页原页图像，核参数表、PE光谱曲线KAPDB0322EB及过压曲线KAPDB0324EA。
- **采用的典型输入**：有效面6×6 mm、50 μm pitch、14400 pixel、增益1.7×10⁶、DCR 2 Mcps（0.5 p.e.阈值）、串扰3%；峰值PDE 40%指定450 nm，不直接作为GAGG全光谱效率。PDE不含串扰/AP，已含74%填充因子，不重复乘该因子。端电容1280 pF、击穿电压53±5 V、工作电压温度系数54 mV/°C和环氧折射率1.55仅作为资料记录；当前没有完整封装、电路或温度反馈模型。
- **近似读图与插值**：只配置400–650 nm、0–3 V，覆盖当前GAGG发光表。3 V光谱点按原图近似读取，说明性读图精度约±1个百分点；2 V按450 nm处0.315/0.400缩放，2 V以下线性连零。波长/过压可分离与低过压行为属项目近似，不是厂家完整二维实测表。当前α有效面光子的满恢复加权PDE为0.34943，属于项目计算。
- **恢复时间代理**：Jaime Rosado, *Performance of SiPMs in the nonlinear region*, NIMA912(2018)39–42，DOI [10.1016/j.nima.2017.10.031](https://doi.org/10.1016/j.nima.2017.10.031)。文件名 `2018_Rosado_NIMA_Performance-of-SiPMs-in-the-nonlinear-region.pdf`；[共享笔记](../../papers/sources/2018_Rosado_NIMA_Performance-of-SiPMs-in-the-nonlinear-region.md)。本轮实际查看原PDF第3页表1，仅采用同系列50 μm的1350CS/3050CS行所列29 ns作代理。该表转引[1]，本轮未重新核查[1]；它不测6050PE，亦未提供可移植温度/误差。AI精读草稿仍待人工审核，未扩大已有用户采用范围；表中25 μm行的pixel数疑点不用于当前配置。
- **剩余未知项**：AP和SPTR暂不附加，零参数表示未建模；低恢复过压PDE、噪声触发曲线、29 ns代理、SPE/前端传递函数及真实光耦合需标定。50 Ω、2/20 ns核、采样/ADC均为项目演示读出参数。拟合恢复代理取20/29/40 ns时，同一α光子历史的离线电荷相对29 ns变化+2.47%/0/−2.76%，是单例敏感性，不是实物误差范围。

本轮运行保留11个新ROOT：500 mV器件诊断6窗、1000 mV量程复跑5窗；每个ROOT仍仅1个20 μs窗。v1的α峰值508.018 mV、4点削顶，v2无削顶；同种子量程对照的光子、雪崩和模拟波形一致。源能/沉积能量仍分别记录，γ展示沿用已知有沉积种子，不能用于效率估计。实现内部一致性与真实器件性能验证分开记录。

### 2026-10-13：当前入口统一为 S13360-6050PE

按用户要求，当前配置目录移除40000 pixel假设器件，笔记主体、原生与离线示例、C++测试配置和器件分析统一使用S13360-6050PE。器件数值沿用已保存的v2配置；不新增未经标定的AP/SPTR或封装参数。原配置、旧笔记与图像归档，历史ROOT、RunInfo和冻结源码不改写。

本次通过 `paper-search` 只读复用共享库索引、Rosado2018完整精读笔记及本节已有选型记录；29 ns继续标为1350CS/3050CS代理，AI草稿待人工审核状态保留，未重新精读PDF或扩展引用证据。新图与分析只读取S13360-6050PE的两批记录；旧器件对照退出当前入口。变更与验证见本地 [统一记录](../jupyter/run_records/sipm_device_only_20261013_143939/validation.json)。

<a id="psd-fit-use-20261013"></a>
### 2026-10-13：α / γ 波形拟合函数的文献核对

用户希望用波形拟合作GAGG + S13360-6050PE粒子鉴别。本轮按 `paper-search` 只读检索本地库，复用来源笔记并定点核原PDF；完整函数、参数表、书目信息/DOI、原文链接、条件和推导已集中迁入 [独立拟合笔记第2节](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#psd-fit-functions-20261013)，核验清单及PDF哈希在 [文献核验记录](../jupyter/run_records/psd_fit_literature_20261013_163822/sources.json)。文献阶段仅核对函数；随后按用户要求完成现有两条独立模拟波形的拟合与模板判别示例，见下方项目计算更新。未重新运行Geant4或更改物理输入。

#### L14：GAGG双指数——分别保留α和γ的衰减常数

- **来源与证据**：延续L04的 [Kobayashi2012笔记](../../papers/sources/2012_Kobayashi_NIMA_Significantly-different-pulse-shapes-for-gamma-and-alpha-rays-in-Gd3Al2Ga3O12-Ce3%2B-scintillating-crystals.md)（DOI 10.1016/j.nima.2012.07.055）与 [Tamagawa2015笔记](../../papers/sources/2015_Tamagawa_NIMA_Alpha-gamma-pulse-shape-discrimination-in-Gd3Al2Ga3O12-GAGG-Ce3%2B-crystal-scintillator-using-shape-indicator.md)（DOI 10.1016/j.nima.2015.05.052）；本轮核两篇PDF2方法及PDF3原图，分别为图1/表2和图3(b)，DOI核各PDF1。
- **可用结论**：衰减段可用 `B + Af exp(-u/tauf) + As exp(-u/taus)`。Kobayashi γ为95.4/351 ns、α145/764 ns；Tamagawa γ91.2/367 ns、α138/599 ns。幅度比例不等于面积比例，后者须乘各自寿命；基线不计入光量。纯衰减段不能在两幅度均自由时独立辨识起点。
- **适用条件与缺口**：前者5 mm立方/H2256 PMT、512平均、1000 ns拟合；后者1×1×3 cm³/R6233K PMT、500 MHz/8 bit、2048 ns；约5.47–5.48 MeV外α与662 keV γ。两者均非S13360读出，无可移植的参数误差/噪声协方差。Tamagawa逐事件鉴别实际用SI，不是逐事件自由双指数拟合。
- **项目采用状态**：作为候选尾部函数和参数初值/失配对照。当前α/γ共享95/351 ns、α快光量0.30的模型仍保留；这一输入与上述实验参数组有已知差别。两篇AI草稿待人工审核，本轮没有新增人工认可或实物验证。

#### L15：SiPM后电压——发光与器件/前端响应卷积

- **来源与证据**：[Wang2019完整笔记](../../papers/sources/2019_Wang_JINST_Improvement-on-pulse-shape-discrimination-performance-of-SiPM-array-coupled-Cs2LiYCl6Ce3%2B-detector.md)，DOI 10.1088/1748-0221/14/11/P11005，本轮实看PDF6式3.2；以及L10的 [Garutti2021完整笔记](../../papers/sources/2021_Garutti_NIMA_Simulation-of-the-response-of-SiPMs%3B-Part-I-Without-saturation-effects.md)，DOI 10.1016/j.nima.2021.165853，本轮实看PDF4式1–3及§2.11。
- **条件与用途**：Wang为CLYC/SensL 64元阵列/OPA657、500 MS/s、12 bit、20 μs，热中子约3.2 MeVee及宇宙线/γ代理；Garutti为低占据SiPM模型，按单雪崩核叠加并与电子响应卷积。它们提供响应分层依据，不提供本项目GAGG α/γ实测常数。
- **项目推导**：双指数发光经单RC核后的每项为 `[exp(-u/tau)-exp(-u/taue)]/(tau-taue)`，按面积归一；完整候选保留光学、SPE、前端及采样孔径。`taue`不等于29 ns恢复代理。当前2/20 ns电流核及2 ns低通应显式卷积或用校准电压模板，不能直接冒充单RC实物模型。
- **边界与状态**：Wang原印R0C0约10.50 μs与“far less than 1 μs”措辞冲突，保留原文疑点且不迁移该数值。有限微胞/恢复非线性可能改变形状，ADC无削顶不代表SiPM线性。两篇AI草稿待人工审核；函数归一、卷积代数和现有模拟样例拟合已检查，未完成器件标定或完整非线性响应反演。

#### L16：两类模板分别拟合——以拟合优度差作PSD候选

- **来源与证据**：S. Marrone等/n_TOF Collaboration，*Pulse shape analysis of liquid scintillators for neutron studies*，NIMA490(2002)299–307；[完整笔记](../../papers/sources/2002_Marrone_NIMA_Pulse-shape-analysis-of-liquid-scintillators-for-neutron-studies.md) / [PDF](../../papers/raw/2002_Marrone_NIMA_Pulse-shape-analysis-of-liquid-scintillators-for-neutron-studies.pdf)。本轮核PDF3/4的式1/2和参数固定流程、PDF6χ²差及PDF7图5；实看3/6/7。DOI原PDF未印，本轮未核实，不由PII猜补。
- **原方法与条件**：NE213/XP2041和C6D6/EMI9823QKA，1 GS/s、8 bit、250 MHz；60Co/AmBe、约200 keVee硬阈值。先标定形状并约束参数，再对同一事件用两类固定模板分别拟合幅度、比较归一χ²。属于液闪n/γ先例，不转用时间常数、阈值或性能数字。
- **项目候选**：在同一采样窗、相同自由参数/噪声规则下，α/γ各拟合面积、观测起点及基线，定义 `D = chi2_gamma - chi2_alpha`；正值倾向α。阈值由独立验证样本确定；协方差随假设变化时需采用完整似然。当前展示的一对波形、首光子真值对齐和不同光量不足以证明分类性能。
- **原文疑点与状态**：式1慢项符号与式2正慢尾不一致，PDF7图5轴差值方向与PDF6正文相反；本项目使用明确的新定义，不宣称逐式复现。来源AI草稿待人工审核，无本函数新增人工采用意见。后续先标定孤立脉冲，再按重建光量分箱报告γ接受率、α误留率及ROC。

<a id="psd-fit-results-20261013"></a>
### 项目计算更新：2026-10-13 孤立 α / γ 波形拟合

以下是当前模拟的计算结果，不是上述文献或真实器件的性能结论。过程、拟合式、静态参数表、残差图和执行输出在 [codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#fit-inputs)；可复用实现为 [alpha_gamma_waveform_fit.py](../analysis/alpha_gamma_waveform_fit.py)，结果与来源信息见 [summary.json](../jupyter/artifacts/alpha_gamma_fit_20261013_165844/summary.json)。

- 采用无ADC削顶的S13360-6050PE ROOT `20261013_09h58m53s`（α）和 `09h59m09s`（γ），各一个已知孤立脉冲。源能5.486/0.661657 MeV，晶体沉积能4.939493/0.202202 MeV；γ是已选定的部分沉积有信号事件。参数从RunInfo与原运行冻结源码读取，未用当前宏推断历史条件。重复种子和仅改变电子噪声的记录不计为独立样本。
- SiPM前对有效面光子到达直方图作bin积分Poisson尾部拟合；SiPM后显式卷积2/20 ns单雪崩核、2 ns低通和0.5 ns采样孔径，拟合analog与ADC。光子诊断使用首光子对齐；电压拟合仅由观测波形估计起点/基线，不传入事件粒子标签、真值时刻或沉积能。
- 固定模拟的95/351 ns寿命，光子侧慢面积占比α/γ为0.70697/0.19464，ADC侧为0.74017/0.20474；ADC残差RMS/峰高分别0.889%/0.833%，单指数对照为2.566%/1.144%。1–2.5 μs拟合窗下固定寿命慢占比变化较小。自由寿命拟合可降低残差，但γ的ADC有效慢常数约217 ns，不能当作351 ns内禀输入的准确恢复。
- 分别用模拟输入定义的α/γ形状拟合同一ADC记录，两例的较小SSE假设都与评价标签一致。当前使用 `D_SSE = SSE_gamma - SSE_alpha`，并非已校准χ²；相邻样点残差相关系数约0.963/0.941，未建立完整协方差。两条不同光量的样例不支持鉴别率、FOM或ROC，也不能据此证明实物适用性。
- 独立状态方程与解析卷积、无噪声参数回收、32次电压优化的多初值收敛、量化与哈希/重复事件检查已通过，见 [数值验证记录](../jupyter/run_records/alpha_gamma_fit_20261013_165844/validation_analysis.json)。后续需要独立校准/验证/测试事件、共同光量区间和实物SPE/前端标定；文献审核状态保持不变。

- **短窗口追加检查**：将ADC拟合终点扩展到阈值后100、200、300、500、750 ns，并保留1–2.5 μs对照，起点仍为阈值前200 ns。固定95/351 ns时，500 ns的慢占比α/γ为0.7365/0.2012，接近2 μs的0.7402/0.2047；100–300 ns自由寿命出现大幅变化、参数触界或分量权重趋零。两类固定模板在这两个样例的全部窗口仍偏向对应标签，不构成最短窗口或鉴别效率验证。100 ns短窗拟合的α自由模型在0–2 μs共同评价区间RMS/峰高约8.48%，显示窗内残差小不保证长尾外推准确。9档窗口共72次拟合已核对样点、残差和2 μs回归，见 [拟合笔记第5.2节](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#short-window-scan-20261013)及 [追加记录](../jupyter/run_records/alpha_gamma_short_windows_20261013_205215/validation_analysis.json)。仍需独立样本决定窗口；已有长记录定位和预触发基线需求没有消失。

<a id="alpha-rate-fitting-20261013"></a>
### 项目计算更新：2026-10-13 纯 α / S13360 计数率与单脉冲恢复

研究与全部图表见新笔记 [codex_0003_Alpha_Pileup_Rate_Fitting.ipynb](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb)。本轮只读复用主题/审核索引、数字脉冲处理、PSD、光产额主题及 L01/L02/L03 来源笔记，并核对 Morad 第553项的实际意见；没有新增论文下载或原文逐页复核。实际采用 L01 的残差发现遗漏成分、L02 的脉冲数与峰数分离、L05 的条件误差与总体效率分母分离。当前逆滤波候选、非负幅度/浮点时间重拟合和停止规则是本项目方案，不称复现作者未完整给出的阈值/χ²检验；L03 的 LASSO 没有作为本次算法实现。

- **数据与模型**：新冻结实验 [alpha_s13360_rate_fit_20261013](../pileup/experiments/alpha_s13360_rate_fit_20261013/README.md)，当前 p–¹¹B 有效单 α 谱、无 γ/其他物理源、S13360-6050PE 原生波长 PDE/共享微胞恢复/串扰/暗计数及 ADC。10⁵/10⁶/10⁷/10⁸ Hz 分别16/4/3/1个独立20 μs窗，中心门[5,15) μs；另1个开发窗排除于结果。原有400窗及保留测试集未用于本次性能分析。SiPM只将计算候选保护1000万→8000万，物理/读出参数未改；实际宏、RunInfo、源码、种子、环境和哈希全部保存。
- **真值隔离与拟合函数**：拟合器只接收总电压、采样时刻、可观察ADC端点码和固定α模板，不接收注入率、真实K/时刻/能量。模板由前次独立孤立α的ADC结果固定95/351 ns、有效慢面积占比0.740174，并卷积现有读出。当前G4发光权重仍0.30/0.70，L04记录的真实GAGG模型差异未消除。实际雪崩按光子祖先归属仅用于评分，保留堆积中的微胞状态，不把独立重数字化的电荷当作本次真实贡献。
- **实际结果**：目标为中心源发射且电荷≥500满恢复雪崩等效数；一对一时间容差±10 ns，恢复还要求电荷误差≤20%。四档分别恢复10/10、35/36、142/204、0/639。10⁷ Hz 时间匹配177/204、漏检27、匹配后电荷未通过35、额外候选13，近邻误合并标志8；零误拆分标志仅限定义的±10 ns邻域。10⁸ Hz这一个中心窗全部ADC到顶，当前整窗拒绝策略使恢复数为0；未削顶模拟analog另作诊断也只通过23/639，不替代ADC性能。模拟analog峰值约2281 mV，当前ADC上限1000 mV。
- **条件误差与限制**：前三档已匹配事件的相对电荷误差中位数约+0.79%、+0.59%、+3.27%，时间误差中位数约−1.73、−0.98、−0.78 ns（相对首个有效面光子，不是源时间）。未校准沉积能量，不能称能量分辨率；电荷/光子数随高率下降含恢复及串扰，也不是厂家PDE测量。小样本探索、已知α模板、固定基线/形状、阈值选择及ADC拒绝政策均限制外推，未给出装置最大计数率或最终测试集性能。
- **验证与开发记录**：构建/SiPM CTest、25份实际RunInfo/宏和ROOT哈希、光子祖先唯一映射、电荷守恒及264个简化控制检查通过，见 [验证记录](../jupyter/run_records/alpha_rate_fit_20261013_211541/validation_analysis.json)。模板+电子噪声双脉冲在本次4–8 ns合并、12–160 ns的每点8个试验通过，不含光子/微胞涨落，不作为实物12 ns能力。最初10⁷ Hz的261–275个候选触发拟合器250保护值，已保存原结果并提高到500；这项数值保护失败不解释为物理不可分。

下一步为独立单脉冲/分能量模板与基线/相关噪声标定、更多独立窗口、增大量程或降低跨阻的新配置对照，并研究削顶区外的局部恢复和非线性电荷反演。论文来源的“AI草稿待审核”及 Morad 的限定确认范围均未改变。

<a id="literature-rate-limits-20261014"></a>
### L17：2026-10-14——拟合解堆积的计数率与恢复质量

**检索与核验范围**：本轮按 `paper-search` 只读检索本地共享库，核对既有精读笔记与审核记录，并定点查原PDF关键文字/图表。完整题录、证据表、装置/噪声条件和待核事项集中在 [纯α堆积笔记第11节](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#literature-rate-limits-20261014)，不重复维护完整报告。以下来源仍为AI草稿待人工审核；`global-review.md` B23归并未新增人工确认。本轮不联网、不修改共享库，不运行论文代码或新增仿真。

- [Zhang2024](../../papers/sources/2024_Zhang_NIMA_Development-of-a-pile-up-pulse-recovery-algorithm-for-the-LaBr3-detector.md)，[DOI](https://doi.org/10.1016/j.nima.2024.169273)：原PDF4–6、图6–7/9/11。同模板、低噪声LaBr₃模拟建议≤6×10⁷ s⁻¹，但本轮粗读图7(a)，10⁷/6×10⁷/10⁸ s⁻¹全能峰效率仅约90%/30%/8%；不是本项目逐事件电荷/时间成功率。实测快中子区≤4×10⁶ s⁻¹的>98%来自模拟推算；γ-flash约6.67×10⁷ s⁻¹由恢复间隔拟合获得，缺同率独立效率真值。
- [Khilkevitch2020](../../papers/sources/2020_Khilkevitch_NIMA_Advanced-algorithms-for-signal-processing-scintillation-gamma-ray-detectors-at-high-counting-rates.md)，[DOI](https://doi.org/10.1016/j.nima.2020.164309)：原PDF8/10/11，§4.5、表2、图17。LaBr₃+PMT实测5.1×10⁶ s⁻¹时fitting的4.44 MeV峰FWHM 0.098 MeV；模拟3×10⁷ s⁻¹时拟合补权后的总率仍偏差约33%。统计补权不恢复每个事件的时间和能量，离线处理不等于在线吞吐。
- [Löher2012](../../papers/sources/2012_Loher_NIMA_High-count-rate-gamma-ray-spectroscopy-with-LaBr3-Ce-scintillation-detectors.md)，[DOI](https://doi.org/10.1016/j.nima.2012.05.051)：原PDF5–6、图7–8/结论；模板矩阵解叠在约0.5–10 Mcps明显改善部分峰计数，最高14 Mcps仍有峰宽恶化，不能称未知脉冲数的全波形拟合或全事件恢复。
- [Shevelev2016](../../papers/sources/2016_Shevelev_NIMA_High-performance-gamma-ray-spectrometer-for-runaway-electron-studies-on-the-FT-2-tokamak.md)，[DOI](https://doi.org/10.1016/j.nima.2016.05.075)：原PDF7、表1/结论；固定形状拟合的FT-2谱仪结论至10⁷ s⁻¹，约2×10⁷ s⁻¹已观察到谱失真。

**本项目采用**：把“能看到恢复后的能谱峰”“给定阈值的逐事件恢复率”“能量偏差/分辨率”和“在线吞吐”分开定义；不移植上述LaBr₃+PMT数值作为GAGG+SiPM α上限。当前95/351 ns模板、微胞恢复和ADC量程，以及源注入率/探测事件率的差别均保留。研究主线、物理参数及拟合器不变；新增记录只说明证据与评价口径。

<a id="pulse-duration-rate-20261014"></a>
### L18：2026-10-14——用信号长度解释计数率与拟合难度

**项目计算与采用口径**：在 [纯α堆积笔记第12节](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#pulse-duration-rate-20261014) 中，沿用独立单α拟合的95/351 ns模板和冻结读出，计算电压FWHM 181.5 ns、10%–90%上升时间25.81 ns、90%/99%面积时长0.728/1.536 μs；这是当前平均模拟模板的结果，不是材料或实物脉宽测量。明确采用泊松占用量 `μ=R_sig W` 和下一邻居重叠概率 `1−exp(−μ)`，另以 `1−exp(−2μ)` 表示两侧任一邻居；`R_sig`为信号事件率，与源注入率及算法恢复率分开。重叠概率不作为拟合失败率，宽度不作为死时间。

**方法条件与实际结果**：已知恰好两个脉冲及其时间、固定基线/正确模板、白噪声时，幅度标准差相对孤立脉冲放大 `1/sqrt(1−ρ²)`；20/50/100 ns间隔对应3.67/1.89/1.35倍，不是盲分辨效率。本次只对已保存盲拟合作真值评价：10⁷ Hz注入组204个阈上目标全有±1.536 μs内有电荷α邻居，仍恢复142个；最近邻<12 ns为3/23、≥50 ns为83/99。3个独立窗的事后分组受幅度、多个邻居和共享微胞状态混杂，未给出通用分辨时间或最大计数率。原漏检/假脉冲、条件电荷误差及10⁸ Hz削顶限制仍保留，未新增能量标定或G4生产。

**文献来源与迁移限制**：复用 [Zhang2024](../../papers/sources/2024_Zhang_NIMA_Development-of-a-pile-up-pulse-recovery-algorithm-for-the-LaBr3-detector.md)，定点查原PDF2正文和PDF3图2：LaBr₃+R6233 PMT的FWHM主带约20个2 ns样点，即约40 ns，且随幅度变化。当前模板/40 ns的FWHM比约4.54，仅用于相同信号率下的占用量比较；不将文献60 Mcps按此比例换算为GAGG工作上限。来源仍为AI草稿待人工审核，本项目解析推导和数值验证不构成对文献的人工批准。

**补核原文讨论**：见 [第12.7节](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#duration-literature-direct-20261014)。[Löher2012](../../papers/sources/2012_Loher_NIMA_High-count-rate-gamma-ray-spectroscopy-with-LaBr3-Ce-scintillation-detectors.md) PDF4 §3.3/图5直接扫描间隔与恢复误差，PDF5式(5)/图8使用`nτ`及60/200 ns时间尺度联系计数率；[Khilkevitch2020](../../papers/sources/2020_Khilkevitch_NIMA_Advanced-algorithms-for-signal-processing-scintillation-gamma-ray-detectors-at-high-counting-rates.md) PDF2明确PMT、ADC带宽和电缆使实测脉冲长于材料衰减；Zhang2024 PDF2/4区分峰/尾堆积并用式(6)指数间隔和图6解释近邻漏失；Shevelev2016 PDF2以旧BGO约320 ns的高率堆积说明更换快LaBr₃的动机。这些支持时间尺度与率联合评价，不提供GAGG通用上限；当前模板数值、近邻分组和ρ误差公式仍明确为本项目结果/推导。保留Löher式(5)对P(0)的文字解释与公式不一致这一限制，审核状态不变。

**验证与待办**：宽度网格/尾部收敛、50万组泊松间隔及每种条件4000次孤立/已知时间双脉冲白噪声最小二乘检查已通过；新Notebook单元已执行，原31单元及输出保留。见 [结果](../jupyter/artifacts/pulse_duration_rate_20261014_105829/summary.json) 和 [验证记录](../jupyter/run_records/pulse_duration_rate_20261014_105829/validation.json)。后续需在独立样本上分别改变慢尾/前沿，明确固定电荷或固定峰高、采样和噪声，并检验幅度比、相关噪声、微胞恢复和ADC量程；当前单一宽度数据不足以测出 `R_max(W)`。

<a id="intermediate-rates-20261014"></a>
### 项目计算更新：2026-10-14 中间源注入率补充

沿用L17/L18的率、效率与时间尺度区分，新增20/40/60/80 MHz各3个独立20 μs窗口，中心门每档共30 μs；物理与0–1000 mV ADC配置不变。实际ADC逐事件恢复为20 MHz：51/443、40 MHz：0/858、60 MHz：0/1223、80 MHz：0/1640；中心ADC端点样点比例依次0.81%、24.82%、65.80%、96.71%。任何端点触及都会按原政策整窗拒绝，故其效率下降不能单独解释为拟合的物理极限。未削顶模拟电压仅对拒绝窗作诊断，不当作实际ADC结果。三窗/档及单一脉宽仍不足以定义最大工作率或`R_max(W)`。

新环境的两份同种子检查物理分支完全一致；14份新ROOT（12独立+2对照）的实际宏/种子/哈希及电荷守恒通过。高密度光子时间+x碰撞改用精确y/能量完成评分归属，原冻结导出器保留；盲拟合不接收真值、源率或已知脉冲数。完整结果、失败分解和输入限制见[纯α笔记第13节](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#intermediate-rates-20261014)及[验证记录](../jupyter/run_records/alpha_intermediate_rates_20261014_114000/validation_results.json)。本轮没有新增文献结论或修改共享库审核状态，后续优先区分量程、局部削顶恢复和固定宽度下的拟合失配。

<a id="full-signal-rates-20261014"></a>
### 项目计算更新：2026-10-14 动态范围足够的完整电压基准

**用户指定的当前假设**：采集动态范围足够，完整信号可被记录；本轮不以ADC位数、量化或削顶限制堆积恢复。统一用已有36个独立窗口的含噪 `analog_mV` 重新盲拟合，保留原SiPM状态、采样/电子学响应、固定模板和评分标准；不传入ADC端点掩码、真值或源率。旧ADC结果和“触顶即整窗拒绝”仅作历史记录，先前提出的量程扫描/局部削顶恢复不作为本轮前置工作。

**实际计算**：八档36窗均完成拟合；新增四档为20 MHz：152/443（34.31%）、40 MHz：140/858（16.32%）、60 MHz：119/1223（9.73%）、80 MHz：68/1640（4.15%）。本次20 MHz包含全部3窗，消除旧analog仅分析拒绝窗的条件选择。其余率、失败分解、条件电荷误差及固定脉宽占用量集中于[纯α笔记第14节](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#full-signal-rates-20261014)。没有新G4运行；原ROOT/宏/冻结实验保留。

**验证与边界**：4项完整电压控制、36窗输入与评分检查通过，见[验证记录](../jupyter/run_records/alpha_full_signal_20261014_143000/validation_results.json)。这是假设动态范围足够且忽略量化的模拟研究，未验证某款ADC；输出电荷尚未转为沉积能谱，小样本与单一宽度不足以定义最大工作率。下一步聚焦近邻间隔、幅度比、模板失配及受控脉宽扫描。本条是项目假设与计算更新，没有新增文献结论，共享库审核状态不变。

<a id="ml-pileup-literature-20261014"></a>
### L19：2026-10-14——机器学习／人工智能的pileup候选方法

**检索与核验**：按`paper-search`只读共享库，优先主题／审核记录并全文补查，选15篇重点来源。各篇题名、作者、DOI、原PDF／笔记链接、页码、训练／测试条件和限制集中于[纯α笔记第15节](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#ml-pileup-literature-20261014)及[证据表](../jupyter/run_records/ml_pileup_literature_20261014_153426/evidence.json)。复用已有完整精读，15个DOI与原PDF首页匹配；定点原页文字／图像核对范围逐篇记录，没有再次全文精读所有论文。共享库保持只读，未联网或训练模型。

**当前采用范围**：将SAEFit分量监督参数回归、1D卷积自编码器和U-Net事件序列展开作为后续候选；Morad的未知到达数／时间检出可参考，Jeon的候选峰高回归需保留检出依赖。Huang2025 Sensors的统计谱／真实率先验与TIM的时序能量分类分开；Lin的识别／计数、Fu／Liu的尾部PSD收益不等同逐事件时间＋电荷恢复。Ma2019已知4／9重数、Kafaee／Wang2022已知双脉冲只用于受控基准。没有据此替换当前算法、移植器件参数或认定GAGG高率性能。

**审核与研究条件**：第545、547、549、553、555项局部意见按原对象与范围保留；其余AI笔记仍待审核，不因本项目引用获得新认可。沿用动态范围足够的完整电压口径，后续同条件比较未知K、时间／电荷成功率、漏检、误合并／误拆分和条件误差，能量标定后另报能谱指标。网络分类99%或谱形接近不直接作为当前恢复率；局部窗容量、长尾历史和GAGG共享微胞状态都要进入训练／测试条件。

<a id="ml-pileup-comparison-20261014"></a>
### L20：2026-10-14——ML与拟合的实际模拟先导比较

**采用方向与实现**：依据L19的Regadío时序事件展开方向，自建含固定电压/导数/逆滤波特征的小型1D膨胀残差CNN，25,707参数、未知K事件概率/时间/电荷输出；不是论文U-Net/LSTM或SAEFit的完整复现。共享库只读，原人工审核范围及待审核状态不变。

**训练与独立评价**：原36个原生共享微胞SiPM窗口先按完整ROOT/种子分为26训练、10验证，随后抽切片。仿真真值仅监督/评分，不输入测试推断。新生成17个独立测试ROOT（0.1 MHz 8窗、1/10 MHz各2窗，其余五档各1窗），每窗20 μs、中央10 μs门；首轮200个保留最终测试未用。2400步CPU训练后由验证集选择第1200步/阈值0.65，保存哈希再测试，无测试调参。沿用完整含噪电压/足够动态范围前提，不使用ADC拒绝或量化限制；拟合保持原独立标定和算法参数。

**本项目计算**：同一目标按q≥500 pe、±10 ns一对一时间匹配且相对电荷误差≤20%计成功；电荷是实际共享微胞雪崩当量，不是沉积能量。

| 源注入率 / MHz | 新测试窗 | 阈上目标数 | 拟合成功 / % | CNN成功 / % |
| --- | ---: | ---: | ---: | ---: |
| 0.1 | 8 | 4 | 4 / 100.00 | 4 / 100.00 |
| 1 | 2 | 14 | 14 / 100.00 | 14 / 100.00 |
| 10 | 2 | 140 | 93 / 66.43 | 91 / 65.00 |
| 20 | 1 | 143 | 58 / 40.56 | 76 / 53.15 |
| 40 | 1 | 281 | 45 / 16.01 | 86 / 30.60 |
| 60 | 1 | 447 | 36 / 8.05 | 91 / 20.36 |
| 80 | 1 | 531 | 30 / 5.65 | 109 / 20.53 |
| 100 | 1 | 654 | 16 / 2.45 | 100 / 15.29 |

**失败与边界**：空/单/阈下和20–200 ns平均模板双脉冲控制正常，但10 ns双脉冲CNN恢复1/2（拟合2/2），40000 pe超训练范围单脉冲CNN0/1且有2个额外预测（拟合1/1）。训练电荷最大约7989 pe，未验证大幅度外推。CNN网格峰距8 ns后允许时间偏移，而拟合还合并最终<8 ns分量，约束不完全相同。高率时间移位对照保留预测数量/电荷，检验偶然匹配；其结果与漏检/误拆分/误合并、条件电荷/时间误差集中于[独立Notebook](../jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb#ml-pileup-comparison-20261014)。不能将高密度时间匹配和电荷先验贡献全部称为真实事件身份恢复。

**验证与下一步**：4项模型检查、同种子物理一致性及17窗输入/哈希/评分/仅电压网络重放通过，见[验证](../jupyter/run_records/alpha_ml_comparison_20261014_162000/validation-results.json)。训练约8.60分钟，离线单窗CPU中位耗时CNN 0.1147 s、拟合 4.723 s；不等同实时吞吐。高率一窗/档、100 MHz无同率验证、已知纯α家族与固定噪声/响应，仅为模拟先导。后续扩大独立曝光、比较神经候选+物理refit、宽度/形状/背景失配及能量标定，不能据此定义最大工作率。

<a id="ml-multimethod-20261014"></a>
### L21：2026-10-14——U-Net、Transformer与CNN＋物理拟合共同测试

**采用方向**：在L19/L20基础上新增小型1D U-Net（四级跳接，无论文LSTM瓶颈）及局部Transformer（32维、4 heads、两层，8 ns token并保留2 ns分支）；不是完整论文复现，没有实施SAEFit。复用[Regadío2021](../../papers/sources/2021_Regadio_NIMA_Unfolding-using-deep-learning-and-its-application-on-pulse-height-analysis-and-pile-up-management.md)、[DOI](https://doi.org/10.1016/j.nima.2021.165403)与[Wang2025](../../papers/sources/2025_Wang_NIMA_Transformer-algorithm-for-pile-up-correction-in-synchrotron-radiation-spectroscopic-detection-experiments.md)、[DOI](https://doi.org/10.1016/j.nima.2025.170215)，沿用已有精读的装置/采样/计数率条件与审核状态，不把其结论移植为GAGG高率能力。共享库只读。

**实现与选参**：保持原CNN和原盲拟合。新网络仍用26训练/10验证完整ROOT、相同特征/损失/随机切片及2400步，由验证集选中U-Net第800步/0.75、Transformer第2400步/0.85。两者输出1024点窗的500点中心，原CNN维持全窗推断，容量和计算时间不完全相同。混合方法使用CNN≥250 pe候选，在实际总电压上联合非负幅度及±6 ns时间调整，最终要求≥500 pe并合并<8 ns分量；`||H A−y||²+λΣ||h_j||²(A_j−A_CNN,j)²`中的λ只由验证集在0/0.1/1中选为0.1。冻结后不再调参，推断不输入真值、源率或已知K。

**本项目计算**：另生成22个新种子共同测试ROOT，低率8窗，其余各2窗，20 μs/窗，中心[5,15) μs。原17个已公开测试仅作回归、未用于选择新模型；原200个保留最终测试仍未使用。保留原生共享微胞状态和电子噪声，用完整电压、足够动态范围前提，不讨论ADC量化/削顶。统一q≥500 pe、±10 ns一对一时间匹配及电荷误差≤20%：

| 源注入率 / MHz | 窗数 | 阈上目标 | 拟合 / % | CNN / % | U-Net / % | Transformer / % | CNN＋物理 / % |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.1 | 8 | 8 | 100.00 | 100.00 | 87.50 | 87.50 | 100.00 |
| 1 | 2 | 10 | 100.00 | 90.00 | 80.00 | 90.00 | 100.00 |
| 10 | 2 | 148 | 57.43 | 60.81 | 57.43 | 62.84 | 62.16 |
| 20 | 2 | 283 | 38.87 | 44.52 | 43.46 | 45.58 | 45.94 |
| 40 | 2 | 594 | 17.00 | 31.82 | 35.35 | 28.96 | 32.15 |
| 60 | 2 | 761 | 10.12 | 23.39 | 25.10 | 20.37 | 23.65 |
| 80 | 2 | 1104 | 4.89 | 16.12 | 19.38 | 13.68 | 16.30 |
| 100 | 2 | 1344 | 2.46 | 17.86 | 18.38 | 15.03 | 13.17 |

**失败、控制与限制**：九组平均模板控制中10 ns双脉冲U-Net/Transformer/拟合2/2，CNN1/2、混合0/2；40000 pe超训练范围单脉冲三网络及混合均0/1、拟合1/1。均不作为实物分辨能力。每档高率仅两窗；100 MHz无同率验证；时间移位显示高密度偶然匹配贡献。条件电荷/时间P16/50/84、漏检、额外预测、误拆分/误合并标志、固定区间示例、残差与离线耗时集中于[独立Notebook第二轮](../jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb#ml-multimethod-20261014)。较低波形RMS不保证逐事件电荷正确，混合不能补出完全没有神经候选的事件。尚无沉积能量标定、背景/形状失配和实物泛化，不给最大工作率或能谱分辨率。

**验证与待办**：5项新算法检查、同种子物理字段一致性、39窗哈希/独立ROOT与种子/仅电压精确重放/统一评分/物理电压检查通过，见[验证记录](../jupyter/run_records/alpha_ml_multimethod_20261014_180000/validation-results.json)。混合优化117次，终止成功117次、达到预算或其它非成功停止0次；终止成功不是事件恢复成功。实际源码、宏、种子、环境及共享目录清单故障恢复记录保留在[新冻结实验](../pileup/experiments/alpha_ml_multimethod_20261014/README.md)。后续扩大独立曝光并验证宽度/形状/噪声/背景失配；候选增峰或幅度范围扩展另建实验，不用本测试调参。

<a id="ml-low-rates-20261015"></a>
### 项目计算更新：2026-10-15——五方法低率0.1/1/10 MHz复查

用户指定只分析已有三档；复用L21的12个新低率ROOT，另12个公开回归窗分别列出，没有G4新运行或重新训练/调参。保持完整电压、共享微胞状态、中央[5,15) μs、q≥500 pe、±10 ns及电荷误差≤20%。真值仅用于推断后的重评分/分组/失败解释。

| 源注入率/MHz | 目标数 | 盲拟合 | CNN | U-Net | Transformer | CNN＋物理拟合 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.1 | 8 | 8/8（100.00%） | 8/8（100.00%） | 7/8（87.50%） | 7/8（87.50%） | 8/8（100.00%） |
| 1 | 10 | 10/10（100.00%） | 9/10（90.00%） | 8/10（80.00%） | 9/10（90.00%） | 10/10（100.00%） |
| 10 | 148 | 85/148（57.43%） | 90/148（60.81%） | 85/148（57.43%） | 93/148（62.84%） | 92/148（62.16%） |

**失败性质**：新0.1 MHz五方法时间匹配均8/8，U-Net/Transformer各1次电荷失败；1 MHz只有Transformer漏检1次，其余方法时间匹配10/10，CNN有1、U-Net有2次电荷失败。此前笼统的低率“漏失更多”已纠正。0.1 MHz孤立655.54 pe脉冲被两网络高估42.14%/23.43%；1 MHz618.19 pe近邻81.16 ns时CNN/U-Net高估24.65%/34.63%，Transformer峰得分0.84670未达冻结0.85；另有7622.61 pe脉冲被U-Net低估25.13%，仍在训练范围内。阈值没有按测试调整。

**与长度关联**：复用L18的当前模板FWHM181.505 ns、W90 728.024 ns。三档实测有电荷α率0.1/0.65/8.7 MHz，`R_charge W90`为0.0728/0.4732/6.3338；全部有电荷α（含保护区、阈下）查两侧近邻，W90内目标为0/8、7/10、148/148。Poisson参考概率不当失败率，分组受幅度/共享微胞混杂，不给因果分辨时间或最高工作率。0.1/1 MHz的8/8与10/10只有少量目标，逐ROOT删除范围不是置信区间，回归/新测试不合并。

**证据和边界**：[同一本Notebook第11–14节](../jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb#ml-low-rates-20261015)保存条件电荷/时间P16/50/84、全部漏检/额外预测/合并拆分、同目标配对、原生失败图及逐窗敏感性；[检查](../jupyter/artifacts/alpha_ml_low_rate_20261015_082901/validation.json)覆盖24低率窗五方法原评分、哈希和计数分解。快照及重跑见[后续实验](../pileup/experiments/alpha_ml_low_rate_20261015/README.md)。没有新增论文结论、改变Regadío待审核或Wang第549项的限定范围；共享库只读，已有文献器件/采样/噪声条件不能直接移植为GAGG低率性能。尚未沉积能量标定或实物泛化；本轮未模拟低于0.1 MHz。

<a id="alpha-gamma-mixed-20261015"></a>
### L22：2026-10-15——α/γ 原生混合波形的拟合先导

**已有依据与审核状态**：只读复用 L14–L16 的 Kobayashi 2012（GAGG/PMT 平均 α/γ 形状）、Marrone 2002（液闪/PMT 形状拟合、实测 γ 波形软件叠加）及相关 SiPM 卷积笔记，完整出处、页码与适用条件见 [既有采用记录](#psd-fit-use-20261013)。共享库的 AI 草稿待人工审核状态保持。没有文献证据可把这些装置下的性能直接视为当前 GAGG/S13360 混合逐事件准确度。

**本项目假设与实现**：两条旧孤立模拟电压固定 α/γ 慢分量权重 0.7414998/0.2047319，快/慢常数95/351 ns，加同一读出核；未知K、时刻、电荷、类型，双模板反卷积给候选，非负最小二乘及近邻离散类型/连续时刻联合搜索。算法在理想空/单/双脉冲控制中修正，正式标签只作评价。有限局部搜索不保证全局最优。另存的矩阵缩减加速试验未达到完整流程的严格数值一致性容差，没有混入正式统计。

**实际输入**：两个独立冻结实验，各自16个新正式20 μs窗，α有效p–¹¹B源注入1/10 MHz、γ注入10 MHz，中央[5,15) μs，每条件8窗。γ分别为原有轫致辐射exp(−E/100 keV)/E（10 keV–5 MeV）与661.657 keV单能源；后者仅在独立快照改写C12源槽三条能量和响应标签，不是C12俘获谱。S13360-6050PE v2，完整窗共享微胞状态，2 ns采样、0.5 ns孔径、0.1 mV电子噪声。使用完整analog电压，遵循动态范围充足前提。

**本次计算结果（不是文献性能）**：完整恢复定义为真值Q≥500 p.e.、一对一类型无关时间匹配≤10 ns、类型正确、电荷误差≤20%。γ注入率均为10 MHz。

| γ 来源 | α 注入率 / MHz | α 完整恢复 | γ 完整恢复 |
| --- | ---: | ---: | ---: |
| 轫致辐射谱 | 1 | 44/64（68.8%） | 2/30（6.7%） |
| 轫致辐射谱 | 10 | 270/571（47.3%） | 3/42（7.1%） |
| 661.657 keV | 1 | 49/58（84.5%） | 29/32（90.6%） |
| 661.657 keV | 10 | 282/605（46.6%） | 12/31（38.7%） |

同时保存检出率、条件分型率、选择纯度、条件电荷/时间误差、未匹配候选、邻接图误拆分/误合并及按独立完整窗bootstrap区间。147次理想控制的孤立18/18成功；双脉冲中的逐脉冲恢复在10 ns为32/36，20–1000 ns各36/36（每档18个双脉冲波形），严格同时受算法8 ns最小间隔限制；不含光子和微胞涨落，不能视为真实分辨时间。五个时间平移对照用于检查高密度偶然匹配，不作为p值或直接效率修正。

**失配诊断与限制**：不计入正式统计的轫致辐射预检查窗中，用真值扣除阈下γ的线性电压，保留实际共享恢复电荷、暗链与噪声，α完整恢复6/10→9/10、γ仍1/4。该单窗辅助诊断支持弱γ干扰α，但不能解释全部分型错误，也不是盲处理或降低γ率的重模拟。电荷是实际单胞等效量，未建立两类沉积能量条件标定，不能称为能量分辨率或能谱已恢复。样本量与模板标定均有限。

**可追溯记录**：[新Notebook](../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb)、[轫致辐射实验](../pileup/experiments/alpha_gamma_brems_fit_20261015/README.md)、[单能γ实验](../pileup/experiments/alpha_gamma_mono662_fit_20261015/README.md)。保留36个新ROOT，两个写入失败原文件同种子重放后排除；正式分母仍32窗。全部有效文件进行各树/分支读取、宏/种子/哈希核验。未修改共享文献库。

**待验证**：扩大独立窗；建立同器件、不同沉积能量和位置/微胞占用下的独立模板与条件电荷能量标定；扫描α/γ率比；比较连续慢分量与弱背景模型。若需新文献精读或入库，仍在papers项目处理。

<a id="alpha-gamma-failure-diagnosis-20261015"></a>
**2026-10-15 失败机制补查与改进方向**：正式32窗中，α源1 MHz时，轫致辐射/单能γ的非零电荷率分别8.600/0.650 MHz，而不是相同的10 MHz探测率；轫致辐射30个过阈值γ中25个已时间匹配，其中22个被判α。其阈下γ占该类总电荷78.1%（α源10 MHz时70.0%），最终只保留大分量的模型存在明显解释缺口。新真值辅助诊断只用两条预检查波形：给定全部非零粒子的时刻与类型、保留小分量并从原电压拟合幅度，轫致辐射窗α/γ完整恢复由盲拟合6/10、1/4变为10/10、3/4。该结果使用真值，不是盲恢复成绩或严格物理上限。放松类型为同一时刻的双模板连续系数后，残差2.821→2.645 mV，α却变为5/10，支持“残差更小不保证分型正确”的本次计算结论。原来的弱γ真值扣除实验仍单独保留。

改进顺序为：分开内部建模阈值与报告阈值，联合描述弱γ并约束复杂度；使用独立孤立样本标定模板分布和信号相关噪声；再比较混合样本CNN初值＋物理联合拟合。保持500 p.e./10 ns/20%评价，同时检查纯度、漏检、假候选、误拆分/合并及电荷偏差，阈值拒判须同时报告效率。此次建议为项目待验证方案，不能从辅助真值实验承诺盲算法增益。

文献依据据既有精读：[Tamagawa等2015完整笔记](../../papers/sources/2015_Tamagawa_NIMA_Alpha-gamma-pulse-shape-discrimination-in-Gd3Al2Ga3O12-GAGG-Ce3%2B-crystal-scintillator-using-shape-indicator.md)，PDF2–4，GAGG/R6233K PMT、500 MHz/2048 ns、选定能窗的平均模板SI及低能误判边界，[出版社/DOI](https://doi.org/10.1016/j.nima.2015.05.052)；[Marrone等2002完整笔记](../../papers/sources/2002_Marrone_NIMA_Pulse-shape-analysis-of-liquid-scintillators-for-neutron-studies.md)，PDF3–4、6–7，NE213/C6D6液闪/PMT标定形状拟合及堆积，[出版社/DOI](https://doi.org/10.1016/S0168-9002(02)01063-X)。本轮核对出版社题录与DOI，未重新全文或PDF精读，AI待审核状态及已有局部采用意见保持。以上文献没有证明当前GAGG/S13360多MHz原生混合流的性能。

详细过程：[现有Notebook第10节](../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-failure-diagnosis-20261015)；[诊断结果](../jupyter/artifacts/alpha_gamma_fit_20261015_093500/oracle_diagnostic/preflight-onset-types.json)；[可复用诊断模块](../pileup/analysis/alpha_gamma_oracle_diagnostic.py)。原有25个Notebook单元及输出保留，新段落记录实际计算与方案边界。

<a id="alpha-gamma-weakfit-pilot-20261015"></a>
**随后实际尝试的盲改进及负结果**：新增[内部门槛与报告门槛分离实现](../pileup/analysis/alpha_gamma_weakfit.py)，弱分量保留在完整电压模型中，输出/评价仍为500 p.e./10 ns/20%。拟合只输入波形和固定标定。两条预检查窗的三次探索中，轫致辐射原始α6/10、γ1/4；内部200 p.e.变为9/10、0/4，内部100 p.e.变为5/10、1/4。单能γ原始α5/7、γ4/6；内部100 p.e.变为1/7、5/6。其γ选择纯度同时从5/5变为5/10，说明只看该类恢复增加会掩盖α泄漏。残差均降低，却未取得同时稳定改善，故不替换原基线、不扩大到32窗宣传新准确度。这些是小样本探索，不是独立性能检验。

当前证据支持下一步同时做弱分量复杂度约束、独立模板/噪声标定与分型置信度，而非单独继续降低门槛。尚未实现正则背景或混合CNN。新控制确认小分量进入模型却不进入过阈值报告，5项既有混合拟合测试与3份保存模型重构检查通过。[Notebook第11节](../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-weakfit-pilot-20261015)保留波形、残差、效率、纯度、假候选及运行时间；[结果索引](../jupyter/artifacts/alpha_gamma_fit_20261015_093500/weakfit_pilot/summary.json)。没有修改文献结论或共享库。


<a id="pileup-fit-optimization-literature-20261015"></a>
### L23：2026-10-15——面向混合α/γ的拟合优化文献

**问题与实际证据**：L22的三次低内部门槛盲试验降低了残差，却未稳定同时改善α/γ恢复与纯度。本次按`paper-search`只读共享库，复用完整精读笔记并定点核5篇原文；没有修改拟合器、启动新仿真或宣称准确度提高。完整题录、笔记/PDF链接、原页位置、仪器/噪声/计数率条件和优先级见[现有Notebook第12节](../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#pileup-fit-optimization-literature-20261015)及[证据记录](../jupyter/run_records/pileup_fit_literature_20261015/evidence.json)。

**采用的方法原则与限制**：

- [Zhang等2024](../../papers/sources/2024_Zhang_NIMA_Development-of-a-pile-up-pulse-recovery-algorithm-for-the-LaBr3-detector.md)，[DOI](https://doi.org/10.1016/j.nima.2024.169273)，本轮PDF3–4、实际看PDF3图2–3/§3.2：幅度相关模板及残差补峰，明确模板失配可制造假脉冲。LaBr3/PMT、500 MHz及同模板合成条件保留；当前已有残差补峰，新增重点是独立模板分布与接受准则，不直接移植10 ns/60 Mcps性能。
- [Lopatin等2012](../../papers/sources/2012_Lopatin_IEEE_Pileup-attenuation-for-spectroscopic-signals-using-a-sparse-reconstruction.md)，5页会议论文，DOI未独立核；本轮PDF3原图/式(4)–(6)：较强LASSO选起点、较弱LASSO估能量，原式没有非负约束且式(5)条件冲突。借鉴两阶段原则；非负/分组候选、离散类型比较和NNLS减偏是本项目拟议改写，不是原文完整复现。HPGe实测谱缺完整采样/率参数，不据此给GAGG上限。
- [Fixsen等2004](../../papers/sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md)，[DOI](https://doi.org/10.1016/j.nima.2003.11.313)，本轮PDF2原图/式(1)：独立标定均值与完整协方差，逆协方差加权。仅TES非线性/非平稳噪声模型验证，不是GAGG混合恢复证据；协方差随参数变化时，完整高斯似然还需log det项，协方差正则化及共享微胞相关性须本项目验证。
- [Luo等2018](../../papers/sources/2018_Luo_NIMA_Pulse-pile-up-identification-and-reconstruction-for-liquid-scintillator-based-neutron-detectors.md)，[DOI](https://doi.org/10.1016/j.nima.2018.03.078)，本轮PDF4文字：液闪n/γ四组合重建仅作联合类型模型先例；当前已实现近邻四组合。原判据为有符号差求和后取绝对值，可能抵消，不直接作为本项目似然；受控双脉冲结果不推广至随机高重数。
- [Novak与Vencelj2009](../../papers/sources/2009_Novak_TNS_Gauss-Seidel-Iterative-Method-as-a-Real-Time-Pile-Up-Solver-of-Scintillation-Pulses.md)，[DOI](https://doi.org/10.1109/TNS.2009.2030725)，本轮PDF2–3文字：已知可靠时刻、固定线性形状下用局部小系统解幅度。可借鉴缓存、局部块及暖启动；不解决漏候选/错类型，求解核吞吐不当完整恢复性能。

**项目下一步（待实现）**：独立原生孤立/空窗的模板与噪声标定 → 加权稀疏/复杂度候选选择 → 保留受约束弱分量并减小幅度收缩 → 类型/数量模型比较与置信度 → 局部求解加速。这个组合是项目方案，不是五篇论文共同验证的算法；保留500 p.e./10 ns/20%及完整电压口径，报告效率、纯度、漏失、误拆分/合并和电荷误差，能量分辨率另需标定。当前32窗已公开分析，新方法须使用独立标定/验证并冻结后另留测试；真值不进入盲推断。

**审核及核验范围**：五篇仍为AI草稿待人工审核；global-review B23-024/B23-052的原式与同模板验证限制保留，其他论文的限定认可不转移。有限网页题录核对没有补足Lopatin的可读出版社DOI证据；未下载新论文。原页文字/图像、共享来源SHA256与验证记录保存在本地jupyter。Trigano2015和Kafaee/Goodarzi2020仅查到引文/队列线索，列为待补原件，未把书目关系当已核算法依据。共享库未写入。


<a id="alpha-gamma-regularized-20261015"></a>
### L24：2026-10-15–16——独立标定、稀疏加权与条件减偏的实际实现

**采用范围**：落实L23中Zhang2024的幅度相关模板方向、Lopatin2012的选择/估计两阶段原则、Fixsen2004的独立均值/噪声标定；近邻四类型模型比较沿用已有实现。所有来源、原页与装置条件见[L23](#pileup-fit-optimization-literature-20261015)及其Notebook链接，五篇的AI草稿待人工审核状态不变。本轮没有新增文献结论或修改共享库，也没有使用Lopatin式(5)的歧义条件。本实现是项目改写，不是原论文完整复现，性能均为项目模拟计算。

**数据与模型**：另建[冻结实验](../pileup/experiments/alpha_gamma_regularized_fit_20261015/README.md)，56个原生孤立标定窗、18个单粒子开发检查窗、4个混合选参窗、8个冻结后新测试窗。保留两个父实验的G4物理、S13360-6050PE v2、2 ns采样/0.5 ns孔径/0.1 mV噪声与完整共享微胞历史；没有按测试改变物理参数。盲拟合只输入完整含噪电压和冻结参数，真值只在标定/评价中使用。保持用户指定的足够动态范围。

**由旧值到新值**：原固定α/γ慢分量0.741500/0.204732，改为在95/351 ns基函数上的独立 `f(Q)=a+b log(Q/1000 p.e.)`；α(a,b)=(0.680950,0.032393)、γ=(0.243212,0.004669)，时延分别−1.9104/−0.3973 ns。56个入射标定事件中仅16个α与13个γ达到Q≥50 p.e.，Q覆盖1545.6–16860.2和256.9–1996.7 p.e.，范围外形状夹住而电荷不夹住。独立残差给出基线方差0.0100869 mV²、信号方差系数0.614075 mV、2 ns相关系数0.899770。中心单粒子不足以标定位置/占用依赖，模板失配也进入残差，不能把全部残差方差称作本征统计噪声。

**算法与验证集选择**：标准化α/γ字典列的非负L1选候选，条件无L1 NNLS减小幅度收缩，保留弱分量，局部时间/类型与近邻四组合重拟合。候选80 p.e.，报告500 p.e.；原2 ns发现、4 ns拟合、2 ns保存模型。候选权重使用80 ns平滑电压，最终协方差由纯电压初步模型产生后固定，是可行加权最小二乘而非完整参数相关似然。λ=3/8只在4个混合验证窗选择，预声明的宏平均完整恢复F1为0.59490/0.47698，选λ=3后于2026-10-15 17:14:52 +08:00冻结，晚于此才生成测试。没有完整后验、拒判或全局最优保证；形状在每轮内固定，条件减偏不保证最终电荷严格无偏。

**同一新测试样本的结果**：α源1/10 MHz，γ源均10 MHz；中央[5,15) μs，真值Q≥500 p.e.、一对一时间匹配≤10 ns、类型正确、电荷误差≤20%。电荷为实际完全恢复单胞等效量，不是沉积能量。

| γ 来源 | α源 / MHz | α完整恢复：旧 → 新 | γ完整恢复：旧 → 新 |
| --- | ---: | --- | --- |
| 轫致辐射 | 1 | 14/16 → 11/16 | 4/14 → 3/14 |
| 轫致辐射 | 10 | 81/153 → 97/153 | 1/11 → 1/11 |
| 661.657 keV | 1 | 15/16 → 13/16 | 1/2 → 2/2 |
| 661.657 keV | 10 | 69/155 → 91/155 | 1/2 → 0/2 |

额外预测53→19、漏检60→61、误拆分标志0→0、误合并标志24→31。高率α改善伴随低率退步、γ不稳定，某些条件的电荷偏差/标准差变差；每条件仅2窗，单能γ每档仅2个目标，不能定义稳定效率或最大率。旧基线只替换带KKT检查的等价线性NNLS缩减，非线性方法不变；单窗与原直接NNLS回查的数量、类别、配对和恢复数一致，但不承诺所有局部搜索逐位一致。原32窗统计与本轮分开。

**保留的失败**：首版权重给噪声下行涨落过高权，单粒子开发检查λ=3只有α2/8、γ1/3恢复且额外13；改后8/8、3/3、额外0，与旧基线相同。开发集已用于修正，不当独立测试。40次同模板双脉冲检查中，2/4 ns网格恢复计数相同；每网格10/20/40/80/200 ns分别为0/8、2/8、8/8、8/8、6/8，含10/20 ns漏候选与200 ns γα局部错解。该控制没有光子/微胞涨落，不能宣称40 ns实物分辨能力。

**追溯与待验证**：已有Notebook[第13节](../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-regularized-20261015)记录全过程、纯度、完整条件误差、固定区间图和逐目标配对。86个ROOT全树/分支、宏/种子/哈希与冻结时序核验，16次测试拟合共同重评分、8个新方法模型重构、一窗只波形重放和13项检查见[验证记录](../jupyter/run_records/alpha_gamma_regularized_20261015/validation.json)。8次L1成功终止只说明满足优化器停止规则，最大投影梯度0.15813（λ=3）仍高于严格1e−5梯度容差；候选选择不是已证实的精确最优解。尚未消融各改动、完成γ位置/占用标定、分型概率校准、沉积能标或实物验证。下一步在新实验验证求解精度、受约束残差增峰、单/双模型竞争与多初值，处理控制中已见的漏候选/局部解；不得以本次已读测试再选参。

<a id="alpha-gamma-model-selection-20261016"></a>
### L25：2026-10-16——α/γ等权目标下的联合模型搜索与独立复验

**采用依据**：沿用[L23](#pileup-fit-optimization-literature-20261015)的Zhang2024残差候选/模板失配警示、Lopatin2012两阶段原则，以及Luo2018联合类别模型先例；原页、装置和采样/噪声限制不变，AI草稿待人工审核状态不变。本次只读复用已有资料，没有新增论文性能结论或修改共享库。以下均为项目自己的实现和计算；复杂度代价、KKT求解和搜索预算不是论文给定参数。

**实现及划分**：用户要求α/γ逐事件恢复同样重要。新增[独立实验](../pileup/experiments/alpha_gamma_model_selection_20261016/README.md)：保持非负L1目标，用活跃二次子问题检查1e−5投影KKT；近邻500 ns内联合两次到达时间、电荷与四类型组合；残差补峰通过固定加权SSE收益与复杂度代价比较。复用18个单粒子和4个混合开发窗，类别/条件等权完整恢复F1为上一版0.59490、增峰代价1为0.62757、代价9为0.62557，冻结代价1、λ=3。09:14:48冻结，09:15:09后生成8个新种子测试窗，每条件两窗，三方法同窗比较。未用新测试择参。

**计算结果**：已公开的同模板/仅电子白噪声控制，10/20/40/80/200 ns × 四类型组合，24/40→40/40脉冲完整恢复；冻结源码重放一致，不是物理10 ns分辨率。新原生低α率轫致辐射条件，相比上一版为alpha 12/16→14/16、gamma 7/12→11/12。高率为brems：α 81/142→78/142，γ 0/6→1/6；mono662：α 85/143→86/143，γ 0/6→1/6。高率未同时可靠恢复两类，不能称为问题已解决。

| γ来源 | α源/MHz | 类别 | 原始拟合 | 上一版稀疏 | 本轮联合 |
| --- | ---: | --- | --- | --- | --- |
| brems | 1 | alpha | 9/16（56.2%） | 12/16（75.0%） | 14/16（87.5%） |
| brems | 1 | gamma | 3/12（25.0%） | 7/12（58.3%） | 11/12（91.7%） |
| brems | 10 | alpha | 61/142（43.0%） | 81/142（57.0%） | 78/142（54.9%） |
| brems | 10 | gamma | 0/6（0.0%） | 0/6（0.0%） | 1/6（16.7%） |
| mono662 | 1 | alpha | 17/19（89.5%） | 17/19（89.5%） | 17/19（89.5%） |
| mono662 | 1 | gamma | 2/2（100.0%） | 2/2（100.0%） | 2/2（100.0%） |
| mono662 | 10 | alpha | 84/143（58.7%） | 85/143（59.4%） | 86/143（60.1%） |
| mono662 | 10 | gamma | 3/6（50.0%） | 0/6（0.0%） | 1/6（16.7%） |

仍用Q≥500 p.e.、类型无关匹配≤10 ns、类型正确且电荷误差≤20%；Q是实际单胞等效电荷，不是沉积能量。源注入率α1/10 MHz、γ10 MHz不等于阈上探测率。G4物理、共享微胞历史和读出未改变，输入完整含噪电压、动态范围足够。

**失败诊断和下一步**：本轮高α率12个阈上γ中，2个完整恢复、1个漏检、8个类型错、1个仅电荷不合格；高率α被误判为γ也从上一版9个增至26个，必须结合纯度判断。在旧开发窗刻意提供真实数量/时刻、用于形状和协方差的真实电荷及协方差用真实类型，另比较固定真实类型/双类型NNLS。高率brems即使固定真实类型也只恢复α53/73、γ1/4；该诊断不是盲算法，也不是严格物理上限。这提示当前平均模板及噪声模型可能受形状失配、涨落或重叠病态性影响，尚未做完整因果分离。开发目标78/161个α低于α标定下限1545.6 p.e.、2/12个γ高于γ标定上限1996.7 p.e.；应优先补齐低Qα、γ及位置/前史占用标定，再验证受约束形状分布模型，不能只靠进一步补峰或降低阈值。

**追溯**：现有Notebook[第14节](../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-model-selection-20261016)保存逐类纯度、电荷条件误差、漏检/误拆分/误合并及图像；[验证记录](../jupyter/run_records/alpha_gamma_model_selection_20261016/validation.json)核对8个ROOT全分支、132份冻结源码、24个拟合的共同评分/模型、一窗盲重放及17项算法检查。本轮线性选择达到KKT容差，但非线性搜索有有限预算，不保证全局最优；每条件两窗，未验证稳定γ效率、沉积能谱或实物泛化。

<a id="alpha-gamma-profile-fit-20261016"></a>
### L26：2026-10-16——全窗电荷重估、相关噪声与受约束类内形状

**文献复用与条件**：只读复用[Fixsen2004精读](../../papers/sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md)的独立平均模板/协方差标定原则（原PDF2 §2–3；TES模型、非平稳噪声，不能移植其性能为GAGG结论），以及[Novak2009精读](../../papers/sources/2009_Novak_TNS_Gauss-Seidel-Iterative-Method-as-a-Real-Time-Pile-Up-Solver-of-Scintillation-Pulses.md)关于固定模板、已知时刻线性系统的条件及有限邻域/模型误差（原PDF1–3，NaI(Tl)离线例与理想模型/FPGA吞吐分开）。两份笔记仍是AI草稿待人工审核。本轮没有新增原文精读或修改共享库；来源哈希见[literature-reuse.json](../jupyter/run_records/alpha_gamma_profile_fit_20261016/literature-reuse.json)。以下算法与数值是项目实现，不是论文复现。

**实现与选择**：新增全窗非负电荷重估与时刻搜索；原AR(1)/对角/均匀及SPE卷积带状协方差比较；电压初始前沿约束；受标定RMS约束的类内慢分量变化。共12组设置，完整失败结果保留。带状协方差的经验乘数2.07796由29条独立孤立标定得到，包含模型失配，不是器件本征参数；未精确描述共享微胞相关性。类内形状均值、时刻与权重固定，两边界模板非负组合，正则化是项目经验规则，不是校准的Bayesian后验。

开发采用原4混合窗，α/γ×四条件等权完整恢复F1：原V3 0.62756657，shape_0.5为0.63009504，其余11候选更低。开发α92/161→93/161，γ5/12不变，额外预测12→11；18个孤立开发窗的阈上α8/8、γ3/3保持，无额外预测。冻结0.5倍原标定RMS、三倍尺度边界，未根据后续测试改参。原G4/SiPM和完整analog_mV保持，Q为实际单胞等效电荷；源注入率与探测率分别理解。

自由时刻开发优化48次中17次正常停止、26次达到60次迭代上限、5次ABNORMAL。负结果只限于当前初始化和预算，不代表该类方法充分收敛后的极限；选中shape_0.5后处理固定时刻。见[停止状态](../jupyter/artifacts/alpha_gamma_profile_fit_20261016/optimizer-status.json)。

**独立复验**：10:29:36冻结，10:30:40后生成8个新原生窗；两γ谱×α源1/10 MHz×每条件2窗，γ源10 MHz，Q≥500 p.e.、时间匹配≤10 ns、类型正确且电荷误差≤20%。两方法同窗、真值在拟合后评分。

本轮同一批8个新原生窗比较，受约束可变形状拟合使α完整恢复从170/290（58.62%）升到173/290（59.66%），γ仍为18/38（47.37%）。额外预测均17个，漏检49→48，误拆分标志均0，误合并标志18→17。高α源10 MHz条件合计：α148/262→151/262，γ5/16不变。改进幅度很小，且每条件只有两窗，不能据此声称稳定提高了实际效率。

保留可变形状作为实验候选，暂不替换上一版比较基线。γ正确分型数21→22，但同时满足类型和电荷要求的预测比例由18/49降为18/51；γ的完整恢复F1从0.41379降到0.40449。因此，本轮没有达到α和γ都可靠逐事件恢复的目标。12组试验也表明，在当前标定与数据条件下，更多拟合自由度、更小残差或更复杂的噪声权重不能自动转化为更好的两类恢复。

**诊断与局限**：使用祖先真值的形状描述量在开发α/γ中位数分别0.7283/0.2321，不能把所有失败归因于两类本来相同；这也不保证叠加电压可稳定分解。显式提供真数量/时刻/Q以及协方差用真实类型后，带状协方差且固定真实类型的高率brems仍仅恢复α52/73、γ2/4；是条件诊断而非盲性能或严格物理上限。下一步优先补齐α约500–1500 p.e.、γ标定上端及不同位置/前史条件的独立响应标定，建立条件化平均模板和形状散布；用含真实光子统计与共享微胞恢复的受控双脉冲扫描间隔和幅度比，检验模板失配、随机涨落与分解病态性各自的作用。当前证据尚不足以把失败归因于单一原因，也不能宣称已达到物理分辨极限。该补充标定与因果扫描本轮尚未实施。

**追溯与验证**：[冻结实验](../pileup/experiments/alpha_gamma_profile_fit_20261016/README.md)、[Notebook第15节](../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-profile-fit-20261016)、[验证](../jupyter/run_records/alpha_gamma_profile_fit_20261016/validation.json)。24项测试、8个新ROOT全分支、144份冻结源码、86个原生拟合模型/评分、40个控制评分与一窗只输入电压的完整重放通过。两方法在20组同模板仅电子白噪声控制均40/40，不是实测10 ns分辨率。旧52个Notebook单元及输出保留。

## 新文献整理需求待办

下列是后续定点核验或文献整理需求；共享库的修改在 `papers` 项目中处理。本轮未使用存在歧义的公式实现算法，因此不以重新全文精读作为启动项目响应检查的前提。

| 发现日期 | 当前问题与已有来源 | 缺少的公式、参数或证据 | 需回查的原文位置／新增整理内容 | 状态 |
| --- | --- | --- | --- | --- |
| 2026-10-15 | L23，两阶段/迭代稀疏重建的期刊方法 | Trigano/Gildin/Sepulcre2015（DOI线索10.1109/LSP.2015.2406911）及Kafaee/Goodarzi2020（TNS67,858–862；DOI未核）的独立原文/完整笔记 | 当前只有他文引文；后者已有共享队列键Kafaee-Goodarzi-2020-TNS67-858。需要时在papers项目核同身份原件并安排精读，不由本项目修改队列 | 书目线索，未读原文；不作为已验证性能证据 |
| 2026-10-15 | L23，模板、噪声与模型选择 | 同GAGG/SiPM的独立模板散布/协方差，弱分量复杂度、非规则增峰检验及类型不确定性 | Zhang PDF3–4、Lopatin PDF3、Fixsen PDF2已定点核；具体实现缺口与L01/L03待办合并追踪。收缩/带状协方差和减偏NNLS由项目独立验证，不冒充原文算法 | L24/L25/L26已记录独立新测试及失败；L26已试受约束类内形状，位置/占用标定与概率校准仍待完成 |
| 2026-10-14 | L19，ML／AI pileup候选 | 未知K事件输出、分量损失／尺度／排序、GAGG训练标签与独立ROOT划分；代码／权重／数据许可与真实泛化 | 优先定点核SAEFit、Regadío U-Net及Morad实现，按需要在papers项目安排新精读；本项目另建同条件空／单／双及多事件基准 | 文献检索完成；L20/L21已实施CNN、U-Net、Transformer、CNN＋物理拟合共同模拟先导，论文完整架构及实物泛化未验证 |
| 2026-09-30 | L01，残差增峰与幅度相关模板 | 式 5 局部坐标；噪声权重、增峰显著性、复杂度及终止规则的完整实现 | Zhang PDF2–4、式 4–5及方法段；若原文未给出，后续按需查看作者代码并记录版本，不把自选参数当文献参数 | 待实现前定点核验；本轮采用方法方向 |
| 2026-09-30 | L02，模型选择与假拆分 | p 值阈值、平滑后相关噪声及拟合区间边界；合成与真实噪声差别 | Adrich PDF2–5，方法、图 6；未给出项由项目独立验证集校准 | 待核；不影响先作单脉冲失配检查 |
| 2026-09-30 | L03，后备稀疏候选方法 | 式 5 冲突条件、分组定义、实际采样和计数率缺口 | Lopatin PDF3 式 5/算法 1；需原文、勘误或作者实现证据，不自行补下标 | 后备路线启用时处理 |
| 2026-09-30 | L04，当前能区的真实 α 响应 | 连续 α 能区的光输出、长尾、位置依赖及与现模型不同的完整参数组 | 先核现有单粒子 ROOT 和训练数据；需要新的低能 α 标定文献时，在 papers 项目单列整理 | 待项目响应检查后确定缺口 |
| 2026-10-10 | L06–L08，通用 SiPM 数字化方案 | 选定型号/过压/温度后的到达光谱 PDE、增益、pixel布局与恢复、SPTR/DCR/串扰/后脉冲、同读出链 SPE 与 ADC 参数 | 已核S13360-6050PE厂家资料；还需实物恢复/SPTR/AP、低过压光谱PDE和SPE/ADC读出链标定；新精读在 papers 项目处理 | 2026-10-13已建立数据表驱动配置，见L13；实物标定待完成 |
| 2026-10-10 | L09–L12，新增P10与现有离线实现 | 串扰候选/观测概率转换、AP分量与条件化、真实SPE和恢复、公共负载参数、前史收敛 | 优先实物标定；必要时核G4SiPM式4–8、Part I表2与Part II式4–6的原实现，不把AI辨析当作者勘误 | 本轮已核方法；参数/代码复现与器件验证未完成 |
| 2026-10-13 | L14–L16，GAGG + S13360波形拟合PSD | 同晶体/温度/偏压/前端的α/γ形状与SPE标定、光量依赖、独立标签和噪声协方差 | 先用独立单脉冲校准/验证数据检查；新增实物拟合文献若需入库，在 papers 项目另行处理 | 函数文献与两条现有模拟波形拟合已完成；独立分类性能、参数不确定度及实物标定待完成 |
| 2026-10-13 | L01/L02/L05，纯α计数率恢复先导试验 | 有SiPM恢复/相关涨落时的候选阈值、局部削顶恢复和电荷→沉积能量标定；最终独立测试的效率不确定度 | 优先本项目独立校准/验证；若需采用论文完整噪声/模型选择算法，再定点核原实现，不把本次自选规则当文献参数 | 新24窗盲拟合与264个简化控制已完成；真实器件与正式性能验证待做 |
| 2026-10-14 | L17，文献计数率与本项目逐事件恢复率的比较 | 统一峰效率/逐事件成功分母，同探测器条件的独立模板、噪声与饱和验证；高率独立真值效率 | 已核四篇关键原页；若采用具体算法，仍需其阈值/噪声/模型选择实现；新精读在papers项目处理 | 本轮完成条件核对；未推定GAGG+SiPM最大计数率 |
| 2026-10-14 | L18，计数率与信号长度 | 不同宽度定义、前沿/慢尾独立变化、固定电荷/峰高及相关噪声下的同指标恢复曲线 | Zhang PDF2–3宽度口径已核；优先项目独立单/双脉冲标定，若需新器件模板文献则在papers项目整理 | 当前模板、泊松占用和原样本近邻分组已分析；宽度因果扫描与最大工作率未验证 |


### 项目验证更新：2026-09-30 α 响应与数据接口

以下是本项目计算，不是上述文献的性能结论。过程、图表和执行输出保存在 [Alpha_Response_Calibration.ipynb](../jupyter/Alpha_Response_Calibration.ipynb)，结果及限制见 [计算摘要](../jupyter/artifacts/alpha_response_20260930_v2/summary.json)。

- 四档各按完整 ROOT 划分为 30 校准、20 验证、50 最终测试；本轮只重新读取 10⁵、10⁶ Hz 的 100 个校准/验证窗。逐文件生产 SHA256、元数据、光子及能量一致性核验通过；总波形和逐 primary 光子贡献在 1/2/5 ns 分箱相等。200 个最终测试 ROOT 未用于本轮分析。
- 以中央发射、剩余记录至少 5 μs、E_dep>0、Q_5us≥20 的条件取得 283 个校准响应与 163 个验证响应。逐事件首光子对齐后，1–5 μs 光输出占比均值 4.0327%；旧模板的源起始 1 μs 记录还包含飞行时间，不能把该均值直接当作每个旧事件的精确补偿系数。
- 线性/二次能标在低沉积能区偏差明显。增加实际冻结 α 累积光产额表乘有效比例的候选，仅在校准集拟合比例 η=0.78348。它与仿真共享物理输入，属于模型内标定；没有更改淬灭参数，也不是 SiPM PDE 测量或低能光产额的独立实验验证。
- 在三种能标共同支持的同一 157 个验证响应上，线性、二次、冻结表候选的能量 RMSE 分别为 80.57、72.31、45.28 keV。它们使用真值分离后的 Q，包含响应涨落，仅是条件能标误差，不能称为 pileup 能量分辨率。全部 163 个验证响应及各方法超范围样本另行保留；不将外推值夹到端点。
- 两个由真值筛出的独立 K=1 验证窗已用总波形做积分/模板检查。未知 K 的识别、假拆分/合并、计数率性能与保护区历史收敛尚未完成。1/2/5 ns 光子守恒也不代表重建性能对分箱已收敛。

当前代码见 [alpha_response.py](../pileup/analysis/alpha_response.py)。与本次结果配套的源码副本、划分清单、输入哈希和环境记录保存在同一结果目录。文献审核状态与项目计算状态分别维护。


<a id="alpha-gamma-context-calibration-20261016"></a>
### 2026-10-16：L27——低电荷、位置及前史标定的实际验证

- **复用依据与状态**：[Fixsen2004精读](../../papers/sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md)，PDF2 §2–3、PDF3 §4（TES模拟，条件平均形状/协方差标定，[DOI](https://doi.org/10.1016/j.nima.2003.11.313)）；[Rosado2019精读](../../papers/sources/2019_Rosado_IEEESensorsJ_Modeling-the-Nonlinear-Response-of-Silicon-Photomultipliers.md)，PDF2–4恢复和非线性响应、PDF5–7两种Hamamatsu器件与LYSO/LED验证，[DOI](https://doi.org/10.1109/JSEN.2019.2938018)。两篇均为“AI精读草稿，待人工审核”，无新增人工确认。仅采用检查条件响应和恢复效应的原则，不移植TES性能或SiPM平均响应为逐事件pileup指标。
- **本项目实现**：独立新增141个标定窗，覆盖低电荷α、三点位置及四种原生双脉冲顺序；新孤立γ只有8个可用例，拟合前完整合并旧56个独立标定窗，所有零信号保留。时间表接口只改变受控源的指定方式，物理和共享微胞/读出参数不改；逐事件祖先标签仅用于标定形状及评价。旧样本集中在中央，合并后的平均不能解释成连续均匀位置的精确平均。
- **方法假设**：新增均值标定；再试电荷/可观测前史依赖的慢分量及常数+1/Q散布约束，0.5/1两尺度。回归和惩罚形式为项目方案，非论文给出的GAGG算法。H是估计起点前电压的代理量，位置/真实微胞占用不输入盲拟合；AR(1)噪声沿用旧独立标定。系数原值/新值和适用范围在[Notebook第16节](../jupyter/codex_0005_Alpha_Gamma_Pileup_Fitting_20261015.ipynb#alpha-gamma-context-calibration-20261016)及[冻结calibration](../pileup/experiments/alpha_gamma_context_calibration_20261016/calibration.json)对照。
- **实际计算**：冻结后的16个新混合窗口中，expanded相对原V3的α完整恢复为374/656（57.01%）→378/656（57.62%），γ为30/64（46.88%）→31/64（48.44%）；额外预测55→44，漏检120→117，误合并标志58→56。α完整预测比例373/563→377/555，γ30/91→31/91，γ预测仍有很大误分类/误恢复负担。各条件并不一致：高α率轫致辐射背景的γ反而4/16→3/16，单能γ背景为2/17→3/17；两种高α率条件合计γ仍仅6/33（18.18%）。按条件分层、成对重采样窗口得到恢复率差的95%百分位区间：α为−0.45至+1.65个百分点，γ为−3.13至+5.80个百分点，均跨零。该区间只包含本次窗口波动，不含标定和物理系统误差，不能据此确认稳定提升。 保留扩展标定及全部失败候选，暂不替换V3默认流程。新增标定补齐了原低电荷α范围，但均值修正尚不足以解决密集混合波形：测试后分组中，500–1000 p.e.的α只恢复38/157；与最近另一Q≥50 p.e.脉冲相距小于50 ns的γ仍为9/36，其中18例被错分型；前史电压≥200 mV的γ仅3/21，其中13例错分型。这些是真值辅助的相关性描述，不能当成因果证明或物理时间分辨极限。两个前史形状候选在开发集的α恢复均从92/161降至83/161，γ仍5/12，因此未进入本轮独立测试。每条件仅4个混合窗，原生双脉冲验证中γ目标也只有2例。 开发按两类和四条件等权F1选择，另外27孤立/16原生双脉冲验证；16新混合窗在冻结后生成。完整漏检/误检/分型/条件Q偏差与散布及窗口bootstrap见[实验记录](../pileup/experiments/alpha_gamma_context_calibration_20261016/README.md)。Q为实现p.e.，不称能量分辨率。
- **验证/待办**：26项Python测试及C++ SiPM检查通过，209个新ROOT、220份拟合复核、一个完整电压盲重跑。下一步优先补齐强前史下的γ条件标定：本次仅4个有效γ前史标定，前史电压上限约172 mV，未覆盖上述≥200 mV失败组。可从独立标定库选取多份有效γ光子输入，与前脉冲在同一微胞阵列按新数字化种子重放，系统覆盖前史、间隔和位置；这属于条件响应检验，不能拿筛选后的样本估计γ探测效率。在独立重复样本上检验条件均值和协方差是否足以描述模板涨落，再在开发集比较类型推断与电荷分配，并用新的冻结后测试判断是否改善。这些是后续方案，本轮尚未实施；本轮16个新测试窗没有用于调参。 仍需处理位置/前史覆盖不足、标定不确定性、当前共享95/351 ns基底与真实器件差异；没有把本轮结果当成物理不可辨识极限。


<a id="alpha-gamma-ml-literature-20261016"></a>
### L28：2026-10-16——α＋γ混合波形的机器学习解堆积候选

**范围与追溯**：使用本地paper-search，只读复用8篇来源笔记，6篇定点核原页；没有联网、新训练或新仿真。
完整书目、DOI、装置/采样/噪声/率条件、原页链接、项目采用范围见[新Notebook](../jupyter/codex_0006_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb)和[evidence.json](../jupyter/run_records/alpha_gamma_ml_literature_20261016/evidence.json)。
本次检索没有核实同GAGG＋共享SiPM条件的α＋γ逐事件完整恢复基准；n/γ文献作为算法类比，不能转移其数值性能。

**采用方向与边界**：Kim2023（DOI 10.1016/j.apradiso.2023.110880，PDF3–6）支持共享CNN的脉高/类型多任务，但以已找到峰为前提；12 ns双脉冲MT分类84.80%，9重降至60.22%，式3/5口径有歧义，不采作能量准确率/召回证据。
Morad2024（10.1109/TNS.2024.3371573，PDF3–7、10）支持数量/到时检测与上下文分类分开；6 Mcps的89.1%是检出后分类，含漏检总准确率72.9%，FA14.1%另列，无逐事件电荷要求。
仅该范围已有[第553项人工限定确认](../../papers/global-review.md#human-review-553-morad2024-confirmed)，其余问题仍待审核；其他7篇保持AI草稿状态。
Deltoro2025（10.1007/s41365-024-01606-y，PDF4–7）支持双解码器分离对照，限250 MHz EJ301的合成nγ/γn、12–160 ns；0.988波形相关不是98.8%事件成功，类别保持约83–85%。
Regadío2021（10.1016/j.nima.2021.165403，PDF2–3、6）支持全波形到事件序列的U-Net方向，没有α/γ分型。
Yang2017（10.1109/TNS.2017.2691729，PDF2–3、5–6）直接做BaF₂ α/γ神经网络甄别但未解堆积，0.21%是符合γ误判α比例；仅用于孤立PSD基准。
Fu2018、Peng2022本轮仅复用笔记：分别保留尾部预测在<60 ns退化、组合分类排除<100 ns的条件，不能作为当前近邻弱γ已可恢复的依据。
Tamagawa2015（10.1016/j.nima.2015.05.052，PDF3）GAGG＋PMT的α138/599 ns与当前共用95/351 ns基底不同，提示域偏移，未更改物理参数。

**项目方案，尚未实施**：优先比较带α/γ事件、时刻、电荷输出的共享1D CNN/U-Net，以及网络/残差候选并集＋允许改类的联合物理拟合；双解码器先验证四种双脉冲与空/单脉冲。固定同一候选的PSD分类作诊断对照。训练需补弱γ、近邻及强前史，保留幅度并用Q匹配对照检查能谱捷径；使用原生共享微胞数字化，按ROOT/种子/父光子身份在增强前分组。α/γ等权开发，同时报告两类完整恢复率、precision/F1、条件电荷偏差/宽度、漏检、错分及误拆分/合并，不以更小波形误差替代事件成功。模板约束保留消融，不能假定混合拟合必优于纯网络。

**现有事实与待办**：L27扩展拟合的高α源10 MHz条件γ仅6/33完整恢复，其中19错分；<50 ns组9/36恢复，≥200 mV前史组3/21，后两组相互重叠、仅作关联。训练覆盖不足与模板失配仍需检验，不能宣称物理不可分。该16窗今后作为已公开回归数据，新的ML最终测试需在网络和阈值冻结后独立评估；源注入率不替代实际探测率，Q不替代沉积能量。本轮Notebook仅汇总旧JSON并生成图表，未给出混合ML性能。共享文献源/原PDF哈希和旧Notebook保持验证见[validation.json](../jupyter/run_records/alpha_gamma_ml_literature_20261016/validation.json)。

**新增文献整理需求（在papers项目处理）**：若精确复现，定点补Kim式3/5与峰检分母，Morad FA分母/NoA先验/基础脉冲划分，Deltoro表2与正文差异/基础事件独立性；同GAGG/SiPM原生α＋γ逐事件标签及模型失配证据仍缺。本项目只记录缺口，不修改库内审核或队列。


<a id="alpha-gamma-ml-20261016"></a>
### L29：2026-10-16—17——α＋γ多任务网络与物理拟合的实际比较

**采用与边界**：落实L28的A/B路线：扩张CNN、U-Net的事件/类型/Q/时间输出，以及神经候选＋可改类的联合物理拟合。文献装置条件、原页与审核状态沿用[L28](#alpha-gamma-ml-literature-20261016)，没有新增人工确认或论文性能外推；这是项目实现，未实现双解码器/幅度单独PSD对照。

**输入和独立性**：原生共享SiPM完整20 μs电压，2 ns采样；52混合＋141孤立/前史训练窗、20混合验证窗，按ROOT/种子先划分。新加12训练窗γ源50 MHz；全部方法冻结后另生成16测试窗（两种γ×α源1/10 MHz×4窗，γ源10 MHz）。两类和四条件等权F1选参数；16旧混合只作回归，43原生控制是开发资料。物理和动态范围假设不变。

**实际结果**：完整恢复定义仍为真值Q≥500 p.e.、时间≤10 ns、类型正确及Q误差≤20%，Q不是沉积能量。

| 方法 | α完整恢复 | γ完整恢复 | α完整precision | γ完整precision | 额外预测 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 扩展标定拟合 | 374/641（58.3%） | 35/75（46.7%） | 70.9% | 37.6% | 41 |
| 多任务CNN | 342/641（53.4%） | 25/75（33.3%） | 67.3% | 41.7% | 20 |
| 多任务U-Net | 375/641（58.5%） | 21/75（28.0%） | 69.3% | 51.2% | 27 |
| U-Net＋物理拟合 | 416/641（64.9%） | 28/75（37.3%） | 69.0% | 38.9% | 64 |

**采用决定**：未达到两类恢复同时改善，不替换已有拟合默认流程。混合相对拟合α+6.55、γ−9.33个百分点；按条件分层的成对窗口bootstrap95%区间分别+2.62至+10.49、−16.88至−1.23个百分点，只含有限窗口波动。弱γ500–1000 p.e.恢复19/49→10/49、错分15/49→26/49；前史≥200 mV的γ2/13→0/13，各分层重叠。高α率γ仍4/31，未解决密集γ恢复。匹配模板双脉冲80事件拟合80/80、CNN46/80、U-Net45/80、混合71/80，提示当前学习和组合仍有缺陷，不能解释成物理极限。

**追溯、验证与待办**：[冻结实验](../pileup/experiments/alpha_gamma_ml_20261016/README.md)、[Notebook第8节](../jupyter/codex_0006_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb#alpha-gamma-ml-implementation-20261016)、[评分/重放复核](../jupyter/run_records/alpha_gamma_ml_20261016/analysis-validation.json)。9项算法检查、28新有效ROOT、420份重评分、四方法电压盲重放；原生控制γ仅3个。保留序列化失败及VirtualBox受阻记录，按原种子在本机缓存完成训练后归档。下一步Q匹配及前史条件原生训练、固定候选PSD/幅度消融与先验校准，开发约束逐类恢复和误报，另用新测试；真实器件、物理域偏移和能量标定仍未验证，不修改共享文献库。

<a id="alpha-gamma-prior-ablation-20261019"></a>
### L30：2026-10-19——固定候选下拆分网络电荷与类别惩罚

**依据与范围**：沿用L28/L29中Kim2023、Morad2024的多任务及检出/分类分开评价方向和原适用条件；本轮只复用既有精读笔记与项目代码，没有重新核原PDF或新增论文性能结论。Morad第553项限定确认、其余AI草稿状态保持。将两个显式惩罚拆分是本项目的受控诊断，文献没有验证当前归一化网络分数为校准后验，也没有验证当前单个Q先验可跨类别共用。

**设计与实现**：20个已公开开发窗（两γ背景×α源1/10 MHz×5窗，γ源10 MHz），791个α/86个γ阈上目标。固定盲候选、初始类别/Q、到达时刻、模板Q及噪声权重，无到时搜索或增峰；仅比较(λQ,λC)=(0,0)/(10,0)/(0,10)/(10,10)。强度沿用旧开发选择，没有新调参。“无显式先验”仍有共同网络候选与初始形状，不能当作纯拟合。盲输入仅完整电压与保存的盲输出，真值在全部拟合后评分。

**实际计算**：无/仅电荷/仅类别/两者的α完整恢复为483/791、517/791、503/791、520/791；γ为42/86、37/86、42/86、38/86；额外预测53/83/59/83。仅电荷项使6个原成功γ失败（5个错分α、1个Q超差），补回1个；γ差−5.81个百分点，成对整窗分层重采样95%范围−10.34至−2.15个百分点。仅类别项α增加20个、γ净变化0（4失4得），γ完整precision37.2%→54.5%，但α完整precision略降、额外预测增加。区间仅描述已公开开发窗的波动，不作为新独立性能或实物结论。

**解释与采用决定**：当前把已解码事件的一个Q中心同时用于α、γ假设，Q惩罚能够影响分型。事后首个错分示例真值γ为554.14 p.e.、网络Q967.02 p.e.；无先验得584.67 p.e.的γ，仅电荷后得948.95 p.e.的α。该例用于解释全窗优化的变化，不把全部失败或物理可辨识性归因于一个因素。暂保留默认流程；下一轮优先检验完整流程中的候选/类别提示与物理Q重估，再做独立的类型/电荷/可观测前史条件Q校准与Q匹配训练，冻结后另留新测试。

**记录与验证**：[实验和复算](../pileup/experiments/alpha_gamma_prior_ablation_20261019/README.md)、[现有Notebook第9节](../jupyter/codex_0006_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb#alpha-gamma-prior-ablation-20261019)、[本地结果](../jupyter/artifacts/alpha_gamma_prior_ablation_20261019/summary.json)。4项数值测试、20份密集输出解码、80组共同候选/目标单调性及预定首窗四组重放纳入验证。全部满足至多6轮的停止容差，但不是全局最优保证。没有新G4、训练、概率校准或能量标定；ROOT和旧结果保留，共享文献库未写入。

<a id="alpha-gamma-full-prior-20261019"></a>
### L31：2026-10-19——仅类别惩罚的完整混合流程及新独立测试

**依据与边界**：沿用L28–L30的既有精读与原审核状态，没有新增原文复核或论文性能结论。拆分惩罚及筛选容差是项目方案；归一化网络事件分数仍不当作已校准后验。当前Q相关模板仍由初始网络/物理电荷决定，去掉显式Q项不等于去掉全部网络Q影响。

**实现与冻结**：仅λQ 10→0、λC保留10，候选/形状初值、固定噪声权重、时刻/类别/电荷搜索及残差补峰预算沿用旧完整流程。20公开开发窗相对物理拟合α471/791→505/791、γ38/86→42/86，额外预测45→57，达到事先设置的逐类恢复、precision与额外预测筛选。工程容差不作为文献阈值或统计等效性结论。其后冻结旧U-Net权重、阈值[0.4,0.55]、标定和拟合配置，另生成16新种子原生窗，两γ背景×α源1/10 MHz×各4窗，γ源10 MHz。物理与完整电压/足够动态范围前提不变。

**独立计算结果**：同窗物理拟合→仅类别混合，α345/643→382/643，γ30/79→29/79；额外预测40→54。完整恢复仍要求Q≥500 p.e.、匹配≤10 ns、类型正确、Q误差≤20%，不是沉积能量分辨率。本批独立数据未通过与开发相同的筛选，不替换默认拟合流程。未通过项：gamma_pooled_recovery, cell_brems_1000000_alpha, cell_brems_10000000_gamma。全部五方法、分条件、错分/漏检与条件Q误差保存在本地summary。- 相对扩展标定物理拟合的alpha恢复数变化+37，恢复率差+5.75个百分点；成对整窗95%范围[+4.10, +7.72]个百分点。 - 相对扩展标定物理拟合的gamma恢复数变化-1，恢复率差-1.27个百分点；成对整窗95%范围[-6.06, +3.70]个百分点。

**解释、验证与待办**：γ损失8个原成功事件、补回7个；弱γ为扩展标定物理拟合 16/52、旧混合：Q＋类别 13/52、新混合：仅类别 13/52。窗口重采样只含有限统计，各分层重叠；若干非线性子问题达到评价预算，不声称全局最优。3项数值检查、旧端点和冻结网络密集图/五方法数值复现、16窗重建文件访问检查、ROOT/RunInfo/原宏与重评分纳入验证。首次缺输出目录在生成ROOT前修复，日志保留；分析外层143中断后仅补跑缺结果的2窗，14个完整结果保持原哈希，不新增样本或调参。下一步在独立标定/开发资料中检验类型/Q/可观测前史条件的模板与不确定度，再设计新实验；此批测试不再用于调参。未验证真实器件、沉积能量与物理域偏移，不修改共享文献库。

**追溯**：[冻结实验](../pileup/experiments/alpha_gamma_full_prior_20261019/README.md)、[Notebook第10节](../jupyter/codex_0006_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb#alpha-gamma-full-prior-20261019)、[本地验证](../jupyter/run_records/alpha_gamma_full_prior_20261019/validation.json)。

**补充分型对照**：同批完整混合由无显式惩罚到仅类别提示，α367/643→382/643、γ32/79→29/79，额外预测49→54；丢失5个原成功γ且均错分α，补回2个。相对物理拟合损失的8个原成功γ也全部错分α。类别提示仍有取舍，下一步优先在训练/独立标定资料中做Q匹配PSD、类别分数校准和可观测前史/Q条件模板检查，不用本批测试选参。γ相对物理拟合的窗口区间跨零，不声称总体显著下降。

<a id="alpha-gamma-recovery-strategy-20261019"></a>
### L32：2026-10-19——针对弱γ错分的文献复核与后续方案

**检索与证据层级**：按`paper-search`只读共享库，复用主题页、完整精读笔记和global-review原始答复。7篇重点来源的题名、作者、DOI、真实笔记/PDF链接、装置/采样/噪声/率条件及逐篇核验范围集中于[Notebook第11节](../jupyter/codex_0006_Alpha_Gamma_Pileup_ML_Literature_20261016.ipynb#alpha-gamma-recovery-strategy-20261019)和[本地证据表](../jupyter/run_records/alpha_gamma_recovery_strategy_20261019/evidence.json)。DOI核PDF首页；实际原页图为Tamagawa3–4、Kim6、Morad10、Fixsen2、Klanner4，其余仅定点文字/既有笔记，不声明本轮重新全文精读。Morad第553项仍仅限检出/分类指标和构造条件，其余AI草稿状态不变。共享库未修改，未联网或追读引用。

**采用范围**：Tamagawa2015（10.1016/j.nima.2015.05.052，PDF3–4）用于总面积归一形状判别对照；99.95% γ保留/0.22% α误留限GAGG/PMT孤立源事件0.5–0.75 MeVee，不是本项目密集堆积成绩。Luo2018（10.1016/j.nima.2018.03.078，PDF4–6）用于多类型假设比较的先例，限BC-501A双脉冲；不采用其会残差抵消的绝对残差总和。Kim2023（10.1016/j.apradiso.2023.110880，PDF6表2/3/5）保留多任务取舍：12 ns分类MT84.80%低于ST87.91%，9重MT总体分类60.22%；原标Acc.H/Recall疑点不消除。Morad2024（10.1109/TNS.2024.3371573，PDF3–7、10）用于检出、上下文分类分开；6 Mcps条件分类89.1%与含漏检72.9%分列，不证明热图分数已校准。

Rosado2019（10.1109/JSEN.2019.2938018，PDF2–5、7）及Klanner2024（10.1016/j.nima.2023.169018，PDF3–4、9）支持SiPM恢复/光时间结构影响响应的机制，验证以指定器件平均电荷为主，不是GAGG逐事件解堆积；不直接移植参数或有疑点的公式。Fixsen2004（10.1016/j.nima.2003.11.313，PDF2–3）只借鉴条件均值/完整协方差标定，限TES模型的局部最小二乘；条件高斯log det及形状分布归一化是本项目拟议统计模型要求，不冒充该文已实现的波形pileup方法。

**项目核查与轻量计算**：重新汇总L31已公开79个γ，物理拟合为成功30/错类型35/正确类型但Q超差4/漏失10；仅类别混合为29/40/2/8。相对物理拟合丢失的8个成功γ全部变错类型；统计区间仍跨零，不宣称总体显著下降。当前混合只用Q条件模板且形状Q初值固定，不使用保存的context_classes；旧物理基线已有Q形状更新，旧context拟合也已试过，不能称本轮才首次想到。旧γ前史参数来自25个有效样本、最大172.20 mV，而公开失败分层有25个≥200 mV目标；这是标定覆盖缺口，不是器件上限。按当前目标，λC=10、假设scoreα/scoreγ=4，会给γ增加27.73白化SSE代价；这是代数量级，不是实测失败分数分布。

**拟议改进，尚未实施**：先固定盲候选，用观测减盲估邻居的残差做独立SI/白化模板/PSD对照；不对已经生成的固定类型模板再分型。用原生Q/可观测前史匹配资料检验幅度捷径，检出和分型输出分开，在独立按ROOT分组的标定集校准融合；同一电压的物理/网络分数不作为独立证据相乘。补弱γ强前史条件响应与不确定度，再检验Q/时刻/形状自洽的双假设比较。自由快慢比例会吸收类间差；若采用条件概率须处理形状分布归一化、条件协方差和形状涨落重复计入，而不是仅增加自由度或给旧经验惩罚机械补项。以上为本项目推导/方案，不是文献性能保证。

**验证与边界**：3个轻量Notebook单元实际执行，旧65单元/图像保留；格式、链接、计数恒等式、来源/输入哈希检查见[本地验证](../jupyter/run_records/alpha_gamma_recovery_strategy_20261019/validation.json)。没有改重建器、训练模型或生成ROOT，默认物理拟合不变。旧公开测试只作失败描述/回归，不择参；新方案冻结后再做独立共同测试，继续要求α/γ完整恢复和误报同时受控。

**新增文献/标定待办**：同GAGG/SiPM原生混合逐事件实测证据、弱γ强前史的条件形状与噪声分解、融合分数校准/重复计入控制仍缺。Kim式3/5、Morad FA分母/NoA先验、Klanner实现异文沿原记录保留；需要新增精读或代码核验时在papers项目另行处理，不从本项目修改共享库。

<a id="free-lifetime-minimum-window-20261019"></a>
### L33：2026-10-19——自由双衰减时间的最短观察窗口诊断

**依据与范围**：沿用L14–L16的完整响应卷积和原S13360两个孤立样例，不新增文献性能结论或物理参数。自由τ是当前电压模型的有效参数；GAGG/PMT平均波形的文献寿命不等于本器件逐事件可辨识精度。原文来源/审核条件保持，原始文献没有规定100/500 ns通用门槛。

**实际计算**：复用原ADC缓存、原约束和四初值，500–750 ns每25 ns扫描，新增22次拟合；两个粒子×500/750 ns共4个历史端点数值回归。T指观测阈值后长度，另用前200 ns且定位/基线依赖更长记录。按两τ均与同事件2.5 μs自由解相差≤5%/10%/20%、无触边且全部更长已测T均满足，两例共同最早已测试T分别为550/550/525 ns。门槛是本轮探索性定义，非置信区间、实物精度或连续意义下最短值；独立粒子仍只有1α+1γ。

**解释与边界**：500 ns自由拟合面积已与长窗相差约3%以内，但γ慢τ仍高约23%；固定τ的100 ns γ慢占比0.3335，相比2 μs的0.2047并不稳定。长窗自由解α80.36/348.05、γ81.41/221.63 ns仍不同于仿真输入95/351 ns，因此本节只回答相对长窗解稳定，未确定真实寿命准确恢复所需T。Poisson前向无新到达概率exp(−RdetT)仅作计数率联系，Rdet非源率、不代表恢复效率，不能由1/T给最大率；后续需独立分Q样本和受控间隔/高率联合拟合。

**追溯与验证**：[原Notebook第8节](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#free-lifetime-minimum-window-20261019)、[新汇总](../jupyter/artifacts/alpha_gamma_window_minimum_20261019/summary.json)、[协议与冻结源码](../jupyter/run_records/alpha_gamma_window_minimum_20261019/protocol.json)、[验证](../jupyter/run_records/alpha_gamma_window_minimum_20261019/validation.json)。原29单元/图像保留，未生成ROOT或替换旧结果；有效寿命与本征发光寿命、观察长度与最大率分开，未将同一事件多个截窗当独立样本。

**固定τ补充对照**：同日仅复算原9档固定95/351 ns结果，仍自由拟合A/B/t0/ws，没有新增拟合。相对各自2500 ns固定τ解，要求本点及所有更长已测点的数值均合格，两类共同最早已测档位为面积≤10%：200 ns、面积≤5%：300 ns、面积和ws各≤5%/10%：500 ns。100 ns的γ面积+11.49%、ws+62.88%，不能说固定τ后100 ns普遍足够。α的1000 ns拟合有基线上界触边；若严格加入上一节的无触边条件，共同起点为1500 ns。必须将这个基线/模型质量问题与数值稳定分开，不能把前述300/500 ns称为已通过全部质量条件的最短窗，也不能把1500 ns解释为物理必需长度。两事件、粗网格及长基线依赖不变，无新增文献或总体精度结论。[Notebook固定τ对照](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#fixed-tau-window-comparison-20261019)、[复算与两种标准](../jupyter/artifacts/alpha_gamma_fixed_tau_window_20261019/summary.json)、[验证](../jupyter/run_records/alpha_gamma_fixed_tau_window_20261019/validation.json)。

<a id="short-window-literature-strategy-20261019"></a>
### L34：2026-10-19——缩短拟合观察窗的方法与信息量诊断

**证据与审核范围**：按`paper-search`只读共享库，复用Belli2008、Radeka1968、Fixsen2004、Nakhostin2019 JINST、Yang2022、Saxena2020六篇完整已有笔记，定点核原页；实际看图Belli3–4、Radeka6–7、Fixsen2、Nakhostin4、Yang4、Saxena2，其他位置分清本轮文字/旧笔记层级。方法与适用条件见[Notebook第9节](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#short-window-literature-strategy-20261019)；完整题名、作者、DOI/未核状态、来源位置和审核记录见[证据表](../jupyter/run_records/alpha_gamma_short_window_strategy_20261019/evidence.json)。两次导航脚本45 s超时后用rg/索引，不据无返回说无文献。六篇仍为AI草稿；第554项认可的是另一篇Nakhostin TNS背景抑制，第658项不转移Fixsen；共享库和审核状态未改。

**采用与限制**：Belli（10.1016/j.nima.2008.06.045，PDF3–4图4–6）支持独立单脉冲截断补尾的误差验证设计；NE213固定7.5/45 ns、200 MS/s条件下7点35 ns/加缓冲50 ns峰距及均值/SD约±8%，不能移为GAGG最短窗。Radeka（DOI本轮未核，PDF6–7/13）用于有限窗口权重及噪声代价原则，串/并联白噪声/冲量信号前提保持；式16–17原印冲突不直接实现。Fixsen（10.1016/j.nima.2003.11.313，PDF2式1–2）用于独立均值/完整协方差标定原则，不作当前SiPM短窗性能证明。Nakhostin（10.1088/1748-0221/14/04/P04012，PDF4/9–10）的NaI离线795 kHz缩尾到20 ns后峰宽16%，100 ns窗含20 ns基线，噪声代价及式2.2/2.3不自洽保留，不照抄为已核代码。Yang（10.1016/j.nima.2022.166376，PDF2–5）的双平均窗尾部预测仍拒绝0.75 μs近邻；Saxena（10.1016/j.nima.2018.09.123，PDF2–5）的HPGe实时自适应成形为1 Mcps输入/370 kcps输出，0.7 μs平顶/2.4 μs复位等待不能照搬。两者只借鉴残尾估计/窗口自适应思想，不据此承诺100 ns逐事件准确恢复。

**项目推导与实际计算**：拟议固定τ后以Af/As线性幅度与少数非线性参数profile搜索、独立基线及条件协方差；相同模型/目标下参数改写不增加信息，不能单独保证精度改善。短窗先截协方差再求逆，不能截长窗逆矩阵；原定位使用整段峰和长前窗，不是严格因果前缀。另实现理想单指数Poisson过程的N/τ Fisher诊断：给定已知起点、任意示例N=10000/τ=351 ns，100 ns总量在τ固定时相对1σ界约2.01%，τ自由时总量约21.25%、τ约24.47%。这是本项目解析说明，非双指数GAGG/ADC下界或实测精度；28项独立积分、7次逆矩阵和2个短窗极限检查通过。

**下一步及未验证项**：先实现完整响应两幅度拟合与只依赖前缀的触发/基线，再检验条件均值/协方差和α/γ形状不确定度，缩尾作有正则的对照；逐事件双τ自由与跨事件共享τ分开。预定扫Tpre=20–200 ns、Tpost=50–1000 ns，按ROOT分标定/开发/新测试，真值只评价；电荷≤5%且事件合格率≥95%/两τ各≤10%只是拟议工程目标。报告两类检测/分类、偏差、分辨率、漏失/误拆/误并及选择效率，校准自适应停止后的覆盖率，不用两条波形的稳定性作“正确”证明。本轮没有新增GAGG拟合、训练或ROOT。[协议（尚未执行）](../jupyter/run_records/alpha_gamma_short_window_strategy_20261019/proposed-protocol.json)、[独立模块](../analysis/short_window_information.py)、[计算结果](../jupyter/artifacts/alpha_gamma_short_window_strategy_20261019/information.json)、[本轮验证](../jupyter/run_records/alpha_gamma_short_window_strategy_20261019/validation.json)。原41 Notebook单元及旧拟合器保留。待补本器件短窗独立原生精度、条件噪声/形状标定、因果总延迟，以及缩尾原式代码/勘误证据；新增精读在papers项目另行处理。


**2026-10-20笔记展示整理**：按用户要求，Notebook 0002聚焦双指数模型，移除单衰减拟合对照及理想信息量诊断的正文／图9-1，合并重复叙述并按研究问题排序；保留其余图号和双指数数值。历史推导、原始结果与验证仍保留在上述独立模块、产物及运行记录中。第9节的第一阶段方案已在下方L35执行；此次整理不新增拟合、文献采用或物理结论。

<a id="short-window-implementation-20261019"></a>
### L35：2026-10-19——因果短前缀拟合的实现与独立原生检验

**文献依据与采用边界**：落实[L34](#short-window-literature-strategy-20261019)的有限窗补尾、线性幅度profile和独立条件均值/协方差原则；Belli NE213固定7.5/45 ns、Radeka理想噪声、Fixsen TES的装置/采样/模型条件与AI草稿审核状态均沿L34，不引入新的文献性能或人工确认。Radeka/Nakhostin原式疑点未消除，未照抄缩尾公式；本轮没有修改共享papers库。α/γ完整模板、条件pilot、50%收缩和停止规则为本项目实现，不冒充论文给出的GAGG性能。

**实现与独立性**：149个旧孤立标定窗、45个旧开发窗，70个Q≥250 pe（52α/18γ）用于三档类别均值/726维协方差；基线均值和电子噪声独立标定，σB=0.2 mV及50%对角收缩为工程设定。四方法是固定95/351 ns两幅度OLS/GLS及独立完整α/γ形状OLS/GLS，输入同一完整原生analog_mV。触发用过去滚动中位数的连续三点2 mV门槛，窗口只读[触发−Tpre,触发+Tpost]；先截协方差再分解，未使用未来峰值/尾部定位。固定类别完整形状比固定双τ且自由比例约束更强，参数改写本身不提高信息量。

配置、标定、4个核心源码哈希冻结后，新增174个独立原生20 μs单事件ROOT（30α/144γ），预定包含γ扩展并保留零信号；单能γ0.1/0.2/0.661657/1.5 MeV、α1.5/1.8/2.5/5.486/8 MeV，三点位置。真值只在该事件44个窗口全部拟合完成后打开，不作为推断输入。Q≥500 pe目标数为α19、γ12；这是实现光子后代雪崩电荷，不是源/沉积/可见能量。

**实际结果**：主方法template_gls。开发选择的前20＋后300 ns候选在新测试通过本批点比例门槛（α19/19（100.0%），γ12/12（100.0%））；目标为两类各≥95%事件Q误差≤5%。测试网格仅作探索描述：前20 ns：后300 ns；前50 ns：后300 ns；前100 ns：后300 ns；前200 ns：后300 ns。小样本点比例和逐条件Wilson区间不能认证总体95%门槛；同事件多个截窗不当作独立样本，未据测试重选算法/停止阈值。前20 ns的关键比较如下。

| 触发后T / ns | 方法 | α电荷误差≤5% | γ电荷误差≤5% |
| ---: | --- | ---: | ---: |
| 100 | physics_ols | 3/19（15.8%） | 0/12（0.0%） |
| 100 | physics_gls | 4/19（21.1%） | 1/12（8.3%） |
| 100 | template_ols | 14/19（73.7%） | 6/12（50.0%） |
| 100 | template_gls | 14/19（73.7%） | 7/12（58.3%） |
| 200 | physics_ols | 10/19（52.6%） | 4/12（33.3%） |
| 200 | physics_gls | 10/19（52.6%） | 7/12（58.3%） |
| 200 | template_ols | 15/19（78.9%） | 3/12（25.0%） |
| 200 | template_gls | 18/19（94.7%） | 12/12（100.0%） |
| 300 | physics_ols | 11/19（57.9%） | 6/12（50.0%） |
| 300 | physics_gls | 18/19（94.7%） | 10/12（83.3%） |
| 300 | template_ols | 17/19（89.5%） | 7/12（58.3%） |
| 300 | template_gls | 19/19（100.0%） | 12/12（100.0%） |
| 500 | physics_ols | 11/19（57.9%） | 3/12（25.0%） |
| 500 | physics_gls | 19/19（100.0%） | 11/12（91.7%） |
| 500 | template_ols | 14/19（73.7%） | 8/12（66.7%） |
| 500 | template_gls | 19/19（100.0%） | 12/12（100.0%） |
| 1000 | physics_ols | 15/19（78.9%） | 6/12（50.0%） |
| 1000 | physics_gls | 19/19（100.0%） | 12/12（100.0%） |
| 1000 | template_ols | 18/19（94.7%） | 10/12（83.3%） |
| 1000 | template_gls | 19/19（100.0%） | 12/12（100.0%） |

**停止与负对照**：开发最坏窗口误差得到k=54.45149，残差阈值2.10599，固定前50 ns；新测试α接受0/19、γ接受0/12，拒绝保留为无条件失败。局部σ和经验放大不能视作已校准顺序置信区间。24个标定均值相加的双脉冲（72次前缀检查）只检查单脉冲假设遇邻居的失效，不是新原生共享微胞双脉冲恢复验证。

**验证与限制**：4项有独立数值参照/未来污染的算法检查通过；174个ROOT的实际RunInfo、宏、种子、配置、哈希与祖先电荷核验，1408个汇总组独立重评分，4次原生未来NaN/平移检查及2次快照源码重放通过。[Notebook第10节](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#short-window-implementation-20261019)保留原50单元和输出；[实验](../pileup/experiments/alpha_gamma_short_window_20261019/README.md)、[测试汇总](../jupyter/artifacts/alpha_gamma_short_window_implementation_20261019/test-summary.json)、[数据验证](../jupyter/run_records/alpha_gamma_short_window_implementation_20261019/native-data-validation.json)。本轮假定至多一个孤立脉冲、不估计数量，尚不能报告盲多脉冲误拆/误并，也不回答free双τ最短窗。旧两条ADC的长窗数值稳定与本轮独立模拟真值精度口径不同，不直接比较缩短倍数。待补有效γ/条件前史/标定不确定度、原生双脉冲、跨事件共享τ与freeτ、实物域偏移及高率联合拟合；不以1/T替代计数率能力。新文献整理仍在papers项目处理。

<a id="advanced-fit-literature-20261020"></a>
### L36：2026-10-20——短窗口GLS之后的条件响应、形状模态与标定统计

**问题与核验层级**：用户询问是否有更好的拟合方法。按`paper-search`只读本地库，复用导航、完整来源笔记及global-review原始答复；两次knowledge_index导航各45秒无输出超时，改用rg，不以无返回判断无文献。五篇正式期刊来源与逐篇实际文本／原页范围、文件哈希、作者、DOI或未核状态、设备／能区／采样／噪声／率条件汇总于[证据表](../jupyter/run_records/alpha_gamma_advanced_fit_literature_20261020/evidence.json)。这是定点证据核对，不是五篇新全文精读或入库；共享库未写入，未联网。过程与候选追加到[笔记0002第9.1–9.2节](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#advanced-fit-literature-20261020)，原代码输出和图号保持。

| 来源 | 实际依据与采用范围 | 限制／审核 |
| --- | --- | --- |
| D. J. Fixsen, S. H. Moseley等，2004，*Pulse estimation in nonlinear detectors with nonstationary noise*；[精读](../../papers/sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md)／[原文](../../papers/raw/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.pdf)；[DOI](https://doi.org/10.1016/j.nima.2003.11.313) | PDF2式1–2（原页实看）、PDF3文字；条件均值和逆协方差随标定能量插值，用于提出自洽条件拟合。 | TES模型0.5–5 eV、每点5000次，采样／率未给；原式无log det。AI待审核，无本篇人工确认。 |
| Chris H. Whitford，2005，*Implementation of optimal filters for pulse height measurement from non-linear detectors with non-stationary noise*；[精读](../../papers/sources/2005_Whitford_NIMA_Implementation-of-optimal-filters-for-pulse-height-measurement-from-non-linear-detectors-with-non-stationary-noise.md)／[原文](../../papers/raw/2005_Whitford_NIMA_Implementation-of-optimal-filters-for-pulse-height-measurement-from-non-linear-detectors-with-non-stationary-noise.pdf)；[DOI](https://doi.org/10.1016/j.nima.2005.09.001) | PDF2式4–5、PDF3图3（两页实看）、PDF4式9–10文字；局部dS/dE与端点自洽插值，并保留时刻／基线代价。 | TES 2–20 keV；20 keV改善19%相对未知相位峰取值，仍差下界14%；40 eV是指定2 keV间距的能标偏差。无堆积或硬件验证，AI待审。 |
| S. E. Busch, J. S. Adams等，2016，*Progress Towards Improved Analysis of TES X-ray Data Using Principal Component Analysis*；[精读](../../papers/sources/2016_Busch_JLTP_Progress-towards-improved-analysis-of-TES-X-ray-data-using-principal-component-analysis.md)／[原文](../../papers/raw/2016_Busch_JLTP_Progress-towards-improved-analysis-of-TES-X-ray-data-using-principal-component-analysis.pdf)；[DOI](https://doi.org/10.1007/s10909-015-1357-z) | PDF3与6原页实看，PDF5／7文字；少量波形变化分量及能量方向标定，作为结构化形状／协方差候选。 | Mo/Au TES、1 MS/s；Mn Gaussian FWHM3.19±0.14→1.56±0.13 eV，Al无改善；仅第658项限定采用，独立划分／权重与4317/4984计数差仍待核。 |
| Roger Barlow, Christine Beeston，1993，*Fitting using finite Monte Carlo samples*；[精读](../../papers/sources/1993_Barlow_CPC_Fitting-using-finite-Monte-Carlo-samples.md)／[原文](../../papers/raw/1993_Barlow_CPC_Fitting-using-finite-Monte-Carlo-samples.pdf)；DOI未核 | PDF2式8–9原页实看，PDF3条件最大化文字；有限MC模板统计进入辅助似然，借鉴到有限标定误差传播。 | 原文为Poisson分箱计数，不能直接用于相关电压；模板均值profile不等于边缘积分。DOI未核、PII不改写为DOI；原式7等疑点保留，AI待审。 |
| Jaime Rosado，2019，*Modeling the Nonlinear Response of Silicon Photomultipliers*；[精读](../../papers/sources/2019_Rosado_IEEESensorsJ_Modeling-the-Nonlinear-Response-of-Silicon-Photomultipliers.md)／[原文](../../papers/raw/2019_Rosado_IEEESensorsJ_Modeling-the-Nonlinear-Response-of-Silicon-Photomultipliers.pdf)；[DOI](https://doi.org/10.1109/JSEN.2019.2938018) | PDF4式19原页实看；PDF2–3／5／7–8文字；SiPM恢复、触发概率与电荷改变，为条件响应提供物理动机。 | S13360-1325CS/1350CS，LYSO3×3×20 mm／LED，300–2100 keV、400 ns电荷积分；仅平均电荷／电流模型，非6050PE波形逆算法，AI待审。 |


**与当前实现的差别**：已核`analysis/alpha_gamma_short_window.py`，现有流程已有线性幅度／基线profile、连续时刻搜索、三Q档类别模板，以及`C=Qpilot*K(tpilot)+sigma_e²I`和50%对角收缩；不是完全固定且与Q无关的协方差。具体缺口是模板和C按单次pilot冻结，拟合Q／t0不再自洽更新。有效标定70条、γ每档6条与726维协方差提示应优先验证统计稳定性及标定范围。L26已有形状试验α170/290→173/290、γ18/38不变；不把增加自由度写成已经改善。

**采用的研究建议／尚未验证**：优先比较自洽条件均值／噪声，随后比较开发集选择的低秩结构与有限标定不确定度。若采用参数相关条件高斯似然，`rᵀ C(theta)⁻¹r + log det C(theta)`是本项目统计建模，不冒充Fixsen原式；旧同前缀共用固定C的log det为常数。PCA旋转不自动白化；形状模态若已在残差C内，不重复添加；有限标定误差会被多事件共享。固定形状Q恢复、带标定先验／共享τ和逐事件完全free双τ分开验证，不能由300 ns小样本100%推出自由寿命最短窗或高率能力。

预定先用开发数据比较100/150/200/300 ns，在公开旧测试上只回归；方法冻结后另设新测试，逐类报告Q误差、分型、漏检、接受度与区间覆盖，保留失败。未实现新拟合器、未重新拟合、未生成ROOT，尚不能排序本项目真实性能。Scoullar2008会议PPR仅复用既有笔记／原修订记录，公开似然及事件数选择细节不足，未作为可直接复现的首选。

**文献整理待办**：如后续采用，先另行确认相关电压的有限模板联合统计模型、类别／Q／时刻条件协方差、PCA能量方向与分量选择、自由双τ可辨识／覆盖。Barlow DOI和原印疑点、Whitford详细TES依赖等仍待核；补充文献整理在papers项目按独立授权处理，不在本项目修改共享库或自动追读。


<a id="advanced-fit-implementation-20261020"></a>
### L37：2026-10-20——自洽条件似然、低秩协方差与标定均值误差的实现

**文献采用与边界**：落实[L36](#advanced-fit-literature-20261020)中Fixsen/Whitford的条件响应思想、Busch的低维变化方向及Barlow–Beeston的有限模板统计原则，原探测器/采样/噪声/能区/率条件和审核状态全部沿L36，没有新文献性能或人工确认。本轮不重读或修改共享库。参数相关高斯目标`rᵀC⁻¹r+log det C`、残差对角的采样解释和32次分层ROOT bootstrap高斯均值误差是本项目方案；不是Fixsen原式或Poisson Barlow–Beeston的原样电压实现。PCA替代相关协方差，不新增无限制形状幅度，不重复计入同一涨落。均值标定误差跨事件共享，当前未完整重采样分档、协方差、器件系统误差。

**实现和分组**：七条流程含原标定GLS、扩充标定GLS、自洽密集条件似然、PCA8/PCA16和两种均值误差组合。只拟合至多一个孤立脉冲的未知类型/Q/连续时刻/基线，完整双分量发光物理保持。基线解析profile及0.2 mV共同软约束；模板随当前log Q插值，协方差随当前Q/时刻更新，先截前缀再分解。新增420个原生20 μs单粒子ROOT，4个预检查计入标定；168标定＋84开发＋168冻结后测试，保留所有零信号γ，不拒绝ADC触顶。α源1.5/2.5/5.486/8 MeV、单能γ源0.1/0.2/0.661657/1.5 MeV、三点位置，保持原物理/SiPM与2 ns analog_mV输入。

加入旧149标定和45开发后，为317/129窗；有效标定Q≥250 pe从70增加到113，γ各档从[6, 6, 6]增加到[14, 14, 13]。旧174个已公开测试不用于择法或本轮评分。开发按两类等权、100/150/200/300 ns等权的Q±5%且类型正确选中`expanded_gls`，前窗20 ns；开发95%点比例规则选择后300 ns。主短窗目标150 ns预先固定，方法/标定/源码冻结后才生成168个新测试；按ROOT/三个种子隔离，真值在全部方法/窗口推断后打开。

**独立结果**：目标为实现光子后代雪崩Q≥500 pe，而非源/沉积/可见能量。本批测试α16、γ15个目标；保留零电荷γ108个。150 ns联合成功由原标定GLS的alpha16/16、gamma11/15，经扩充标定GLS的alpha16/16、gamma9/15，到所选方法的alpha16/16、gamma9/15。完整关键窗口如下，同事件截窗不当作独立样本。

| 触发后T / ns | 同批方法 | α Q±5%且类型正确 | γ Q±5%且类型正确 |
| ---: | --- | ---: | ---: |
| 100 | legacy_gls | 14/16 | 10/15 |
| 100 | expanded_gls | 12/16 | 10/15 |
| 150 | legacy_gls | 16/16 | 11/15 |
| 150 | expanded_gls | 16/16 | 9/15 |
| 200 | legacy_gls | 15/16 | 15/15 |
| 200 | expanded_gls | 15/16 | 15/15 |
| 300 | legacy_gls | 16/16 | 15/15 |
| 300 | expanded_gls | 16/16 | 15/15 |
| 500 | legacy_gls | 16/16 | 15/15 |
| 500 | expanded_gls | 16/16 | 15/15 |

150 ns配对比较：alpha对legacy_gls差+0.00个百分点、配对ROOT重采样95%范围[+0.00,+0.00]；alpha对expanded_gls差+0.00个百分点、配对ROOT重采样95%范围[+0.00,+0.00]；gamma对legacy_gls差-13.33个百分点、配对ROOT重采样95%范围[-40.00,+13.33]；gamma对expanded_gls差+0.00个百分点、配对ROOT重采样95%范围[+0.00,+0.00]。范围仅表示指定模拟设计下有限ROOT统计，不是含标定/器件系统误差的总体保证。所选方法在新测试网格中两类≥95%且所有更长档通过的探索性最早点为后300 ns；未据此重新选择方法或窗口，小样本100%不能认证总体≥95%。电荷偏差/SD为匹配条件统计，联合成功和覆盖保留漏检为失败。

**区间与失败控制**：所选方法局部名义95%区间：alpha后150 ns覆盖14/16、有效区间16/16；gamma后150 ns覆盖10/15、有效区间15/15；alpha后300 ns覆盖7/16、有效区间16/16；gamma后300 ns覆盖12/15、有效区间15/15。区间条件于当前类别和模型，不声称已覆盖类别切换/模型失配或顺序停止。所选方法在150 ns、108个零信号窗上有0个触发、0个Q≥500 pe误报。优化器状态、起点/Q触边、无效区间和分Q结果全部保留。另有24条标定均值相加＋白噪声双脉冲、48个前缀负对照，只评分第一个脉冲，未模拟共享微胞非线性或恢复第二个脉冲，不是原生盲解堆积。

**验证和结论适用性**：9项独立数值参照检查、4个旧开发端点回归通过；420个ROOT的RunInfo/实际宏/三个种子/配置/原文件哈希/祖先电荷核验，280组独立重评分，4组原生未来NaN/时移/缓存重放和2组冻结源码重放通过。原短窗源码及全部历史结果哈希保持，Notebook原51单元的代码/输出/图像保留，新增图11-1/11-2/11-3。只验证同G4/SiPM模型、有限能区/位置的孤立固定形状Q；未给free双τ最短窗、实物精度、原生高率误拆/误并或最高源率，不用1/T替代计数率能力。下一步仍需条件噪声/标定系统误差、跨事件共享τ、原生双脉冲及混合高率的独立开发与测试。

**复现与记录**：[实验及冻结源码](../pileup/experiments/alpha_gamma_advanced_fit_20261020/README.md)、[Notebook第11节](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#advanced-fit-implementation-20261020)、[完整报告](../jupyter/artifacts/alpha_gamma_advanced_fit_20261020/report.json)、[验证](../jupyter/run_records/alpha_gamma_advanced_fit_20261020/native-data-validation.json)。文献追加整理仍在papers项目按独立授权处理。


<a id="shared-tau-implementation-20261020"></a>
### L38：2026-10-20——跨事件共享双τ与独立标定约束的实现

**文献与项目方案边界**：延续L36–L37中Fixsen2004/Whitford2005条件响应与统计标定的有限借鉴范围；原TES设备、模拟/采样/率条件、缺口及AI待审核状态保持，不新增人工确认或文献性能。本轮跨事件共享τ、log τ高斯约束、逐ROOT线性化bootstrap与慢分量比例分类器是本项目统计方案，不声称原论文已经给出此GAGG方法或普适短窗长度。共享papers库未修改，也未新增原文核验。

**方法改变**：保持双分量发光和原生SiPM物理不变；同样的单位面积双指数与读出/孔径卷积，快慢幅度非负，基线/时刻逐事件拟合。比较固定95/351 ns、独立标定联合共享τ、逐事件带log τ约束、逐事件free双τ，完整类别模板GLS另作外部对照。四种物理方法共用同一前缀的电压pilot冻结协方差，先截矩阵后求解；variable projection消去两幅度/基线。free仍有fast25–250 ns、slow250.5–2000 ns和触发相对时刻−20至＋2 ns范围，不叫无界拟合。原有模型/结果不变，新图12-1至12-4。

**独立标定与选择**：沿用317标定、129开发，不使用之前174或168个公开测试作本轮择参/评分。113条Q≥250 pe标定在触发后1000 ns联合估计共享τ＝91.743114/364.822935 ns，每事件各自幅度/基线/时刻；Q与类型真值只用于独立标定样本和后置分类器，事件推断不读真值。该τ属于有效电压响应模型，不直接等于材料发光寿命。2000次逐ROOT profile-score/information线性化bootstrap近似标定log τ协方差，并加0.02 log标准差下限；不是完整bootstrap refit或全系统误差。约束宽度1/3/10倍仅由开发按两类等权、100/150/200/300 ns联合Q±5%且类型正确选择，冻结为3倍。分类器按各流程独立标定的三Q档慢比例均值/SD，先拟合后分型，未把真类型输入事件拟合。

**新测试与实际结果**：冻结源码/标定/物理快照后生成168个原生20 μs单粒子ROOT（24α/144γ），源α1.5/2.5/5.486/8 MeV、γ0.1/0.2/0.661657/1.5 MeV，三点位置、5 μs发射、2 ns analog电压；保持所有零信号γ和完整动态范围，不拒绝削顶。仍假设至多一个孤立脉冲，不估计数量。下表以实现光子后代雪崩Q≥500 pe为分母，不是能量分辨率；同事件各截窗不作独立样本。

| 触发后T/ns | 方法 | α Q±5% | γ Q±5% | α Q＋类型 | γ Q＋类型 |
| ---: | --- | ---: | ---: | ---: | ---: |
| 100 | fixed_tau | 3/17 | 3/12 | 3/17 | 3/12 |
| 100 | shared_tau | 5/17 | 2/12 | 5/17 | 2/12 |
| 100 | constrained_tau | 5/17 | 1/12 | 5/17 | 1/12 |
| 100 | free_tau | 3/17 | 2/12 | 3/17 | 2/12 |
| 100 | template_gls | 11/17 | 8/12 | 11/17 | 8/12 |
| 150 | fixed_tau | 5/17 | 6/12 | 5/17 | 6/12 |
| 150 | shared_tau | 7/17 | 5/12 | 7/17 | 5/12 |
| 150 | constrained_tau | 9/17 | 3/12 | 9/17 | 3/12 |
| 150 | free_tau | 5/17 | 1/12 | 5/17 | 1/12 |
| 150 | template_gls | 15/17 | 8/12 | 15/17 | 8/12 |
| 200 | fixed_tau | 8/17 | 7/12 | 8/17 | 7/12 |
| 200 | shared_tau | 9/17 | 4/12 | 9/17 | 4/12 |
| 200 | constrained_tau | 13/17 | 4/12 | 13/17 | 4/12 |
| 200 | free_tau | 1/17 | 0/12 | 1/17 | 0/12 |
| 200 | template_gls | 15/17 | 10/12 | 15/17 | 10/12 |
| 300 | fixed_tau | 16/17 | 9/12 | 16/17 | 9/12 |
| 300 | shared_tau | 15/17 | 6/12 | 15/17 | 6/12 |
| 300 | constrained_tau | 15/17 | 5/12 | 15/17 | 5/12 |
| 300 | free_tau | 2/17 | 3/12 | 2/17 | 3/12 |
| 300 | template_gls | 16/17 | 11/12 | 16/17 | 11/12 |
| 500 | fixed_tau | 17/17 | 11/12 | 17/17 | 11/12 |
| 500 | shared_tau | 17/17 | 11/12 | 17/17 | 11/12 |
| 500 | constrained_tau | 17/17 | 11/12 | 17/17 | 11/12 |
| 500 | free_tau | 9/17 | 8/12 | 9/17 | 5/12 |
| 500 | template_gls | 17/17 | 12/12 | 17/17 | 12/12 |
| 1000 | fixed_tau | 17/17 | 12/12 | 17/17 | 12/12 |
| 1000 | shared_tau | 17/17 | 12/12 | 17/17 | 12/12 |
| 1000 | constrained_tau | 17/17 | 12/12 | 17/17 | 12/12 |
| 1000 | free_tau | 16/17 | 10/12 | 16/17 | 9/12 |
| 1000 | template_gls | 17/17 | 12/12 | 17/17 | 12/12 |

**寿命与最短窗的解释**：固定/shared的τ来自设定/独立标定，不记为单事件寿命恢复。原生free τ与G495/351 ns的偏差混合有限统计、不可辨识、优化触边及前向响应失配；不能仅据长窗收敛认证材料寿命正确。每类各≥95%并且所有更长档合格的探索性最早网格、分Q结果、偏差/SD/RMSE、类型错误、漏失、触边及局部区间覆盖全部保存；不据测试重新选择窗口/方法，小样本点比例不认证总体95%。

**控制与验证**：另有24个匹配均值模型标定及48个独立测试，单位面积双指数卷积＋0.1 mV白高斯电子噪声，无光子/微胞涨落；Q为均值面积参数，观察参考时刻预设，不混作原生实现Q或真实触发延迟。这一对照分开评价Q±5%和两个逐事件τ±10%。另24条均值相加双脉冲、48个前缀只检查单脉冲模型失效，不恢复第二个事件，也无共享微胞非线性。9项独立数值检查、168个ROOT的RunInfo/三个种子/实际宏/配置/祖先电荷、360组独立重评分、4组原生未来NaN/时间平移和2次冻结源码重放通过。

**采用范围和待办**：本次如实保存成功与失败，不把共享τ当作每事件free测量或认为约束必然改善。后续需拆开未观察尾部的随机涨落、模型失配与类别退化，校验真实SiPM前向响应/占用依赖及标定系统误差，再研究原生混合高率；逐事件能量标定、实物、盲多脉冲误拆/误并和最高源率未验证，不用1/T替代计数率能力。新文献精读/入库仍在papers项目按独立授权处理。

**复现**：[Notebook第12节](../jupyter/codex_0002_GAGG_SiPM_Alpha_Gamma_Fitting.ipynb#shared-tau-implementation-20261020)、[实验及冻结源码](../pileup/experiments/alpha_gamma_shared_tau_20261020/README.md)、[完整结果](../jupyter/artifacts/alpha_gamma_shared_tau_20261020/report.json)、[独立验证](../jupyter/run_records/alpha_gamma_shared_tau_20261020/native-data-validation.json)。


**2026-10-20笔记整体复核与整理**：Notebook0002将当前结论、三批独立样本和评价定义前置，区分两事件长窗稳定性、原生实现电荷精度、有效电压τ及匹配均值对照；把第9节旧待实施计划更新为第10–12节已完成／尚未完成状态，精简重复解释，修正条件均值和约束逆矩阵公式排版。全部32个代码单元、保存输出与17个固定图号保持。此次仅整理展示，不新增拟合、仿真、文献采用或科学性能结论。


<a id="pure-alpha-fit-methods-20261020"></a>
### L39：2026-10-20——纯α盲拟合的算法核对与方法升级依据

**问题与记录**：用户针对0003询问如何在拟合方法上提高事件恢复、是否均为blind。完整复核与候选方案追加至[Notebook0003第16节](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#pure-alpha-fit-methods-20261020)，不修改原71个单元、计算或图像。主结果为已知α/独立固定双τ模板下的未知K/t/Q重建；图12-2(a)使用已知K=2及到达时间，须明确属于条件分析。

**实际实现**：`alpha_pileup_fit.py`已有逆滤波候选、全窗非负幅度联立、初值±12 ns内循环调时、≤64候选时的联合非线性优化、一轮残差补峰、8 ns合并和500 p.e.输出阈值；未加权残差。第14节改用完整含噪电压，仍沿用基线，真值在推断后才打开评分。本轮不把工程门槛或低残差当物理分辨率/正确K的证明。

**只读证据与审核范围**：复用L01–L03、L23–L26及现有完整来源，定点提取Adrich PDF2–3、Zhang PDF3–4、Lopatin PDF3–4、Fixsen PDF2–3；实际看Adrich PDF2和Lopatin PDF3，其他两篇本轮为定点文字回查。B23归并无新增科学确认，四篇AI待审状态不变；没有联网或修改共享库。真实来源、哈希和本轮证据层级见[结构化记录](../jupyter/run_records/pure_alpha_fit_literature_20261020_170045/sources.json)。

- **单/双模型比较**：[Adrich与Zychor2021笔记](../../papers/sources/2021_Adrich_NIMA_Least-squares-fitting-algorithm-for-peak-pile-up-correction-in-gamma-ray-spectroscopy.md)、[原文](../../papers/raw/2021_Adrich_NIMA_Least-squares-fitting-algorithm-for-peak-pile-up-correction-in-gamma-ray-spectroscopy.pdf)、[DOI](https://doi.org/10.1016/j.nima.2020.164962)，PDF2式1–2/流程；LaBr₃/PMT、约80 kcps实测，采样/噪声/接受阈值缺项。可借鉴主动试拆隐藏双脉冲；p值不是模型后验，不能移植间隔效率。
- **条件模板与残差增峰**：[Zhang等2024笔记](../../papers/sources/2024_Zhang_NIMA_Development-of-a-pile-up-pulse-recovery-algorithm-for-the-LaBr3-detector.md)、[原文](../../papers/raw/2024_Zhang_NIMA_Development-of-a-pile-up-pulse-recovery-algorithm-for-the-LaBr3-detector.pdf)、[DOI](https://doi.org/10.1016/j.nima.2024.169273)，PDF3–4/图2–4；LaBr₃/R6233、500 MHz，同模板σ=7道三幅度模拟，率10⁶–10⁸ s⁻¹。GAGG须独立标定；式5局部坐标和接受准则缺口保留。
- **两阶段稀疏**：[Lopatin等2012笔记](../../papers/sources/2012_Lopatin_IEEE_Pileup-attenuation-for-spectroscopic-signals-using-a-sparse-reconstruction.md)、[原文](../../papers/raw/2012_Lopatin_IEEE_Pileup-attenuation-for-spectroscopic-signals-using-a-sparse-reconstruction.pdf)，PDF3式4–6/算法1；会议论文、DOI未核，HPGe/Cs-Am实测率/采样缺项。选点与估能分开；项目可用无L1非负幅度重拟合，但不称原文逐式复现。式5同下标冲突继续保留。
- **完整协方差加权**：[Fixsen等2004笔记](../../papers/sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md)、[原文](../../papers/raw/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.pdf)、[DOI](https://doi.org/10.1016/j.nima.2003.11.313)，PDF2式1–2；TES、0.5–5 eV模型、每能量5000次，采样/噪声配置不完整且未验证逐事件pileup。只采用标定均值/协方差的思想，不当PCA或包含logdet的完整似然。

**待验证**：优先在0003纯α中对照局部K增删/拆并、多初值联合重拟合与独立空/单样本校准的复杂度代价；再消融稀疏候选、阈下干扰分量、相关噪声和条件形状。保持双τ结构和无真值盲输入。混合任务已有相关实现，但不能替代纯α独立同条件验证。当前36窗均已检视，仅作回归；冻结后留新测试。保持完整动态范围及500 p.e./10 ns/20%口径，同时报告漏检、额外预测、误拆分/合并、电荷/时间误差。本轮没有修改拟合器、生成新ROOT或取得性能提升证据；不将建议写为已验证结果。


<a id="pure-alpha-local-model-selection-20261020"></a>
### L40：2026-10-20——纯α局部脉冲数搜索的实际实现与独立测试

落实[L39](#pure-alpha-fit-methods-20261020)的主动拆峰、局部K比较与多初值联合时间/非负电荷重拟合；保持95/351 ns和原慢权重，完整含噪电压输入，无K/t/Q或源率真值输入。以独立48个单α窗最大拆峰SSE收益/Q标定代价系数0.07574128，24个开发单α及三个已公开1/10/40 MHz窗选择factor=1，于2026-10-20T17:34:07.113436+08:00冻结后生成23个新原生测试窗。四篇文献条件/审核边界沿用L39；局部搜索/代价不是原论文逐式复现或校准的统计检验。

八档源率0.1/1/10/20/40/60/80/100 MHz分别8/3/2/2/2/2/2/2窗，20 μs/窗、中央10 μs。原生共享微胞、2 ns采样、0.1 mV噪声不变；Q≥500 p.e.、±10 ns和±20%完整恢复。逐率结果、额外预测、漏检及条件电荷误差见[Notebook0003第17节](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#pure-alpha-local-model-selection-20261020)及[完整结果](../jupyter/artifacts/alpha_local_model_selection_20261020/test-summary.json)。非均衡曝光合计4215个目标，原→新完整恢复476→691、额外预测246→284、漏检2413→2071；不能替代逐率结果或窄区间推断。

另有174条匹配均值软件控制，仅电子白噪声，不当实物时间分辨率。保留开发首版源码及26个完成结果；局部列更新在冻结前通过五项数值检查。新[实验快照](../pileup/experiments/alpha_local_model_selection_20261020/README.md)含实际C++/配置/二进制身份、Python源码、宏和种子；[验证记录](../jupyter/run_records/alpha_local_model_selection_20261020/validation.json)核对冻结时序、ROOT/宏/RunInfo、46组重评分/模型及原Notebook输出。未改变物理、做GLS/条件形状消融或完成沉积能标；每个高率仅2窗，局部搜索不保证全局最优，低残差不等于正确事件数。

100 MHz首次因处理候选数保护上限10000000而中止，保留两份不完整ROOT/日志；在独立alpha_local_capacity_20261020快照仅将计算上限提高至30000000，以相同1 MHz种子的16数组逐项一致验证后重试原计划种子。此工程容量调整不改变物理/读出概率或冻结拟合参数，验证重放不计入测试分母；已知失败不混入恢复率。

一份100 MHz ROOT的主光子时间分支无法解压；以WaveformEvent的冗余时间及唯一的x/y/能量浮点位精确一对一连接恢复评分字段，要求覆盖全部光子。另一完整100 MHz文件的16数组等价验证及恢复前后4个总波形数组相同见[恢复验证](../jupyter/run_records/alpha_local_model_selection_20261020/redundant-time-validation.json)。原ROOT保留且该原分支仍不可直接读取，评分重导出使用保存的适配器；盲拟合输入、源码及参数不变。


<a id="new-pileup-literature-p11-p12-20261020"></a>
### L41：2026-10-20——新增P11/P12文献对盲堆积恢复的适用性复核

**请求与范围**：用户要求查看新增文献。索引P11新增11篇、P12新增15篇；另有补充材料/手册并入既有来源。按当前0003盲纯α拟合问题完整回读10篇来源笔记，定点回查7篇14页PDF文字，实际看Cleland PDF8、Mohammadian PDF3、Kim PDF6、Zhao PDF3；其余16篇只作索引筛选，未声称重读26篇原文。记录追加到[Notebook0003第18节](../jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb#new-pileup-literature-p11-p12-20261020)；[完整证据](../jupyter/run_records/new_pileup_literature_20261020_211656/review.json)、[路径/哈希及页码](../jupyter/run_records/new_pileup_literature_20261020_211656/collection.json)。共享库只读，无联网/下载/入库操作。

**审核范围**：10篇均为AI草稿待人工审核，无新增科学确认。已回查global-review第668–670项原答复“1—3均认可”：Molina2026仅拒绝工程证据；Ma2022/Kim2026仅条件化谱计数变化；El-Tokhy2021仅方法参考且定量保留。它们是另4篇的限定采用，不扩为整篇、事件恢复证明或对本轮10篇的认可。Molina拒绝标志不作为解堆积输出。

**证据与本项目判断（数值均为作者报告/已有笔记，不是本轮新计算）**：

- **W. E. Cleland等1994**：*Signal processing considerations for liquid ionization calorimeters in a high rate environment*；W. E. Cleland、E. G. Stern；[笔记](../../papers/sources/1994_Cleland_NIMA_Signal-processing-considerations-for-liquid-ionization-calorimeters-in-a-high-rate-environment.md)、[原文](../../papers/raw/1994_Cleland_NIMA_Signal-processing-considerations-for-liquid-ionization-calorimeters-in-a-high-rate-environment.pdf)、[DOI](https://doi.org/10.1016/0168-9002(94)91332-3)。PDF7–9，式21–44、图4。SSC液氩量能器设计；约400 ns漂移、16 ns束团间隔，已知脉形；串联热噪声与堆积背景协方差；理论与Monte Carlo，非GAGG实测。 用完整协方差逆矩阵加权；一阶时间线性化下，最优滤波与加权χ²拟合相等。 主要抑制目标束团之外的随机堆积背景，没有恢复未知K个事件；最优性只属于给定模板/协方差及线性无偏约束。 项目判断：优先方法依据；多事件扩展和实际收益待本项目验证。 DOI仅据本地题录，未独立外核。

- **M. Nakhostin等2010**：*A digital method for separation and reconstruction of pile-up events in germanium detectors*；M. Nakhostin、Zs. Podolyak、P. H. Regan、P. M. Walker；[笔记](../../papers/sources/2010_Nakhostin_RSI_A-digital-method-for-separation-and-reconstruction-of-pile-up-events-in-germanium-detectors.md)、[原文](../../papers/raw/2010_Nakhostin_RSI_A-digital-method-for-separation-and-reconstruction-of-pile-up-events-in-germanium-detectors.pdf)、[DOI](https://doi.org/10.1063/1.3499241)。PDF3–5，式1–2、图2–7；本轮定点文字回查PDF4。同轴Ge/Cs-137；75 MHz、14 bit，约14 μs记录；慢成形2.4 μs、快速定时支路40 ns。实测单脉冲人工叠加及约100 kHz下硬件拒绝器漏过的事件；噪声幅度未报。 两个起点由CFD估计后固定，双分量拟合只估幅度；典型600/200 keV对严格大于250 ns的扫描点幅度误差小于1%。 起点是由波形估计，不是真值；幅度拟合条件于K=2及已估时刻。不是未知多事件盲恢复率，250 ns点不满足上述精度，也不是GAGG最短窗。 项目判断：定时/电荷误差分解的诊断依据；不直接移植成形常数或250 ns。

- **Hong-Quan Huang等2017**：*Estimation method for parameters of overlapping nuclear pulse signal*；Hong-Quan Huang、Xiao-Feng Yang、Wei-Cheng Ding、Fang Fang；[笔记](../../papers/sources/2017_Huang_NST_Estimation-method-for-parameters-of-overlapping-nuclear-pulse-signal.md)、[原文](../../papers/raw/2017_Huang_NST_Estimation-method-for-parameters-of-overlapping-nuclear-pulse-signal.pdf)、[DOI](https://doi.org/10.1007/s41365-016-0161-z)。PDF3–5式10–16/流程，PDF5–7表1–4，PDF8结论；本轮回查PDF3、8。Sallen–Key数值二/四脉冲；200 MHz、5 ns；共同输入τ=100 ns，白噪声标准差5 counts；二脉冲另测共同τ在95–105 ns；预设到时，无真实随机到达率。 把整段重叠波形的幅度、到时、共享形状参数编码为染色体，以无权残差搜索。 实例按已知二/四分量编码；逐减M并取最小残差的选阶没有复杂度惩罚/盲计数验证；没有全局最优证明、完整随机流能谱或实时结果。 项目判断：值得做小规模全局搜索对照，回答是否主要困于局部极小值；不直接全窗高维GA。

- **Mohammad-Reza Mohammadian-Behbahani等2023**：*Non-iterative pulse tail extrapolation algorithms for correcting nuclear pulse pile-up*；Mohammad-Reza Mohammadian-Behbahani；[笔记](../../papers/sources/2023_Mohammadian-Behbahani_NET_Non-iterative-pulse-tail-extrapolation-algorithms-for-correcting-nuclear-pulse-pile-up.md)、[原文](../../papers/raw/2023_Mohammadian-Behbahani_NET_Non-iterative-pulse-tail-extrapolation-algorithms-for-correcting-nuclear-pulse-pile-up.pdf)、[DOI](https://doi.org/10.1016/j.net.2023.08.011)。PDF3式3–5/图4，PDF6表2–4。2×2英寸NaI(Tl)/PMT成形器输出，24 MHz、典型宽约15 μs；Cs源距80/20 cm两组约1 s离线波形；真实输入率/噪声强度未报；保留前峰至下一起点之间干净尾段的后半。 在策略1高堆积数据中NIM 17.95 s、默认MATLAB单指数NLS 1020.10 s，约57倍；两者校正谱相近。 约57倍是该实现/单指数尾窗下的离线比，不是恢复准确度或对当前NNLS的速度比。式3反向累加与式4–5的符号/有限尾端常数不自洽，未获作者澄清；近上升沿重叠缺干净尾段。 项目判断：作为速度与窗口失配的证据；目前不优先直接实现原印NIM。

- **Alberto Regadío等2024**：*Three topologies of deep neural networks for pulse height extraction*；Alberto Regadío、J. Ignacio García Tejedor、Luis Esteban、Sebastián Sánchez-Prieto；[笔记](../../papers/sources/2024_Regadio_ARI_Three-topologies-of-deep-neural-networks-for-pulse-height-extraction.md)、[原文](../../papers/raw/2024_Regadio_ARI_Three-topologies-of-deep-neural-networks-for-pulse-height-extraction.pdf)、[DOI](https://doi.org/10.1016/j.apradiso.2023.111137)。PDF3式1–3，PDF7表1；PDF8–10图14–19据既有笔记。固定CR–RC差指数1/0.1 μs，20 MHz、1024点；归一均匀幅度、白/布朗噪声及裁剪模拟；率与时间轴原印不一致，不能换算为可信MHz上限。 输出幅度冲激序列；MMSE取整段标签统一平移k=0…5后的最小MSE。表1中GRU 9569参数、CNN 3274305参数，MMSE分别0.48/0.30×10^-4；MSE指标下排序不同。 不是相同参数/计算预算比较；损失不等于事件效率/真实能谱分辨率。5×50 ns=250 ns是允许标签延迟的换算，非最小峰距；不能在评分时用真值逐事件重新对齐。 项目判断：补充ML基准与损失消融的依据；不是已经证明可恢复当前40–100 MHz。

- **W. Kim等2023**：*Optimizing deep learning-based piled-up pulse height correction method for high radiation-field application*；W. Kim、K. Ko、S. Lee、J. Park、G. Song、K. Lim、G. Cho；[笔记](../../papers/sources/2023_Kim_JINST_Optimizing-deep-learning-based-piled-up-pulse-height-correction-method-for-high-radiation-field-application.md)、[原文](../../papers/raw/2023_Kim_JINST_Optimizing-deep-learning-based-piled-up-pulse-height-correction-method-for-high-radiation-field-application.pdf)、[DOI](https://doi.org/10.1088/1748-0221/18/12/C12002)。PDF6式2.1–2.2/窗口，PDF8表2；最终窗口/表3见PDF9，据既有笔记。1.5×1.5英寸LaBr3，500 MS/s、2 ns/点；百万实测单脉冲合成1–10分量；训练PTP 16–80 ns；噪声随实测单脉冲携入，未单列SNR/真实连续输入率。 表2最短4 bins=8 ns，检峰计数比50.10%→52.50%；12 bins=24 ns为100%，16 bins=32 ns降至96.44%，并非单调。最终CNN12_16/DNN16_16均16点=32 ns，但这是候选峰附近脉高回归窗口。 RR是检出峰数/真峰总数，缺一对一匹配及假阳性统计；式2.2实际为相对误差形式、无1−，与约99%且越高越好的accuracy表述冲突。漏检峰不能由单脉高回归自动找回。 项目判断：优先用于评价设计及候选消融；不以99%作为我们的事件恢复目标已获文献保证。

- **Guang Zhao等2024**：*Peak finding algorithm for cluster counting with domain adaptation*；Guang Zhao、Linghui Wu、Francesco Grancagnolo、Nicola De Filippis、Mingyi Dong、Shengsen Sun；[笔记](../../papers/sources/2024_Zhao_CPC_Peak-finding-algorithm-for-cluster-counting-with-domain-adaptation.md)、[原文](../../papers/raw/2024_Zhao_CPC_Peak-finding-algorithm-for-cluster-counting-with-domain-adaptation.pdf)、[DOI](https://doi.org/10.1016/j.cpc.2024.109208)。PDF3式1–3/图3，PDF4图4–6/表1。He/isobutane漂移室电子峰，1 GS/s、12 bit、15-bin归一片段；伪域噪声10%→20%、上升2→4 ns，另有180 GeV/c μ子实测；无GAGG/电荷回归。 伪数据目标域AUC由直接迁移0.878升至半监督0.912，ideal为0.926；利用高置信CWT正例及低幅弱负例对齐域。 AUC不是逐事件恢复效率；实测无逐峰真值。逐片段标准差归一丢失幅度尺度；式2–3源/目标索引待核；样本划分和超参缺失。 项目判断：保留迁移方向，后置于同模型盲恢复与可验证失配测试。

- **Qiuyi Wang等2026**：*Adaptive nuclear pulse processing method with dual-model architecture*；Qiuyi Wang、Hongquan Huang、Zhiyong Lan、Xingke Ma；[笔记](../../papers/sources/2026_Wang_RDTM_Adaptive-nuclear-pulse-processing-method-with-dual-model-architecture.md)、[原文](../../papers/raw/2026_Wang_RDTM_Adaptive-nuclear-pulse-processing-method-with-dual-model-architecture.pdf)、[DOI](https://doi.org/10.1007/s41605-025-00606-9)。PDF10–16，式24–25、表2–4，图16–17见PDF17–19；本轮仅完整既有笔记。单个三指数响应的理想归零截断/S–K成形模拟，5 ns/点；辨参1500点、修复500点；A=300 mV附近、白噪声σ=5（单位不全）；未给实测率。 CNN-LSTM辨参后按漂移选择BiLSTM修复分支；15%/25%漂移算例匹配模型较固定模型好。 third-order是三个响应指数，不是三事件；Δ≥30%丢弃。15%/25%算例预设已识别参数，不是端到端全漂移保证；无未知K事件分离。 项目判断：模型失配辅助资料；不作为当前恢复率提升证据。

- **唐琳等2023**：*Truncated pulse height estimator based on an improved UNet model / 基于改进UNet模型的截断脉冲高度估计器*；唐琳、周爽、廖先莉、刘泽、李波；[笔记](../../papers/sources/2023_Tang_NuclTech_Truncated-pulse-height-estimator-based-on-an-improved-UNet-model.md)、[原文](../../papers/raw/2023_Tang_NuclTech_Truncated-pulse-height-estimator-based-on-an-improved-UNet-model.pdf)、[DOI](https://doi.org/10.11889/j.0253-3219.2023.hjs.46.110505)。PDF3–7，式1–6、图3–7、表1–2；本轮仅完整既有笔记。FAST SDD/X射线、开关复位/CR输出；正文20 MHz/50 ns与另处50 MHz冲突；合成128点，最短示例20点；噪声加于阶跃，实测8组离线。 UNet-LSTM示例20脉冲平均相对误差2.31%；全部200验证序列2.18%；实测8组均值2.36%为作者报告。 截断比例不是截断深度；合成结果来自验证集，实测参考高度来源不清；没有随机堆积检出/误拆分/恢复能谱或最短窗扫描。 项目判断：短窗/截断问题的边界资料；不是目前双τ盲解堆积方法的直接替代。

- **P. Windischhofer等2023**：*Passive quenching, signal shapes, and space charge effects in SPADs and SiPMs*；P. Windischhofer、W. Riegler；[笔记](../../papers/sources/2023_Windischhofer_NIMA_Passive-quenching,-signal-shapes,-and-space-charge-effects-in-SPADs-and-SiPMs.md)、[原文](../../papers/raw/2023_Windischhofer_NIMA_Passive-quenching,-signal-shapes,-and-space-charge-effects-in-SPADs-and-SiPMs.pdf)、[DOI](https://doi.org/10.1016/j.nima.2022.167627)。PDF3–8，式6–29、图4–9；本轮仅完整既有笔记。硅SPAD一维/1.5维解析与MC；厚0.5 μm、直径10 μm、约16 fF/200 kΩ算例；没有多胞/前端完整响应或pile-up率验证。 近击穿线性化、大雪崩、绝热且快速放电忽略恢复等条件下得到两倍电压摆幅；皮秒淬灭与RqCd恢复时间不同。 不是本项目S13360的独立参数测量，也不含完整共享微胞占用/串扰/前端；不能据此把单胞电荷乘2或发光寿命改成皮秒。 项目判断：保留物理模型边界，当前参数不变。

**下一步候选，尚未实现**：先在0003的现有双τ/局部K搜索中对照完整相关残差加权，固定同候选/开发校准后比较；若C依赖拟合参数，应保留高斯近似的logdet且不能称普通NNLS已精确解完整似然。用少量K=1/2/3受控局部的遗传/群体搜索判断局部极小值影响，多种子/相近计算预算；不是全局最优证书。GRU/允许延迟的损失可作补充消融，保留事件定时和候选召回；现有CNN/U-Net/Transformer/物理重拟合不重复算新实现。域适应、漂移修复后置；非迭代单指数尾法不替换GAGG双发光衰减。保持完整动态范围与无真值主推断，真K/时刻只作标注oracle。

**待办/缺口**：NIM积分方向/离散边界、Kim accuracy定义/峰匹配、Regadío率与时间轴、Zhao索引/弱标签与划分、真实响应/相关涨落标定；需新开发/冻结后独立测试，已公开23窗不再作未见测试。保留Q≥500 p.e./±10 ns/±20%、漏检/额外预测/误拆误并、电荷偏差与宽度及时间平移负对照，Q不当沉积能量。新精读/库维护留在papers项目按独立授权处理。本轮没有训练、拟合、生成ROOT、改变物理参数或取得新的恢复性能；原Notebook88单元与图输出保持。

<a id="alpha-correlated-global-gru-20261021"></a>
## L42. 相关残差加权、局部种群搜索与轻量GRU的实际对照（2026-10-21）

**研究问题与采用依据**：实施L41提出的三个方向，沿用其定点原文证据及审核边界，不重复整理精读报告。Cleland/Stern（1994，PDF7–9式21–44、图4；DOI `10.1016/0168-9002(94)91332-3`仅本地题录）启发全协方差权重，其液氩/已知模板结论不等于未知K恢复；Huang等（2017，PDF3、8及原笔记表1–4，DOI `10.1007/s41365-016-0161-z`）启发有预算种群搜索，但原文没有全局最优证书或完整随机流未知K验证；Regadío等（2024，PDF3/7、原笔记PDF8–10，DOI `10.1016/j.apradiso.2023.111137`）启发GRU对照，本轮未复现CR–RC模型或允许标签延迟的MMSE。三篇的完整sources路径和适用条件沿用L41；均保持AI草稿/待人工审核，既有其他论文的审核认可不转移。

**本项目假设和实现**：固定95/351 ns双发光衰减、慢分量0.740174；纯p–¹¹B有效单α谱，无γ。源注入率0.1–100 MHz，不是探测/反应率。原生共享微胞20 μs、2 ns完整模拟电压，中央[5,15) μs评分，动态范围充分，无ADC整窗拒绝。

采用经验$C_{ij}=s_i s_j\rho^{|i-j|}+\sigma_e^2\delta_{ij}$、$s_i^2=v_d+\eta\max(\mu_{pilot,i},0)$。C由原始观测pilot整窗固定，不读真值/源率，不随候选K/时间/Q改变，logdet因此为常数。innovation白化后非负幅度变量投影，局部边缘协方差提出方案、全记录固定C接受。K增删代价在单粒子标定重新定标，并非χ²显著性；与旧SSE差异不能只归于C一项。对角WLS为非对角相关项的消融。

48个已知单α标定中33窗有观测候选，得σe=0.1 mV、vd=0 mV²、η=0.7976401192 mV、ρ=0.9819124013。vd落在回归非负边界，原生2 MHz暗计数仍开；残差包含模板失配，非真实微观SiPM协方差。白化后相邻相关0.403、二阶矩1.353，保留模型不充分的证据。24个单α及1/10/40 MHz三个开发窗按预定误报/F1规则在0.5/1/2中均选0.5；没有在新测试上重选C或K代价。

局部DE＋SLSQP每窗最多12次DE，预算分布于采集时间；盲流水线实际提出一/二分量，已知K=2/3诊断单列，未声称未知K三分量全局恢复。GRU 7435参数，原CNN的26训练/10验证ROOT、五特征与损失、2400步/batch16；第2400步、阈值0.85由旧验证选择，训练327.5秒。架构/种子不同，非等训练墙钟；离线滑动窗口，不宣称实时硬件延迟。

原划分的100 MHz仅1个训练ROOT、无同率验证，最高验证率为80 MHz；100 MHz结果按外推报告，没有利用新测试重选阈值。

**新计算结果**：共同冻结`2026-10-21T08:24:47.916277+08:00`后生成23个新ROOT，与标定/开发、神经训练/验证及已公开测试身份/种子隔离。0.1/1/10/20/40/60/80/100 MHz各8/3/2/2/2/2/2/2窗。目标Q≥500 pe，时间≤10 ns且Q相对误差≤20%才完整恢复。Q为实际祖先归属雪崩等效电荷，不是沉积能量。

| 方法 | 完整恢复/目标 | 漏检（无时间匹配） | 额外输出 | 合并/拆分标志 |
|---|---:|---:|---:|---:|
| 原始SSE | 486/4333 | 2480 | 253 | 871/0 |
| 局部K-SSE | 726/4333 | 2118 | 282 | 978/7 |
| 对角WLS | 776/4333 | 2043 | 269 | 1010/4 |
| 相关GLS | 636/4333 | 2243 | 114 | 962/1 |
| GLS＋局部DE | 635/4333 | 2244 | 116 | 962/1 |
| 旧CNN | 1038/4333 | 1448 | 604 | 960/113 |
| 轻量GRU | 982/4333 | 1837 | 199 | 1006/38 |


同一批23个新原生窗中，旧局部K、对角WLS、相关GLS、GLS＋局部DE、旧CNN和轻量GRU的完整恢复数依次为726、776、636、635、1038、982（共同目标4333个）。这是非均衡曝光总计，逐率结果及额外输出须同时阅读。 条件Q偏差/标准差、时间误差、逐真值ID新增/损失事件、墙钟在同一本笔记中；并发任务影响墙钟，不作严格硬件吞吐比较。100 MHz实际/五种循环错位平均完整匹配为：局部K-SSE 4.73% / 2.79%；相关GLS 3.49% / 2.31%；GLS＋局部DE 3.35% / 2.33%；旧CNN 17.60% / 10.09%；轻量GRU 15.93% / 7.04%。错位保留输出数量/电荷，不是精确FDR、无偏效率校正或真实性能上限。

216个匹配模板盲软件窗含白/独立相关噪声的空/单/双脉冲，双间隔4/8/12/20/40/80/200 ns、比0.125–4，各2次；没有原生光子/共享微胞涨落，不用于宣称探测器分离极限。36个已知K条件诊断中DE在33/36降低目标；K=2全部分量完整恢复1/18→13/18，评估次数944→4163；K=3全部分量完整恢复1/18→7/18，评估次数1515→6174。DE额外评估，且K已知，不进入盲效率表。

**追加相同工作量/多随机种子对照**：在查看上述粗初值诊断后，另存新协议/源码SHA，再生成36个新的合成噪声案例，每例3个优化种子。已知K及固定触发区域，两策略均含同三组粗初值，严格各500次实际白化NNLS变量投影评估；一方随机多初值SLSQP，另一方将公共初值后的剩余预算70%给DE、再抛光和重启填满。K=2：随机多初值SLSQP 36/54，DE混合44/54；K=3：随机多初值SLSQP 20/54，DE混合26/54。每K为18噪声×3优化种子，54次不当54个独立噪声样本；不等于严格相同墙钟，非法时间向量的开销另记。补充只比较已知K条件精度，没有改变原生拟合/阈值、训练、种子或测试评价。结果在optimizer-budget-controls.json，源码另存supplemental-search-source。单元检查以独立调用计数核验实际profile预算、返回目标及固定种子复现。

**记录与验证**：本地 `jupyter/codex_0003_Alpha_Pileup_Rate_Fitting.ipynb` 第19节、固定图19-1至19-5，原94单元/图输出不变；原生索引同步本地Run_Catalog。实现位于 `pileup/analysis/alpha_correlated_*.py`，实际源/宏/二进制/标定/权重快照和无需Git切换的重跑入口在 `pileup/experiments/alpha_correlated_global_20261020/`。结果与独立验证分别位于本地 `jupyter/artifacts/alpha_correlated_global_20261020/` 和 `jupyter/run_records/alpha_correlated_global_20261020/`。

**I/O故障单独记录**：19个G4任务全部正常结束后，60 MHz第二个ROOT的全局光子x分支解压失败。由另一份PrimaryPhoton记录按共同时间/y/能量完整值唯一一对一映射16,534,845条记录恢复重复x字段，原ROOT不改写，4个电压数组逐位不变。新增适配只用于评分导出，缺唯一键/双份字段均坏/覆盖不足即拒绝；6项恢复/拒绝测试通过。checked hash加速导出在正常0.1/10/40 MHz的48个数组上与原导出器逐位相同。I/O源码/验证另存io-source及io-recovery.json，不追写原推断冻结；原种子、物理和拟合均不改变，剩余四个预定80/100 MHz由2进程继续、生成后顺序导出。23个预定窗口全部保留，不因算法表现换种子或删测试。

通过6项针对性数值检查、白化快慢核及完整拟合逐位等价；重新计算161份保存评分、69个模型/加权代价和祖先电荷守恒，核验新旧源码SHA、冻结时序/间隔和数据隔离，复读0.1/10/100 MHz实际ROOT，禁用NumPy文件加载的冷盲拟合及冻结源码新路径复跑。原ROOT保留，读取使用SHA相同普通本机字节副本，无临时链接。

**限制与待办**：高率每档仅2短窗、低率目标少；相关模型不足、固定均值模板/局部K预算、电荷选择效应和高密度偶然匹配仍存在。原pilot初始候选>500则退出、局部新增最多160，属于计算预算，不是物理上限；本批初始候选保护触发0/23窗，状态逐窗保留。不宣布最高可用率/通用最小间隔/沉积能谱分辨率；free-τ最短窗及α/γ混合不由本节回答。后续需在新开发协议研究更符合原生涨落的C及未知K候选/联合搜索，再冻结并另作独立测试；这轮测试已公开，不能继续用来选参。NIM、Kim accuracy、Regadío率轴等L41原文缺口仍待文献库独立维护，本轮未改papers。

<a id="alpha-ml-refresh-20261021"></a>
### L43：2026-10-21——纯α机器学习与当前WLS重新比较

**采用依据与边界**：沿用L21的Regadío2021/Wang2025架构方向及L40–L42的局部K、相关噪声与GRU来源。共享库只读，原器件、采样、噪声、计数率和人工审核范围均保留；没有新精读或将其他装置的性能移植为GAGG能力。新混合设计是本项目针对旧候选缺失、电荷先验和边缘外推的待验证假设，不是论文已证明的组合算法。

**原值→本轮实现**：旧CNN＋物理采用CNN候选、神经电荷惩罚λ=0.1、±6 ns调整和8 ns最终合并；本轮保留旧流程作对照，另实现观测pilot＋CNN候选并集、取消神经电荷幅度惩罚、采用当前WLS的±24 ns局部调整／3 ns分量间隔及残差增峰、拆分、合并、删除。CNN的262点边缘上下文只约束神经候选，完整电压、保护区与物理pilot边缘分量仍参加拟合。250 p.e.候选门沿用旧混合，报告门500 p.e.。原观测pilot固定经验方差，K代价和预算继承独立标定／开发选择；未根据本轮测试再选参。这些改变作为整体实验，没有单项消融。

**数据与实际运行**：CNN、U-Net、Transformer、GRU的权重和阈值不变，没有重训。原26训练／10验证ROOT保持分离；10验证窗只诊断。已有23个公开窗口、4333个目标作共同复评，不是第二批独立重复。源码／配置／权重冻结后生成16个新确认ROOT，八档源注入率0.1/1/10/20/40/60/80/100 MHz为8/2/1/1/1/1/1/1窗。纯p–¹¹B有效单α谱、固定95/351 ns、共享微胞、20 μs记录、2 ns采样、中央[5,15) μs评分；完整模拟电压与充分动态范围前提保持，未改物理参数或ADC政策。

**本项目新确认结果**：共同目标为实际雪崩等效电荷Q≥500 p.e.，时间±10 ns最大基数一对一匹配且Q相对误差≤20%。Q不是沉积能量。

| 方法 | 完整恢复/目标 | 无时间匹配漏检 | 已匹配电荷失败 | 额外输出 |
| --- | ---: | ---: | ---: | ---: |
| 原SSE | 258/2098 | 1178 | 662 | 104 |
| 当前WLS | 419/2098 | 943 | 736 | 112 |
| CNN | 518/2098 | 673 | 907 | 270 |
| U-Net | 524/2098 | 791 | 783 | 214 |
| Transformer | 459/2098 | 823 | 816 | 107 |
| GRU | 472/2098 | 868 | 758 | 83 |
| 旧CNN＋物理 | 513/2098 | 695 | 890 | 216 |
| 新CNN＋残差WLS | 446/2098 | 865 | 787 | 133 |


**判断**：完整恢复增加伴随额外输出增加，属于取舍，未达到两项同时改善；不据此自动替换默认方法。相对旧混合，完整恢复513→446，额外输出216→133。总计受高率事件主导；逐率、漏检、条件Q偏差与宽度、误合并／误拆分标志和循环时间错位共同报告。额外输出可能源于阈下信号或错误拆分，不直接等同“此处完全无信号”。

**实际验证与溯源**：461份保存评分重算、265个模型重构，16个实际ROOT的RunInfo／宏／三个种子及冻结后时序核对；原权重和推断源码不变，训练／验证／公开／新确认身份隔离。原pilot方差与新混合的logdet最大差0。66个匹配均值模板软件控制、8幅全窗全部分量求和、3项新模块检查、禁止评分真值读取的盲预检查及冻结源码8方法冷复跑通过。候选总预算触发0窗，pilot失败0窗；预算不是物理上限。原Notebook14个代码单元与十幅图保留，新过程位于本地`jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb`第15–21节，ROOT和实际宏另记Run_Catalog。

来源：[新冻结实验](../pileup/experiments/alpha_ml_refresh_20261021/README.md)、[全部实际表](../jupyter/artifacts/alpha_ml_refresh_20261021/report-tables.json)、[检查](../jupyter/artifacts/alpha_ml_refresh_20261021/validation.json)。推断副本可独立重跑到新输出目录，无须Git回退；本地模型／数据依赖仍须保留。

**限制与后续**：高率每档仅一窗、低率目标少，100 MHz没有同率神经验证；移位重复不是新ROOT，条件电荷宽度不是能量分辨率。没有背景／形状失配、free-τ最短窗、实物或通用最高可用率结论。后续训练、候选或损失修改使用开发／验证数据另定协议，并留新测试；不在这16窗上反复挑参数。新文献整理仍由papers库另行维护，本项目不修改该库。

<a id="alpha-ml-literature-20261021"></a>
### L44：2026-10-21——更适合逐事件解堆积的机器学习方法复核

**问题、检索与核验层级**：围绕0004的漏检、电荷分配失败和额外输出，只读共享库。依AGENTS/index、vocabulary/topic-index/review-index及三份优先主题页定位；复用只读knowledge_index检索、rg补查来源题名与正文。回读下列五篇完整来源笔记，定点提取S1 PDF4–7、S2 PDF3、S3 PDF10–11；实际看S1 PDF4/6和S2 PDF3原页。S4/S5本轮未回查PDF。没有联网、下载、训练、仿真或改写papers，未递归精读引文；其他数值明确据既有精读笔记。页码从本地PDF第1页计。

**已有人工意见**：回查global-review第555项原答复，只认可Lin2026的表征／识别／计数用途，并保留面积代理标签；不是逐事件t/Q恢复认可。另核第668项，Molina2026仅限NaI(Tl)–SiPM堆积拒绝工程证据，不作为分解算法。S1/S2/S4/S5仍AI草稿待审，没有新增用户科学确认。

**文献依据、条件与实际采用范围**：

- **S1：Bashir Sadeghi; Dustin Scriven; Aaron Chester; Ronald Fox; Sean Liddick; Giordano Cerizza，2026（2025在线）**。*A machine learning framework for accurate and robust analysis of radiation detector pulses*。[笔记](../../papers/sources/2026_Sadeghi_NIMA_A-machine-learning-framework-for-accurate-and-robust-analysis-of-radiation-detector-pulses.md)、[原文](../../papers/raw/2026_Sadeghi_NIMA_A-machine-learning-framework-for-accurate-and-robust-analysis-of-radiation-detector-pulses.pdf)、[DOI](https://doi.org/10.1016/j.nima.2025.170971)。PDF4式9–12、图4；PDF5式13/数据；PDF6表2；PDF7§4.6.5。条件：CeBr3 50×50×3 mm³/H13700 PMT；250 MSPS、4 ns、250点/1000 ns；实测单脉冲及同一单脉冲缩放/平移合成双脉冲；未给独立连续输入率和完整噪声强度。结论：联合分类、参数化响应和自编码重建；训练使用已知两个分量。表2分量ASNR 20.20/20.32 dB vs LMFit 13.67/12.75 dB，整体ASNR 21.94 vs 22.28 dB；Δt RMSE 0.71 vs 27.06 SP，1 SP=4 ns。真实175双脉冲无独立分量真值。 项目采用：采用分量监督与整体/分量指标分开的设计依据；未知K、多重数原生共享微胞推广为本项目待验证方案。待核：真实逐分量真值、未知K≥3、同波形副本相关、LMFit设置、噪声与输入率、式12符号异文及公开代码版本。

- **S2：Alberto Regadío; Luis Esteban; Sebastián Sánchez-Prieto，2021**。*Unfolding using deep learning and its application on pulse height analysis and pile-up management*。[笔记](../../papers/sources/2021_Regadio_NIMA_Unfolding-using-deep-learning-and-its-application-on-pulse-height-analysis-and-pile-up-management.md)、[原文](../../papers/raw/2021_Regadio_NIMA_Unfolding-using-deep-learning-and-its-application-on-pulse-height-analysis-and-pile-up-management.pdf)、[DOI](https://doi.org/10.1016/j.nima.2021.165403)。PDF3图2–3；采样/训练条件据既有笔记PDF2，模拟检出据PDF6表1。条件：模拟CR–RC 1/0.1 μs、50 MHz/20 ns、成形前白噪声；训练到达率说明10⁴–2×10⁵ s⁻¹，表1为10⁴–4×10⁴ s⁻¹；另有NaI/前放/2 μs成形实测谱，无真实高率检出分母。结论：一维U-Net含中央LSTM，输出幅度承载的短事件序列；当前本项目U-Net没有该LSTM。模拟97.5%仅表1噪声0.10/4×10⁴ s⁻¹条件，不是实测MHz恢复率。 项目采用：多尺度卷积与时序瓶颈的结构候选；U-Net–GRU为本项目替代设计，不是原文复现。待核：dense线性/图中ReLU三层异文，20 ns/图5采样异文、检出误报及真实高率验证缺口。

- **S3：Congyu Lin; Xiaoying Zheng; Tom Trigano; Dima Bykhovsky; Yongxin Zhu; Li Tian，2026**。*Spectroscopic Pulse Embeddings by Contrastive Learning from Unlabeled Data for Pile-Up Analysis*。[笔记](../../papers/sources/2026_Lin_Sensors_Spectroscopic-Pulse-Embeddings-by-Contrastive-Learning-from-Unlabeled-Data-for-Pile-Up-Analysis.md)、[原文](../../papers/raw/2026_Lin_Sensors_Spectroscopic-Pulse-Embeddings-by-Contrastive-Learning-from-Unlabeled-Data-for-Pile-Up-Analysis.pdf)、[DOI](https://doi.org/10.3390/s26072138)。本轮PDF10–11图7/式11–15；结果与限制据既有笔记PDF13、17、20–24、29。条件：SSRF Vortex-ME4 SDD/复位CSA，16 bit/250 MHz，64点聚集段、SNR 26.2 dB，80/450 kcps与1 Mcps；真实未标注预训练+模拟标签；NSCL CeBr3/PMT数据另训，实测单脉冲合成。结论：真实波形对比预训练后做堆积分类及聚集段多重数；没有逐事件t/Q回归。既有笔记表5 Cu 1 Mcps总体98.65%，表4峰堆积子集85.73%，不作事件恢复率。 项目采用：实测少标签/模拟到实测的后续训练策略，后置于当前同模型分量恢复。待核：面积代理标签和批次隔离、InfoNCE掩码/分母、64点边界、图表异文、真实逐事件真值。

- **S4：J. M. Deltoro et al.（其余17位作者见原笔记首页），2025**。*Reconstruction of pile-up events using a one-dimensional convolutional autoencoder for the NEDA detector array*。[笔记](../../papers/sources/2025_Deltoro_NST_Reconstruction-of-pile-up-events-using-a-one-dimensional-convolutional-autoencoder-for-the-NEDA-detector-array.md)、[原文](../../papers/raw/2025_Deltoro_NST_Reconstruction-of-pile-up-events-using-a-one-dimensional-convolutional-autoencoder-for-the-NEDA-detector-array.pdf)、[DOI](https://doi.org/10.1007/s41365-024-01606-y)。据既有精读笔记PDF4图3、PDF5–6数据与式2、PDF7表1–2、PDF7–8图10；本轮未回查PDF。条件：NEDA EJ301/PMT；250 MHz、14 bit；真实单脉冲合成n–γ/γ–n对，12–160 ns、150–7000 ADC；CC标签，无TOF；未给连续随机到达率。结论：共享编码器+两个按到达顺序输出的波形解码器；相关系数0.988与γ/n标签保持84.7%/83.3%并存。 项目采用：分量波形监督的旁证及评价反例，不选作高率未知K主方案。待核：固定双异类脉冲、训练/测试母波形隔离、PSD标签、表文成功率及分母、真实高率多重堆积。

- **S5：Junhyeok Kim; Byoungil Jeon; Jisung Hwang; Gyohyeok Song; Myungkook Moon; Gyuseong Cho，2023**。*Pulse height estimation and pulse shape discrimination in pile-up neutron and gamma ray signals from an organic scintillation detector using multi-task learning*。[笔记](../../papers/sources/2023_Kim_ARI_Pulse-height-estimation-and-pulse-shape-discrimination-in-pile-up-neutron-and-gamma-ray-signals-from-an-organic-scintillation-detector-using-multi-task-learning.md)、[原文](../../papers/raw/2023_Kim_ARI_Pulse-height-estimation-and-pulse-shape-discrimination-in-pile-up-neutron-and-gamma-ray-signals-from-an-organic-scintillation-detector-using-multi-task-learning.pdf)、[DOI](https://doi.org/10.1016/j.apradiso.2023.110880)。据既有精读笔记PDF3–4图4–5/损失、PDF6表2–3、PDF7表5；本轮未回查PDF。条件：2英寸stilbene/H7195 PMT，250 MS/s、12 bit、248点；已寻峰附近16点+梯度；合成12–80 ns、1–9重，真实源最高28943 cps；传统CC标签。结论：共享特征的高度/PSD多任务回归依赖先寻峰。MT非所有指标更好，九重分类60.22% vs ST 63.84%；Acc.H/Recall公式与文字存在异文。 项目采用：多任务不保证检出与电荷同时改善；当前三输出头已是多任务，不能把新增多任务名称当新方法。待核：漏检寻峰分母、标签/母波形隔离、任务权重、accuracy及recall定义、高率真实未知K评价。

**本项目判断与拟议路线，不是新计算结果**：新确认2098目标中CNN完整恢复518、已匹配Q失败907、漏检673；原训练已使用事件／logQ／时间三头，因此再以“多任务”为名不构成新方法。优先检验SAEFit式分量监督和训练中的参数化物理解码；以现有GRU作结构对照，分别消融分量损失、U-Net–LSTM/GRU瓶颈与训练Q覆盖。原论文的固定单／双输出不直接搬进未知K主表；事件序列或有限槽位／存在性推广为项目假设，容量、近邻及重复输出单列。本轮未找到这些推广在当前GAGG条件的直接文献证明。

原生祖先分量只作训练／评价，推理不读K、到时、Q真值或源率；按完整ROOT／种子先划分。继续固定95/351 ns双发光衰减，充分动态范围前提不变；共享微胞分量不能当独立低率模板随意叠加，均值物理解码失配须检查。Q为雪崩等效电荷，未转沉积能量；阈下α和暗计数仍进入波形模型。Lin式对比预训练后置于实测少标签／域迁移，不把其分类准确率当本项目事件恢复率。

**待验证和文献整理待办**：新的开发／验证协议分别比较损失与结构，多训练种子；控制间隔／幅比／Q外推，再测八档原生流。在相近额外输出预算下比较完整恢复，同时报告Q偏差／宽度、漏检、误拆分／误合并和循环错位诊断。高率按独立ROOT统计；已公开23窗和新16窗现均已查看，之后只作开发诊断／回归，另留冻结后测试。补100 MHz训练／验证覆盖及跨率对照。SAEFit式12符号和代码版本、原噪声／输入率与未知K扩展，Regadío架构／时间异文，以及Lin对比掩码和代理标签等缺口，若需澄清应在papers库另立维护任务，本项目不改该库。未把缺口写成已完成性能改进。

记录保存在本地 `jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb` 第22节；[来源、页码与哈希](../jupyter/run_records/paper_search_ml_20261021_141317/evidence.json)。旧33个代码单元和23幅图保留，本轮仅增加研究文字与目录入口，不改变推断和评分。

<a id="alpha-component-ml-20261021"></a>
### L45：2026-10-21——分量监督物理解码与U-Net–GRU的实际实现和独立比较

**采用依据与范围**：沿用L44的Sadeghi/SAEFit（DOI `10.1016/j.nima.2025.170971`，既有定点PDF4–7式9–13/表2）与Regadío2021（DOI `10.1016/j.nima.2021.165403`，PDF3图2–3）设计启发。前者为CeBr3/PMT单/合成双脉冲，后者主要为CR–RC模拟及NaI谱学；原噪声/采样/率和人工审核边界仍按L44。共享库只读，本轮没有新精读或更新其审核。此处为项目改编，不把论文分类率或分量ASNR移植为GAGG效率。

**原值→新实现**：保留原GRU、原事件BCE＋2logQ Huber＋0.5时间Huber，另实现固定95/351 ns双发光衰减及已标定电子学的可微物理解码。新增项为实际原生祖先分量误差（权重1）与这些分量的条件求和误差（权重0.1），在训练真值最近bin计算，支持[-16,1200] ns、2 ns采样；它不是完整含暗电压的无监督重建，也没有更改最终密集事件表示/8 ns候选峰选择。相同bin多事件仅从新增项排除，原目标和测试不排除。实际雪崩保留共享微胞历史，只重放线性电子学，不把独立低率分量随意叠加。

新增U-Net–GRU为多尺度与时序结构的整体对照，未单独区分两者贡献。四方案各两个种子、2400步/batch16；同结构/同种子配对初始权重、采样、数据和优化器。7435/159315参数不当作等计算量比较。49训练/26验证按完整ROOT划分：旧26/10保留角色，已公开23测试转训练、16确认转验证；不再称它们新测试。每400步及第1步以验证逐率宏电荷合格F1在预声明阈值网格选检查点；另按配对GRU的验证额外输出预算选阈值，新测试不保证相等。

**新数据与物理条件**：全部权重/阈值/源码/测试计划冻结后生成60个新ROOT，0.1/1/10/20/40/60/80/100 MHz各40/8/2/2/2/2/2/2窗；低率扩展在冻结前完成，初始23种子保留。纯p–¹¹B有效单α谱，无γ；共享微胞、暗链、噪声、完整20 μs模拟电压、中央[5,15) μs评分，充分动态范围。Q为实际雪崩等效电荷，不是沉积能量；源注入率不当作探测/反应率。推理只读观测电压与冻结参数，无真值、K或源率输入。

**本项目实际结果**：同一新批共同目标4253个，时间±10 ns最大基数一对一匹配后Q误差≤20%才完整恢复；两个训练种子共享同批ROOT，不当作两倍物理曝光。

| 方案 | 完整恢复（两种子） | 额外输出（两种子） |
| --- | ---: | ---: |
| GRU原目标 | 1041/1100 | 525/498 |
| GRU分量目标 | 856/832 | 512/521 |
| U-Net–GRU原目标 | 1124/1164 | 462/455 |
| U-Net–GRU分量目标 | 922/830 | 475/489 |
| 当前WLS | 812 | 300 |
| 旧冻结GRU | 1036 | 217 |


**本轮取舍**：保留U-Net–GRU原目标作为候选，暂不采用当前分量损失；两个结构/两个种子都未获益。相对同批新GRU原目标，U-Net–GRU原目标完整恢复增加83/64、额外输出减少63/43，属于小幅改善。旧冻结U-Net恢复1148、额外466，新结构尚未明显超越它；旧数据划分不同，不单独归因于GRU瓶颈。

**探索性诊断**：新增分量方案在时间匹配者上的Q误差中位数为约−21.1%至−13.3%；事后在真值最近bin读取电荷头没有消除下偏。该oracle检查不修改盲输出或主评分，既没有提供恢复率提升，也未独立确定失配、损失几何或训练预算的因果。实际源和结果另存`alpha_component_posthoc.py`与本地`anchor-head-diagnostic.json`。下一轮优先检查未知K事件表示与事件级Q/计数约束；若再试分量项，先核对孤立/受控电荷无偏，另定开发与新测试。

合计受高率目标主导；Notebook另报逐率、条件Q偏差/宽度、漏检、电荷失败、拆并标志、同事件得失和额外输出。Q宽度不是能量FWHM，较小总波形残差不保证正确分量，额外输出也不一定对应物理上完全无信号。

**实际检查与I/O记录**：1320份新测试评分重算、1200份神经输出复解码、120个物理总模型重构，60个实际RunInfo/宏/三个种子及冻结时序核对；75分量缓存、配对配置、禁用数组加载的盲推断、冷复跑及8个完整原生分量图通过。另逐项核对60个原ROOT采样时间/电压与导出数组；独立进程冷复跑22份数值拟合一致，仅WLS的`refinement_seconds`计时字段不参与数值相等断言，验证器适配器保留原冻结文件。62软件控制含白/AR1噪声、空/单/双及高幅外推，1364份评分与124个物理模型重算；不提供原生最小可分间隔。原控制保存函数在数值工作完成后遇到数组JSON转换错误，独立保存适配器以原种子重生成同62例，物理数组另存NPZ，冻结推断源码/参数不改。推断可在G4其余窗口生成时顺序处理已完成记录；流式驱动只增加I/O编排，不改变计划或模型。

**绘图恢复**：全窗图首轮导出以143中断，无Python回溯，原因未确认；部分输出保留。单独适配器仅固定图例位置，重导出相同案例/分量/采样/图号，已完成首轮数组与重导出逐项一致。数值推断与评分不改。

**限制与待办**：高率每档2短窗、同窗相关、原2 ns/8 ns表示、均值响应失配及高密度偶然匹配仍限制结论。缓存字段名`mean_decoder_floor`实际为真实Q/时间下的参考误差，不是最小残差下界。新增损失只在训练标注中心约束，未直接优化最终盲列表；一次权重/预算下的结果不否定所有分量监督。后续若改权重、表示或响应自由度，需另定开发/验证协议并留新测试，不在此批上反复挑选。对比预训练/实测迁移、混合α/γ、free-τ、最短窗和沉积能标未由本轮回答。

过程位于本地`jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb`第23–28节。来源：[冻结实验](../pileup/experiments/alpha_component_ml_20261021/README.md)、[实际结果](../jupyter/artifacts/alpha_component_ml_20261021/report-tables.json)、[验证](../jupyter/artifacts/alpha_component_ml_20261021/validation.json)。实际源码/宏/运行器副本保留，无需Git回退；本地输入、模型与环境仍须保留。


<a id="ml-training-revision-20261021"></a>
### L46：2026-10-21——针对分量监督负结果修订机器学习训练设计

**问题与证据层级**：用户要求用paper-search尝试解决ML解堆积不佳。复用本地topic/review/vocabulary索引、三份优先主题页及完整来源笔记，定点回查四篇原PDF。实际提取Sadeghi PDF3–6、Deltoro PDF4–6、Kim PDF5/7–9、Regadío PDF3/7；实际看原页分别4、4、8、3。页码按本地PDF第1页，Kim PDF8为印刷页6。索引search-topic延迟后完成、命中7个显式关联主题；search-review进程退出143、没有输出且原因未确认，以rg和global-review原答复替代。不据工具无输出断言无相关文献。本轮不联网、不递归精读新引文、不改papers，不训练或新模拟。

**文献结论、条件与采用范围**：

- **R1：Bashir Sadeghi; Dustin Scriven; Aaron Chester; Ronald Fox; Sean Liddick; Giordano Cerizza，2026**，*A machine learning framework for accurate and robust analysis of radiation detector pulses*。Nuclear Instruments and Methods in Physics Research A 1082, 170971。[笔记](../../papers/sources/2026_Sadeghi_NIMA_A-machine-learning-framework-for-accurate-and-robust-analysis-of-radiation-detector-pulses.md)、[原文](../../papers/raw/2026_Sadeghi_NIMA_A-machine-learning-framework-for-accurate-and-robust-analysis-of-radiation-detector-pulses.pdf)、[DOI](https://doi.org/10.1016/j.nima.2025.170971)。PDF3式3–8；PDF4图4/式9–12；PDF5数据及式13–15；PDF6表2。**条件**：CeBr3 50×50×3 mm³/H13700 PMT；250 MSPS/4 ns；250点/1 µs；实测单脉冲及其缩放、平移副本合成双脉冲；没有独立连续随机源率或完整噪声强度。 **原文支持**：先拟合共享响应参数，再训练逐事件参数修正；单/双分类、整体与两个已知分量及自编码重建联合训练。先按母单脉冲70/10/20划分再合成。 **项目差异**：当前仅在训练真值最近bin的Q/时间头上加固定响应分量项；没有论文的单/双分类、共享响应初始化、逐事件形状修正和自编码分支。论文的两阶段指共享参数→神经修正，不等于已经验证单→双→多的课程训练。 **限制/待核**：固定单/双任务；同一母波形副本相关；175条真实双脉冲没有独立分量真值。ASNR和分类准确率不等于本项目Q±20%、时间±10 ns完整恢复率；式12单脉冲共享参数符号与式11不一致，保留待核。 **审核**：AI精读草稿待人工审核；没有针对本篇科学结论的新人工确认。

- **R2：J. M. Deltoro; G. Jaworski; A. Goasduff; V. González; A. Gadea; M. Palacz; J. J. Valiente-Dobón; J. Nyberg; S. Casans; A. E. Navarro-Antón; E. Sanchis; G. de Angelis; A. Boujrad; S. Coudert; T. Dupasquier; S. Ertürk; O. Stezowski; R. Wadsworth，2025**，*Reconstruction of pile-up events using a one-dimensional convolutional autoencoder for the NEDA detector array*。Nuclear Science and Techniques 36, 32。[笔记](../../papers/sources/2025_Deltoro_NST_Reconstruction-of-pile-up-events-using-a-one-dimensional-convolutional-autoencoder-for-the-NEDA-detector-array.md)、[原文](../../papers/raw/2025_Deltoro_NST_Reconstruction-of-pile-up-events-using-a-one-dimensional-convolutional-autoencoder-for-the-NEDA-detector-array.pdf)、[DOI](https://doi.org/10.1007/s41365-024-01606-y)。PDF4图3/§3.3；PDF5–6数据构造；其他成绩据既有完整笔记PDF7–8。**条件**：NEDA EJ301/PMT；250 MHz/4 ns、14 bit；实测单脉冲人工合成n–γ/γ–n；间隔12–160 ns、幅度150–7000 ADC；CC标签，无TOF真值；无连续随机源率。 **原文支持**：共享编码器、两个独立解码器，分别输出先/后两个实际组成波形，以合成前两个波形为训练标签。 **项目差异**：这是直接分量波形输出；当前网络实际只输出事件概率、Q和时间，分量是固定物理解码器派生的辅助项。 **限制/待核**：只验证固定两个异类脉冲；没有同类αα、未知K≥3或高率原生共享SiPM检验。母波形是否跨训练/测试复用未交代清楚；相关系数0.988并非逐事件能量/时间准确率。论文没有验证PIT/集合匹配。 **审核**：AI精读草稿待人工审核；入库授权不等于科学结论认可。

- **R3：W. Kim; K. Ko; S. Lee; J. Park; G. Song; K. Lim; G. Cho，2023**，*Optimizing deep learning-based piled-up pulse height correction method for high radiation-field application*。JINST 18 C12002；正式出版会议论文。[笔记](../../papers/sources/2023_Kim_JINST_Optimizing-deep-learning-based-piled-up-pulse-height-correction-method-for-high-radiation-field-application.md)、[原文](../../papers/raw/2023_Kim_JINST_Optimizing-deep-learning-based-piled-up-pulse-height-correction-method-for-high-radiation-field-application.pdf)、[DOI](https://doi.org/10.1088/1748-0221/18/12/C12002)。PDF5数据；PDF7寻峰与网络；PDF8图4/表2（印刷页6）；PDF9表3及窗口选择。**条件**：LaBr3 1.5×1.5英寸/PMT、宽带放大；500 MS/s/2 ns、496点/992 ns；实测单脉冲合成人工1–10重堆积，训练间隔16–80 ns；噪声未完整量化，无实测连续Mcps恢复验证。 **原文支持**：检峰和条件脉高回归分别评价；表2中8 ns间隔导数检峰为105003/200000=52.50%，16 ns为150684/200000=75.34%，因此回归精度不能代替检出效率。 **项目差异**：本项目输入已有导数特征；要改的是输出容量/候选选择与评价，不能把再加导数称作新方案。 **限制/待核**：不直接提供逐事件总效率。式2.2写平均绝对相对误差，而表3约99%按准确度解释，尚有定义异文；不采用99%作本项目恢复率或最小窗依据。 **审核**：AI精读草稿待人工审核；已入库，没有新增科学采用确认。

- **R4：Alberto Regadío; J. Ignacio García Tejedor; Luis Esteban; Sebastián Sánchez-Prieto，2024**，*Three topologies of deep neural networks for pulse height extraction*。Applied Radiation and Isotopes 204, 111137。[笔记](../../papers/sources/2024_Regadio_ARI_Three-topologies-of-deep-neural-networks-for-pulse-height-extraction.md)、[原文](../../papers/raw/2024_Regadio_ARI_Three-topologies-of-deep-neural-networks-for-pulse-height-extraction.pdf)、[DOI](https://doi.org/10.1016/j.apradiso.2023.111137)。PDF3式1–3；PDF7表1。**条件**：固定CR–RC 1/0.1 µs、20 MHz/50 ns、1024点；均匀0–1幅度；成形前人工白/Brown噪声；非GAGG/SiPM实测。率/图轴单位有原印异文。 **原文支持**：MMSE对整条目标事件序列选择一个0–5采样的公共延迟后计算MSE，允许训练输出整体延迟。 **项目差异**：这是可单独消融的对齐损失；不是每个事件任意移动、不是一对一集合匹配，也不是改变本项目物理τ。 **限制/待核**：缺逐事件Q/时间及误检联合评价；不能靠事后真值移动得到测试改善。当前oracle最近真值bin仍有Q偏差，因此该损失只列次优先级。 **审核**：AI精读草稿待人工审核；已入库，没有针对本篇科学结论的人工确认。

**实际项目诊断，非论文数值**：读取L45同批保存结果4253目标；U-Net–GRU原目标两种子已匹配Q失败占1725/2849=60.55%、1620/2784=58.19%，分量方案为1895/2817=67.27%、2011/2841=70.78%。因此候选数与Q校准应分开诊断，不把所有失败归于漏检。均值模板软件single_4000白噪声控制种子1原/分量方案Q=4373/5012 pe，WLS=4000.12 pe；10 ns等幅双控制完整2/0/2。控制分布与原生训练不同；原生分量模型总体偏负，而该控制偏正，不能据此单独确诊错误或用统一比例修正。原有oracle最近bin诊断仍有Q下偏，非盲提升。

**采用决定与拟议路线**：不重复把同一个真值bin辅助项称作完整SAEFit。优先A：训练集简单案例过拟合检查＋独立同分布控制，再做原生孤立Q/t训练验证；B：同结构同种子消融原目标、逐真Q²归一化分量项和训练确定的固定尺度分量项，检查Q覆盖、绝对尺度及收敛；C/D：直接两个实际分量输出，与两个存在性/Q/t物理解码分支比较，加入空/单/双使0/1/2未知，系统扫描间隔/幅比/Q/噪声，与WLS及当前密集头对照；E/F：单/双同时改善后，再扩展多槽＋存在性＋一对一分配和原生高率。课程训练、固定尺度损失、空窗扩展、PIT/集合匹配为本项目待验证提案，不是上述论文已证明的GAGG成果。SAEFit“两阶段”是共享响应→逐事件神经修正，不能引用成单→双→多课程的直接证据。MMSE仅列对齐消融，不借测试真值移动输出。

**数据与评价**：先按母ROOT/种子分组，再裁剪/合成；均值线性叠加与原生共享微胞分量分开，暗链/阈下α保留；推断不读真K、Q、时刻或源率。暂保持95/351 ns双发光衰减及充分动态范围；Q为实际雪崩等效pe，非沉积能量。每2 ns一槽和8 ns候选筛选可能造成竞争，但偏移后输出时间差不一定≥8 ns。计数、条件Q/时间、完整效率、误拆/误并、额外输出同时评价，不能以总残差/相关度替代。已查看60窗今后只作开发/回归，新方案冻结后另留测试；不重复使用此批挑选后再叫独立确认。高率每档仅两ROOT的限制保留。

**既有人工采用意见的边界**：核对[第545项原答复](../../papers/global-review.md#human-review-545-huang2025-confirmed)，Huang2025只认可给定谱库、已知率和模拟条件下统计能谱重建；[第555项原答复](../../papers/global-review.md#human-review-555-lin2026-confirmed)，Lin2026只认可表征/识别/多重数且保留面积代理标签。第668/669项的堆积拒绝/谱学采用边界亦保留；没有将这些意见转移为逐事件Q/t方法认可。

**待验证及库维护待办**：先验证A/B，再决定C/D实现；不得把本次方案写成已改善。SAEFit式12符号/作者代码版本/真实分量真值，Deltoro/Kim母波形隔离与损失超参，Kim准确度公式，Regadío率/图轴与实现异文仍待核。若需新精读/补充资料，另在papers库处理，本项目不修改共享库。单/双任务的文献结果不能推成GAGG 80–100 MHz未知K恢复保证；free-τ、混合α/γ、沉积能标不由本次回答。

本地记录：jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb第29节；[证据、条件与输入哈希](../jupyter/run_records/paper_search_ml_training_20261021_164120/evidence.json)。保留此前58个代码单元、39幅PNG和全部评分；新增只读取保存结果的轻量复核表，无新图/训练/推断。原Notebook/文档在同运行目录备份。


<a id="alpha-pair-curriculum-20261021"></a>
### L47：2026-10-21–22——单脉冲校准、直接分量头与共同微胞双脉冲的实际训练

**依据与边界**：按L46/Notebook0004第29节实现。Sadeghi2026（10.1016/j.nima.2025.170971，PDF3–6式3–12/表2）、Deltoro2025（10.1007/s41365-024-01606-y，PDF4图3、PDF5–6数据）的原探测器/采样/噪声/固定双脉冲边界仍按L46；AI待审状态没有新增科学确认。本轮不重新检索/全文阅读或改写papers。课程训练、观测积分先验、两槽存在性和训练集合分配为项目改编，不称论文原样复现。

**原值→实现**：由旧2ns密集头/8ns选峰/真值训练anchor辅助项，改做完整2.048µs局部帧的已知单事件校准与未知0/1/2分离。单事件Q用观测积分和预测到时的均值尾部修正，加有界残差；比较无分量、真实分量范数归一化、训练固定尺度三目标，各两种子/1200步。真实分量范数含形状涨落，近似Q²而非旧固定均值响应范数的逐行复现；共同积分结构的影响不能单独归因于归一化。两槽头分别为正观测电压的分量掩码和固定双τ物理解码，先比较无序最小分配，再在开发失败后、未生成新测试时追加有序目标＋独立计数＋全局条件解码，共四方案各两种子/2400步。此联合改法不分离各因素因果；推断不提供K/t/Q真值或源率。

**数据与训练**：85/35旧孤立α母ROOT按原标定/开发分组，合格59/26。双分量前200步空/单预热，之后50%共同微胞光子重数字化、25%线性分量控制、25%均值控制。512训练/128验证联合数字化案例来自各自组内母光学输运；共同推进恢复/PDE/XT/暗链及噪声，电子学只在线性分量重放时拆分。线性缩放/插值不冒充原生微胞堆积。固定95/351ns和充分动态范围，Q为雪崩等效pe，非沉积能量。

**冻结后测试与规则**：14模型、阈值、传统拟合惩罚、源码冻结后生成24个新原生G4孤立ROOT；独立均值、线性分量和256联合重数字化控制分域报告，不累加成独立ROOT曝光。新传统对照盲选0/1/2、NNLS消去Q并搜索时刻；验证空/单额外输出帧预算≤5%。网络refit保留网络报告K，局部±12ns再拟合。已知K拟合只作单列诊断。普通时间10ns/Q20%之外，新增min(10ns,真间隔/4)严格判据，真间隔只进入评分，以免近邻中点输出被误称分开。

**原生裁帧补充**：原始真时刻居中选帧会提供时间先验，保留为诊断。初次测试推断完成后、尚未审阅准确度时发现这一问题；冻结权重/阈值不变，另按实际宏中的源触发时刻随机64–256ns预触发裁帧，选帧不读真实到时。主要单脉冲表与native比较采用此jitter-native视图；它与旧视图共享同24母ROOT，不是新增独立曝光，仍是给定外部触发的局部任务，不是完整20µs盲定位。设计时序和864行补充评分核验单独保留，未冒称最初冻结前已预声明。

**实际结果**：新已知单事件合格分母18，single_logq两种子分别为18/18、18/18；single_relative两种子分别为18/18、18/18；single_fixed两种子分别为18/18、18/18。共同微胞控制盲拟合完整179/267、额外35；pair_direct：完整118/122, 额外12/6; pair_physical：完整140/121, 额外14/11; pair_ordered_direct：完整166/168, 额外16/18; pair_ordered_physical：完整163/168, 额外15/14。有序物理头＋局部重拟合，两种子分别完整恢复184/267和191/267，额外输出分别14和13；81个两事件均合格的帧中，严格双事件成功分别37和41，盲拟合为32。两训练种子分别保留，范围不是置信区间。 部分流程在这批共享微胞控制上同时增加完整恢复、保持或减少额外输出；这只说明本局部任务的改善，不能推成高率未知K结论。

**实际验证与限制**：重算全部普通/严格评分、重解码原始槽预测，验证分量非负/和守恒、九个数字化案例精确重放、四独立进程禁止额外NPZ读取的冻结盲推断、24原生ROOT哈希/实际宏/三种子/冻结时序；3项单元检查及训练集工程过拟合通过。图30-1/30-2/31-1/32-1/32-2/33-1与原39图保留。两槽容量不支持任意高率完整20µs未知K；光子重数字化是离线同模型概率实现、非C++随机序列逐位等价。多案例共享24母ROOT，种子差异非置信区间；严格近邻失败、弱/零响应和额外输出均保留。共享状态4ns合格双事件有序物理头＋refit两种子均0/7，未得普适间隔/最高率、free双τ、混合α/γ或沉积能标结论。本次是局部无加权拟合对照，不能等同旧全窗WLS/GLS。先解决单/双失败，再决定未知K≥3/随机流扩展，未将超容量样本删后宣称高率成功。计数与单事件Q/形状涨落是下一轮开发检查项，改法须用另一批独立母ROOT测试。

本地Notebook0004第30–34节及Run_Catalog记录实际过程；[独立实验](../pileup/experiments/alpha_pair_curriculum_20261021/README.md)、[结果](../jupyter/artifacts/alpha_pair_curriculum_20261021/report-tables.json)、[验证](../jupyter/artifacts/alpha_pair_curriculum_20261021/validation.json)。实际源码/宏/运行器副本保留，复现不要求Git回退；模型、输入和环境仍需本地保留。


<a id="pileup-resolution-literature-20261022"></a>
### L48：2026-10-22——解堆积最小间隔：不可辨识、条件恢复率与最短测试范围

**问题与采用结论**：用户询问“多少时间间隔就一定不能解pileup”。本轮本地只读检索8篇既有完整来源笔记，定向回查15个原PDF页、目视核对6页图表；没有新增全文精读、联网、训练或G4运行。未找到能直接用于本GAGG/SiPM配置、对所有方法成立的正数ns不可分定理，不将检索范围说成文献穷尽。采用条件双事件恢复曲线及固定误报预算，区分寻峰、Q/t联合恢复、统计能谱校正与最短测试范围。

**检索与审核**：从topic-index/review-index及vocabulary导航，优先读数字脉冲处理、PSD及光产额/能量分辨率主题页，再回来源与global-review。两条只读knowledge_index查询各45秒超时，已改用rg和手工索引；超时/无命中不表示文献不存在。global-review新增47篇归并M2/M4/M5和逐篇B23-018/019/023/044/052/054仍为AI待审；Kim JINST与Sadeghi来源笔记也未获本轮新增科学确认。其他文章的既有局部意见不转移到本题。未在本轮库内匹配Georgiev2024 heuristic deconvolution，不把先前联网预印本当本轮已入库证据。

以下每篇的条件与局限据实际来源保存；页码从本地PDF第1页计。

- **R1：P. Adrich; I. Zychor（2021）**，*Least-squares fitting algorithm for peak pile-up correction in gamma-ray spectroscopy*。NIMA 990, 164962；[DOI 10.1016/j.nima.2020.164962](https://doi.org/10.1016/j.nima.2020.164962)。[精读笔记](../../papers/sources/2021_Adrich_NIMA_Least-squares-fitting-algorithm-for-peak-pile-up-correction-in-gamma-ray-spectroscopy.md)、[本地原文](../../papers/raw/2021_Adrich_NIMA_Least-squares-fitting-algorithm-for-peak-pile-up-correction-in-gamma-ray-spectroscopy.pdf)。**本轮局核**：PDF3装置与插值；PDF5图6及成功定义；其余装置/构造细节据既有精读笔记。**条件**：1×1英寸LaBr3:Ce/R6231-100 PMT/DT5730；低率实测单脉冲相加并平移，真实高率测量约80 kcps。20–80%上升8.5 ns；论文未明确实际采样设置，噪声由实测波形带入而未独立扫描。 **原文支持**：峰值时刻间隔约8.5 ns对应50%双脉冲成功，约17 ns超过90%。成功要求两脉冲都恢复且各能量在预期值±3σ内；σ=1.55%（661 keV），不是任意幅度/时间联合正确率。 **项目采用**：采用随间隔报告成功率和明确成功定义的方式，不把LaBr3数值移植为GAGG阈值。 **缺口**：图6是效率曲线，不是小于8.5 ns必然失败的界限；未给各间隔样本数、完整幅比扫描和同一误拆预算。p值不能作双事件后验概率。 **审核**：既有AI精读草稿待人工审核；本轮没有新增科学确认，不扩大其他论文的局部采用意见。

- **R2：J. Z. Zhang; H. Y. Wu; W. Jiang; Z. H. Wang; R. R. Fan; X. C. Ruan（2024）**，*Development of a pile-up pulse recovery algorithm for the LaBr3 detector*。NIMA 1063, 169273；[DOI 10.1016/j.nima.2024.169273](https://doi.org/10.1016/j.nima.2024.169273)。[精读笔记](../../papers/sources/2024_Zhang_NIMA_Development-of-a-pile-up-pulse-recovery-algorithm-for-the-LaBr3-detector.md)、[本地原文](../../papers/raw/2024_Zhang_NIMA_Development-of-a-pile-up-pulse-recovery-algorithm-for-the-LaBr3-detector.pdf)。**本轮局核**：PDF3模板与拟合；PDF4图6(a)、模拟幅度/噪声/率及10 ns解释。**条件**：2×2英寸LaBr3:Ce/R6233 PMT，500 MS/s（2 ns）、14 bit；幅度相关模板。合成与恢复共用模板；幅度2000/5000/8000 ADC，白噪声σ=7 ADC，随机事件率10^6–10^8 s^-1。装置细节据既有笔记。 **原文支持**：作者据恢复间隔谱图6(a)给约10 ns的准确定位下限；小于该尺度的脉冲容易合并，并形成和峰。 **项目采用**：把约10 ns作为一个有条件的算法实例，用短间隔缺失和和峰诊断合并。 **缺口**：这是给定模板、噪声、幅度和算法的结果；没有证明所有算法在10 ns以下不可分，也没有给统一逐对Q/t成功率。 **审核**：既有AI精读草稿待人工审核；本轮没有新增科学确认，不扩大其他论文的局部采用意见。

- **R3：F. Belli; B. Esposito; D. Marocco; M. Riva; Y. Kaschuck; G. Bonheure; JET EFDA contributors（2008）**，*A method for digital processing of pile-up events in organic scintillators*。NIMA 595, 512–519；[DOI 10.1016/j.nima.2008.06.045](https://doi.org/10.1016/j.nima.2008.06.045)。[精读笔记](../../papers/sources/2008_Belli_NIMA_A-method-for-digital-processing-of-pile-up-events-in-organic-scintillators.md)、[本地原文](../../papers/raw/2008_Belli_NIMA_A-method-for-digital-processing-of-pile-up-events-in-organic-scintillators.pdf)。**本轮局核**：PDF3寻峰与截断补尾；PDF4图5、图6及50 ns/87%说明。**条件**：NE213/EMI9134B PMT，200 MS/s（5 ns）、14 bit；固定双指数峰后补尾，实测n/γ，最高应用率约230 kHz。未给通用噪声/SNR扫描。 **原文支持**：寻峰最小间距人为设为3样点=15 ns以减少噪声假峰；7拟合点加3点缓冲对应约50 ns峰距，积分误差均值/标准差约在±8%量级，约87%的已识别堆积满足此间距。 **项目采用**：严格分开第二事件检出、积分质量和选择效率。 **缺口**：15 ns是程序规则；50 ns是该补尾法的谱学质量实例。±8%不是每事件保证；87%分母不含未发现的<15 ns堆积。 **审核**：既有AI精读草稿待人工审核；本轮没有新增科学确认，不扩大其他论文的局部采用意见。

- **R4：B. Löher; D. Savran; E. Fiori; M. Miklavec; N. Pietralla; M. Vencelj（2012）**，*High count rate gamma-ray spectroscopy with LaBr3:Ce scintillation detectors*。NIMA 686, 1–6；[DOI 10.1016/j.nima.2012.05.051](https://doi.org/10.1016/j.nima.2012.05.051)。[精读笔记](../../papers/sources/2012_Loher_NIMA_High-count-rate-gamma-ray-spectroscopy-with-LaBr3-Ce-scintillation-detectors.md)、[本地原文](../../papers/raw/2012_Loher_NIMA_High-count-rate-gamma-ray-spectroscopy-with-LaBr3-Ce-scintillation-detectors.pdf)。**本轮局核**：PDF4 §3.3、图4–5；约30/60 ns分区由图5读取。**条件**：LaBr3:Ce/PMT实验500 MS/s；受控双脉冲为60 ns指数尾、高斯平滑约10 ns上升、幅比1:2、SNR10/20，0–200 ns按4 ns扫描，每点5000白噪声实现。 **原文支持**：本触发/矩阵校正方法在约30 ns以下合并；约30–60 ns为误差放大的过渡区，60–200 ns接近本征幅度不确定度。 **项目采用**：采用幅比、SNR和间隔联合扫描；能量波动与恢复效率分开。 **缺口**：作者明确限定current analysis；固定幅比和SNR，图中误差棒为分布FWHM而非均值标准误，不能作为探测器的普适死区。 **审核**：既有AI精读草稿待人工审核；本轮没有新增科学确认，不扩大其他论文的局部采用意见。

- **R5：W. Kim; K. Ko; S. Lee; J. Park; G. Song; K. Lim; G. Cho（2023）**，*Optimizing deep learning-based piled-up pulse height correction method for high radiation-field application*。JINST 18, C12002（正式出版会议论文）；[DOI 10.1088/1748-0221/18/12/C12002](https://doi.org/10.1088/1748-0221/18/12/C12002)。[精读笔记](../../papers/sources/2023_Kim_JINST_Optimizing-deep-learning-based-piled-up-pulse-height-correction-method-for-high-radiation-field-application.md)、[本地原文](../../papers/raw/2023_Kim_JINST_Optimizing-deep-learning-based-piled-up-pulse-height-correction-method-for-high-radiation-field-application.pdf)。**本轮局核**：PDF5数据构造；PDF8表2（印刷页6）。**条件**：1.5×1.5英寸LaBr3/PMT与宽带放大，500 MS/s（2 ns），992 ns帧；实测单脉冲合成1–10脉冲，训练PTP16–80 ns；噪声未量化为完整SNR扫描，无连续高率实测恢复验证。 **原文支持**：表2每行200000个真峰；8/16/24 ns处导数寻峰输出105003/150684/200000个，峰计数比52.50%/75.34%/100%。 **项目采用**：候选检出和后续网络回归分别评价，不用条件回归精度替代全链效率。 **缺口**：这是找到的峰数/真峰数，不是逐对正确分离率，也不是同时满足脉高与时间误差的成功率；100%计数也不能排除漏检与误检抵消。 **审核**：既有AI精读草稿待人工审核；本轮没有新增科学确认，不扩大其他论文的局部采用意见。

- **R6：J. M. Deltoro et al.（2025）**，*Reconstruction of pile-up events using a one-dimensional convolutional autoencoder for the NEDA detector array*。Nuclear Science and Techniques 36, 32；[DOI 10.1007/s41365-024-01606-y](https://doi.org/10.1007/s41365-024-01606-y)。[精读笔记](../../papers/sources/2025_Deltoro_NST_Reconstruction-of-pile-up-events-using-a-one-dimensional-convolutional-autoencoder-for-the-NEDA-detector-array.md)、[本地原文](../../papers/raw/2025_Deltoro_NST_Reconstruction-of-pile-up-events-using-a-one-dimensional-convolutional-autoencoder-for-the-NEDA-detector-array.pdf)。**本轮局核**：PDF6数据构造；原审核限制据来源笔记及global-review原始归并。**条件**：NEDA EJ301/PMT，250 MS/s（4 ns）、14 bit；只合成n–γ/γ–n双脉冲，150–7000 ADC，实测噪声随单波形带入；不是GAGG α/γ。 **原文支持**：从实测单脉冲平移3–40样点，即12–160 ns，生成训练和测试样本。 **项目采用**：引用时明确最短测试间隔；不能把未测试区写成失败区。 **缺口**：12 ns是数据覆盖下限，未验证更短间隔，不是不可分定理。波形相关度、粒子标签保持、事件次序成功率不同；原表/正文异文及母波形隔离仍待核。 **审核**：既有AI精读草稿待人工审核；本轮没有新增科学确认，不扩大其他论文的局部采用意见。

- **R7：B. Sadeghi; D. Scriven; A. Chester; R. Fox; S. Liddick; G. Cerizza（2026）**，*A machine learning framework for accurate and robust analysis of radiation detector pulses*。NIMA 1082, 170971（2025年上线）；[DOI 10.1016/j.nima.2025.170971](https://doi.org/10.1016/j.nima.2025.170971)。[精读笔记](../../papers/sources/2026_Sadeghi_NIMA_A-machine-learning-framework-for-accurate-and-robust-analysis-of-radiation-detector-pulses.md)、[本地原文](../../papers/raw/2026_Sadeghi_NIMA_A-machine-learning-framework-for-accurate-and-robust-analysis-of-radiation-detector-pulses.pdf)。**本轮局核**：PDF5合成策略；PDF7 §4.6分箱与真实分量限制；PDF10图10。**条件**：CeBr3/H13700 PMT，250 MS/s（4 ns）、1 μs帧；母单脉冲先分组，再用同一波形缩放并加平移副本，比例0.1–10；形状/噪声相关。无GAGG高率随机流验证。 **原文支持**：最短时间箱为[2,10]个采样周期，即8–40 ns，以箱上限绘图。分类与分量指标随间隔/幅比变化，图10显示SAEFit短间隔分类优于作者所用LMFit。 **项目采用**：学习按间隔与幅比分层、对分量监督，但不引用其最短箱作为GAGG最小可分间隔。 **缺口**：8–40 ns箱平均不能当8 ns点成功率；分类准确率或总体时间差RMSE不是逐对Q/t恢复保证。175条真实双脉冲没有分别可观测的分量真值，LMFit实现细节不足以代表所有拟合法。 **审核**：既有AI精读草稿待人工审核；本轮没有新增科学确认，不扩大其他论文的局部采用意见。

- **R8：W. Guo; R. P. Gardner; C. W. Mayo（2005）**，*A study of the real-time deconvolution of digitized waveforms with pulse pile up for digital radiation spectroscopy*。NIMA 544, 668–678；[DOI 10.1016/j.nima.2004.12.036](https://doi.org/10.1016/j.nima.2004.12.036)。[精读笔记](../../papers/sources/2005_Guo_NIMA_A-study-of-the-real-time-deconvolution-of-digitized-waveforms-with-pulse-pile-up-for-digital-radiation-spectroscopy.md)、[本地原文](../../papers/raw/2005_Guo_NIMA_A-study-of-the-real-time-deconvolution-of-digitized-waveforms-with-pulse-pile-up-for-digital-radiation-spectroscopy.pdf)。**本轮局核**：PDF2引言；PDF10讨论；本地PDF页眉保留Article in press。**条件**：LSO/NaI(Tl)及前放模型；闪烁测量10 MS/s/8 bit、脉冲发生器另用250 MS/s/8 bit，不能混用；比较DPD/LLD/LSF。4 Mcps谱例为模拟，不是在线吞吐实测。 **原文支持**：作者区分可确定性恢复的事件与过近/真正同时到达的堆积；后者在附加谱形、时间分布等假设下可作统计谱校正。 **项目采用**：区分统计谱校正与逐事件反演；严格重合的数学例子在本项目另列假设推导。 **缺口**：没有给对所有探测器和方法有效的正数ns界限；统计能谱校正不等于恢复每个真实事件。本文自己的“truly coincident”讨论不代替GAGG可辨识性证明。 **审核**：既有AI精读草稿待人工审核；本轮没有新增科学确认，不扩大其他论文的局部采用意见。


**项目推导，非论文定量结论**：在已知同型固定线性模板、自由幅度且无额外通道/先验时，Δ=0有Q₁h+Q₂h=(Q₁+Q₂)h，单条波形不能唯一确定个体Q或未知K。α/γ若模板线性独立需另论，不能将同型退化推广到所有组合；Δ>0也不保证可稳定恢复。已知K/时刻、固定正定噪声协方差下，无偏无约束GLS的分量方差较孤立情形放大1/(1−ρ²)，仅作条件诊断，不冒充未知K、自由τ或随机形状下的全局界限。本轮未计算GAGG数值信息界限。2 ns采样、95/351 ns衰减、上升时间、记录窗长度和双事件间隔分别定义。

**已有结果复核，非新实验**：L47联合微胞控制、有序物理头＋局部refit两种子，严格双事件成功为4 ns 0/7与0/7，8 ns 1/10与2/10，10 ns 0/11与0/11，12 ns 3/10与3/10，20 ns 2/3与2/3，40 ns 12/17与12/17，80 ns 9/11与10/11，200 ns 10/12与12/12。Q≥500 pe、各Q误差≤20%、各时刻误差≤min(10 ns,真实首光子间隔/4)；额外输出另列。测试复用24母ROOT，种子与重复控制不当独立曝光；各间隔Q/幅比/母样本不同，点间不单调不是连续分辨曲线。小样本0成功不能证明所有方法成功概率为零。均值模板4 ns控制为13/28、8/28，已知K拟合诊断27/28；三者不移用作真实涨落或盲重建保证。当前只能称4–12 ns恢复困难，不能宣称4/8/10/12 ns绝对不可分。

**下一步，尚未执行**：在共同响应、固定τ起点下，分层扫描Δ/Q/幅比/噪声/采样相位，并按独立母ROOT保留测试；先做已知K/时刻的协方差/似然诊断，再比较未知K拟合与ML，真值只进入诊断/评分。冻结空/单事件误报预算（如≤5%额外输出帧），报告双事件Q/t联合效率、漏检/误拆/误并、Q偏差和分辨率、不确定度及Δ50/Δ90；若证据不足不强给阈值。保持足够动态范围，不新增ADC削顶规则。源率与探测率分开，Poisson相邻间隔概率1−exp(−Rdet·d)不是漏检率，1/d不是最高可测率。

**待核/库维护待办**：各论文跨幅比/SNR的独立测试与误报预算；R1采样设置与逐点分母；R5计数比/准确度定义；R6母波形隔离与表文异文；R7最短箱内逐点结果、真实分量真值和LMFit实现；适合本GAGG随机光子响应的双事件可辨识性/Fisher或似然界限。新精读、核作者代码或入库另在papers中处理，本项目不修改共享库。

本地记录：jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb第35节。[逐条证据/条件/哈希](../jupyter/run_records/paper_search_pileup_resolution_20261022_085713/evidence.json)、[原PDF局核清单](../jupyter/run_records/paper_search_pileup_resolution_20261022_085713/pdf-check-manifest.json)。原Notebook/文档在同目录备份；本轮只增加研究文字及现有数据复核表，保留全部原计算与图像，不创建新图号或新Notebook编号。


<a id="alpha-pair-resolution-20261022"></a>
### L49：2026-10-22——条件双事件分辨扫描的实现与独立母样本检验

**采用来源与范围**：落实L48的条件恢复率、幅比/噪声分层和固定误报预算建议；来源位置与适用条件沿用L48 R1–R8，不移植LaBr3/有机闪烁体的ns数值。本轮为项目模拟计算，没有新增文献精读或科学审核；来源原有AI待审状态保持。区分逐事件Q/t恢复、输出计数与统计谱学性能，也区分算法失败和可辨识性结论。

**冻结协议与实际数据**：测试前冻结K=0/1/2变投影拟合、两个已训练有序物理网络及其局部重拟合、全部阈值/惩罚；未训练或据新测试选阈值。固定95/351 ns双τ、2 ns采样、2.048 μs局部电压、充分动态范围。Δ按有效面首光子时间定义，Q为实际雪崩等效pe，不是源或沉积能量。17个间隔为0/2/4/6/8/10/12/16/20/30/40/60/80/120/200/400/800 ns，每条件每间隔40对；Δ=0只记总Q/计数。mean七组共5880帧；240个新单α母ROOT（源2.5/5.486/8 MeV各80，中央入射）离线联合推进微胞产生shared共3720帧；同一中能光学母样本另作linear及shared-noise各840帧。联合数字化不缩放/抽稀光子，也不冒充新原生双事件ROOT。跨间隔/条件复用母样本，总11280帧不等于独立曝光数；新母样本与已用144个旧ROOT分离。

**指标**：各Q误差≤20%、各时刻误差≤min(10 ns,Δ/4)、两真事件一对一恢复且无额外输出；真Q≥500 pe。空窗和相关各单事件组分别≤5%额外输出帧。Δ90为该点及之后全部已测正间隔均≥90%的最小网格点，不插值；仅盲方法且误报预算通过时报告合格Δ90。逐点Wilson95区间不当作全曲线同时保证。漏配、额外输出、邻近匹配定义的误拆/误并标志及条件Q误差分布全部保存。两个网络是同一架构的冻结种子；本轮不是五个新架构。已知K候选选择及已知K/t GLS仅作真值辅助诊断，不是全局性能上界。

**实算结果**：均值参考Q4000/4000 pe、0.1 mV白噪声下，合格Δ90为盲拟合16 ns、两套ML＋拟合均6 ns。真实光子/共同微胞的中能组在80 ns实际Q中位9546/9398 pe；40 ns三者均29/40，60–800 ns各已测点均40/40，合格Δ90均60 ns。低/低组盲拟合为20 ns（36/40；逐点下界支持的阈值为60 ns）；高/高组混合方法均80 ns；先弱后强均120 ns。低/低混合流程因400/800 ns外推退化未达到全扫描Δ90；先强后弱混合流程也未达到。高电荷原生单事件盲拟合46/80出现额外输出，故含高电荷条件即使双事件曲线达到90%，也不报告预算合格阈值；两网络及混合方法的各空/单控制额外输出均0/80，全部单事件80/80恢复。

**重要反例与下一步取舍**：中能8 ns原始网络分别17/40、35/40严格成功，同一输出经固定模板重拟合均0/40；12 ns为29/40、32/40降至5/40、1/40。不是分类/波形相似度，而是同一Q/t联合评分。不能说8 ns普适不可分，也不能把35/40或在测试集逐点择优当作稳定90%新方法。中能20 ns混合种子1的25个失败均通过时间匹配而Q不合格；真实K/t固定模板GLS也是40个Q失败，说明当前模板/涨落处理应优先检查，但该诊断不证明信息不存在。同一光学母样本20 ns的线性叠加盲拟合40/40、共同微胞11/40；微胞实现不同，尚非纯占用因果消融。同一共同实现加白噪声至1 mV，盲拟合仍11/40，混合种子1为14/40（原15/40）。拟采用Q/状态相关响应及形状协方差校准、在开发集选择是否接受重拟合，并固定高Q单事件误拆预算后另用新母样本验证；尚未实施或宣称改善。

**验证与限制**：240个有效ROOT的实际宏、RunInfo、三类种子、传感器/读出配置核对；观测文件只含电压，推断读标签拦截检查通过；所有分量加回暗计数后闭合最大误差约3.87×10⁻¹² mV；共78400条评分、1568个曲线点通过计数/预算检查。六组标量/优化数字化对照的雪崩、末态与最终电压完全一致，再用原标量实现重放强弱顺序8 ns例一致。初次入口spawn分派失败及1个ROOT短写均保留记录；后者以原种子/配置重跑到新时间戳，不挑种子。实际运行源码、宏、权重和环境均留快照。均值Fisher导数近似的短间隔敏感性单列，不改变推断或评分。最低原生单事件Q约2554 pe，500–1000 pe真实弱事件未覆盖；1000 pe仅均值控制。结果限于纯α、固定双τ、已知基线、局部K≤2；未给free-τ、α/γ、实物或随机高率结论，1/Δ90不是最高计数率。

本地过程追加至jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb第36–40节，新增固定图号0004-37-1/37-2/38-1/38-2/38-3/38-4/39-1；Run_Catalog登记240个有效ROOT与保留的失败尝试。[独立实验与重放](../pileup/experiments/alpha_pair_resolution_20261022/README.md)、[完整评分](../jupyter/artifacts/alpha_pair_resolution_20261022/resolution-report.json)、[数据与评分验证](../jupyter/run_records/alpha_pair_resolution_20261022/science-validation.json)。文献库本轮只读，L48待核事项保留。

<a id="alpha-pair-conditional-20261022"></a>
### L50：2026-10-22——条件响应标定、相关拟合和重拟合接受规则

**依据与采用范围**：落实L49提出的响应失配检查。条件均值／非平稳残差的思想参考Fixsen等，*Pulse estimation in nonlinear detectors with nonstationary noise*，[DOI 10.1016/j.nima.2003.11.313](https://doi.org/10.1016/j.nima.2003.11.313)；来源位置为共享库 `sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md`，关键方法对应原文第2页式1–2及第3页示例。该来源是TES模拟，已有笔记仍为AI草稿待审，不是GAGG或堆积恢复性能证明。本轮复用已有阅读记录，没有新增精读或改变审核状态。若将来推广参数相关完整协方差似然，仍须处理logdet及正则化；不能把加权残差最小直接称为最大似然。

**实际实现**：用已暴露的上一轮240个母ROOT开发，明确撤销其本轮“独立测试”角色。按母样本块分成训练／选择各120个母ROOT，无跨组复用。训练利用祖先分量、实际Q和首光子时刻，拟合9项带岭惩罚及平滑的经验响应修正；特征含自身Q、估计邻居Q及相对时刻的固定尺度包络。邻居状态是可从候选Q/t计算的代理，不输入真实微胞状态。未观测尾部不作为零目标，修正保持因果起始。生成器仍固定双τ95/351 ns，未改晶体／SiPM物理参数，未逐事件free-τ。

形状涨落复用既有异方差AR(1)白化器，训练残差标定其参数；每帧C由固定的观测平滑先验形成，各K/Q/t候选共享C，logdet在比较中为常数。这是项目的经验近似，9项特征和AR(1)均不是直接照搬Fixsen的完整能量相关协方差。新均值／新GLS在盲起点±16 ns联合优化时刻和log Q，max_nfev=60，Q边界50–50000 pe；原固定模板使用±12 ns变投影。因此旧／新是流程比较；新均值与新GLS的优化设置相同，可作加权消融。

**比较与选择**：保留原拟合、两个原始网络、原网络＋固定模板、仅波形起点的条件GLS、增加网络起点的条件GLS、新均值／GLS混合及选择性接受流程。两网络为已有同架构的两个冻结种子，本轮没有重训。仅拟合对照继承原固定模板选出的K=0/1/2候选，再作局部优化；尚未对新均值模型完成新的多起点／全局搜索。新K惩罚在开发验证集选为0.003；网络起点版本平均整对成功0.770，仅波形起点0.552。两接受规则均选GLS、间隔下限0、最大log Q改变量10，验证整对成功分别0.696／0.721。这近似接受所有收敛重拟合；验证未支持更严格的间隔／改变量门槛，不能声称门控本身带来上述收益。以上均为开发选择指标，不是独立结果。

**已完成的独立协议**：生成前冻结代码、模型、校准及阈值，随后新生成源能量2.5/5.486/8 MeV各80个中央入射单α母ROOT。五组电荷组合、4/8/12/20/30/40/60/80/120/200 ns、每点40对，以及80空窗和240原生单事件；双事件为新光子的共同微胞离线数字化，未缩放光子或相加独立电压。主指标沿用各Q误差≤20%、各t误差≤min(10 ns,Δ/4)、两事件一对一恢复且无额外输出；Q≥500 pe，空／单各组额外输出点估计≤5%。跨间隔／条件复用母样本，不把帧数当作独立曝光数。

<!-- conditional-independent-results:start -->
**独立实算结果**：中／中组40 ns实际Q中位约9437/9299 pe；8 ns原始网络两种子19/40、31/40成功，固定模板重拟合均0/40，新均值＋GLS均40/40。中／中及高／高组在8–200 ns各已测点均40/40，40/40的逐点Wilson95区间为约91.2%–100%，不是普适或全曲线保证。中／中8 ns的仅均值修正种子1为1/40，仅波形起点的新GLS为0/40，加入网络起点则40/40；条件响应、相关加权和起点搜索共同影响结果，不能归因于接受规则。先弱后强8 ns新GLS两种子均36/40，40–200 ns各点均40/40。完整[逐点评分](../jupyter/artifacts/alpha_pair_conditional_20261022/test-summary.csv)及[全部方法与配对比较](../jupyter/artifacts/alpha_pair_conditional_20261022/test-report.json)保留13流程、30160条逐帧评分；边界样本的百分位bootstrap可能退化，不当作精确概率保证。

**保留的退化与误报控制**：低／低30 ns旧混合40/40，新GLS两种子均34/40；低电荷原生单事件新GLS电荷相对误差中位约−14.24%，整帧76/80，均值修正而无相关加权为−0.68%、80/80，说明该经验协方差并非各Q都合适。先强后弱200 ns两网络分别40/40、36/40只输出一个事件，保持K的重拟合仅0/40、4/40成功；重新选择K的条件拟合36/40，旧盲拟合39/40。后者的高电荷单事件却有34/80额外输出，因此不能只按双事件结果判优。新GLS和条件K选择的空／三档单事件各组额外输出均0/80，逐点Wilson95上限约4.58%。不根据测试结果逐条件拼接方法，也不把种子差异当置信区间。[冻结输出诊断](../jupyter/artifacts/alpha_pair_conditional_20261022/test-diagnostics.json)只解释结果，没有调整推断参数。

**已完成验证**：240个ROOT完成RunInfo／实际宏／种子检查，所有源种子与25份既有规范运行清单分离；生产在冻结之后、推断标签访问拦截通过，观测文件只含电压。分量闭合最大误差3.18×10⁻¹² mV；预先固定的新4 ns样本用原标量数字化器重放，与优化实现逐项完全一致。临时PyTorch环境缺失后，按旧wheel哈希离线恢复2.14.1+cpu；两网络各64份旧观测的Q、t和分类概率与原输出完全一致。三项针对性测试及2320帧／13流程的计数、评分、预算检查通过，见[科学验证](../jupyter/run_records/alpha_pair_conditional_20261022/science-validation.json)。这些结果限于纯α、固定双τ和局部K≤2，不是free-τ、α/γ、实物或随机高率证据。
<!-- conditional-independent-results:end -->

**记录与后续**：独立实验为 [alpha_pair_conditional_20261022](../pileup/experiments/alpha_pair_conditional_20261022/README.md)。过程已追加至本地codex_0004第41–44节，保留旧计算、输出和图号；新增图号0004-41-1、42-1、42-2、43-1。首次开发的平滑前沿问题及修正前结果保留，全部冻结推断采用修正后重跑版本。下一步先在开发数据修正低Q相关权重偏差，引入残差新增候选及K比较，检查纯拟合搜索充分性；再检验新的源能量／位置、α/γ及K≥3随机流。测试后不回调本轮门槛，不提供最高计数率结论；文献库保持只读。

<a id="alpha-pair-residual-20261023"></a>
### L51：2026-10-23——相关权重正则化和残差新增事件候选

**复用依据**：据已有精读笔记，Fixsen等2004的条件均值／非平稳权重方法（DOI [10.1016/j.nima.2003.11.313](https://doi.org/10.1016/j.nima.2003.11.313)，共享库 `sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md`，原PDF2式1–2）支持分别检验均值与权重的适用性；其对象为TES模型模拟，没有给GAGG协方差正则化处方。Zhang等2024（DOI [10.1016/j.nima.2024.169273](https://doi.org/10.1016/j.nima.2024.169273)，`sources/2024_Zhang_NIMA_Development-of-a-pile-up-pulse-recovery-algorithm-for-the-LaBr3-detector.md`，原PDF2–4、式1–5及图3–4）在LaBr₃／R6233 PMT、500 MHz／14 bit条件下使用幅度相关模板和残差增峰，提醒单事件模型失配可能产生幽灵脉冲；增峰阈值和显著改善规则未充分公开。两篇均保持AI草稿待审；本轮读取完整既有笔记及主题／审核入口，没有重新阅读PDF、增加科学审核或修改共享库。`global-review.md` B23-052所列同模板模拟／独立验证及PMT饱和边界继续保留。

**项目实施**：沿用L50的均值响应和两套网络权重，不重训网络；在观测决定且每帧固定的AR(1)精度中比较原ρ≈0.9783及缩减ρ=0.9／0.8，并保留早期GLS／OLS混合探查。该缩减是经验正则化，不是器件参数新测量。残差字典扫描先投影已有事件Q／t切向，提出最多4个新增起点，另筛4个离散分裂起点；用原非线性条件响应联合拟合，先20次函数评估再精修最好的两个至60次，Jacobian调用另计，局部半径24 ns。保留不使用网络起点的版本、保持网络K的版本及原权重扩充搜索消融。上述投影、ρ网格和门槛是项目方案，不冒充论文原算法或全局最优。

**开发／测试划分**：原响应训练母样本作小样本诊断；L50已公开的240个独立母ROOT本轮改作开发选择，每条件／间隔前20对、全部80空窗和240单事件，共1320帧。开发门槛同时限制各控制额外输出≤5%、各单事件成功≥95%、各单事件电荷偏差中位绝对值≤5%；合格方案按五类组合等权整对成功选择。之后冻结参数，另生成240个新单α母ROOT，联合数字化2320帧。严格指标、充分动态范围、固定双τ、Q单位、跨点复用母样本等边界沿用L50；同一测试不得事后逐条件择优或回调阈值。

<!-- residual-independent-results:start -->
**独立结果**：完整开发选择为 `rho90`，主方案每事件惩罚 `0.001`。冻结后新增240个光学母ROOT、2320帧，比较14条流程；两个网络权重未改变。50个双事件扫描点相对同批原GLS（seed1）为改善22点、相同25点、下降3点；跨点复用母样本，这不是50次独立统计检验。

低Q原生单事件的匹配后电荷相对误差中位数为 -13.59% → -3.82%，严格整帧成功 75/80 → 72/80。中／中8 ns整对成功 40/40 → 40/40。上述原GLS对照保持网络给出的事件数。

与原来同样允许未知K选择的条件拟合相比，先强后弱40 ns为 0/40 → 39/40，200 ns为 37/40 → 40/40；原来保持网络K的GLS在200 ns为 0/40。同时保留这两类基线，避免把放开网络事件数带来的收益全部归于新搜索。箭头比较同一批新测试，不拼接上一批数字。

主方案额外输出为空窗0/80、低Q8/80、中Q0/80、高Q0/80，5%点估计预算未全部通过；该预算不是总体误报率的置信保证。完整结果同时保留漏配、额外输出、误拆／误并、匹配Q误差标准差与分位数、逐点Wilson区间及同点配对比较。

主方案未通过全部预设控制要求，完整新流程暂不替换默认方案。相对旧未知K条件拟合，50个扫描点为改善20点、相同25点、下降5点。加入网络起点相对纯波形起点的逐帧整对成败为0胜／0负；本批没有观察到网络起点提高整对成功率。


下降最多的已测点：low_high、8 ns：37/40 → 30/40；low_high、12 ns：36/40 → 34/40；low_high、40 ns：40/40 → 38/40。这些退步未删除或用其他方法替换。
相对旧未知K条件拟合的全部下降点：low_high、8 ns：37/40 → 30/40；high_high、4 ns：37/40 → 32/40；low_high、4 ns：4/40 → 2/40；low_high、12 ns：36/40 → 34/40；low_high、40 ns：40/40 → 38/40。

验证完成：6项针对性测试；128帧旧网络输出精确重放；新样本首对以原标量数字化代码精确复算；240个RunInfo、实际宏和种子核验；推断禁读真值；新旧冻结文件校验。codex_0004第45–48节新增输出已执行，历史单元与图像保留，Run_Catalog新增240条ROOT记录。范围仍限纯α、固定双τ和局部K≤2，不回答free-τ、α/γ、实物能谱或最高计数率。
<!-- residual-independent-results:end -->

**记录与待核**：实验为 [alpha_pair_residual_20261023](../pileup/experiments/alpha_pair_residual_20261023/README.md)，本地过程已追加至codex_0004第45–48节，固定图号0004-45-1、46-1、46-2、47-1。文献方面仍待核Zhang的具体噪声权重／增峰准则及Fixsen插值／正则化实现边界；本项目后续需检查不同位置／连续电荷范围、α/γ和K≥3。本轮不提供free-τ或最高计数率结论。


<a id="alpha-pair-null-gate-20261023"></a>
### L52：2026-10-23——真实弱α、模型权重与单事件经验增峰筛选

**文献与审核边界**：据共享库已有精读笔记，Zhang等2024（DOI [10.1016/j.nima.2024.169273](https://doi.org/10.1016/j.nima.2024.169273)，`sources/2024_Zhang_NIMA_Development-of-a-pile-up-pulse-recovery-algorithm-for-the-LaBr3-detector.md`，笔记定位原PDF3–4及图3–4）提示模板失配会诱发幽灵脉冲，增峰改善标准未充分公开。其条件是LaBr₃／R6233 PMT、500 MHz／14 bit，科学内容保持AI草稿待审；`global-review.md` B23-052的同模板模拟／独立验证及硬件适用边界保留。本轮只读主题、审核入口和既有笔记，未重读PDF、未联网或修改共享库。当前模型权重、Q分档和经验秩均为项目实现，不冒充该文原算法。

**开发与实际采用**：先导16个单αROOT确定源能量1.625 MeV对应实际电荷约713 pe，随后新增80个弱母ROOT，与已公开较强母光子共同数字化出1400双事件及80原生单事件。未缩放光子或相加独立电压。初始增峰筛选使低Q单事件额外输出12/160→3/160，但弱单事件Q偏差仍约−10.10%，因此未满足全部要求。

20个开发弱单事件的加权对照中，普通最小二乘Q偏差约+1%，观测先验GLS约−10%，模型先验GLS约−1%。随后以单事件OLS估计的平滑模型形成每帧固定GLS精度矩阵，ρ仍0.9，均值模板与物理参数不变；所有K/Q/t候选共享矩阵。该诊断支持检查权重，不单独证明唯一的微观偏差来源。

完整模型权重开发集3120帧，包括两批已公开数据各每点前10对及全部控制，以及新弱事件数据。使用560个单事件按观测单拟合Q分档，比较(C1−C2)/C1和(C1−C2)/C0两个搜索后评分；以保守并列和加一修正的上尾秩p≤0.025替代K1→K2固定惩罚，零／单检测保留旧判据。选择约束为空／单额外输出≤5%、各单成功≥95%、各单Q偏差中位绝对值≤5%，再提高强／弱组合等权平均恢复。选择relative_residual；弱单事件Q偏差约−1.06%，弱双事件平均成功约49.7%。控制也参与经验秩标定，不能把开发通过当独立误报认证。

**冻结后的独立检验**：新生成1.625／2.5／5.486／8 MeV各80个母ROOT，十种组合、4/8/12/20/40/80/200 ns，每点40对，另有80空窗与320原生单事件，共3200帧。t指SiPM有效面首光子时刻；Q为实际雪崩等效pe，不能直接叫源／沉积／可见能量。完整模拟与分析保持固定双τ、充分动态范围和局部未知K≤2。两个真实Q均≥500 pe才进入整对分母，同时保留全部帧与阈下数量。逐点要求各Q误差≤20%、各t误差≤min(10 ns,Δ/4)且无额外输出。跨条件／间隔复用母样本，不累计独立曝光。

**独立实算结果**

320个新母ROOT、3200帧、五种流程；本批独立空／单控制全部通过预定要求。

| 控制 | 原额外输出 | 新额外输出 | 原严格成功 | 新严格成功 | 原Q误差中位 / % | 新Q误差中位 / % |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| all | 0/80 | 0/80 | 80/80 | 80/80 | — | — |
| weak | 0/80 | 0/80 | 80/80 | 80/80 | -9.74 | -0.89 |
| low | 6/80 | 0/80 | 74/80 | 80/80 | -3.87 | -1.07 |
| mid | 1/80 | 2/80 | 79/80 | 78/80 | -1.85 | -1.16 |
| high | 0/80 | 0/80 | 80/80 | 80/80 | -1.53 | -1.14 |

| 双事件条件 | Δ / ns | 原观测权重 | 新模型权重＋经验筛选 | 已知K诊断 |
| --- | ---: | ---: | ---: | ---: |
| high_low | 40 | 39/40 | 39/40 | 39/40 |
| high_low | 200 | 40/40 | 40/40 | 40/40 |
| low_weak | 8 | 10/40 | 4/40 | 10/40 |
| low_weak | 40 | 31/40 | 27/40 | 27/40 |
| low_weak | 200 | 37/40 | 40/40 | 40/40 |
| high_weak | 40 | 0/40 | 0/40 | 4/40 |
| high_weak | 200 | 21/40 | 19/40 | 19/40 |
| weak_high | 8 | 0/40 | 0/40 | 12/40 |
| weak_weak | 8 | 16/40 | 21/40 | 21/40 |
| weak_weak | 40 | 33/40 | 34/40 | 34/40 |

70个双事件扫描点中，相对同批原观测权重流程，改善13点、相同32点、下降25点。跨点复用母样本，这不是70次独立统计检验。

包含弱事件的35点等权描述性平均由51.50%变为49.93%；原五种较强组合35点由85.50%变为85.71%。本轮支持单事件电荷偏差及误拆控制改善，尚未显示弱双事件总体恢复改善，不以开发集收益替代独立结果。

控制通过仅指预定的样本点估计要求，不等于所有双事件条件通过或总体误报保证；新流程mid单事件仍有2/80额外输出，完整表保留Wilson95区间。预先固定的两个完整波形示例均未严格成功：40 ns例漏掉弱事件，200 ns例的弱事件时刻偏移约17 ns。

已知K诊断由评分端使用真实脉冲数，不是盲方法或理论上限。全部逐点评分保留在本地test/summary.csv和test/report.json；匹配电荷误差不能代替能量分辨率，成功率分母和阈下数均单独保存。

四个Q标定档实际上来自四种固定源能量，未在档内连续覆盖；当前经验秩不能外推为整个Q区间、不同位置或实物模型失配下的误报保证。主方案仍作为研究候选，不直接替换随机高率全窗默认流程。

验证：416个本轮完成ROOT的哈希、实际宏和种子；227个冻结推断／物理文件及两轮旧快照；开发／测试首对标量与优化数字化器完全一致；10项针对性测试；16000条评分及375条汇总一致性；笔记旧内容、编号与图像保留。

本地codex_0004新增第49–52节，固定图号0004-49-1、50-1、50-2、50-3、51-1；Run_Catalog记录416个完成ROOT及4次中断尝试。

图49-1经展示核对后使用与算法相同的加一／并列上尾秩；原展示图保留在本地figure-history。重放production.py后运行correct_figure.py得到最终图；推断、选择及测试评分未改变。

测试后失败分解：high_weak、40 ns的真实K诊断中，24/40未在10 ns容限内匹配两个时刻，12/40时刻已匹配但电荷误差超标，4/40完整成功。按先时刻后电荷互斥分类，前者也可能有电荷误差；这是诊断而非唯一因果解释，没有据此改动冻结门槛。

**记录与待办**：实验为[alpha_pair_null_gate_20261023](../pileup/experiments/alpha_pair_null_gate_20261023/README.md)，过程记录在本地codex_0004第49–52节。下一步优先核验强脉冲背景中的弱分量均值／协方差及时间—电荷耦合，再扩展连续Q、位置、α/γ与K≥3。文献仍待核Zhang具体增峰规则；经验秩的非规则模型选择适用性如需正式统计论证，另在共享库中安排文献整理，不从本项目改库。本轮不提供free-τ或最高计数率结论。

<a id="alpha-pair-weak-diagnosis-20261026"></a>
### L53：2026-10-26——用真值辅助对照定位强→弱事件的拟合失败

**研究依据与边界**：沿L50–L52的条件响应／模型失配问题继续实验。这里的结果来自项目模拟，不是文献给出的GAGG分辨极限。共享库Fixsen2004笔记（`sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md`，DOI [10.1016/j.nima.2003.11.313](https://doi.org/10.1016/j.nima.2003.11.313)，PDF2式1–2、PDF3图1–3）区分非线性条件均值和非平稳协方差，使用TES模型、每能量5000次随机重复；仍为AI草稿待人工审核。本文不含本项目的低秩截断、正则化、未知K候选搜索或GAGG弱双事件验证。本轮复用既有笔记，未重新精读PDF，未修改共享库。

**实际诊断**：上一轮320个母ROOT的测试已公开，现作为开发资料。预先选十种组合×8/40/200 ns×40对和全部400个空／单控制，共1600帧；无新ROOT或网络训练。固定95/351 ns双τ、2 ns采样、完整2.048 μs模拟电压、充分动态范围。Q为实际雪崩等效pe，t为有效SiPM面首光子时刻。两真实Q≥500 pe、各Q误差≤20%、各t误差≤min(10 ns,Δ/4)且无额外输出才算整对成功。跨点复用母样本，不累计独立曝光。

| high_weak方法 | 8 ns | 40 ns | 200 ns |
| --- | ---: | ---: | ---: |
| 原模型权重＋经验筛选，盲 | 0/40 | 0/40 | 19/40 |
| 指定真实K，原候选搜索 | 0/40 | 4/40 | 19/40 |
| 固定真实t，非线性拟合Q，GLS | 3/40 | 13/40 | 28/40 |
| 从真实K/Q/t启动联合GLS | 1/40 | 8/40 | 26/40 |
| 精确减去真实另一分量再拟合 | 40/40 | 40/40 | 40/40 |

最后一行另使用真实邻居Q/t与弱事件真实初值，返回的强事件直接为真值；绝非盲性能或理论上限。固定时刻时的Q初值由观测总量分配，无真实Q初值。对应OLS诊断完整保留，不能将不同权重或减去分量后的代价混比。

GLS代价重放最大相对差9.73×10⁻¹⁶；40 ns组真值启动后32个失败解均比真值处代价更低。优化噪声后残差下降本身正常，不单独证明模型有偏；结合参数错误、固定时刻与精确分量扣除对照，证据支持检查强分量涨落／条件响应及Q–t耦合，不能只增加起点后就期待消除全部失败。

后验收敛核验：40和200 ns的真值启动GLS均40/40正常停止，最大评估次数25和17，小于120上限；8 ns有1/40达到上限，未选择性重跑。正常停止仅为局部数值条件，不证明全局最优，但40／200 ns的失败不能直接解释为迭代预算耗尽。

**采用与待办**：此处采用诊断顺序和证据边界，继续检验经验协方差的盲实现。尚未证明8 ns可盲恢复、普适最小间隔或最高计数率；free-τ、混合α/γ及K≥3均未验证。运行源码和冻结清单在[诊断实验](../pileup/experiments/alpha_pair_weak_diagnostics_20261026/README.md)，本地Notebook0004第53–54节保留固定图0004-53-1、54-1、54-2，全部逐帧指标与失败记录保留。

<a id="alpha-pair-empirical-covariance-20261026"></a>
### L54：2026-10-26——孤立脉冲经验协方差的四种盲实现及开发否定性结果

**文献条件**：沿用L53列出的Fixsen2004条件均值／协方差思想及AI待审边界；TES模拟不能提供GAGG最小间隔或逐事件恢复性能。这里的低秩分解、对角正则化、份额、Q分档和经验增峰秩均为项目实现。未修改物理参数或共享文献库。

**实际实现**：560个已有原生单事件按replicate奇偶分成280个协方差训练母样本、280个门槛标定母样本，与320个已公开开发验证母样本无种子交叉。训练每Q档40/80/80/80例，真实Q/t仅用于对齐训练残差；从原20 μs记录读取−64到2112 ns的完整观测，不把未观测区当零。中心化残差均值仅保存，不加回条件响应。

以经验方差和残差SVD模式构成正定的对角＋低秩协方差，比较纯对角、16模式×0.5、16模式×0.9、32模式×0.9。推断只以观测电压的单事件OLS初估选择Q档、平移权重，每帧所有K/Q/t候选共享固定C；保持原搜索预算。四方案各拟合1880帧，共7520次盲帧分析；其中280帧只标定门槛，1600帧才用于开发验证。worker禁止读ROOT、truth和case标签。

每个方案使用独立280个单事件、同一(C1−C2)/C1经验上尾秩、p≤0.025及1500/5000/12500 pe分档；额外将旧AR权重用同样280例标定，控制门槛样本量差异。所有开发验证资料已在上一轮公开，不能再称为独立测试。

**预定选择**：各空／单控制额外输出≤5%，各单严格成功≥95%、Q偏差中位绝对值≤5%；弱组合平均须超过两种旧盲基线，较强组合相对原观测权重最多下降2个百分点。两类平均各含五种组合×8/40/200 ns，共15点，和L52七间隔35点平均覆盖不同；跨点母样本复用，不累计独立曝光。

| 同点集方法 | 弱组合平均 / % | 较强组合平均 / % | 空／单控制 |
| --- | ---: | ---: | --- |
| 原观测权重＋惩罚 | 61.00 | 93.17 | 弱Q偏差−9.74%；low额外6/80 |
| 原模型权重＋560例经验门槛 | 58.83 | 93.33 | 通过点估计要求 |
| 旧AR权重＋共同280例门槛 | 58.83 | 93.33 | 通过点估计要求 |
| 经验对角 | 59.67 | 89.50 | mid额外6/80、成功74/80 |
| 16模式，份额0.5 | 57.83 | 86.83 | 通过点估计要求 |
| 16模式，份额0.9 | 56.17 | 87.33 | 通过点估计要求 |
| 32模式，份额0.9 | 56.67 | 87.50 | high额外6/80、成功74/80 |

high_weak、8及40 ns的四候选均0/40；200 ns依次20/40、9/40、8/40、8/40，原观测／模型权重分别21/40、19/40。局部收益并存：weak_high、8 ns由旧AR的0/40变为25/40、24/40、24/40、23/40，但其他方向损失使四方案均未满足整体条件。

**采用结论**：保留四种实现和全部负结果，chosen=null，不升级默认流程，未生成新独立测试或新ROOT。完整报告同时保存漏检、额外输出、误拆／误并和条件匹配Q/t误差；未在10 ns容限匹配的输出也计为额外，不能一概叫作凭空产生的事件。Q为实际雪崩等效pe，误差分布不是能量分辨率。

16项针对性测试通过；重放9200条真值辅助评分、8000条盲开发评分与175条汇总；训练／门槛／验证母样本无种子交叉；冻结源码与原输入哈希保持。图和笔记记录在本地Notebook0004第55–57节，固定图0004-55-1、56-1、56-2，实验为[alpha_pair_empirical_covariance_20261026](../pileup/experiments/alpha_pair_empirical_covariance_20261026/README.md)。本轮也保存120组同门槛样本量的配对变化，作为描述性开发证据。

**待验证**：四档离散Q、每档有限校准样本、线性时间对齐以及单事件初估权重的限制仍在。下一步检查双事件条件均值、联合协方差及强弱分量互相关；候选依赖C的高斯近似需同时处理log det C。模拟真分量可用于训练／诊断，但盲推断仍只能使用观测。候选通过开发约束后再冻结并检验新母样本，随后研究随机高率、K≥3和α/γ；当前结果不提供free-τ或普适分辨极限。

<a id="alpha-pair-joint-likelihood-20261026"></a>
### L55：2026-10-26——联合分量、面积守恒与含logdet的候选依赖拟合

**文献采用边界**：复用共享库Fixsen2004笔记，路径 `sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md`，DOI [10.1016/j.nima.2003.11.313](https://doi.org/10.1016/j.nima.2003.11.313)，PDF2式1–2及PDF3；条件为TES模拟、每能量5000次，AI稿待人工审核。采用条件均值／非平稳协方差的思路，不把原文加权最小二乘等同本项目含logdet的高斯近似，也不把TES结果移作GAGG分辨界限。此次复用既有笔记，未新精读或修改共享库。

**数据与核验**：既有6120帧含5400联合微胞双事件、560原生单事件、160空窗；按母事件replicate整除2再模2分折，各280独立母事件、3060帧。留出折的280原生单事件标定门槛，与另320个已公开开发评价母事件无交叉。不同间隔和条件复用母事件，不累计独立曝光。原生总signal_mV包含暗计数；正式分量按avalanche_owner重放，560例粒子＋暗信号闭合误差最大1.48×10⁻¹² mV，实际metadata配置逐例核对。误含暗计数的初次汇总保留并明确停用，没有进入模型标定。

**实际诊断**：从实际2/20/2 ns读出电路构造固定AR3变换，舍弃前3瞬态样点，分析1021维边缘数据；不修改95/351 ns闪烁双τ。留出high_weak的8/40/200 ns处，互协方差对总方差贡献分别约−0.684%、−0.210%、−0.014%，暂不支持其为主要失败来源；残差主要在相邻点相关。这是项目模拟的条件统计，并非普适独立性结论。

**模型修正**：Q为实际雪崩电荷，线性读出的完整面积应为Q×13.618501389 mV·ns/pe。旧经验条件修正会改变这个积分；新完整方案重新标定均值，并将每个候选分量按完整面积归一化。留出分量RMSE由旧均值4.486 mV，经重标定4.470 mV，到守恒均值4.327 mV。积分残差方差相对未条件化短程协方差降到约0.004–0.005，支持检验总电荷条件化；这些诊断不能代替盲逐事件恢复。

**项目实现**：协方差由候选Q/t和邻居均值决定，包括短程方差／相邻相关／分量交叉项、实际白电子噪声经滤波的带状矩阵及暗噪声。完整方案再用固定0.95份额的秩≤2条件化修正；训练估计约0.995，取0.95保留余量。目标为rᵀC⁻¹r＋logdet C，用带状Cholesky及小矩阵求值，已和独立稠密联合矩阵核对。这里的滤波、面积约束、协方差形式及未知K搜索均为项目推导，不归于Fixsen原文。

比较 `old_mean_nll`（旧均值＋新似然，无积分条件化）和 `joint_conserved_nll`（重标定守恒均值＋新似然＋积分条件化）。NLL可以为负，增峰用max(NLL₁−NLL₂,0)，不用相对残差比；拟合Q按1500/5000/12500 pe分档，以280留出原生单事件的加一／并列经验上尾秩p≤0.025选第二事件。首事件需Q≥500 pe且NLL低于空模型。保留原观测／模型GLS基线，另将旧模型GLS按同一280门槛母事件重标定。推断只读电压，worker禁止读取ROOT、truth和case标签。

每种新流程分析1880帧，开发评价为十种组合×8/40/200 ns×40对及400空／单控制，共1600帧。两类平均各15点；保持L54的恢复、控制与强组合损失约束。比较同时改变代价、候选排序和门槛，完整方案还改均值与积分约束，不能视为单因素因果对照。

实验和快照为[alpha_pair_joint_likelihood_20261026](../pileup/experiments/alpha_pair_joint_likelihood_20261026/README.md)，本地Notebook0004第58–62节，固定图0004-58-1、59-1、60-1、60-2、61-1。初次小样本的K=0空分量轴软件故障发生在写出预测前，原源码保留；v2只修正FIR空轴处理，物理、参数和搜索预算不变。

**仍待验证**：有限离散Q／位置覆盖、短程高斯近似、边界历史、候选搜索和母事件数量。真K诊断只能区分门槛与已有双事件参数失败，不是盲性能或理论上限；尚不能据此回答free-τ最短窗、随机高率K≥3、混合α/γ或普适最小间隔。新文献精读需求仍应在共享库单独安排，本项目仅只读复用。

**开发结果与选择**：两种新方案均通过预定控制。`old_mean_nll`弱／较强组合15点平均为65.00%／92.33%，`joint_conserved_nll`为65.17%／94.17%；原观测GLS为61.00%／93.17%，原模型GLS原560例门槛为58.83%／93.33%，按本轮相同280母事件重标定后为58.33%／93.33%。预定规则选中面积守恒完整方案。其空／weak／low／mid／high额外输出为0/0/0/2/0（每组80例），四档单事件Q偏差中位为−0.069%／+0.024%／−0.044%／−0.051%。

high_weak在8/40/200 ns的盲成功为0/0/22（每点40对），原模型GLS为0/0/19。仅在评分端指定真实K、直接选已有双事件候选后，新方案为0/3/22；不重新拟合，也不使用真实Q/t作为初值。40 ns仍不能仅靠放宽门槛解决。4项数值测试、4800条评分及105条汇总重放通过；完整方案有40/1600最终候选未获优化器成功状态，保留其结果，不以局部数值停止证明全局最优。开发通过后，另冻结新母事件确认实验，不将上述收益当独立性能结论。

<a id="alpha-pair-joint-confirmation-20261026"></a>
### L56：2026-10-26——联合似然的冻结后独立光学母事件确认

**来源与采用边界**：沿L55的Fixsen2004条件均值／非平稳协方差思路，保留原TES模拟、每能量5000次、AI笔记待审等条件。以下均为本项目新模拟与统计，不是文献给出的GAGG性能；没有新增精读或修改共享库。

**固定设计与执行**：L55按预定条件选中 `joint_conserved_nll` 后冻结均值、协方差、完整面积约束、搜索和门槛，新增1.625/2.5/5.486/8 MeV单α源能量各80个母ROOT，共320个；源能量不是沉积／可见能量。源、雪崩、读出和联合数字化随机流在生产前核对旧记录，物理源码及二进制保持。每ROOT一个20 μs窗，用保留光子共同推进微胞构造十种组合×8/40/200 ns×40对，加80空窗及320原生单事件，共1600帧。首对优化数字化与原标量实现逐数组一致，不叠加独立电压、不缩放光子。

初次并行进程收到SIGTERM，原因未确定；2个已核验ROOT保留，4次中断宏／日志／原文件保留，以原预定种子和单进程完成剩余样本。未根据物理结果选择重跑。320个完成ROOT的实际宏、种子、哈希与RunInfo对应核验，全部登记Run_Catalog。

仍为固定95/351 ns、充分动态范围模拟电压、2 ns采样、2.048 μs局部未知K≤2；Q为实现雪崩等效pe，t为有效面首光子。所有1200对的两真实Q均≥500 pe，最低正Q约501.86 pe。各Q误差≤20%、各t误差≤min(10 ns,Δ/4)、无额外输出为整对成功，另外完整保留漏检、额外输出、误拆／误并及匹配Q/t精度。

| 冻结方法 | 弱组合15点平均 / % | 较强组合15点平均 / % | 空／单控制通过 |
| --- | ---: | ---: | --- |
| 原观测GLS | 60.50 | 92.83 | 否：weak Q偏差−9.04% |
| 原模型GLS，原560例门槛 | 58.83 | 92.67 | 是 |
| 原模型GLS，共同280例门槛 | 58.50 | 92.67 | 是 |
| 面积守恒联合NLL | 64.50 | 93.33 | 是 |

**不确定性**：不同条件和间隔重用母事件，未按1600帧独立重抽。用replicate整除2的40个母事件块，前20／后20块分层，各层重抽20块，跨方法、能量、条件和间隔共享权重，2000次配对bootstrap同时重算阈上分母；该分析计划在本批拟合前登记。相对原观测GLS，弱组合差+4.00个百分点，95区间[+1.00,+7.17]；较强组合差+0.50，区间[−1.17,+2.33]。相对原模型GLS，弱组合差+5.67，区间[+2.67,+8.67]。区间条件于当前四档能量、固定响应和分块设计，不含所有模型失配或训练不确定性；逐点另存Wilson95区间。

新方案空／weak／low／mid／high控制额外输出为0/0/0/1/0（每组80例）；对应单事件Q偏差中位+0.076%／−0.007%／−0.021%／−0.043%。这一电荷精度仍不是沉积能量分辨率。

**关键取舍**：weak_high表示先弱α后强α，8 ns由原观测／模型GLS的0/40提高到27/40，是主要收益来源；high_weak在8及40 ns仍0/40，200 ns新20/40，低于原观测23/40和原模型26/40。仅在评分端选已有K=2候选，新high_weak的8/40/200 ns成功为0/3/20，各40对；普通10 ns下时刻匹配两事件分别13/15/27例，但这不等于Q及严格近邻时间都恢复正确。不能按真实强弱顺序逐帧择优，更不能把平均改善写成所有条件改善。

**采用结论与验证**：本批通过预定整体确认要求，支持将该方法保留为当前纯α局部双事件研究候选；不自动替换随机高率全窗默认流程。320个ROOT／实际宏／种子、冻结文件、固定门槛、母事件独立性、6400条评分及140条汇总核验通过。新方法41/1600最终候选未获优化器成功状态，原观测7、原模型2；全部计入结果，未选择性重跑，不能以有限局部搜索宣称全局最优。

实验为[alpha_pair_joint_confirmation_20261026](../pileup/experiments/alpha_pair_joint_confirmation_20261026/README.md)，本地Notebook0004第63–65节保留固定图0004-63-1、63-2、64-1，实际源码、模型、门槛、宏、逐帧指标和失败记录均保留。下一步集中研究强尾中弱事件的时间—电荷可辨识性和涨落模型，再扩展连续Q／位置／器件失配、K≥3与α/γ。尚无普适最小间隔、free-τ最短窗或最高注入率结论。

<a id="alpha-pair-weak-search-20261026"></a>
### L57：2026-10-26——弱分量起点、局部剖面与同预算搜索对照

**文献与方法边界**：沿L55–L56复用Fixsen等2004既有精读笔记，来源为 `sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md`（DOI [10.1016/j.nima.2003.11.313](https://doi.org/10.1016/j.nima.2003.11.313)，PDF2式1–2、PDF3；TES模拟每能量5000次，AI草稿待人工审核）。本轮未重新精读PDF或修改共享库。只沿用区分条件均值与非平稳噪声的思路；以下剖面、搜索规则及GAGG结果属于本项目，文献不提供其恢复率或普适分辨下界。

**固定方案**：保持面积守恒响应、候选依赖协方差及logdet、95/351 ns双τ和2 ns×1024点完整模拟电压。Q为实现雪崩等效pe，t为有效SiPM面首光子。旧联合似然开发批与已公开的确认批各1600帧均作为本轮开发资料；280个门槛母样本与训练280、两批评价各320个母样本四组互不交叉。每种新搜索3480帧，280帧只标定门槛。

`legacy_extra`使用原15%/50%/85%份额，`weak_extra`另加1%/3%/5%/10%及其对称份额。两种方向、4/8/12/20/40/80/160/320 ns间隔保持；每个间隔／方向取初始NLL最好的份额，每方向选四个起点，均至多八次20迭代筛选、两次70迭代精修，原最佳双事件候选保留。弱份额方案初始代价计算更多，但额外非线性优化预算相同。K=0/1及原始检出保护不变，worker禁止读真值，两个方案各自在相同280原生单事件上重新标定上尾秩门槛。

**预定开发选择**：latest批high_weak、40 ns至少增加10个百分点；older批40 ns和两批200 ns不得退步。两批weak_high、8 ns最多下降5个百分点，弱组合平均不下降、较强组合平均最多下降2个百分点；空／单额外输出≤5%，各单成功≥95%，Q偏差中位绝对值≤5%。整对仍须两真实Q≥500 pe、各Q误差≤20%、各t误差≤min(10 ns,间隔/4)、无额外输出。满足要求的候选才进入新母样本确认；不按真实类型或强弱顺序逐帧择优。

**已完成诊断**：latest的high_weak、weak_high各120对进行真实Q/t初始化和固定真实时刻拟合Q。另预定high_weak每间隔第0、1对，固定弱Q/t、重新优化强Q/t，得到六张13×9局部图；40 ns第0对细化至25×17。扫描范围和一个强事件初值含真值，局部最低值不冒充全局最优或盲性能，不按Wilks近似标置信区间。1127个网格点代价重放最大差5.46×10⁻¹²，共用网格点差为0；16个局部 nuisance 优化点未获成功停止标志，保留并注明。

| latest high_weak诊断 | 8 ns | 40 ns | 200 ns |
| --- | ---: | ---: | ---: |
| 原已有双事件候选，评分端指定K | 0/40 | 3/40 | 20/40 |
| 真实Q/t初始化联合拟合 | 12/40 | 9/40 | 29/40 |
| 固定真实时刻，仅拟合Q | 8/40 | 12/40 | 31/40 |

事后仅对已有候选做排序审计：把真实初始化及固定时刻得到的解加入原候选，再按同一NLL选最小值，high_weak仍为0/3/19（各40对）。这是额外的真值辅助诊断，没有新增拟合或改动盲实验；说明这些样本中参数更接近真值的候选未必代价最低。噪声数据的拟合值低于真值点本身正常，不能仅据此宣称模型有偏或达到物理分辨极限。

**完整盲对照结果**：两种新搜索各3480帧已全部完成，并分别重新标定门槛。

| 方法 | older弱组合平均 / % | latest弱组合平均 / % | older较强组合平均 / % | latest较强组合平均 / % | 两批空／单控制 |
| --- | ---: | ---: | ---: | ---: | --- |
| 原联合NLL | 65.17 | 64.50 | 94.17 | 93.33 | 通过 |
| 增加原份额起点 | 65.00 | 64.33 | 94.00 | 93.17 | 通过 |
| 补充弱份额起点 | 64.67 | 64.17 | 93.83 | 93.33 | 通过 |

latest批high_weak在8/40/200 ns三方法均为0/0/20（各40对）；older批200 ns为22/22/21，其余两个间隔仍为0。weak_high、8 ns保持older的24/40与latest的27/40。弱份额方案相对原NLL的弱组合平均差，older为−0.50个百分点、母事件配对95区间[−1.00,0.00]，latest为−0.33、区间[−0.83,0.00]；均为开发描述性区间，不是独立确认。

**候选与门槛诊断**：latest批high_weak、40 ns原K=2候选有37/40在参数精度上已经失败；其余3个正确候选的拟合Q都≥500 pe，但增峰上尾秩p约为0.383、0.321、0.222，未通过0.025门槛。这3例不是Q下限造成的拒绝，单独降低电荷下限不能补回；放松证据门槛必须另验误报。补充弱起点后仍3/40个参数正确候选，盲恢复仍0/40。该事后审计不改变预定方法或选择规则。

**采用与验证**：两种新搜索均未通过预定要求，本轮不升级方法，也未生成新母ROOT。7项数值检查、9600条评分、210条汇总、247条诊断重放通过，原输入／冻结源码／母事件独立性及各搜索门槛均核验。候选代价重放差为0；全部未收敛结果保留：原／增加原份额／补充弱份额最终选中候选未收敛数为older的40/43/40、latest的41/44/43，各1600帧。额外计算中位约2.46/2.00秒每帧，仅是缓存原搜索之后的增量，含复核且并发阶段不同，不作严格速度排名。

首次完整并行两任务收到SIGTERM，原因未确定；240条盲预测和214条诊断检查点原样保留，降低并发继续预定未完成任务，未按结果选择重跑。实验为[alpha_pair_weak_search_20261026](../pileup/experiments/alpha_pair_weak_search_20261026/README.md)，本地Notebook0004第66–70节保留六幅固定图号和全部结果；图号、原340个历史单元和历史计算均保留。报告v2仅补子图字母与解释，首版报告源码和输出备份保留，物理及推断快照不变。

**下一步待验证**：在相同Q、间隔和评价标准下，比较由当前拟合均值／协方差生成的匹配模型数据与原生光子／共同微胞数据，区分该统计模型内的信息／搜索限制与响应／噪声近似问题。此项后续执行记录见L58；不从现有局部剖面推出实物探测器绝对分辨极限，也不直接推广到free-τ、α/γ或随机高率K≥3。

<a id="alpha-pair-model-closure-20261026"></a>
### L58：2026-10-26——原生与匹配高斯模型的恢复对照

**文献及审核边界**：只读复用共享库 `sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md`，DOI [10.1016/j.nima.2003.11.313](https://doi.org/10.1016/j.nima.2003.11.313)，PDF2式1–2、PDF3；TES模拟、每能量5000次，AI草稿待人工审核。本轮从主题／审核索引与相关方法页定位后复用已有精读，没有重新检查原PDF或修改共享库。只采用条件均值与非平稳协方差思路，原文加权最小二乘不等同本项目带logdet的高斯似然；下面精确采样、交叉门槛及GAGG结果均为项目实现。

**预定设计**：从已公开latest批选760原生帧，含mid_mid、high_weak、weak_high×8/40/200 ns×40对及400空／单控制。在各帧数值相同的Q/t上生成两次独立条件高斯噪声，共1520评价帧；另280门槛母事件参数点各一帧，共1800次新盲拟合。原生评价320母事件与训练280、门槛280互不交叉，条件／间隔复用母事件；重复抽样不增加光学母事件。没有新ROOT或独立实物测试。

保持95/351 ns双τ、2 ns×1024点完整模拟电压、充足动态范围、局部未知K≤2。原生Q为实现雪崩等效pe，高斯Q为数值相同的条件参数，不是新实现雪崩电荷或沉积能量；t为有效面首光子。用固定FIR后1021维协方差C=B−VVᵀ的精确平方根采样，逆变换时前三个残差样点为零。模型未给出该边界的联合分布；另做16个预定FIR零空间前缀扰动，而不凭空补窗前模型。

主比较使用同一原生280单事件门槛，次比较使用独立280高斯单事件标定门槛并交叉应用。两者的增峰统计、拟合Q分档、加一／并列上尾秩p≤0.025、Q≥500保护均相同；未用评价真值调门槛。均值、候选依赖C及logdet、搜索预算不变。推断worker只接收电压，禁止读取ROOT、truth和case标签。

**实际盲结果**：严格成功仍为两真实Q≥500、各Q误差≤20%、各t误差≤min(10 ns,间隔/4)、无额外输出。共用原生门槛得到：

| 组合与输入 | 8 ns | 40 ns | 200 ns |
| --- | ---: | ---: | ---: |
| mid_mid，原生 | 40/40 | 40/40 | 40/40 |
| mid_mid，高斯 | 80/80 | 80/80 | 80/80 |
| high_weak，原生 | 0/40 | 0/40 | 20/40 |
| high_weak，高斯 | 0/80 | 0/80 | 42/80 |
| weak_high，原生 | 27/40 | 34/40 | 40/40 |
| weak_high，高斯 | 56/80 | 75/80 | 80/80 |

换用高斯门槛后，上表严格成功数全部不变，但原生high单控制额外输出0/80→7/80；高斯high控制2/160。原生门槛下，控制仅mid有原生1/80、高斯2/160额外输出，其余组为零，四档单Q偏差中位绝对值均小于0.10%。误拆／误并、漏检、额外输出和匹配后Q/t分布保留；额外输出包含时间未匹配情形，不都等于完全无信号时的假脉冲。

**诊断**：原生240对复用此前诊断，高斯480对做同预算真值辅助拟合。high_weak的已有K=2／真实初始化／固定真实时刻，在40 ns分别为原生3/9/12（各40对）、高斯2/24/29（各80对）。高斯8 ns为0/12/15，200 ns为42/56/56。40 ns高斯有21对在真实初始化后由错变对，但同一似然代价没有更低；仅扩大搜索、降低NLL不保证恢复准确。这些结果使用真实K或Q/t，不是盲性能，也不是优化器／信息论上界。

**噪声和不确定性**：真参数白化的high_weak、40 ns平方范数／1021中位为原生0.9967、高斯0.9919，偏度中位却为0.251与−0.0049；原生仍有正偏、较重正尾和少量滞后相关。总平方范数接近并不证明模型完全准确，但模型失配不是匹配高斯短间隔失败的充分解释。χ²₁₀₂₁仅作真参数高斯参照，不是拟合后约化χ²。

40个母事件块、前20／后20分层、2000次共同权重bootstrap保留跨条件复用。high_weak、200 ns高斯相对原生差+2.5个百分点，95区间[−17.5,+18.75]，不能宣称两域等价。全成／全败点差值重抽区间[0,0]只是样本退化；逐帧Wilson仅补充描述，不消除母事件依赖或模型不确定性。

**采用与验证**：不更换默认算法或门槛。5项数值测试及四参数点各4096次采样检查通过，1800个盲预测、16个边界预测、4560条盲评分／56汇总、2880条辅助评分／48汇总、2280帧真参数残差核验。盲候选代价重放差0，辅助候选最大相对差1.99×10⁻¹⁴；原生组与上一轮一致。16个前缀扰动均未改变输出K或严格成功。原生门槛下最终选中候选未收敛为原生15/760、高斯32/1520，全部保留。

实验与冻结快照为[alpha_pair_model_closure_20261026](../pileup/experiments/alpha_pair_model_closure_20261026/README.md)，本地Notebook `jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb` 第71–76节，固定图0004-72-1、72-2、73-1、74-1、75-1。完整0–2048 ns示例预先指定40/200 ns的第0对；原生弱雪崩分量和高斯条件弱均值区分标注。已有364单元及图像保持。

**下一步待验证**：拆开第二事件检出与两个Q/t的精度，量化弱分量时间—电荷联合不确定性，再比较有针对性的估计／搜索修正。高斯匹配不佳不能证明所有算法都不可能分离，真值辅助改善也不能冒充盲性能；尚不回答普适最小间隔、free-τ最短窗、α/γ或随机高率K≥3。新文献精读需求仍在共享库单独处理，本项目不修改该库。

<a id="alpha-pair-information-20261027"></a>
### L59：2026-10-27——恢复阶段、电荷分配与局部信息诊断

**文献采用及审核状态**：本轮只读共享库索引、审核入口、相关方法页及已有sources笔记，没有联网、重查原PDF或修改文献库。

- S. N. Liddick、I. G. Darby、R. K. Grzywacz，2012，*Algorithms for pulse shape analysis using silicon detectors*，NIMA669，70–78，DOI [10.1016/j.nima.2011.12.034](https://doi.org/10.1016/j.nima.2011.12.034)。来源`papers/sources/2012_Liddick_NIMA_Algorithms-for-pulse-shape-analysis-using-silicon-detectors.md`，既有笔记定位PDF3–8、图3–13。约65 µm DSSD、40 MSPS／12 bit、25 µs记录，低率α链，超过8000个由实测单脉冲复制／缩放／平移构造的双脉冲；不是独立原生GAGG双事件。采用双事件识别效率与参数估计精度分开评价的原则，不转移100 ns／250 ns等依赖幅比的成绩，也不把定时散布当最小间隔。AI草稿待人工审核；`global-review.md` B23-022明确没有新增人工科学确认，χ²／logχ²门选异文和实现缺口仍保留。
- Fixsen等2004，*Pulse estimation in nonlinear detectors with nonstationary noise*，DOI [10.1016/j.nima.2003.11.313](https://doi.org/10.1016/j.nima.2003.11.313)。来源`papers/sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md`，PDF2式1–2、PDF3。TES模拟、0.5–5 eV、每能量5000次；采用条件均值／协方差思想。原文WLS不是本项目完整高斯logdet或下面的信息公式，也没有证明GAGG堆积恢复。AI草稿待人工审核，其他论文的人工审核不转移到本篇。

**项目推导，不冒充文献结论**：对条件高斯模型直接求期望score协方差，得到`Fij=mu_i.T C^-1 mu_j + tr(C^-1 C_i C^-1 C_j)/2`。以`t1/(10ns),t2/(10ns),ln(Q1/pe),ln(Q2/pe)`为坐标，必须包含候选相关C的导数。`F^-1`仅作局部正则／正确模型／适当无偏性条件下的参照；本轮没有补齐该公式的一手文献核查，已有6项独立数值检查支持实现，不将公式归给Fixsen。需正式写作时另在文献库安排相应原始统计论证的核查。

**设计与边界**：保持纯α、95/351 ns固定双τ、2 ns×1024完整电压和充足动态范围、FIR后1021维边缘分布。原生Q为实现雪崩等效pe，高斯Q为相同数值的条件参数；t为有效SiPM面首光子。复用L58的760原生／1520匹配高斯评价帧，拆解1080条双事件结果并检查2160条已有真值辅助拟合误差；没有新盲拟合、网络训练或ROOT。所有公开母样本属于开发诊断，跨条件复用和两次高斯噪声抽样不作为新增独立母事件。

444点局部信息使用真实Q/t和已知K=2，绝不输入盲重建：360个实际参数点；80个名义点，强Q=18000 pe、弱Q=500/700/1500/3000 pe、间隔2/4/8/12/20/40/80/120/200/400 ns、两种顺序、首时刻160.37 ns；4个强→弱40 ns、18000/700 pe相位点。已知强Q/t仍保留强事件随机噪声，不等于减去真实强波形。弱Q较大的lnQ尺度不能直接写成同数值的相对标准差；容差±20%对应[ln0.8,ln1.2]，不将曲线过线认定分辨极限。

**实际分层结果**：先强后弱40 ns，原生输出K=2为7/40，两个时刻均在严格容差内2/40，两Q/t严格正确0/40；匹配高斯为9/80、4/80、0/80。8 ns原生为3/40、0/40、0/40，高斯5/80、1/80、0/80；200 ns原生39/40、26/40、20/40，高斯77/80、61/80、42/80。K=2输出不等于正确弱事件检出，严格成功还要求各Q误差≤20%且没有额外输出。普通10 ns匹配、漏检、额外输出、误拆／误并及全部空／单控制原值保留。

**实际局部信息**：下表为360个实际参数点各组的弱lnQ尺度中位，是模型内参照，不是实际拟合相对电荷标准差。

| 组合、间隔 | 全部参数未知 | 两真实时刻已知 | 强事件Q/t已知 |
| --- | ---: | ---: | ---: |
| 先强后弱，8 ns | 0.91495 | 0.88009 | 0.04974 |
| 先强后弱，40 ns | 0.47785 | 0.47425 | 0.04871 |
| 先强后弱，200 ns | 0.19703 | 0.19645 | 0.04429 |
| 先弱后强，8 ns | 0.20906 | 0.19382 | 0.04732 |
| 先弱后强，40 ns | 0.08294 | 0.08195 | 0.04090 |
| 先弱后强，200 ns | 0.03450 | 0.03447 | 0.02782 |

40 ns先强后弱的弱t局部尺度中位3.643 ns，强弱Q局部相关−0.99480；探索性矩阵分拆中只固定强Q得到弱lnQ尺度0.04871，只固定强t仍0.47509。该比较支持电荷分配歧义，不能据此宣称盲算法已经得到这个改善；漏检和定位损失仍存在。名义40 ns四相位弱lnQ尺度0.43765–0.43916、弱t3.105–3.792 ns，仅限该条件，不排除任意相位或插值影响。

**实际误差对照**：原生240对、高斯480对的已有K=2／真实初始化／固定真实时刻三类拟合全部保留，包括失败、低Q和未收敛。弱Q/t误差按真实参数处的逆信息归一化，χ²₂的5.991和标准正态的±1.960仅是局部95%参照区域，未校准成实际盲区间。40 ns高斯先强后弱已有K=2为22/80在二维区域内，真实初始化74/80但严格成功24/80；固定真实时刻73/80在一维区域内但严格成功29/80。参照区域比±20%要求宽，区域内不能冒充逐事件准确恢复。

固定时刻同条件弱Q误差16/50/84%分位：高斯−46.0%／+10.4%／+59.9%，原生−24.3%／+21.5%／+43.1%。真实初始化直接给Q/t起点，且限制t在真值±48 ns，误差集中甚至近零不是盲估计器突破信息下界。模型失配、偏差、非线性／边界、有限搜索以及母事件依赖均限制区域内比例的解释，没有按2160条独立样本给“覆盖率”置信区间。

**完成的验证与记录**：6项数值测试通过，含140000次独立有限似然差分score抽样；三个先导点通过后4进程计算。444点均通过预定差分减半检查，426／9／9点分别用0.0005／0.00025／0.000125步长；最大矩阵变化0.1582%、尺度变化1.9592%，白化相对误差≤1.32×10⁻¹⁵，信息条件数最大约1.75×10⁶。未用伪逆／正则化，无未解决或被剔除点；1080条分层成功与原缓存逐条一致，2280帧主评分和2160条辅助结果输入哈希核验。实际源码、参数与报告分别冻结，实验为[alpha_pair_information_20261027](../pileup/experiments/alpha_pair_information_20261027/README.md)。本地Notebook `jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb` 第77–81节、固定图0004-77-1、78-1、79-1、80-1、80-2，保留旧387单元。

**下一步及待办**：围绕总Q与分配比例的盲估计检验似然谷和局部搜索；换坐标不增加信息，学习或电荷先验只能由训练母样本确定，避免借四档真实能量限制测试答案。连续Q／幅比外推与相同误拆预算比较后，冻结方案并另取独立母样本；本轮不升级默认方法。尚不回答普适最小间隔、free-τ最短窗、α/γ、随机高率K≥3或实物性能。高斯局部信息与有限样本／非规则模型选择界限的一手统计论证，以及Liddick门选异文，需在共享库单独安排核查，本项目不修改该库。

<a id="alpha-pair-charge-allocation-20261027"></a>
### L60：2026-10-27——总电荷／分配坐标与共同剖面搜索的实际比较

**证据与审核边界**：本轮只读共享库入口、主题／审核索引、相关方法页及Fixsen2004完整既有笔记，并复用L59已核的Liddick2012记录，没有联网、重查原PDF或修改文献库。Fixsen等，*Pulse estimation in nonlinear detectors with nonstationary noise*，DOI [10.1016/j.nima.2003.11.313](https://doi.org/10.1016/j.nima.2003.11.313)，来源`papers/sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md`，PDF2式1–2、PDF3；TES模拟0.5–5 eV、每能量5000次，采用条件均值／协方差原则，原文WLS不等于当前完整高斯NLL。Liddick、Darby、Grzywacz，2012，*Algorithms for pulse shape analysis using silicon detectors*，DOI [10.1016/j.nima.2011.12.034](https://doi.org/10.1016/j.nima.2011.12.034)，来源`papers/sources/2012_Liddick_NIMA_Algorithms-for-pulse-shape-analysis-using-silicon-detectors.md`，PDF3–8、图3–13；约65 µm DSSD、40 MSPS／12 bit、25 µs记录及实测单脉冲复制／缩放／平移合成双脉冲，只采用检测与估计分开评价。两篇科学内容仍为AI草稿待人工审核；B23-022没有新增人工科学确认，门选异文保留。当前新坐标、缩放及数值profile是项目方案，不声称是文献原算法或已获GAGG证实。

**预定设计**：280个原生独立单事件门槛帧＋760已公开原生评价＋760匹配高斯draw0评价，共1800帧。每种输入三组合×8/40/200 ns×40对，加空和四档单事件各80帧；draw1未运行。训练280／门槛280／评价320母事件互不交叉，条件间仍复用，不算新独立确认。纯α、固定95/351 ns双τ、2 ns×1024完整电压、充分动态范围、未知K≤2。原生Q为实现雪崩等效pe，高斯Q为相同数值条件参数；t为有效面首光子，不作沉积能量。没有新ROOT、训练或电荷先验，worker阻止打开truth／case／ROOT。

保持均值、协方差、logdet及K0/K1；所有变体保留原最佳K2。component_local与allocation_local都在原候选上最多70次局部迭代；后者用`t1/10,t2/10,100lnS,10r`，S=Q1+Q2，Q1=L+(U−L)r，L=max(50,S−50000)，U=min(50000,S−50)，保持各Q为50–50000 pe。MLE目标不加Jacobian，时刻范围仍起点±24 ns；数值缩放和有限差分的物理步长也改变，不宣称增加信息。

两个profile变体共享完全相同起点：观测单事件时刻前／后4/8/12/20/40/80/160/320 ns，第一Q份额0.01/0.03/0.05/0.10/0.15/0.50/0.85/0.90/0.95/0.97/0.99；每点完整可行域中数值优化lnS，最多24次、容差10⁻⁴，保留原值。各方向选4个间隔的最佳份额，共最多8个起点，各20次筛选后两个候选各70次精修；复用local结果。有限标量与局部搜索不是全局极小值证明。各新搜索用相同280单事件重新校准NLL增峰参考分布，Q分档与加一／并列上尾秩p≤0.025保持；原数值门槛另作敏感性对照。

**实际结果**：各Q误差≤20%、各t误差≤min(10 ns,间隔/4)、无额外输出，真Q均≥500，每点40对：

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

原生40 ns的allocation_profile有31/40候选代价下降，强制K2成功3/40，实际盲成功0/40，分别记录。全部条件和输入的漏检、额外输出、误拆／误并、Q偏差／sigma68误差宽度、t误差、K2与时间匹配、阈下候选及未收敛结果保留。sigma68为匹配条件下相对Q误差宽度，没有沉积能标，不作能谱FWHM。

- `baseline`：原生空／weak／low／mid／high额外输出帧依次为0/0/0/1/0（每组80帧）。
- `component_local`：原生空／weak／low／mid／high额外输出帧依次为0/0/0/1/0（每组80帧）。
- `allocation_local`：原生空／weak／low／mid／high额外输出帧依次为0/0/0/1/0（每组80帧）。
- `component_profile`：原生空／weak／low／mid／high额外输出帧依次为0/0/0/1/0（每组80帧）。
- `allocation_profile`：原生空／weak／low／mid／high额外输出帧依次为0/0/0/1/0（每组80帧）。

**选择要求与决定**：原生空／单额外输出≤5%，单事件严格成功≥95%、Q偏差中位绝对值≤5%；先强后弱40 ns至少增10个百分点、200 ns不减；先弱后强8 ns最多降5个百分点，中／中每点≥95%，六个强弱组平均不下降。四种变体均未同时通过预定原生开发约束，不选择新默认流程。

| 变体 | 开发要求全部通过 | 未通过项 |
| --- | --- | --- |
| component_local | False | high_weak40_gain |
| allocation_local | False | high_weak40_gain |
| component_profile | False | high_weak40_gain |
| allocation_profile | False | high_weak40_gain |

**不确定性与成本**：40母事件块按前20／后20分层，2000次共同权重重抽保留方法／条件／间隔／输入复用。中／中覆盖40块，强弱覆盖前20块、每块2对；不是逐帧独立区间。全成／全败的[0,0]为样本退化，不是总体保证。原生每帧增量：component_local 0.094 s、73次函数调用；allocation_local 0.095 s、71次函数调用；component_profile 3.979 s、3039次函数调用；allocation_profile 4.049 s、3080次函数调用；profile包含独立使用时的共享剖面成本，本次实际只计算一次。时间排除缓存基线，并发逐帧wall time不作实时吞吐指标。

**完成的验证**：7项数值检查和32帧先导通过；1800帧／7200变体候选组、7600主评分／13680含门槛敏感性评分、252汇总、3600强制K2诊断核验。原基线逐条一致，K0/K1不变，K2最低代价不劣化，最大重放相对差0；266981个标量profile点中950个未报收敛，全部保留。源码、模型、配置和报告冻结，未改默认方法。实验为[alpha_pair_charge_allocation_20261027](../pileup/experiments/alpha_pair_charge_allocation_20261027/README.md)。本地`jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb`第82–86节，图0004-82-1、83-1、83-2、84-1、85-1；旧410单元保留。预定40／200 ns第0对显示完整0–2048 ns实际原生分量、条件均值与最终盲输出，选K=1时弱输出画零。

**待验证**：结果限于已公开三组合／三间隔、固定双τ、局部K≤2。连续Q／幅比外推、位置和噪声失配、α/γ、free-τ、高率K≥3及实物均未验证。若同一NLL搜索收益不足，优先检验估计目标、响应／噪声模型和训练得到的先验；不将四档真实能量用作测试先验。符合开发要求的方案仍需扩大条件后冻结、新母样本确认。有限样本估计与非规则模型选择的一手统计论证、Liddick门选异文等继续列为共享库单独核查待办，不从本项目修改文献库。

<a id="alpha-pair-risk-training-20261027"></a>
### L61：2026-10-27——连续Q监督目标与物理拟合的实际比较

**既有证据与审核状态**：只读共享文献库的主题/审核入口、数字脉冲/PSD/光产额页及两篇既有sources记录，本轮未联网、未重查PDF、未修改库。Morad等，2024，*Model-Based Deep Learning Algorithm for Detection and Classification at High Event Rates*，DOI [10.1109/TNS.2024.3371573](https://doi.org/10.1109/TNS.2024.3371573)，来源`papers/sources/2024_Morad_TNS_Model-Based-Deep-Learning-Algorithm-for-Detection-and-Classification-at-High-Event-Rates.md`，PDF3–5式4–9、算法1–2，PDF5–7及PDF10表II。25 mm CLYC、PMT800 V、1 GS/s、17 mV阈值，Poisson线性叠加+AWGN、K≤3；采用检测/到达时刻/分类分开评价，不将6 Mcps检出81.8%、条件分类89.1%或联合72.9%误作本项目电荷恢复率。原文误报14.1%的分母仍有解释缺口。人工审核第553项仅限定确认联合指标区别及适用条件，未全面确认方法；原笔记关于NoA posterior/likelihood及再乘Poisson项的疑问保持。当前K2 logit不当似然，不重复乘源率先验。

Kim等，2023，*Restoring original signals from pile-up using deep learning*，NIMA1055/168492，DOI [10.1016/j.nima.2023.168492](https://doi.org/10.1016/j.nima.2023.168492)，来源`papers/sources/2023_Kim_NIMA_Restoring-original-signals-from-pile-up-using-deep-learning.md`，PDF2–6、表1–2。TexAT Micromegas，原512×40 ns、100点峰对齐输入，约40000组合样本60/20/20，CNN波形MSE；弱Type2幅度误差15.2%与强Type1约5%说明弱目标单列的必要性，不移用为GAGG性能。时间表1.245×40=49.8 ns与正文48 ns差异保留；其第二信号常作为噪声，不等于本项目两个α都要恢复。科学内容仍为AI精读稿待人工审核；关联书目更新不提升审核级别。

**本项目方案及实际执行**：NumPy MLP1408→128→96→773，观测单事件初拟合创新作为输入；只在训练使用模拟K/Q/t。三组相同初始化/小批次，绝对((Qhat−Q)/5000)²、ln(Qhat/Q)²、相对(Qhat/Q−1)²。总损失BCE+3Q+0.4时间CE，3000步、batch64、lr0.0007，初始化270310010、抽样270310011；不是上述文献的原网络或性能复现。9600训练/960验证来自连续Q设计分布，非四档测试真值先验。初版总Q头修正过窄，先导/训练输入核查后、正式训练和评价前另建v2，将exp(±0.1)改成exp(±1)；数据及生成函数一致并有哈希、初版保留。

每500步按验证强制K2严格成功、Q绝对ln误差中位、count BCE选择。选择结果：absolute第3000步，验证强制K2成功30.33%，总训练49.5s；logarithmic第2500步，验证强制K2成功38.00%，总训练48.4s；relative第3000步，验证强制K2成功41.00%，总训练53.6s。混合方案选择`relative`，两种物理精修共用原K2/网络两起点、各70次迭代及±24 ns范围；带惩罚者增加Σ[(lnQ−lnQnet)/σ]²，σ取验证全部双事件68%绝对lnQ误差、下限0.2，当前为[0.20801263325687447, 0.21970885078116778]。这是同一波形导出的数据依赖正则，不是独立贝叶斯先验。网络分数及各输出用同一280原生单事件分档校准加一上尾秩≤.025，K0/K1沿用基线，每候选Q≥500；所有失败保留，worker禁读truth/case/ROOT。

纯α、95/351 ns固定双τ、2 ns×1024完整模拟电压、充分动态范围、未知局部K≤2。原生Q为实现雪崩等效pe，高斯Q为条件生成参数，t为有效面首光子，不是沉积能量。训练/验证/评价种子基数270310001/2/3，每帧base×100000+index，噪声SeedSequence([seed,993])。冻结权重后新高斯1000帧含600对和400控制；已有原生760帧是公开迁移诊断，另280门槛。高斯沿用同一已标定模型，不是新G4/物理确认。

严格整对正确为每Q±20%、每t±min(10 ns,gap/4)、无额外输出。实际结果：

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

| 方法 | 开发要求全部通过 | 未通过项 |
| --- | --- | --- |
| nn_absolute | False | native_extra、native_single_success、continuous_pair_gain、native_pair_no_drop |
| nn_logarithmic | False | continuous_pair_gain、native_pair_no_drop |
| nn_relative | False | native_extra、native_single_success、continuous_pair_gain、native_pair_no_drop |
| mlp_nll | False | continuous_pair_gain、native_pair_no_drop |
| mlp_regularized | False | native_extra、native_single_success、continuous_pair_gain、native_pair_no_drop |

通过预定开发要求的方案：无；不自动替换默认方法。 合计条件配比由本轮设计，不代表源率。完整误差偏差/σ68、漏检/额外输出/误拆/误并、强制K2、母事件块配对区间均保留；σ68是时间匹配条件下相对Q误差宽度，非能谱FWHM。6项数值检查、24帧先导、模型/代价重放及10560主评分/5760诊断通过；不依残差下降或训练loss宣称恢复改善。

**记录与待验证**：`jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb`第87–92节，六幅图0004-87-1、88-1、89-1、90-1、90-2、91-1；冻结实验`pileup/experiments/alpha_pair_risk_training_v2_20261027/`，本地结果/记录在同名artifacts/run_records。旧431单元及此前研究保留。仅一个训练种子，尚需训练随机性、连续能量新光学母事件、位置/噪声失配、完整随机高率和α/γ验证。后续先按时间候选、Q偏差和计数门槛分解失败；如采用学习型先验，须单独校准覆盖/偏差，不使用旧测试四档Q作答案。需要新统计/后验校准文献整理时交由文献库处理，本项目不修改共享库。

<a id="alpha-pair-physical-gate-20261027"></a>
### L62：2026-10-27——物理判决、学习型分数和网络候选的分离对照

**依据与审核边界**：只读复用共享库已有精读及主题／审核记录，本轮未联网、未重查PDF、未修改共享库。Liddick、Darby、Grzywacz，2012，*Algorithms for pulse shape analysis using silicon detectors*，DOI 10.1016/j.nima.2011.12.034，来源为papers/sources/2012_Liddick_NIMA_Algorithms-for-pulse-shape-analysis-using-silicon-detectors.md，PDF3–8、图3–13：约65 μm DSSD、40 MSPS／12 bit、25 μs记录，以选出的实测单脉冲缩放／移位构造8000余双脉冲。仅借鉴“单／双判别与参数精度分别评价”和慢上升／噪声引起的假双控制；等幅100 ns与幅比5的250 ns不是同一条件，不移用为GAGG性能。该科学笔记仍为AI稿待人工审核，参数数目3／4和门选记号疑问保留，不移用其6.8／8.5数值门槛。

Morad等，2024，*Model-Based Deep Learning Algorithm for Detection and Classification at High Event Rates*，DOI 10.1109/TNS.2024.3371573，来源为papers/sources/2024_Morad_TNS_Model-Based-Deep-Learning-Algorithm-for-Detection-and-Classification-at-High-Event-Rates.md，PDF3–5方法、PDF10表II：25 mm CLYC、PMT800 V、1 GS/s、17 mV阈值，Poisson线性叠加＋AWGN、K≤3。仅采用检测、条件分类和联合成功须区分的原则，不把其分类率当逐事件电荷恢复率。人工审核第553项只是限定确认指标区分和适用条件，未全面确认方法；NoA posterior／likelihood、Poisson项和误报分母的疑问保持。当前NLL差、秩校准与候选消融是本项目实现，不是两篇论文的算法复现。

**项目实施**：在相同缓存拟合候选上分离事件数判决：原生先强后弱200 ns的学习型门槛拒绝13个严格正确K2，改用物理似然判决后严格恢复7/40→20/40；40 ns仍0/40。追加局部精修／网络起点＋物理判决，原生均241/360，高斯分别500/600、501/600。门槛修复通过开发控制，但网络相对局部精修的增益未达两域各5个百分点要求；保持原默认流程。

| 流程 | 原生360对 | 已公开连续Q高斯600对 |
| --- | ---: | ---: |
| baseline | 241/360（66.94%） | 500/600（83.33%） |
| learned_gate | 226/360（62.78%） | 482/600（80.33%） |
| local_physical | 241/360（66.94%） | 500/600（83.33%） |
| mlp_physical | 241/360（66.94%） | 501/600（83.50%） |
| mlp_fixed_gate（门槛敏感性） | 241/360（66.94%） | 501/600（83.50%） |

纯α、固定95／351 ns、2 ns×1024、充分动态范围、局部未知K≤2；Q为等效pe、未作沉积能标。复用relative网络的未正则化物理NLL候选；local只取原K2起点的70步精修，mlp额外含网络起点的70步精修，两者都保留初值。没有新优化或训练，也不是等计算预算的网络结构比较。K0／K1及单事件存在规则不变；物理增峰分数max(NLL1−NLL2,0)，每候选Q≥500、正证据、同280原生单事件按拟合单Q分档1500／5000／12500 pe作加一上尾秩≤.025。搜索改变后重新校准同一种规则；沿用原数值门槛单列敏感性，不混进四主流程。物理判决不使用网络logit或学习型Q惩罚；最低代价候选阈下时不偷换高代价过门候选。

**故障分解**：原生200 ns三个库各20/40正确K2，学习型门槛拒绝13个，物理门槛全部接受；配对恢复差+32.5个百分点，母事件块95%区间+17.5至+47.5。40 ns三个库均只有3/40正确K2且都被秩门槛拒绝，最终0/40，显示候选参数与判决两方面问题，不能据此推导普适分辨极限。网络相对局部对照的新增成功仅高斯先强后弱128–256 ns的1/40，组内差+2.5个百分点、区间0至+7.5。每条件2000次bootstrap，种子270327062；原生按母事件块，高斯按独立抽样。区间为已公开开发样本描述，跨条件复用及模型误差不作独立帧推断。

**控制与代价**：原生空／weak／low／mid／high额外输出帧，local和mlp均0／0／0／1／0；高斯local为0／0／0／2／1、mlp为0／0／0／2／0，每组80。高斯基线mid为1/80，故不可称所有误拆都降低。每组额外≤5%、单事件严格成功≥95%、匹配Q偏差中位绝对值≤5%的开发要求均通过。每Q±20%、每t±min(10 ns,gap/4)、无额外为严格整对成功；另保留10 ns匹配下的漏检／额外／误拆／误并及Q偏差、sigma68和时间误差。sigma68非能谱FWHM。网络降低NLL>10⁻⁶的帧为原生92/686、高斯136/920，两域下降量中位为0，不能据数值改善宣称真实恢复。

**记录与验证**：280校准＋760原生＋1000高斯共2040帧全部代价重放，最大相对差0；6项数值测试、32帧先导、7040主评分＋1760门槛敏感性、2880强制K2诊断，旧基线和学习型流程精确重现。采用上一轮缓存、全部未收敛结果保留，推断禁读truth／case／ROOT。单worker、BLAS1的11.48 s只为缓存重放，不作端到端吞吐。实验为[alpha_pair_physical_gate_20261027](../pileup/experiments/alpha_pair_physical_gate_20261027/README.md)，本地jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb第93–97节，图0004-93-1／93-2／94-1／95-1／96-1。旧455单元、输入、模型与此前源码保留。

**下一步及待核查**：优先量化弱Q分配歧义及区间覆盖，在空／单控制下展示恢复率—额外输出曲线，再评估多解或不同估计目标；不再把单纯更多起点当作已证实方向。学习型分数需另行校准而非当作物理似然。新增统计／区间校准的一手文献整理、Liddick门选异文及Morad统计量疑问仍交共享库处理。当前没有新母事件、噪声／位置变化、free-τ、混合α/γ或随机高率K≥3验证，全部评价已公开，不更新默认流程。

<a id="alpha-pair-evidence-profile-20261027"></a>
### L63：2026-10-27——检测门槛取舍与条件电荷相容范围

**证据及审核状态**：本轮通过paper-search只读复用本地文献库的主题、审核及完整来源笔记，没有联网、重查PDF或修改共享库。Liddick、Darby、Grzywacz，2012，Algorithms for pulse shape analysis using silicon detectors，NIMA669,70–78，DOI 10.1016/j.nima.2011.12.034，PDF3–8／图3–13；约65 μm DSSD、40 MSPS／12 bit、25 μs记录，选出的实测单脉冲缩放移位构造8000余人工双脉冲。借鉴检测、参数精度及假双控制分别报告；不移用其数值门槛、100 ns等幅结果或250 ns幅比5结果。来源papers/sources/2012_Liddick_NIMA_Algorithms-for-pulse-shape-analysis-using-silicon-detectors.md，仍为AI精读稿待人工审核，门选记号和参数数目疑问保留。

Baker与Cousins，1984，Clarification of the use of chi-square and likelihood functions in fits to histograms，NIM221,437–442，来源papers/sources/1984_Baker_NIM_Clarification-of-the-use-of-chi-square-and-likelihood-functions-in-fits-to-histograms.md；PDF1、PDF4§5–7、PDF5§8区分点估计、区间覆盖与拟合优度，提醒有限样本需按具体模型研究抽样分布。其Poisson／multinomial直方图模型不是本项目相关电压的统计模型，不移植计数似然或声称本文保证GAGG区间覆盖。该笔记暂无人工科学确认，DOI未由已有PDF核实，不补造。式3模型概率记号等原文疑问保留。本轮的秩扫描、触发局部参数域、剖面与1／4／9层级均为项目诊断设计，非上述文献复现。

**项目实际结果**：网络起点＋物理判决在两个域所有空／单组额外输出帧≤5%的开发约束下，原生40 ns严格成功为0/40；即使允许所有合资格候选通过也只有3/40。放宽判决只能接纳已有正确候选，不能修正其余候选的电荷／时刻错误。ΔJ≤4的160条条件剖面中，13条触及50 pe下界、13条存在多个分离区段；0条真实Q2评分探针把参考最小代价再降低超过0.01。这些数值及有限搜索误差需与包含率一起解释，不将该层级称为95%置信区间。

| 方法 | 开发控制允许的最大α | 原生40 ns：默认→放宽 | 原生200 ns：默认→放宽 | 原生总成功 | 高斯总成功 |
| --- | ---: | ---: | ---: | ---: | ---: |
| baseline | 0.05000 | 0→0/40 | 20→20/40 | 241/360 | 501/600 |
| learned_gate | 0.04878 | 0→0/40 | 7→12/40 | 232/360 | 491/600 |
| local_physical | 0.05000 | 0→0/40 | 20→20/40 | 241/360 | 501/600 |
| mlp_physical | 0.05000 | 0→0/40 | 20→20/40 | 241/360 | 502/600 |

表中α为回看两个域十个空／单控制组、每组80帧均≤5%额外输出时可取最大值，是开发约束边界，不作为新冻结默认。127点为同280标定帧的上尾秩全部可取值及预定0／.025／.05／.1／.25／.5／1，分档样本40／80／80／80，加一分母41／81。除α外候选、Q≥500、存在规则和物理正证据均不变；α=1不是无条件强制K2。物理及学习型门槛分别保留原校准。严格恢复每Q±20%、每t±min(10 ns,gap/4)、无额外；普通10 ns匹配的漏检／额外／误拆／误并及Q偏差／sigma68、时间误差、未收敛完整保存。sigma68不是沉积能谱FWHM。

**条件剖面**：原生先强后弱40／200 ns及高斯24–64／128–256 ns各40条，共160条，全部预定、不按成功挑选。强制K=2，Q2为按时间排序的后槽；每固定Q2重新拟合Q1和两时刻，使用原相关高斯J=rᵀC⁻¹r+ln|C|。两时刻均在观测单拟合触发附近[−96,+512] ns、裁至记录内，明确不同于旧候选全记录域；旧候选域外与真时刻域外分别计数。基础33个对数Q点50–50000 pe，加500、原Q2及混合Q2的五个比例点，正反扫描、每点两个观测起点各50迭代；保留初值、未收敛及固定其他参数切片，未证明全局最优。分离相容区段保留，lnQ插值端点仅为近似。

| 输入／间隔 | ΔJ≤4包含真实Q2 | 范围宽度／真实Q2中位 | 触及50 pe下界 | 分离区段 | 真Q探针降低参考最小值 |
| --- | ---: | ---: | ---: | ---: | ---: |
| native／40 | 36/40 | 2.693 | 7/40 | 12/40 | 0/40 |
| native／200 | 39/40 | 0.751 | 0/40 | 0/40 | 0/40 |
| continuous／middle | 37/40 | 1.246 | 6/40 | 1/40 | 0/40 |
| continuous／long | 35/40 | 0.563 | 0/40 | 0/40 | 0/40 |

ΔJ≤1／4／9为预定代价层级，不对应本项目已验证的68／95／99%覆盖。主剖面冻结后，另以真实Q2固定的评分探针重新拟合其余参数；不提供真t或真Q1，不反馈到主剖面、门槛或盲输出。真实Q处实际代价用于经验包含率；如果低于主参考最小值>0.01，标明格点／搜索遗漏，不能隐藏。区间宽度仍需加密和独立重复抽样检验。报告全部40条与原物理门槛选K2子集，避免忽略选择影响；Q真值只用于评分／画图归一。

**条件与验证**：纯α、固定95／351 ns、2 ns×1024完整模拟电压、充分动态范围、Q为等效pe而非沉积能量。门槛沿用280原生单事件；原生760、高斯1000评价均已公开。本轮无新训练／随机抽样／母ROOT。6项数值检查、4条先导、160剖面与160评分探针，主曲线代价重放、原α=.025四流程7040条逐帧评分精确复现；127×7040只是同一端点缓存的894080次选择，不是新增独立样本。所有未收敛保留，主剖面worker禁读truth／case／ROOT；真实Q2仅在另行标记的评分阶段传入。最多8核，实际主剖面7worker／BLAS1。

实验为[alpha_pair_evidence_profile_20261027](../pileup/experiments/alpha_pair_evidence_profile_20261027/README.md)；本地jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb第98–102节、图0004-98-1／98-2／99-1／100-1／101-1。五幅图含四条预定完整0–2048 ns波形；图101-1为已知K2相容模型，不能当作最终盲恢复。旧475单元与所有既有输入／结果保留。

**后续与文献待办**：先检查格点／搜索稳定性，再在冻结后的独立模拟重复抽样中校准所需覆盖并检验原生迁移；必要时给出多解与不可靠标志，禁止真值挑解。正式区间的有限样本／非规则模型选择一手统计证据、Baker文献DOI及相关边界理论可由共享库另行整理，本项目不修改库。新光学母事件、训练随机性、位置／噪声失配、α/γ、free-τ和随机高率K≥3尚未验证。

<a id="alpha-pair-profile-stability-20261027"></a>
### L64：2026-10-27——条件电荷剖面的数值稳定性审核

**文献依据与审核范围**：本轮通过paper-search只读复用共享库的主题、审核索引及两篇完整精读笔记，没有重新读取PDF或联网。Baker与Cousins（1984），Clarification of the use of chi-square and likelihood functions in fits to histograms，NIM221,437–442，已有笔记定位PDF1、PDF4§5–7、PDF5§8：点估计、区间覆盖及拟合优度应分开，有限样本需按具体模型研究抽样分布。其对象为Poisson／multinomial直方图，不是本项目的相关电压；不移用其计数似然或渐近覆盖结论。来源papers/sources/1984_Baker_NIM_Clarification-of-the-use-of-chi-square-and-likelihood-functions-in-fits-to-histograms.md，仍为AI草稿待人工审核，暂无科学采用确认；DOI未由已有笔记核实，式3概率记号疑问保留。

Liddick、Darby与Grzywacz（2012），Algorithms for pulse shape analysis using silicon detectors，NIMA669,70–78，DOI10.1016/j.nima.2011.12.034，已有笔记定位PDF3–8／图3–13，约65μm DSSD、40MSPS／12bit、25μs记录，8000余人工叠加样本。仅沿用检测、参数精度与假双控制分开报告的思路；不移植100ns、幅比5或ADC限制。来源papers/sources/2012_Liddick_NIMA_Algorithms-for-pulse-shape-analysis-using-silicon-detectors.md，仍为AI草稿待审，门选记号及自由参数数目缺口保留。增加起点、网格加密及本轮0.05／5%数值筛查门槛是项目预定诊断设计，非这两篇论文提供的定理或性能保证。

**项目实际计算**：32条预定波形中，32条在重点ΔJ≤4层级通过末次加密的数值筛查；30条在ΔJ≤1／4／9三层级同时通过。增加搜索和加密网格没有改善这批样本的已知K2严格恢复数，原生40ns组3/8条的分离区段仍存在。0条真实Q2评分探针把最终参考代价降低超过0.01；该探针只用于评分，不反馈主剖面。仍有数值变化或搜索遗漏，不能把当前范围当作已稳定的置信区间。本轮完成数值审核，独立覆盖率校准留待处理这些变化后开展。

| 条件 | 最后加密三层级均稳定 | ΔJ≤4宽度中位：旧→新／真实Q2 | 多区段：旧→新 | 已知K2严格正确：旧→新 |
| --- | ---: | ---: | ---: | ---: |
| native／40 | 8/8 | 12.809→14.237 | 3→3/8 | 1→1/8 |
| native／200 | 8/8 | 0.928→0.950 | 0→0/8 | 1→1/8 |
| continuous／middle | 6/8 | 1.423→1.394 | 0→0/8 | 1→1/8 |
| continuous／long | 8/8 | 0.736→0.764 | 0→0/8 | 6→6/8 |

末次加密在ΔJ≤1／4分别为32/32和32/32通过；未通过的是continuous第160条（ΔJ≤9，端点变化6.61%）、continuous第190条（ΔJ≤9，端点变化5.40%）。这两条的区段数、边界标志及参考代价不变，只是端点变化超过预定5%。

**必须按同一批样本比较**：原生40ns组本轮8条的范围跨度／真实Q2中位数为12.81→14.24；上一轮全部40条的旧中位数为2.69。前两者才是加密前后比较，不能说从2.69恶化到14.24。系统抽取8条不保证统计代表性；跨度也包含分离区段之间的空隙。四组已知K2的严格恢复数仍分别为1／1／1／6（各8条），增加搜索没有改善这批样本的逐事件准确度。

本轮额外27180次局部优化中，1248次保留未收敛标志；这不是失败波形条数。新真实Q2评分探针有1条未收敛，取旧／新较低代价后实际选用的探针有0条未收敛。有限搜索仍可能遗漏更低解。

另做了**评分后的一致性检查**：图中区段来自网格的lnQ插值，真Q包含计数来自实际拟合探针，两者不能定义混用。本轮32条×3层级共96个包含判定有0个不一致；插值ΔJ与探针ΔJ的最大绝对差为0.369。包含判定相同不等于整条连续曲线已收敛。这项检查不调整主曲线，也不是新的独立覆盖率检验。

每组原有40条，预先按第0／5／10／15／20／25／30／35条系统取8条，不按本轮拟合好坏选样。旧→新均为同一32条；不能直接与上一轮全部160条的中位数比较，也不认证未重算的128条。四个条件为原生先强后弱40／200ns及条件高斯24–64／128–256ns；两个域的Q／t分布不相同，不能把差异全归模型失配。

**数值方案**：保留原条件均值、候选依赖相关协方差、残差二次型加logdet，以及Q为50–50000pe、两时刻在观测单事件触发[−96,+512]ns裁至[0,2046]的参数域。先在原网格每方向用6个分离的观测起点、最多150迭代；再在每单元加几何中点，最后在1／4／9跨界单元及最小值两侧继续加几何中点，新节点每方向4起点／150迭代。起点来自观测缓存、旧剖面、邻点及预定触发偏移／间隔，按代价与数值坐标距离选择；没有传入真t／真Q。所有旧点的更优代价保留，阶段最优点的固定其余参数切片也纳入全节点包络，所以网格步骤不只是插值、也可能改善旧点。各阶段计算量不同，未证明全局最优。

末次加密的每层级仅在参考最小值改善≤0.05、相容区段数及边界标志不变、每个配对端点乘性变化≤1.05时记数值筛查通过。该规则先于计算冻结，是项目容差而非文献定理或严格数值误差界；没有穷举起点或再次全域四分网格。区段间隙保留，lnQ线性插值端点仍为近似。真值只用于评分和画图归一。

主结果哈希冻结后，另固定真实后槽Q2做评分探针，6观测起点／150迭代，未知Q1和两时刻；保留旧／新评分探针的较低实际代价。探针不反馈到主曲线、默认候选或门槛；低于最终参考>0.01单独记搜索遗漏，不能隐藏或解读成更好的覆盖。包含率不称置信水平。完整已知K2评分同时保留Q偏差／sigma68、时间误差、漏检、额外、误拆／误并和严格恢复；sigma68非能谱FWHM，Q未标定为沉积能量。

**复现与验证**：纯α固定95／351ns、完整1024×2ns电压、充分动态范围不变。本轮没有新ROOT、随机样本或训练；最多8核，先导4worker、主计算7worker且BLAS1。6项数值检查覆盖解析相关二次型、网格、分离区段与保留旧点；主曲线所有实际节点及评分探针重放代价，worker阻断truth／cases／ROOT访问。共享盘全仓库Git差异扫描的准备中断已记录；保留备份和协议，改为限定路径的版本差异记录后完成快照，不影响数值配置。

实验[pileup/experiments/alpha_pair_profile_stability_20261027](../pileup/experiments/alpha_pair_profile_stability_20261027/README.md)保存实际源码、模型／参数与协议哈希；本地Notebook为jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb第103–107节，图0004-103-1、104-1／104-2、105-1、106-1。原497单元、既有160条剖面和盲判决完整保留。

**待办**：先针对变化区段检查数值搜索及响应插值，冻结规则后再做独立重复抽样覆盖校准与原生迁移。小样本、条件K2、边界与选择不能忽略；尚未研究free-τ、α/γ、高率K≥3或普适分辨极限。正式有限样本／非规则剖面区间的一手证据及Baker DOI仍可由共享库另行整理，本项目不修改库。

<a id="alpha-pair-contour-coverage-20261027"></a>
### L65：2026-10-27——边界括区细化与固定模型覆盖率先导

**依据与审核范围**：本轮只读复用共享库的主题／审核索引及已有完整精读笔记，没有重读PDF或联网。Baker、Cousins（1984），Clarification of the use of chi-square and likelihood functions in fits to histograms，NIM221,437–442：已有笔记定位PDF1、PDF4§5–7、PDF5§8，区分点估计、区间覆盖与拟合优度，有限样本按具体模型研究抽样分布。其研究对象是Poisson／multinomial直方图，不是相关电压；本文没有提供本项目数值阈值或边界算法。来源papers/sources/1984_Baker_NIM_Clarification-of-the-use-of-chi-square-and-likelihood-functions-in-fits-to-histograms.md，AI草稿待人工审核、无科学采用确认，DOI未核实，式3概率记号疑问保留。

Liddick、Darby、Grzywacz（2012），Algorithms for pulse shape analysis using silicon detectors，NIMA669,70–78，DOI10.1016/j.nima.2011.12.034：已有笔记PDF3–8／图3–13，约65μm DSSD、40MSPS／12bit、25μs记录、8000余人工叠加样本。仅沿用检测、参数精度和假双背景分开评价；不移植100ns、幅比5或ADC限制。来源papers/sources/2012_Liddick_NIMA_Algorithms-for-pulse-shape-analysis-using-silicon-detectors.md，仍为AI草稿待审，门选记号和自由参数数目缺口保留。数值括区1%／代价0.01、两参数点及抽样数量是项目设计，不是文献性能保证。

**实际结果**：5条定向数值诊断中，5/5通过新增的边界括区与更强局部重拟合检查。随后冻结方案，生成32条新条件高斯双脉冲和12条空／单控制；32/32条双脉冲通过数值筛查，全部32条都保留在主要包含率分母中。40ns组ΔJ≤4范围包含真实Q2为14/16，已知K2严格整对恢复为3/16。200ns组ΔJ≤4范围包含真实Q2为16/16，已知K2严格整对恢复为10/16。这检验的是冻结模型内、固定参数点的重复抽样表现；每组16条不足以认证95%覆盖，更不能推广为原生SiPM的盲恢复或普适最短间隔。

两条的变化来源并不一样：第160条仅同网格加强局部搜索，ΔJ≤9端点就移动6.72%，原节点最大代价下降3.021；后续加点与复核再改变约0.071%。第190条同网格端点变化可忽略，后续加点与复核改变约3.45%。因此前者明确存在原节点的局部搜索不足；后者主要在后续细化中变化，但该步骤也重新优化了未知参数，不能全归因于画图插值。

5条最终括区均小于1%，更强的末次局部搜索通过筛查；最小代价的总改善仍很小，5条的严格整对恢复结果均未改变。这是定向数值诊断的结论，不是总体恢复率或全局最优证明。

| 间隔/ns | ΔJ上限 | 范围包含真Q2 | 该比例的精确二项95%区间 | 跨度/真Q2中位 | 分离区段帧数 |
| --- | --- | --- | --- | --- | --- |
| 40 | 1 | 9/16 | 29.9%–80.2% | 0.780 | 2/16 |
| 40 | 4 | 14/16 | 61.7%–98.4% | 1.605 | 5/16 |
| 40 | 9 | 16/16 | 79.4%–100.0% | 4.049 | 6/16 |
| 200 | 1 | 10/16 | 35.4%–84.8% | 0.345 | 0/16 |
| 200 | 4 | 16/16 | 79.4%–100.0% | 0.690 | 0/16 |
| 200 | 9 | 16/16 | 79.4%–100.0% | 1.035 | 0/16 |

**范围宽度也要看**：ΔJ≤4时，40／200ns两组的总跨度／真Q2中位分别为1.605和0.690；40ns有5/16条出现分离区段，3/16条碰到50pe下界。宽范围可以包含真值，同时仍无法精确分配两个事件的电荷。

**事件数选择会改变分母**：40ns原盲基线只在6/16条输出两个事件，其中ΔJ≤4区段包含真值为4/6；全体已知K2样本的14/16不能直接当作最终盲输出的覆盖率。200ns原盲基线严格整对恢复12/16，条件剖面最优点为10/16；这批小样本未显示逐事件恢复的普遍改善，也未检验差异显著性。

**数值一致性**：32条×3层级的区段与真实Q2探针包含判定有0个不一致；0个探针低于主参考代价超过0.01，0个探针未收敛。插值ΔJ与实际探针的最大差仍为0.508，不能把成员判定一致扩大为整条连续剖面已全局精确。

**方法及区别**：两条上轮ΔJ≤9端点变化超过5%的高斯波形160／190，加原生440／480和高斯200为定向数值诊断。先同网格6起点／250迭代重拟合跨界与最小值邻点，再几何二分，新节点4起点／150迭代；括区≤1%后以8起点／300迭代核验，要求端点≤1%、参考代价改善≤0.01、区段数及边界状态不变，最多12轮。均正反扫描、保留旧更优点及阶段最优固定其余参数切片。有限搜索不保证全局最优，括区1%不代表电荷精度1%。波形模板插值规则、均值和协方差未改。

数值门通过后冻结完整新样本算法：33点基础logQ网格加观测候选，4起点／120迭代；所有单元中点4起点／150迭代，随后上述边界程序。Q=[18000,800]pe、t1=200.3ns、gap40／200ns各16次；另空／弱单／强单各4次。NumPy SeedSequence[202610271801,10000+index,994]独立噪声；固定模型FIR空间相关高斯、前三残差条件为零、逆变换成1024×2ns float64电压。纯α固定95／351ns、充分动态范围、Q域50–50000pe、观测K1局部时域和已知K2保持。冻结旧relative-MLP只给起点，没有学习型Q惩罚或新训练。无新ROOT／光学母事件；不能称原生物理独立确认。

主要包含率以实际输出区段集合包含真Q计数，所有32条含数值失败都留在分母。真实Q2探针只在主结果锁定哈希后评分：固定Q2、8观测起点／300迭代拟合未知Q1/t，不反馈主曲线。插值与探针不一致、探针低于参考、未收敛和域边界单列。表中精确二项95%区间只描述包含比例的抽样不确定性，不把ΔJ=1／4／9命名为68／95／99%电荷区间；每组16条只属先导，不用它反调阈值。盲空／单控制沿用原baseline及原门槛，未拿新剖面替换盲候选。完整指标同时记录Q偏差／sigma68、时误、漏检、额外、误拆／误并与严格恢复；Q非沉积能量，sigma68非FWHM，受控gap非源率。

**运行与复现**：实际源码、模型和协议在[pileup/experiments/alpha_pair_contour_coverage_20261027](../pileup/experiments/alpha_pair_contour_coverage_20261027/README.md)；主数值推断阻断truth／cases／ROOT。首轮仅因NumPy计数JSON序列化失败，旧源码/清单/空输出保留，修复类型并加回归后按同配置重算，在新抽样前完成。使用最多8核、BLAS1。Notebook路径jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb，第108–113节、图0004-108-1至112-1；旧518单元及结果保留。

**待验证**：处理残余数值问题后扩大独立重复次数和Q／gap参数点；另立校准与确认样本再选择覆盖阈值；原生迁移、free-τ、α/γ与高率K≥3分别检验。有限样本／非规则剖面统计的一手证据及Baker DOI仍需共享库另行补充，本项目不修改文献库。

<a id="alpha-pair-selection-20261027"></a>
### L66：2026-10-27——事件数选择、弱电荷偏差与独立选门确认

**依据及审核**：本轮按paper-search只读检索共享库的主题／审核索引、数字脉冲处理／PSD／光产额主题页，复用Liddick、Darby、Grzywacz（2012）Algorithms for pulse shape analysis using silicon detectors完整既有精读笔记，未重读PDF或联网。来源papers/sources/2012_Liddick_NIMA_Algorithms-for-pulse-shape-analysis-using-silicon-detectors.md；DOI10.1016/j.nima.2011.12.034；据笔记PDF3–8、图3–13。对象为约65μm DSSD、40MSPS／12bit、25μs波形及8000余人工叠加双脉冲，采用“识别效率、已识别样本参数精度及假双背景分开报告”的方法思路；不移植100ns、幅比5或ADC精度结论。笔记仍为AI草稿待人工审核，门选记号与自由参数个数疑问保留，没有新的科学采用确认。

本轮分母分解、固定单脉冲校准点、秩门槛和样本数是项目的诊断设计，不是该论文提供的选择后置信区间算法。ΔJ=1／4／9仍不赋予名义置信水平。选择后区间校准和具有未知参数时的误报控制仍需专门的一手统计依据；共享库未修改。

**已完成**：旧32条区间作探索性分母分解，40ns、ΔJ≤4的P(C)=14/16、P(C|S)=4/6、P(C且S)=4/16；不修改旧区间，不称新覆盖校准。主算法waveform_candidates、均值／协方差／logdet、固定双τ与全模拟电压保持，仅比较选门。冻结源码／模型／种子后生成256条校准单脉冲，800／3000／8000／18800pe各64，按观测K1电荷分档。新秩门槛冻结后才生成448条确认：空／四种单／40ns双／200ns双各64，双Q=18000／800pe、t1=200.3ns。

预先规定原／新校准各门槛参数0.025／0.05／0.10，主比较是0.025；没有根据确认真值反调参数或选择获胜门槛。分数D=max(J1−J2,0)，K2任一Q<500置零；存在性、两Q≥500、D>0及秩条件共同决定接受。有限固定Q点和观测分档不保证未知参数上统一控制误报；秩参数不是电荷置信度。校准样本数不同也改变秩离散度，例如n64的0.025仅接受1/65，n80可到2/81；因此比较同时包含噪声来源、有限抽样和秩网格的差异，不能全归因于模型匹配。逐档n和实际切点已保存。

| 间隔/ns | 门槛校准来源 | 输出双事件 | 严格整对恢复 | 时间漏检事件 | 额外输出事件 |
| --- | --- | --- | --- | --- | --- |
| 40 | 原门槛 | 17/64 | 0/64 | 54/128 | 7 |
| 40 | 新条件高斯单脉冲 | 33/64 | 1/64 | 48/128 | 17 |
| 200 | 原门槛 | 63/64 | 46/64 | 8/128 | 7 |
| 200 | 新条件高斯单脉冲 | 63/64 | 46/64 | 8/128 | 7 |

40ns：同为0.025门槛参数，改用新单脉冲校准后，输出两个事件由17/64变为33/64，严格整对恢复由0/64变为1/64。即使直接给这批真实双脉冲选K=2，当前候选中严格成功也只有10/64；这是同一有限候选的诊断上限，不是物理可辨识极限。

40ns被选中候选的有序Q2误差中位由+61.1%变为+39.5%，这是选择群体改变，不是同一条拟合参数变准；包含时间未匹配候选，不称能量分辨率。时间未匹配的额外输出由7增为17个，多检出没有转化成相应的准确恢复。额外输出表示10ns内找不到匹配真事件，不必然表示整窗原本无信号。

200ns：同为0.025门槛参数，改用新单脉冲校准后，输出两个事件由63/64变为63/64，严格整对恢复由46/64变为46/64。即使直接给这批真实双脉冲选K=2，当前候选中严格成功也只有46/64；这是同一有限候选的诊断上限，不是物理可辨识极限。

四档单脉冲共256条，原／新主门槛误报双事件分别为1/1条；含任意时间未匹配输出的帧分别为1/1条。各Q档结果另列，不能把这个等权混合比例当作任意输入谱的误报率。

**解释边界**：包含率、选出双事件和严格恢复三者分开。强制K2诊断只指同一有限候选的成功上限，不是信息论或物理极限。单脉冲假双与时间未匹配输出另报。Q为等效pe非沉积能量，sigma68非FWHM；误差分位含按时匹配版本和全体有序后槽版本，后者保留错误时刻／未选候选。确认误差条为固定校准后计数的精确二项95%区间，不含有限校准不确定性、不等于电荷区间。受控时距非源率。

**复现**：实验[pileup/experiments/alpha_pair_selection_20261027](../pileup/experiments/alpha_pair_selection_20261027/README.md)含实际源／模型哈希、设计和六门槛。新SeedSequence[202610272001,split,index,995]；同一冻结JointNoise高斯、FIR头三残差为零、1024×2ns float64完整电压，充分动态范围。无ROOT、原生光学母事件或训练，纯α固定95／351ns、局部K≤2。推断worker阻断truth／cases／ROOT；评分在预测锁定哈希后进行。6worker校准保存200条后约20分钟无新写入，原因未定；保留原200条及中断记录，2worker仅续算56条，确认也用2worker；计算参数和数据不变。原Notebook543单元及已有结果保留，新第114–119节，图0004-114-1、115-1、116-1／116-2、117-1、118-1。

**下一步／文献待办**：把电荷／时刻正确性作为候选改进目标，在冻结完整选择流程后另做选择后区间校准及独立确认；扩大Q／gap及原生SiPM迁移。选择后区间、未知参数下误报控制的一手统计依据需共享库另行整理。free-τ、α/γ、K≥3和源率扫描未验证。原默认不变，本轮无commit／push。


## L67 — 2026-10-28：条件协方差与弱脉冲均值估计方程

本轮只读复用共享库的主题／审核索引、数字脉冲处理／PSD／光产额主题页和Fixsen等（2004）Pulse estimation in nonlinear detectors with nonstationary noise既有精读笔记，来源papers/sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md，DOI10.1016/j.nima.2003.11.313。据笔记PDF2式(1)–(2)、PDF3图1–3及第4–5节，研究对象为非线性TES及非平稳噪声，0.5–5eV，每个能量点5000次模型噪声抽样；未给出可移植到本项目的采样率或独立真实探测器pile-up恢复验证。笔记仍为AI草稿待人工审核；符号及RMS/Variance术语疑问保留，不把2026-10-08用户提供正式PDF视为科学采用确认，也不沿用其他论文的审核状态。

采用的是条件均值与完整协方差建模的思路。原文插值均值与权重，并不等于本项目的“内层固定C、外层更新C”算法；本轮IRGLS、阻尼、搜索预算和推进标准均为项目试验设计。全高斯似然在当前模型下仍有正确的统计含义，不能把去掉协方差导数或logdet说成理论纠错，也不宣称普遍无偏／最优。未重读PDF、未联网、未修改共享库。

40ns的严格整对恢复，原NLL／扩大搜索NLL／固定一次权重GLS／迭代GLS分别为2/16、2/16、2/16、2/16。所有失败与未收敛样本均保留在分母中。

同一组在固定真实时刻／固定真实电荷／从全部真实参数启动后自由拟合时，严格成功分别为4/16、7/16、6/16；这三项使用额外真值信息，只是误差来源诊断。

200ns的严格整对恢复，原NLL／扩大搜索NLL／固定一次权重GLS／迭代GLS分别为10/16、10/16、11/16、11/16。所有失败与未收敛样本均保留在分母中。

同一组在固定真实时刻／固定真实电荷／从全部真实参数启动后自由拟合时，严格成功分别为13/16、14/16、13/16；这三项使用额外真值信息，只是误差来源诊断。

扩大搜索使4/32条的J降低超过0.05；从真值出发仍有4/32条找到比扩大搜索更低超过0.05的J。有限起点和有限预算不能证明全局最优；J更低也不等于Q和t更接近真值。0.05只用于描述数值差异，不是统计显著性阈值。

共有7/32条出现：真值启动后的自由拟合严格成功，扩大搜索候选失败，但成功解的J反而高出0.05以上。其中40ns为4条、200ns为3条。这些反例说明同一代价的排名与本次逐事件误差并不一致；增加搜索起点无法保证选中更准确的那个解，也不能据此判定物理上一定不可分辨。

没有方法同时通过预先冻结的四项推进条件，因此本轮不增加112条确认样本、不改变默认重建流程。这是44条已暴露模型样本上的开发筛查结论，不是物理不可分辨极限，也没有证明其他优化器无效。

项目位置：`pileup/analysis/alpha_pair_mean_score.py`、`pileup/experiments/alpha_pair_mean_score_20261028/`；本地Notebook：`jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb`第120–125节。采用层级为方法启发，不是原论文算法复现或物理性能迁移。待办：若继续使用均值估计方程，需另找参数相关协方差下的有限样本偏差／效率、收敛及区间校准一手统计依据；新精读需求在共享库单独处理。本轮没有名义区间覆盖率、原生迁移或新独立确认结论。


<a id="research-next-steps-20261028"></a>
## L68 — 2026-10-28：后续研究顺序、可恢复范围与全流程评价

**本轮性质**：回应下一步研究建议；只读本地文献库，据既有精读笔记，未重读PDF、未联网、未启动数值实验。建议先冻结当前盲重建并独立验证参数歧义/区间与质量标记，再扫描间隔、幅比、相位和观测窗口，最后做原生SiPM、α/γ与连续计数率验证。没有把建议当作已完成算法或性能改善。

| 来源及真实文件 | 笔记记录的原文位置 | 本项目采用内容与保留条件 |
| --- | --- | --- |
| M. D. Haselman、J. Pasko、S. Hauck、T. K. Lewellen、R. S. Miyaoka (2012), *FPGA-Based Pulse Pile-Up Correction With Energy and Timing Recovery*, IEEE TNS59(5),1823–1830；[DOI](https://doi.org/10.1109/TNS.2012.2207403)；[精读笔记](../../papers/sources/2012_Haselman_TNS_FPGA-Based-Pulse-Pile-Up-Correction-With-Energy-and-Timing-Recovery.md)；[PDF](../../papers/raw/2012_Haselman_TNS_FPGA-Based-Pulse-Pile-Up-Correction-With-Energy-and-Timing-Recovery.pdf) | PDF3–5图2–7，PDF6–7图9–11/表I | 采样相位、部分积分和幅度/定时耦合；计数、能量和定时分别评价。LFS-3/MAPD-N实测孤立脉冲合成流，65–1000MHz采样且带宽随配置优化，近峰仍拒绝；不移植最短间隔或高率性能。global-review N29-010为AI待审归并，表I分母和流水误差待核。 |
| P. Žugec等／n_TOF Collaboration (2016), *Pulse processing routines for neutron time-of-flight data*, NIMA812,134–144；[DOI](https://doi.org/10.1016/j.nima.2015.12.054)；[精读笔记](../../papers/sources/2016_Zugec_NIMA_Pulse-processing-routines-for-neutron-time-of-flight-data.md)；[PDF](../../papers/raw/2016_Zugec_NIMA_Pulse-processing-routines-for-neutron-time-of-flight-data.pdf) | PDF7–9式13–23、图9–10 | 记录模板拟合及残差质量，保留相位和顺序扣尾的条件。n_TOF多类探测器、典型100MHz–1GHz FADC、离线长记录；固定模板要求形状稳定。D不是校准概率，示例不提供普适逐事件准确率。AI待审；笔记提及“旧555意见”不等于本篇确认，审核原文第555项对象实际为Lin2026。 |
| Shoaib Usman、Amol Patil (2018), *Radiation detector deadtime and pile up: A review of the status of science*, NET50,1006–1016；[DOI](https://doi.org/10.1016/j.net.2018.06.014)；[精读笔记](../../papers/sources/2018_Usman_NET_Radiation-detector-deadtime-and-pile-up-a-review-of-the-status-of-science.md)；[PDF](../../papers/raw/2018_Usman_NET_Radiation-detector-deadtime-and-pile-up-a-review-of-the-status-of-science.pdf) | PDF1–3、5–9，图1–2、6–10 | 区分计数链损失、拒绝、计数补偿与事件/谱形恢复。各模型依赖输入统计和链路，不用固定材料衰减常数充当全系统死时间。global-review B23-036为AI待审归并，公式及活时间异文保留，本轮不采用未闭合公式。 |

**已有结果的约束**：L67的40ns四种已知K估计均2/16，200ns为10/16、10/16、11/16、11/16；是已暴露开发数据。L66的40ns独立模型确认，选双从17/64增至33/64，严格恢复仅0/64至1/64，额外输出7至17。L65的14/16区段包含不能认证95%覆盖。三项不能混为同一评估样本或据此证明物理不可辨识。

**实施建议**：完整候选和事件数选择冻结后，校准/确认分离；真值只用于模拟和评分，质量标记只用观测与拟合输出。区段宽度、分离、域边界、数值稳定性与实际覆盖并列，固定K、选择后条件结果和全部真事件中的有效输出效率分列。ΔJ=1/4/9不直接命名为置信水平。相位、幅比、光量和短窗口扫描逐步展开；观测触发裁窗与真值对齐诊断分开。原生数据按光学母事件独立划分。α/γ同等评价；先0.1/1/10MHz再按结果扩至更高率，未知K≥3、窗前历史、边缘及能谱选择效率纳入。Q为等效pe，未标定前不等同沉积能量。

**证据缺口与待办**：这三篇仅支撑评价与适用条件，没有提供本项目参数相关协方差、边界/多峰、未知K及选择后区间的统计保证。有限样本校准或参数自助法须另有专门一手统计依据及独立验证；共享库的新精读需求在该库单独处理。未知参数下的统一误报/覆盖、free-τ、原生模型失配和高率混合性能仍待验证。提高被保留样本准确率不等于增加正确恢复事件。

**记录与检索范围**：本地Notebook `jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb` 第126节；本轮记录 `jupyter/run_records/research_next_steps_20261028/`。只读knowledge_index的search-topic/search-review在共享盘长时间无输出后被终止（退出130），未重建索引；实际依据为vocabulary/topic-index/review-index导航及rg补查sources、methods、concepts、synthesis、global-review，随后完整回读上列三篇。选择后覆盖相关检索无直接命中不表示全库不存在。本轮未发现这三篇专属的新增科学采用确认，保持原AI待审范围，未修改库或审核状态。


## L69 — 2026-10-28：盲选后的联合参数代价区域与独立歧义验证

本轮承接第126节及文献采用记录L68，复用Haselman等（2012，DOI10.1109/TNS.2012.2207403，笔记定位PDF3–7图2–11/表I）的能量、定时、计数分开评价；Žugec等（2016，DOI10.1016/j.nima.2015.12.054，PDF7–9式13–23/图9–10）的模板拟合质量记录；Usman、Patil（2018，DOI10.1016/j.net.2018.06.014，PDF1–3/5–9）的计数补偿与事件恢复区分。文献的LFS-3/MAPD-N合成流、n_TOF多类探测器和综述模型条件均不移植为GAGG数值保证。沿用既有AI待审状态；本轮没有重新检索或阅读PDF，没有修改共享库。

联合代价集合、经验上分位标定、八个容差外约束搜索及数值预算是项目试验设计，不是上述论文算法的复现。三篇均未提供本项目选择后、有限样本、边界/多解条件下的置信区间保证。相应一手统计证据仍列为共享库独立整理待办。

用校准集中盲选为双事件的36条真实双脉冲，经验95%上分位得到联合代价阈值c=11.864；40/200ns分别贡献6/30条。这是当前选门和设计下的混合样本标定，不是统一95%置信保证。

40ns：原盲流程输出双事件14/64，严格整对恢复0/64；选出双事件后的联合区域包含真实Q1/Q2/t1/t2为13/14，按全部真实双脉冲计算则为13/64。时间漏检57个，额外输出7个。

200ns：原盲流程输出双事件63/64，严格整对恢复38/64；选出双事件后的联合区域包含真实Q1/Q2/t1/t2为58/63，按全部真实双脉冲计算则为58/64。时间漏检15个，额外输出14个。

找到达到或超出容差边界的替代解：79条盲选双事件，其中真实双脉冲77条、控制样本误报双事件2条；严格正确为38/79。分母为0时准确率未定义。

有限搜索未找到替代解：0条盲选双事件，其中真实双脉冲0条、控制样本误报双事件0条；严格正确为0/0。分母为0时准确率未定义。

质量标记没有改变原拟合或增加正确恢复事件；找到替代参数说明该代价区域内存在歧义，有限搜索未找到不证明全局唯一或物理可辨识。

**方法及边界**：参考代价为原K2候选，不声称全局最小值。校准144条，在盲选K2的真实双事件上取ceil(0.95n)排序分数并与0取大，域外真值保留为无穷大，至少20个有效选中对且阈值有限才推进；冻结后生成独立288条。分别报告选双后的四参数联合包含及全体真实双脉冲中的输出且包含。95%为经验目标，没有名义覆盖或未知参数统一保证；不把旧Q2剖面阈值等同当前四维区域。逐固定参数组给精确二项抽样区间，固定分层合并只报计数，区间不含校准及物理模型误差。

**歧义标记**：以Qhat/1.2、Qhat/0.8及时间±10ns定义八个闭合容差外约束，保持时刻排序和原参数域。观测起点2次/100迭代后4次/180迭代复核，保留旧更优点，最小外部代价变化≤0.05为有限数值稳定。实际找到区域内替代点说明存在歧义；没找到不证明全局唯一，收敛/边界不满足时数值未定。原输出不改变；辅助诊断失败仍保留原拟合。没有生成各参数完整投影区间、宽度和分离拓扑。所有失败、漏检、额外输出和质量筛选损失保留；Q是等效pe，sigma68不是能谱FWHM。

**复现与验证**：固定95/351ns、1024×2ns完整电压和充足动态范围；条件高斯FIR首三残差为零。SeedSequence[202610280701,split,index,997]，校准/确认分别split1/2。实际源码、模型、门、协议、运行器快照：pileup/experiments/alpha_pair_blind_uncertainty_v2_20261028/；本地Notebook jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb第127–132节；本地运行记录jupyter/run_records/alpha_pair_blind_uncertainty_v2_20261028/。首轮依赖快照导入失败发生在数据生成前，旧目录完整保留，v2补全absolute analysis依赖。未修改旧生产数据或共享库，未commit/push。

**待办**：各参数完整投影、扩大幅比/相位/窗口、原生独立母事件及模型失配、α/γ、free-τ、连续K≥3和计数率仍需分阶段验证。经验标定和有限搜索不解决已有一手统计依据缺口，也不能证明物理不可分辨。


## L70 — 2026-10-28：原点重启、多起点与Powell的独立优化消融

据共享库既有精读笔记，Haselman等（2012，DOI10.1109/TNS.2012.2207403，PDF3–7图2–11及表I）支持分开评价电荷、定时与保留计数；原条件为LFS-3/MAPD-N孤立实测脉冲构成的合成流，65–1000MHz采样、100kcps–1Mcps，不是本项目GAGG性能保证。Žugec等（2016，DOI10.1016/j.nima.2015.12.054，PDF7–9式13–23、图9–10）提供模板拟合、亚样点相位及质量检查的依据；n_TOF多种探测器、典型100MHz–1GHz采样、形状稳定模板和顺序扣尾，不等同本轮联合高斯优化。

本轮只读核对主题、审核入口、两篇来源笔记及global-review中N29-010所在原始AI审核整理；两篇仍为AI待审，没有新增人工认可，未重读PDF、联网或修改共享库。延长迭代、多起点和Powell组合是项目数值消融设计，不是这两篇的算法复现；不把较低代价解释成真值误差必然更小或全局最优。

| 采用依据 | 本轮核查层级与边界 |
| --- | --- |
| Haselman等，FPGA-Based Pulse Pile-Up Correction With Energy and Timing Recovery，TNS59(5),1823–1830(2012)；[DOI](https://doi.org/10.1109/TNS.2012.2207403)；[既有笔记](../../papers/sources/2012_Haselman_TNS_FPGA-Based-Pulse-Pile-Up-Correction-With-Energy-and-Timing-Recovery.md)；[原PDF](../../papers/raw/2012_Haselman_TNS_FPGA-Based-Pulse-Pile-Up-Correction-With-Energy-and-Timing-Recovery.pdf) | 仅读既有完整来源笔记；PDF3–7图2–11/表I页码沿用笔记。global-review的N29-010所属段明确新增判断仍AI待审；约30%额外保留和能量/定时变化有原合成流条件，不作当前优化器效果依据。 |
| Žugec等，Pulse processing routines for neutron time-of-flight data，NIMA812,134–144(2016)；[DOI](https://doi.org/10.1016/j.nima.2015.12.054)；[既有笔记](../../papers/sources/2016_Zugec_NIMA_Pulse-processing-routines-for-neutron-time-of-flight-data.md)；[原PDF](../../papers/raw/2016_Zugec_NIMA_Pulse-processing-routines-for-neutron-time-of-flight-data.pdf) | 仅读既有完整来源笔记，定位PDF7–9式13–23/图9–10；形状质量和重构曲线不代替逐事件效率。AI待审，旧555意见不扩展为本篇人工认可。 |

**项目设计**：不改变物理响应、候选相关协方差和完整Gaussian quadratic+logdet。原未知K选门固定，仅原K2参数重拟合；K0/K1逐项保留，不在优化后重新按Q≥500筛选。宽域重启同时改变参数化/边界和迭代预算，不能把效果全归因于迭代次数。观测起点、六起点L-BFGS-B及Powell配置在新测试生成前冻结；此前替代解与真值不输入优化器。

**实际结果**：

独立288条确认采用同一波形的配对比较，原事件数门和K判决保持不变；主方案为预先固定的多起点后Powell，另外两种仅作过程消融。

40ns：整对严格恢复0/64→0/64；改善0、退化0，配对精确检验p=1。时间漏检50→51，额外输出4→5。

200ns：整对严格恢复36/64→34/64；改善1、退化3，配对精确检验p=0.625。时间漏检18→20，额外输出15→17。

80条被盲选为K2的确认波形中，20条代价降低超过0.05。控制组原误报双事件1/160，K误报数保持不变是设计约束，不能作为新筛选器有效的证据。

较低代价与逐事件准确度须分别解释；有限搜索未找到更优解不证明物理不可辨识。本轮不替换默认算法，也不外推为原生SiPM、α/γ、free-τ或高计数率性能。

| 间隔/ns | 原流程 | 宽域重启 | 多起点 | 再加Powell（主方案） | 改善/退化 |
| --- | --- | --- | --- | --- | --- |
| 40 | 0/64 | 0/64 | 0/64 | 0/64 | 0/0 |
| 200 | 36/64 | 36/64 | 34/64 | 34/64 | 1/3 |

**评价与限制**：全体64条/间隔的严格整对恢复、电荷误差中位数/sigma68、定时误差、漏检、额外输出及误拆/误合并同时保留。比较同一波形的改善与退化，不拿不同随机批次的恢复率直接评算法。主方案预定为polished；其他阶段探索性，逐组精确区间/配对检验不提供模型误差或普适保证。事件数不变不能证明新误报控制有效；未选中双事件的漏检无法被本轮后处理挽回。

**复现**：纯α固定95/351ns、完整1024×2ns电压、充足动态范围，Q为等效pe而非沉积能量；FIR首三残差零，沿用条件Gaussian模型。15条旧开发样本明确暴露；新288条使用SeedSequence[202610281501,3,index,997]且与旧432条哈希不重复。源码/模型/协议在pileup/experiments/alpha_pair_optimizer_polish_20261028/，本地结果jupyter/artifacts/alpha_pair_optimizer_polish_20261028/、记录jupyter/run_records/alpha_pair_optimizer_polish_20261028/，Notebook为jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb第133–137节。未修改共享库、原实验或默认流程，未commit/push。

**待验证**：本次有限搜索不是全局优化或物理不可辨识证据。后续需分别检查选门漏检和参数估计，再扩展幅比/间隔/采样相位/窗口；原生独立母事件、α/γ、free-τ和连续高率仍需独立验证。若提出全局性或统计覆盖保证，需另行在共享库整理对应一手依据。


<a id="manuscript-readiness-20261028"></a>
## L71 — 2026-10-28：论文准备度、可支撑主张与最小补充证据

**研究问题与核验层级**：用户询问项目距离写文章还缺什么。只读本地知识库，据既有精读笔记核对下列五篇的方法、结果及限制；本轮未重读原PDF、未联网、未启动实验。不是穷尽综述或首创性证明。完整跨项目评估写入本地 `jupyter/Research_Plan_20260921.ipynb` 第14节，保留历史计算/图像；本节只维护采用依据、当前判断及待办。

| 来源、题目与真实文件 | 既有笔记记录的原文位置 | 本项目采用、原条件与审核边界 |
| --- | --- | --- |
| Masaaki Kobayashi、Yoichi Tamagawa、Shougo Tomita、Akihiro Yamamoto、Izumi Ogawa、Yoshiyuki Usuki (2012), *Significantly different pulse shapes for γ- and α-rays in Gd3Al2Ga3O12:Ce3+ scintillating crystals*；[DOI](https://doi.org/10.1016/j.nima.2012.07.055)；[精读笔记](../../papers/sources/2012_Kobayashi_NIMA_Significantly-different-pulse-shapes-for-gamma-and-alpha-rays-in-Gd3Al2Ga3O12-Ce3%2B-scintillating-crystals.md)；[PDF](../../papers/raw/2012_Kobayashi_NIMA_Significantly-different-pulse-shapes-for-gamma-and-alpha-rays-in-Gd3Al2Ga3O12-Ce3%2B-scintillating-crystals.pdf) | PDF2–3，图1–2、表2 | 5 mm GAGG、H2256 PMT、5 GHz采样、512条平均；0.662 MeV γ与约5.47 MeV α。平均形状及1 μs积分条件不能替代本项目SiPM逐事件标定。 补齐材料/读出响应及τ失配验证；AI待审，B23-021整理不等于人工科学确认。 |
| D. J. Fixsen、S. H. Moseley、B. Cabrera、E. Figueroa-Feliciano (2004), *Pulse estimation in nonlinear detectors with nonstationary noise*；[DOI](https://doi.org/10.1016/j.nima.2003.11.313)；[精读笔记](../../papers/sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md)；[PDF](../../papers/raw/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.pdf) | PDF2式1–2及§3，PDF3图1–3 | TES非线性/非平稳噪声模型，每能量5000次模拟；依赖逐能量均值及完整协方差标定，没有GAGG堆积效率或真实器件高率结果。 说明条件均值/协方差标定的已有方法基础；本项目联合NLL、logdet与选门须另验。AI待审，Busch的第658项不移植到本篇。 |
| P. Adrich、I. Zychor (2021), *Least-squares fitting algorithm for peak pile-up correction in gamma-ray spectroscopy*；[DOI](https://doi.org/10.1016/j.nima.2020.164962)；[精读笔记](../../papers/sources/2021_Adrich_NIMA_Least-squares-fitting-algorithm-for-peak-pile-up-correction-in-gamma-ray-spectroscopy.md)；[PDF](../../papers/raw/2021_Adrich_NIMA_Least-squares-fitting-algorithm-for-peak-pile-up-correction-in-gamma-ray-spectroscopy.pdf) | PDF2–5，式1–2、表1、图5–6 | 1英寸LaBr3:Ce/PMT；低率实测单脉冲模板、高率约80 kcps测量及人工叠加双脉冲。候选接受、求和峰校正和有真值的整对成功分母不同；采样设置、噪声/p值及误拆分细节不全。 同任务模板基线、整对参数准确度与谱学结果联合评价；不移植其ns阈值。AI待审，B23-044为待审归并。 |
| Itai Morad、Max Ghelman、Dimitry Ginzburg、Alon Osovizky、Nir Shlezinger (2024), *Model-Based Deep Learning Algorithm for Detection and Classification at High Event Rates*；[DOI](https://doi.org/10.1109/TNS.2024.3371573)；[精读笔记](../../papers/sources/2024_Morad_TNS_Model-Based-Deep-Learning-Algorithm-for-Detection-and-Classification-at-High-Event-Rates.md)；[PDF](../../papers/raw/2024_Morad_TNS_Model-Based-Deep-Learning-Algorithm-for-Detection-and-Classification-at-High-Event-Rates.pdf) | PDF2–3、7–10，图6–13、表II | 25 mm CLYC/PMT、1 GS/s单脉冲库；Poisson到时、线性叠加及合成噪声的高率测试。输出检出/分类，没有逐事件能量谱恢复或持续8 Mcps硬件吞吐验证。 分开检出、条件分类、全部真事件准确率及虚警；第553项仅限定确认此口径与构造流条件，整篇其余内容仍待审。 |
| J. M. Deltoro等 (2025), *Reconstruction of pile-up events using a one-dimensional convolutional autoencoder for the NEDA detector array*；[DOI](https://doi.org/10.1007/s41365-024-01606-y)；[精读笔记](../../papers/sources/2025_Deltoro_NST_Reconstruction-of-pile-up-events-using-a-one-dimensional-convolutional-autoencoder-for-the-NEDA-detector-array.md)；[PDF](../../papers/raw/2025_Deltoro_NST_Reconstruction-of-pile-up-events-using-a-one-dimensional-convolutional-autoencoder-for-the-NEDA-detector-array.pdf) | PDF4–8，图3–10、表1–2 | EJ301/NEDA，250 MHz、14 bit；实测单脉冲人工相加，12–160 ns的异类双脉冲。CC标签无独立TOF真值；底层波形划分、表2分母及表文冲突保留。 分量监督可作对照；波形相关度不能替代Q/t、粒子类型和计数恢复。未验证同类/K≥3或真实连续高率；AI待审，B23-054不构成整体确认。 |

**人工意见核对**：回读[Morad第553项原答复](../../papers/global-review.md#human-review-553-morad2024-confirmed)，仅确认高率联合检出/分类的分母及构造流条件；不扩大为能量恢复或在线实现。global-review的B23-021/044/054归并明确是AI待审；Fixsen未获本篇科学确认，Busch第658项不迁移。本轮没有新增人工科学审核或修改共享库。

**项目证据与判断**：已有完整G4/SiPM链、母事件划分、盲推断隔离、8流程先导、原生双脉冲独立确认和数值负结果，可以开始写引言/方法与部分结果。复算L56实际JSON的15个弱组合点：原观测GLS均值60.50%、原模型GLS58.83%、联合NLL64.50%；沿用L56按母事件成块的配对95%区间，联合相对观测GLS为+4.00个百分点[+1.00,+7.17]。这不是全部计数率收益；其先强后弱8/40 ns仍各0/40。L70独立288条条件高斯确认，优化后40 ns严格整对0/64→0/64，200 ns36/64→34/64，配对p=0.625；不能据较低代价宣称净改善。两批条件、数据域和对照不同，不混为一个总体效率。

**主对照选择补核**：原观测GLS未通过弱单事件电荷偏差控制；原模型GLS与联合NLL均通过该批控制。主比较更宜采用58.83%→64.50%，配对差+5.67个百分点、母事件块95%区间[+2.67,+8.67]；已读取实际 `paired-confirmation.json` 核验，仍限该条件设计。

已有2026-09-30光子层283/163响应及共同157例能标比较（冻结表RMSE45.28 keV），但使用真值分离的光子Q。尚需把**当前SiPM盲恢复的Q**接到独立验证的E_dep能标和最终计数谱。当前95/351 ns、充分动态范围、数据表/代理器件参数是明确模型条件；不自动等于实物GAGG/S13360性能。Kobayashi的不同α响应用于提出失配检查，不据另一晶体/PMT测量静默替换现有物理模型。

**优先补充**：①冻结论文主张/代表方法/主要端点/误报预算与独立划分；②当前SiPM端能标、低率参照及响应失配；③复用旧17间隔作为开发基础，完成最新版的Δt×光量/幅比×先后顺序条件图并在关键区用新母事件确认；④先0.1/1/10 MHz全窗未知K（含K≥3、前史/边缘）验证，再按拟写主张扩率，报告全部事件效率、漏检、误拆/误并、额外输出、能量偏差/分辨率、计数谱和计算成本；⑤共同波形/误报条件下配对比较，统计按母事件/窗口分块，公开过的确认数据不再冒充后续未见测试。

**随主张升级的要求**：若首篇声称α/γ混合，两类逐事件恢复同等重要，须补四类双组合、多重堆积、两种γ谱及比例/低电荷/强前史验证；不以α改善抵消γ退化。若声称实物或在线性能，分别补装置标定/真实高率测量或端到端持续吞吐；纯模拟研究可明确限定条件。fixed/shared双τ可以作为主方案；逐事件free-τ最短窗不必成为所有堆积结论的前提。短窗主张须使用观测触发、记录第二事件是否在窗内，不以1/T或1/Δ换算最高率。保持充分动态范围前提，不新增整窗触顶拒绝要求。

**写作与复现待办**：整理响应/能标、条件恢复、公平比较、完整波形、率依赖和绝对计数能谱六组主文图表；现有编号保留，未生成图不冒充结果。明确数据/模型可用性（ROOT、权重和本地笔记当前不随GitHub分发）；不必把所有历史算法列为主方案。统计区间/选择后的可靠性保证及方法新颖性若升级为主文主张，还需定向查证对应一手统计文献与直接竞争方法；新的共享库整理在papers中另行处理。有限搜索失败不是普适信息极限，现有负结果可限定有效范围。

**本轮实际验证与记录**：相关README/实现/协议/结果JSON与Notebook已读取；对上述两组关键数字复算，检查Notebook格式、新内容链接、原有代码及输出保留。记录位于 `jupyter/run_records/manuscript_readiness_20261028/`。search-topic/search-review返回原文件检索结果，脚本试写缓存被只读库拒绝，未发生库内修改；补用rg并回读原答复。没有新ROOT、训练、重拟合、参数调整、Git提交或上传。


<a id="sipm-energy-transfer-20261028"></a>
## L72 — 2026-10-28：当前 SiPM 电荷能标与独立能量端点评价

**采用依据及边界**：承接L71的论文评价建议，复用其中三篇既有来源记录：Kobayashi等（2012，DOI10.1016/j.nima.2012.07.055，PDF2–3图1–2/表2）的GAGG/PMT响应条件，不能替代本项目SiPM逐事件标定；Fixsen等（2004，DOI10.1016/j.nima.2003.11.313，PDF2式1–2/§3与PDF3图1–3）的条件均值/协方差标定思想，原对象为TES，不提供当前堆积能标保证；Adrich、Zychor（2021，DOI10.1016/j.nima.2020.164962，PDF2–5式1–2/表1/图5–6）的参数恢复、接受比例和谱学评价区分，原对象为LaBr3/PMT。真实来源文件及原条件完整登记于L71；本轮没有新检索、原PDF回查或共享库修改，沿用AI待审状态。

**项目试验设计**：四候选为Q=aE、Q=aL(E)、归一化光产额的Q=ax/(1+bx)、源能组均值节点单调PCHIP；源能分组只用于标定与评分。使用对应冻结材料表的累计α光产额，已含淬灭，不重复Birks修正。目标为E_dep/MeV，不是E_source、入射能量或keVee；Q为初级粒子后代的实际雪崩电荷/q0，排除暗计数祖先。候选形式、均方损失、拒绝惩罚、域边界及通过阈值是本项目设计，不是上述文献的算法复现或数值结论。

560个既有原生单事件母ROOT划分280拟合、280选择；四源能组等权，组内固定平均Q缩放，不用逐事件观测Q倒数加权。选择损失为四组相对能量MSE平均，域外事件罚1；反演限于拟合集实际E_dep范围，域外显式无效而不裁到端点。选定yield_scale后不混入选择集重拟合。主对照为冻结模型GLS与联合NLL，原候选/选门不变，能量推断不读E、真值类型/数量/时刻或ROOT。

评价母事件320个，与560个标定母事件无重叠；来自既有原生确认批次，已公开Q/t结果，本轮新增能量端点，不称作新未见输运确认。共1600帧：320单、80空、1200双，10个组合×8/40/200ns×40对。保持纯α、固定95/351ns、完整1024×2ns模拟电压、充分动态范围及未知K≤2。真值Q≥500pe；主端点为E误差≤10%、时间≤min(10ns,Δt/4)、整对无额外输出，另存5/20%能量容差及原检测指标。

**主要结果：能标未通过预设单事件控制**。真实分量Q的四组有限估计均80/80，但E_dep相对偏差中位数为+13.52%、+5.91%、−1.57%、−5.95%，相对误差sigma68为6.47%、2.24%、0.77%、0.51%；控制要求每组有效率≥95%且绝对偏差≤5%。这些宽度是条件误差分位宽度，不是能谱FWHM。

选择组边界诊断：yield_saturation在四组有限值偏差中位数约−0.11%、−0.10%、−0.25%、−0.12%，但最高源能组12/80电荷超过反演域，其总选择损失0.038214高于yield_scale的0.005401；PCHIP损失0.035492，线性能量比例0.279162。由真实E的观察边界限定平均Q曲线，可能把正常观测涨落拒绝；结果说明当前覆盖/拒绝规则可使偏差更大的候选获选。不能仅删掉域外事件宣称饱和候选已通过，也不能据评价集静默改选。

| 主比较端点：能量±10%严格整对恢复 | 模型GLS | 联合NLL | 配对差 / 百分点；母事件块95%区间 |
| --- | ---: | ---: | --- |
| 含weak的15个组合/间隔均值 | 21.33% | 24.83% | +3.50；[+0.67,+6.50] |
| 其余15个组合/间隔均值 | 47.33% | 49.50% | +2.17；[−1.17,+5.67] |

全部1200对中两方法严格成功412/1200、446/1200；2400个目标事件中漏检152、161，额外输出22、27，域外能量输出47、55，严格能量事件恢复1447、1476。这里的额外输出包括时间无法匹配的预测，不单指完全无信号时的触发。误拆/误合并沿用原评价的操作性标志，不能将这些计数全解释成某种物理失败原因。True Q诊断使用真K/t，仅换算真实分量Q，整对711/1200；它不属于盲性能，也不是严格上限。

配对区间用40个复用母事件块，分两层各20块，2000次bootstrap，种子202610283501；保留跨组合/间隔相关，不含能标选择、参数与物理模型误差。三间隔条件图及每档800个真事件的绝对计数谱已生成；计数谱是四能点受控混合，不是p–¹¹B连续谱或注入率扫描。能量±10%和历史Q±20%端点不同，不能将百分比变化当作算法变化。单事件能标控制失败后，条件能量图仅供诊断，未推进为最终能谱结果或扩大最终条件网格。

**复现与验证**：8项数值/时间匹配测试通过；560/320母事件不重叠、器件配置、原ROOT身份及当前派生文件哈希核验；能量推断的三次truth/labels/ROOT读取探针被阻止；3200条旧盲输出Q/t与检测数一致；4800条评分含1600条真Q诊断。原ROOT字节本轮未全部重哈希，沿用原始记录与导出元数据。首轮因ROOT单精度源能量的严格相等检查中止，v2改用2e−7相对/1e−9MeV绝对容差并保留两种原值，旧目录完整保留。没有新ROOT、训练、波形重拟合或物理参数改变。

实验实际源码/协议/配置及模型：`pileup/experiments/sipm_energy_calibration_v2_20261028/`；本地结果、记录分别为`jupyter/artifacts/sipm_energy_calibration_v2_20261028/`、`jupyter/run_records/sipm_energy_calibration_v2_20261028/`；新Notebook `jupyter/codex_0007_SiPM_Energy_Calibration.ipynb`，固定图号0007-2-1、0007-3-1、0007-4-1/2、0007-5-1。研究计划追加第15节，保留旧计算/输出。

**下一步**：先在目标能区上下边界外补能标保护点，并补窄能带之间的中间点，分开拟合/选择/新母事件确认；通过后扩最新版Δt×光量×顺序图，再验证0.1/1/10MHz完整窗口未知K和绝对计数谱。共享微胞前史与能标误差需分开检查，本轮未做单因素恢复开关消融，不能将全部真Q偏差归因于微胞恢复。α/γ两类同等恢复目标保持；混合、free-τ、模型失配与实物条件仍未由本轮验证。


<a id="sipm-energy-guardband-20261029"></a>
## L73 — 2026-10-29：SiPM能标保护能点、插值与新母事件确认

**依据与审核边界**：承接L71–L72，按共享库topic-index/review-index及响应主题导航，只读复用Rosado2018完整既有精读笔记：[Performance of SiPMs in the nonlinear region](../../papers/sources/2018_Rosado_NIMA_Performance-of-SiPMs-in-the-nonlinear-region.md)，NIMA912,39–42，[DOI](https://doi.org/10.1016/j.nima.2017.10.031)。笔记定位PDF2–4式5–10、图1/3：非线性平均积分电荷标定后再评价能量分布，电荷峰压窄不等于能量分辨改善。原条件为四款S13360与LYSO/LFS/CsI(Tl)、约300–2100keV γ光电峰、完整脉冲积分；采样率、积分门长、温度和样本数不全，强饱和已显示模型失败。保留AI待审状态；本轮未回查原PDF、联网或修改共享库。该原文条件不移植成当前GAGG/6050PE α结论。

**项目设计**：仍比较Q=aE、Q=aL(E)、Q=ax/(1+bx)、单调PCHIP。x为归一化累计α光产额；有理饱和候选为项目经验式，不是Rosado完整模型复现。实际材料表已含淬灭，不重复Birks。E_source、E_dep/MeV与Q/等效pe分别记录；Q只取初级粒子后代的真实雪崩电荷，排除暗计数祖先。

上一轮选择损失被硬反演边界影响。新增源1.5/9MeV保护点及2/3.5/4.5/6.5MeV中间点；先导16例沉积范围确认：下端0.134–0.156MeV，上端8.60–8.65MeV。目标按内部源能组预定，保留组内所有事件，不按真实E或观测Q删边缘。反演仍限于拟合集实际E范围、不外推/裁端点；保护组不进入主选模损失，避免把新标定极值再次作为主目标。

新ROOT共496：先导16、拟合96、选择64、冻结后确认320。加旧560个母事件后，拟合376、选择344；按源能组等权、组内固定平均Q缩放。选择损失为八个内部组相对E误差MSE等权平均，域外罚1；选定后不混入选择集重拟合。确认十档源能量各32，其中1.8/7MeV未参加拟合或选模。能标与实际分析源码冻结后才开始确认生产。

**新单事件结果**：选择`empirical_pchip`，十组预设控制通过；有限逆变换320/320，能量±10%内318/320，最大组偏差中位数绝对值3.48%。控制为每组有效率点估计≥95%、绝对偏差中位数≤5%；32例/组不提供总体95%保证，全有效32/32的精确双侧95%二项下限约89.1%。这些是已知原生分量Q到E_dep的独立能标确认，不是未知事件数的盲检出率。

| 源能量 / MeV | 有限估计 | E_dep偏差中位数 | 相对误差sigma68 | 能量±10%内 |
| --- | ---: | ---: | ---: | ---: |
| 1.625 | 32/32 | +3.48% | 3.17% | 31/32 |
| 1.8 | 32/32 | -0.45% | 4.82% | 31/32 |
| 2 | 32/32 | -0.68% | 3.51% | 32/32 |
| 2.5 | 32/32 | -0.09% | 2.40% | 32/32 |
| 3.5 | 32/32 | -0.71% | 1.04% | 32/32 |
| 4.5 | 32/32 | +0.06% | 0.95% | 32/32 |
| 5.486 | 32/32 | +0.13% | 0.85% | 32/32 |
| 6.5 | 32/32 | -0.09% | 0.76% | 32/32 |
| 7 | 32/32 | -0.21% | 0.49% | 32/32 |
| 8 | 32/32 | +0.19% | 0.44% | 32/32 |

误差sigma68=(p84−p16)/2，不是能谱FWHM；组内中位数bootstrap区间另存。旧固定能标也在同一320例上评价，保留其失效，不从新确认重新选择候选。新增保护和中间样本、内部目标设计共同变化，没有单因素消融证明各自贡献。PCHIP与有理饱和候选的选择损失0.00065313/0.00065852很接近，不据一次选择认定普遍优劣。

最低源能组1.625MeV仍有+3.48%中位偏差，条件bootstrap95%区间[+1.45%, +3.80%]。通过预设±5%判据不表示无偏；保留这一残余，不在确认集上改模型。

依冻结曲线Q=f(E)作局部小误差传播，ΔE/E≈[Q/(E f′(E))]ΔQ/Q；这与Rosado式9的局部斜率解释一致，未采用其具体器件数值或FWHM作本项目结果。弱源1.625MeV的标定参考点相对换算因子约1.437，能量±10%局部对应Q约±6.96%；源8MeV参考点因子约0.864，对应Q约±11.57%。该近似不是大误差、偏斜分布或堆积模型失配的严格保证；实际成功率仍逐事件非线性反演评分。

**历史盲拟合重放**：原1600帧已公开过Q/t及能量结果；当前仅更新能量换算，原3200条两方法Q/t/事件数、漏检、额外输出与误拆/误并标志一致。能标母事件与这批320个历史光学母事件仍分离，但不能称作新堆积确认。真值Q诊断使用已知K/t，不属于盲恢复或严格上限。

含弱α的15个组合/间隔等权均值，能量±10%的严格整对恢复由旧能标下模型GLS 21.33%、联合NLL 24.83%，变为新能标下31.67%、38.17%。新能标下两方法配对差+6.50个百分点，母事件块95%区间[+3.00,+10.17]。重采样沿用40母事件块分两层各20、2000次配对；区间不含能标参数/选择和物理模型不确定度。它说明能量换算影响评价，波形检出数量未改变；不与历史Q±20%端点混比。

按8/40/200ns输出条件图和每档800个真事件的绝对能量计数谱，不将恢复谱归一为真值总数。计数谱是四能点受控混合，不是p–¹¹B连续谱或注入率扫描；域外能量、漏检和额外输出保留。

含弱组合的True Q诊断严格整对均值为64.17%；某些组合低于盲拟合，进一步说明它不是严格性能上限。当前数据提示需要检验静态能标与分量电荷分配误差是否抵消，但尚不能仅据该图确定成因。能量端点改善不证明实际电荷分解也改善。
**复现与完成验证**：496个ROOT均为一个20μs原生单α窗，中央入射，固定95/351ns双τ，2ns采样，最多4任务且启动间隔至少1s。实际宏与RunInfo、配置、三类种子、ROOT/普通临时副本哈希核对；无临时链接、无旧数据覆盖或物理参数改动。源种子690000000/690010000/690020000/690030000分别起始先导/拟合/选择/确认，SiPM及读出各加100000/200000，检查已有运行清单冲突。母事件按ROOT、SHA256及种子划分；能量推断三次truth/labels/ROOT探针拒绝，预测冻结后才连接新确认E_dep。

实际共创建499份ROOT：496份通过完整流程纳入研究，另3份因处理错误或中断保留、不纳入统计。按原种子重试没有增加独立母事件；原尝试的ROOT结构/RunInfo另行核对，底层处理错误原因未定位。原文件、实际宏和日志保留，实验retry-audit.json及Run_Catalog登记重试关系，不根据科学结果选择重试。

实际源码、配置、宏、模型与冻结记录在`pileup/experiments/sipm_energy_guardband_20261029/`；本地结果`jupyter/artifacts/sipm_energy_guardband_20261029/`，记录`jupyter/run_records/sipm_energy_guardband_20261029/`。codex_0007追加第8–14节，图号0007-8-1、9-1、10-1、11-1、12-1；原18个单元及输出保留。Run_Catalog登记所有496个实际ROOT/宏/种子，Research_Plan追加第16节。未更换默认拟合、未提交或推送。

**后续**：新单事件控制通过时，先检查共享微胞前史下静态能标的条件偏差与实际盲Q误差，扩最新版间隔×光量×顺序图，再推进0.1/1/10MHz完整窗口未知K和绝对计数谱；未通过则保留失败并另开冻结诊断。中央位置、固定器件模型的能标仍需位置/参数失配、稀疏能带间更密扫描及实物校准；α/γ同等恢复、free-τ和K≥3未由本轮验证。

上述前史诊断另只读参考[Rosado2019既有精读笔记](../../papers/sources/2019_Rosado_IEEESensorsJ_Modeling-the-Nonlinear-Response-of-Silicon-Photomultipliers.md)，IEEE Sensors Journal19(24),12031–12039，[DOI](https://doi.org/10.1109/JSEN.2019.2938018)，PDF3–4式14/19与PDF5–8验证。文献平均积分电荷依赖光子时间分布、像素恢复与相关噪声，验证为S13360-1325CS/1350CS、LYSO约300–2100keV γ的400ns积分，以及连续绿色LED；未完整报告采样率/样本数。其近似限于适中非线性和相关噪声，1350连续光平台约20%差异仍保留。将其用于提出“孤立能标可能不适用于不同脉冲前史”的检查是项目推论，不是GAGG逐事件性能证据。式19及恢复积分尚未在本轮实现；原文上升/下降双指数差不是当前95/351ns两发光衰减混合。保留AI待审状态，未回查原PDF或修改共享库；后续若实现原式，先定点核对公式及其假设。
