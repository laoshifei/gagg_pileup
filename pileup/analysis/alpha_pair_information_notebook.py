"""Append executed local-information diagnostics to Notebook 0004."""
from __future__ import annotations

import base64
import hashlib
import json
import os
import re
import shutil
import tempfile
from pathlib import Path

import nbformat

from .alpha_pair_information_study import PROJECT, OUTPUT, RECORD, load, save, digest


def deliver():
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output

    assert load(RECORD/'science-validation.json')['passed']
    assert load(RECORD/'figure-delivery.json')['passed']
    report=load(OUTPUT/'report.json')
    path=PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    previous=path.read_bytes()
    notebook=nbformat.reads(previous.decode(),as_version=4)
    before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4)
    assert len(before.cells)==387
    assert notebook.cells[:len(before.cells)]==before.cells
    if len(notebook.cells)>len(before.cells):
        assert digest(path)==load(RECORD/'notebook-delivery.json')['notebook_sha256']
        assert all(c.get('id','').startswith('information-20261027-') for c in notebook.cells[len(before.cells):])
        notebook.cells=notebook.cells[:len(before.cells)]
    shutil.copyfile(path,RECORD/'notebook.before-delivery.ipynb')
    original=json.loads(json.dumps(notebook.cells))
    cells=[]

    def md(text):cells.append(nbformat.v4.new_markdown_cell(text))
    def code(text):cells.append(nbformat.v4.new_code_cell(text))
    def figure(number,title):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{title}**')
        code(f'display(Image(filename=str(information_output / "figures/{number}.png")))')

    md('''<a id="information-stages-20261027"></a>
## 77. 从输出两个脉冲，到分别估准两个事件

第76节提出的下一步在这里实际执行。将已有盲结果依次拆为：输出K=2；两个时刻均在10 ns内；两个时刻均在严格容差min(10 ns,Δt/4)内；再要求两个Q均在真实值±20%内且没有额外输出。各项均以全部双事件为分母，输出K=2本身不等于真实弱事件已经检出。

复用上一轮760帧原生数据、1520帧匹配高斯评价数据及其冻结拟合；其中双事件分别360、720帧，三种组合×8/40/200 ns，原生每点40对、高斯每点80对。这里**没有新拟合波形、训练网络或生成ROOT**，没有改变增峰门槛。全部未收敛结果保留；空／单控制和误拆／误并／漏检／额外输出仍沿用上一轮原值。已公开母样本用于诊断，不能当作新独立测试。

本轮新增444个已知K=2、已知真实Q/t的数学模型信息计算：360个实际参数点，80个预定名义网格点，4个采样相位检查点。原生Q是实现雪崩等效pe，高斯Q是数值相同的条件参数，不是沉积能量或新实现的雪崩电荷；t是有效SiPM面首光子时刻。固定95/351 ns发光双τ、2 ns×1024点、完整0–2048 ns模拟电压和充分动态范围；FIR后实际使用1021维边缘分布。

文献依据：Liddick等2012（DOI：10.1016/j.nima.2011.12.034，既有笔记定位PDF3–8、图3–13）分别讨论双脉冲识别与参数精度。其对象为约65 µm DSSD、40 MSPS／12 bit、25 µs记录，超过8000个由实测单脉冲复制、缩放、平移构造的双脉冲，不是GAGG共同微胞原生样本。本文的定时散布也不是最小可分间隔。本项目只采用分开评价的原则，没有移用其数值分辨率或ADC条件。既有AI精读稿待人工审核；global-review的B23-022明确没有新增人工科学确认。本轮复用精读笔记，没有重查原PDF。
''')
    code('''from pathlib import Path
import json
import numpy as np
import pandas as pd
from IPython.display import display, Image
information_root = Path.cwd().resolve()
while not (information_root / "pileup/experiments/alpha_pair_information_20261027").is_dir():
    if information_root.parent == information_root:
        raise RuntimeError("请从项目根目录或jupyter目录运行")
    information_root = information_root.parent
information_output = information_root / "jupyter/artifacts/alpha_pair_information_20261027"
information_record = information_root / "jupyter/run_records/alpha_pair_information_20261027"
information_report = json.loads((information_output / "report.json").read_text())
information_stages = pd.DataFrame(information_report["stages"])
print("444个参数点；2280帧缓存盲评分；其中1080帧双事件；2160条真值辅助拟合误差。")
print("新盲拟合：", information_report["new_blind_fits"], "；新ROOT：", information_report["new_roots"])
display(information_stages[["domain", "condition", "gap_ns", "frames", "K2_selected", "both_times10", "both_times_strict", "strict_success"]])
''')
    figure('0004-77-1','四层结果，各层分母均为全部双事件；(a–c)原生中／中、先强后弱、先弱后强，(d–f)匹配高斯对应条件。')
    code('''display(information_stages[["domain", "condition", "gap_ns", "frames", "missed", "extra", "split", "merge"]])
information_cached = pd.DataFrame(information_report["primary_cached_summary"])
display(information_cached[information_cached.kind != "pair"][["method", "condition", "eligible_frames", "success", "extra_frames", "q_error_p16_p50_p84"]])
print("额外输出可包括未与真实时刻匹配的输出；不是全部发生在无信号位置。误拆、误并为既有匹配图的操作性标志。")
''')

    md(r'''<a id="information-fisher-20261027"></a>
## 78. 用局部信息分开时间耦合与电荷分配歧义

在冻结条件高斯模型内，令$z=(t_1/(10\,\mathrm{ns}),t_2/(10\,\mathrm{ns}),\ln(Q_1/\mathrm{pe}),\ln(Q_2/\mathrm{pe}))$，$\mu=\mu(z)$，$C=C(z)$。下式由本项目对高斯对数似然求导并取期望得到：

$$F_{ij}=\mu_i^\mathsf{T}C^{-1}\mu_j+\frac12\operatorname{tr}(C^{-1}C_iC^{-1}C_j).$$

第二项来自参数相关协方差，不能仅用均值模板的梯度代替。Fixsen等2004（DOI：10.1016/j.nima.2003.11.313，既有笔记PDF2式1–2、PDF3；TES模拟，每能量5000次；AI稿待人工审核）支持条件均值／协方差处理，但其原文是加权最小二乘；这里的完整高斯信息公式不是该文直接给出的GAGG结论。

定义$V=F^{-1}$，以$\sqrt{V_{ii}}$描述局部参数尺度。这是固定正确模型、已知K=2、局部可微和适当无偏性等条件下的信息参照。它不是当前有限样本拟合的实际标准差，不自动覆盖有偏、有边界或错误极小值的估计。特别是较大的$\sigma_{\ln Q}$不能直接写成同数值的相对电荷标准差；±20%的电荷容差对应非对称对数区间$[\ln0.8,\ln1.2]$。

对比三种条件：四个参数均未知；两个真实时刻固定后只估Q；强事件真实Q/t固定后只估弱事件Q/t。前者弱参数取$F^{-1}$的子块，后两者先取剩余参数的信息子块再求逆，不能混淆。固定强Q/t**仍保留强事件的随机噪声**，不是精确减去真实强波形。图中的相关系数来自这个局部逆矩阵，不是跨样本真实Q的相关系数。

在预定“强Q/t已知”对照之后，另作探索性的矩阵分拆：只固定强Q或只固定强t，以定位哪类干扰参数起主要作用。这不改变拟合算法，也没有将真值交给盲推断。

数值计算采用中心差分，先比较缩放坐标步长0.001和0.0005，矩阵与局部尺度相对变化均≤2%才通过，否则按预定规则再减半。主步长的时刻扰动为0.005 ns。响应含经验插值和分段限制，差分收敛只检验当前点的局部稳定性；不代表全局光滑。非正定或条件数超过10¹²的点记录未解决，不用伪逆或正则化掩盖。
''')
    figure('0004-78-1','实际360个参数点的局部信息结果；(a–c)先强后弱，(d–f)先弱后强。各列为弱lnQ尺度、弱t尺度、参数相关；线为中位，阴影为参数点间16–84%分位，不是置信带。')
    code('''information_local = pd.DataFrame(information_report["information"])
display(information_local[["condition", "gap_ns", "frames", "weak_sigma_time_ns", "weak_sigma_log_charge", "fixed_times_sigma_log_charge", "known_strong_sigma_log_charge", "correlation_weakQ_strongQ"]])
display(information_local[["condition", "gap_ns", "known_strong_Q_sigma_log_charge", "known_strong_t_sigma_log_charge", "correlation_weakQ_weakt"]])
print("各列表顺序为16%、50%、84%分位；只固定强Q／只固定强t是探索性派生对照。")
''')
    strong40=next(r for r in report['information'] if r['condition']=='high_weak' and r['gap_ns']==40.)
    weak40=next(r for r in report['information'] if r['condition']=='weak_high' and r['gap_ns']==40.)
    md(f'''40 ns先强后弱：弱lnQ局部尺度中位为 **{strong40['weak_sigma_log_charge'][1]:.4f}**；固定两个真实时刻后为 **{strong40['fixed_times_sigma_log_charge'][1]:.4f}**，几乎不变；固定强事件Q/t后为 **{strong40['known_strong_sigma_log_charge'][1]:.4f}**。只固定强Q得到{strong40['known_strong_Q_sigma_log_charge'][1]:.4f}，只固定强t仍为{strong40['known_strong_t_sigma_log_charge'][1]:.4f}。弱Q与强Q的局部相关系数中位{strong40['correlation_weakQ_strongQ'][1]:.4f}。

这支持当前模型内“两个事件之间的电荷分配歧义”是主要限制之一，单独改善起始时刻未必解决。它不是把强Q已知的结果变成可实现算法。反向的先弱后强在40 ns时，弱lnQ尺度中位仅{weak40['weak_sigma_log_charge'][1]:.4f}，体现顺序相关性；不是只有间隔一个变量控制恢复。
''')

    md('''<a id="information-grid-20261027"></a>
## 79. 同时改变间隔、电荷比和采样相位

预定名义网格固定强Q=18000 pe、弱Q=500/700/1500/3000 pe，两种先后顺序，Δt=2/4/8/12/20/40/80/120/200/400 ns，首事件160.37 ns，共80点。此处Q是模型参数，不是新模拟的实现电荷；曲线描述局部已知K信息，不是新测试的恢复率，也没有将不同Q的结果合成一个最短间隔。

图中ln(1.2)水平线只表示20%上侧误差的对数尺度。即使局部尺度低于这条线，也不能宣称单次事件必定在±20%内，更不能把交点叫作最小可分时间。固定的一种相位也不足以代表任意采样位置。

另对强→弱、40 ns、18000/700 pe，将首时刻设为160 ns加0.13/0.37/0.73/1.13 ns，保留四个相位点；其中0.37与网格重复，用于一致性检查。只在这一条件检查相位，不据此排除所有插值或相位效应。
''')
    figure('0004-79-1','名义网格的弱lnQ局部尺度；(a)先强后弱，(b)先弱后强。强Q固定18000 pe，颜色为弱Q，包含2–400 ns全部预定点；不是盲恢复曲线。')
    code('''display(pd.DataFrame(information_report["phase_check"])[["phase_ns", "weak_sigma_time_ns", "weak_sigma_log_charge", "fixed_times_sigma_log_charge", "known_strong_sigma_log_charge", "derivative_step", "derivative_converged"]])
information_all = pd.DataFrame(json.loads((information_output / "local-information-summary.json").read_text())["rows"])
display(information_all[information_all.role == "nominal_grid"][["condition", "gap_ns", "weak_charge_pe", "weak_sigma_time_ns", "weak_sigma_log_charge", "derivative_converged"]])
''')

    md(r'''<a id="information-errors-20261027"></a>
## 80. 实际拟合误差是否服从这个局部尺度？

复用原生240对与匹配高斯480对的三类真值辅助拟合，共2160条：已有K=2候选、真实Q/t初始化联合拟合、固定真实t仅拟合Q。每条都进入统计，包括低Q、失败和优化器未收敛输出，没有只选成功样本。它们的约束、搜索预算不同，不能单凭这些对照认定全局极小值或最佳算法。

对于前两类，以真实参数处的弱$(t/10\,\mathrm{ns},\ln Q)$子协方差，将误差归一化，计算$d^2=e^\mathsf{T}V_w^{-1}e$。若局部线性、无偏高斯估计近似成立，参照分布是$\chi^2_2$。对于固定真实时刻，用仅Q的信息逆矩阵归一化弱lnQ误差，参照$|N(0,1)|$。图中竖线分别为5.991和1.960。

真实初始化给了真实Q/t起点，并将时刻搜索限制在真值±48 ns内；有限优化可能停留在起点附近。其误差特别集中甚至接近零，不代表一个合法盲估计器突破了信息下界，也不能据此校准真实算法的不确定性。

**这些只是局部高斯参照区域，不是经过校准的95%置信区间。**计算使用真实Q/t和真实K，未建立盲算法中的覆盖保证；原生模型失配、非线性、边界、偏差和有限搜索都可能使误差偏离参照。高斯两次抽样也不是两个新光学母事件。表中区域内比例只作描述，不将2160条或跨条件复用的样本当作独立观测建立覆盖率置信区间。
''')
    figure('0004-80-1','弱Q/t归一化误差累计分布；(a–c)先强后弱的8/40/200 ns，(d–f)先弱后强。颜色区分原生／高斯，线型区分已有K=2／真实初始化，黑虚线为局部高斯参照。全部误差保留。')
    figure('0004-80-2','固定两个真实时刻后的弱lnQ归一化误差；(a–c)先强后弱的8/40/200 ns，(d–f)先弱后强。竖线1.960不是已校准的实际区间保证。')
    code('''information_errors = pd.DataFrame(information_report["errors"])
display(information_errors[["domain", "condition", "gap_ns", "method", "frames", "strict_success", "inside_reference95", "not_converged", "weak_charge_fractional_error", "weak_time_error_ns"]])
print("Q误差分位包含全部辅助拟合、按先后顺序对应真实分量，不是只在成功匹配后统计。")
print("inside_reference95表示落入真参数处局部高斯参照区域，不能直接称作盲算法的95%覆盖率。")
''')
    md('''40 ns先强后弱的匹配高斯样本中，已有K=2候选仅22/80落入二维局部95%参照区域，真实初始化为74/80；后者严格整对成功仍只有24/80。固定真实时刻后73/80落入一维参照区域，但严格整对成功只有29/80。这两套比例并不矛盾：参照区域可以比±20%的电荷要求宽很多，区域内并不等于准确恢复。

同条件固定真实时刻后的弱Q相对误差16/50/84%分位，高斯为−46.0%／+10.4%／+59.9%，原生为−24.3%／+21.5%／+43.1%。原生偏差与高斯宽度不同，不能仅用一条局部尺度把实际误差完全解释完。已有K=2候选还存在很大的时刻离群，表和图均没有截掉这些失败。
''')

    md('''<a id="information-conclusions-20261027"></a>
## 81. 本轮结论与下一步

已完成从输出数、时间匹配到电荷精度的分层，以及444点局部信息计算和2160条实际辅助拟合误差对照。40 ns先强后弱盲结果，原生输出K=2有7/40，两个时刻正确2/40，严格整对成功0/40；高斯对应9/80、4/80、0/80。因此漏检、定位和电荷精度都有损失。当前模型内，即使已知K且处于真参数附近，强弱电荷的反相关仍是主要困难之一；固定真实时刻几乎不能缩小弱Q局部尺度，固定强Q则可以。这个结果说明不能只围绕寻峰时间继续加搜索预算。

实际误差与局部高斯参照的偏离保留在第80节，不把矩阵求逆等同于逐事件可信误差条。当前模型仍有原生非高斯残差，已知K的有限搜索也存在局部解和偏差，因此本轮没有证明所有方法均无法恢复，更没有给出固定的“多少ns必定不能分”。

下一步应针对电荷分配做可盲实现的对照：将联合估计写为总Q与分配比例，检查较长的似然谷；若使用电荷先验或学习型估计，先验只能由训练母样本确定，必须与不加先验方案在同一误拆预算下比较，并检验连续Q、极端幅比和模型失配。换坐标本身不增加信息，利用离散四档真实能量作为先验则会虚增成绩，不能这样做。先在现有数据上开发，冻结方案后再生成独立母样本测试。

本轮没有改动默认算法、训练网络或物理模型；不提供free-τ最短窗、混合α/γ、随机高率K≥3或实物性能结论。源码及条件冻结，旧387个单元及输出保留；图77–80的编号固定。

完整结果保存在本地[分层统计](artifacts/alpha_pair_information_20261027/stage-summary.csv)、[局部信息逐点值](artifacts/alpha_pair_information_20261027/local-information.csv)、[实际误差汇总](artifacts/alpha_pair_information_20261027/error-summary.csv)及[报告](artifacts/alpha_pair_information_20261027/report.json)。
''')
    code('''information_checks = json.loads((information_record / "science-validation.json").read_text())
print("通过：", information_checks["passed"], "；参数点：", information_checks["tasks"], "；实际误差：", information_checks["normalized_diagnostic_scores"])
print("未解决参数点：", information_checks["unresolved"], "；差分未收敛：", information_checks["derivative_unconverged"])
print("最大白化相对误差：", information_checks["maximum_whitening_error"])
print("差分减半后：矩阵最大相对变化", information_checks["maximum_derivative_matrix_change"], "；局部尺度最大变化", information_checks["maximum_derivative_sigma_change"])
print("评分与缓存一致：", information_checks["strict_scores_match_cache"], "；输入与冻结源码保持：", information_checks["inputs_unchanged"] and information_checks["freeze_unchanged"])
''')
    for i,cell in enumerate(cells):cell['id']=f'information-20261027-{i:02d}'
    InteractiveShell.instance()
    namespace={};executed=0;cwd=os.getcwd()
    try:
        os.chdir(PROJECT/'jupyter')
        for cell in cells:
            if cell.cell_type!='code':continue
            with capture_output(stdout=True,stderr=True,display=True) as capture:
                exec(compile(cell.source,'<information-notebook>','exec'),namespace)
            assert not capture.stderr,capture.stderr
            cell.outputs=[]
            if capture.stdout:cell.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=capture.stdout))
            for out in capture.outputs:cell.outputs.append(nbformat.v4.new_output('display_data',data=out.data,metadata=out.metadata))
            executed+=1;cell.execution_count=executed
    finally:os.chdir(cwd)
    notebook.cells.extend(cells)
    assert notebook.cells[:len(original)]==original
    anchors=re.findall(r'id="codex-figure-([^\"]+)"','\n'.join(c.source for c in notebook.cells if c.cell_type=='markdown'))
    assert len(anchors)==len(set(anchors))
    for cell in cells:
        if cell.cell_type=='markdown':
            assert not re.search(r'\]\([^)]*\.(?:pdf|ipynb|md)(?:[#?)][^)]*)?\)',cell.source,re.I)
    images=[out['data']['image/png'] for c in cells if c.cell_type=='code' for out in c.outputs if 'image/png' in out.get('data',{})]
    expected=load(RECORD/'figure-delivery.json')['figures']
    assert len(images)==len(expected)==5
    assert {hashlib.sha256(base64.b64decode(i)).hexdigest() for i in images}==set(expected.values())
    nbformat.validate(notebook)
    encoded=(nbformat.writes(notebook)+'\n').encode()
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=path.stem+'.',suffix='.partial',delete=False) as stream:
        stream.write(encoded);stream.flush();os.fsync(stream.fileno());temporary=Path(stream.name)
    assert temporary.read_bytes()==encoded
    nbformat.validate(nbformat.reads(encoded.decode(),as_version=4))
    assert path.read_bytes()==previous,'Notebook changed during delivery; retain user edits'
    temporary.replace(path)
    assert path.read_bytes()==encoded
    save(RECORD/'notebook-delivery.json',dict(passed=True,previous_cells=len(original),new_cells=len(cells),
        executed_cells=executed,history_preserved=True,figures=expected,notebook_sha256=digest(path),source_sha256=digest(__file__)))
    print('INFORMATION_NOTEBOOK_DELIVERED',len(cells),executed,flush=True)


if __name__=='__main__':deliver()
