"""Executed local Notebook delivery and preservation checks for this study."""
from __future__ import annotations

import ast
import base64
import hashlib
import os
import re
import subprocess
import sys
import tempfile
import nbformat

from .alpha_pair_contour_study import PROJECT,EXPERIMENT,OUTPUT,RECORD,INPUT,OLD_RECORD,PREVIOUS,NOTEBOOK,NAME,check_freeze
from .alpha_pair_information_study import load,save,digest,rows

MARKER='alpha-pair-contour-coverage-20261027'
LITERATURE='''本轮只读复用共享库的主题／审核索引及已有完整精读笔记，没有重读PDF或联网。Baker、Cousins（1984），Clarification of the use of chi-square and likelihood functions in fits to histograms，NIM221,437–442：已有笔记定位PDF1、PDF4§5–7、PDF5§8，区分点估计、区间覆盖与拟合优度，有限样本按具体模型研究抽样分布。其研究对象是Poisson／multinomial直方图，不是相关电压；本文没有提供本项目数值阈值或边界算法。来源papers/sources/1984_Baker_NIM_Clarification-of-the-use-of-chi-square-and-likelihood-functions-in-fits-to-histograms.md，AI草稿待人工审核、无科学采用确认，DOI未核实，式3概率记号疑问保留。

Liddick、Darby、Grzywacz（2012），Algorithms for pulse shape analysis using silicon detectors，NIMA669,70–78，DOI10.1016/j.nima.2011.12.034：已有笔记PDF3–8／图3–13，约65μm DSSD、40MSPS／12bit、25μs记录、8000余人工叠加样本。仅沿用检测、参数精度和假双背景分开评价；不移植100ns、幅比5或ADC限制。来源papers/sources/2012_Liddick_NIMA_Algorithms-for-pulse-shape-analysis-using-silicon-detectors.md，仍为AI草稿待审，门选记号和自由参数数目缺口保留。数值括区1%／代价0.01、两参数点及抽样数量是项目设计，不是文献性能保证。'''


def summary(report):
    text=f"5条定向数值诊断中，{report['refinement_stable']}/5通过新增的边界括区与更强局部重拟合检查。"
    if 'coverage' not in report:
        return text+'独立抽样尚未开展；不能把未完成的统计检验写成已校准区间。'
    text+=f"随后冻结方案，生成32条新条件高斯双脉冲和12条空／单控制；{report['stable_frames']}/32条双脉冲通过数值筛查，全部32条都保留在主要包含率分母中。"
    for band in ('40','200'):
        c=next(c for c in report['coverage'] if c['band']==band and c['level']==4.)
        m=next(m for m in report['conditional_metrics'] if m['band']==band and m['kind']=='pair' and m['method']=='conditional_profile_K2')
        text+=f"{band}ns组ΔJ≤4范围包含真实Q2为{c['included']}/16，已知K2严格整对恢复为{m['success']}/16。"
    text+='这检验的是冻结模型内、固定参数点的重复抽样表现；每组16条不足以认证95%覆盖，更不能推广为原生SiPM的盲恢复或普适最短间隔。'
    return text


def coverage_table(report):
    if 'coverage' not in report:return '本轮尚无独立覆盖率抽样。'
    lines=['| 间隔/ns | ΔJ上限 | 范围包含真Q2 | 该比例的精确二项95%区间 | 跨度/真Q2中位 | 分离区段帧数 |',
           '| --- | --- | --- | --- | --- | --- |']
    for c in report['coverage']:
        lo,hi=c['binomial95']
        lines.append(f"| {c['band']} | {c['level']:g} | {c['included']}/16 | {100*lo:.1f}%–{100*hi:.1f}% | {c['width_over_trueQ'][1]:.3f} | {c['disconnected']}/16 |")
    return '\n'.join(lines)


def boundary_explanation(report):
    a=next(r for r in report['refinement'] if r['domain']=='continuous' and r['index']==160 and r['level']==9.)
    b=next(r for r in report['refinement'] if r['domain']=='continuous' and r['index']==190 and r['level']==9.)
    return f'''两条的变化来源并不一样：第160条仅同网格加强局部搜索，ΔJ≤9端点就移动{100*a['same_grid_endpoint_fraction']:.2f}%，原节点最大代价下降{a['same_grid_node_gain']:.3f}；后续加点与复核再改变约{100*a['grid_endpoint_fraction']:.3f}%。第190条同网格端点变化可忽略，后续加点与复核改变约{100*b['grid_endpoint_fraction']:.2f}%。因此前者明确存在原节点的局部搜索不足；后者主要在后续细化中变化，但该步骤也重新优化了未知参数，不能全归因于画图插值。

5条最终括区均小于1%，更强的末次局部搜索通过筛查；最小代价的总改善仍很小，5条的严格整对恢复结果均未改变。这是定向数值诊断的结论，不是总体恢复率或全局最优证明。'''


def result_explanation(report):
    if 'coverage' not in report:return ''
    a=next(r for r in report['coverage'] if r['band']=='40' and r['level']==4.)
    b=next(r for r in report['coverage'] if r['band']=='200' and r['level']==4.)
    baseline=next(r for r in report['conditional_metrics'] if r['band']=='200' and r['method']=='original_blind_baseline')
    known=next(r for r in report['conditional_metrics'] if r['band']=='200' and r['method']=='conditional_profile_K2')
    maximum=max(abs(r['interpolation_minus_probe']) for r in report['frames'])
    return f'''**范围宽度也要看**：ΔJ≤4时，40／200ns两组的总跨度／真Q2中位分别为{a['width_over_trueQ'][1]:.3f}和{b['width_over_trueQ'][1]:.3f}；40ns有{a['disconnected']}/16条出现分离区段，{a['lower_boundary']}/16条碰到50pe下界。宽范围可以包含真值，同时仍无法精确分配两个事件的电荷。

**事件数选择会改变分母**：40ns原盲基线只在{a['baseline_selected_K2_frames']}/16条输出两个事件，其中ΔJ≤4区段包含真值为{a['baseline_selected_K2_included']}/{a['baseline_selected_K2_frames']}；全体已知K2样本的{a['included']}/16不能直接当作最终盲输出的覆盖率。200ns原盲基线严格整对恢复{baseline['success']}/16，条件剖面最优点为{known['success']}/16；这批小样本未显示逐事件恢复的普遍改善，也未检验差异显著性。

**数值一致性**：32条×3层级的区段与真实Q2探针包含判定有{report['inclusion_disagreements']}个不一致；{report['probe_below_reference']}个探针低于主参考代价超过0.01，{report['failed_truth_probes']}个探针未收敛。插值ΔJ与实际探针的最大差仍为{maximum:.3f}，不能把成员判定一致扩大为整条连续剖面已全局精确。'''


def notebook():
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    assert load(RECORD/'science-validation.json')['passed'] and load(RECORD/'figure-visual-review.json')['passed']
    original=NOTEBOOK.read_bytes();book=nbformat.reads(original.decode(),as_version=4);before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4)
    assert len(before.cells)==len(book.cells)==518 and before.cells==book.cells
    report=load(OUTPUT/'report.json');cells=[]
    def md(s):cells.append(nbformat.v4.new_markdown_cell(s))
    def code(s):cells.append(nbformat.v4.new_code_cell(s))
    def fig(number,title):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{title}**')
        code(f'display(Image(filename=str(contour_output / "figures/{number}.png")))')
    md('''<a id="contour-boundary-refinement-20261027"></a>
## 108. 把两条仍变化的边界算细

承接第103–107节：上次重点ΔJ≤4的32条均通过末次加密筛查，但条件高斯中间隔第160／190条的ΔJ≤9端点仍变化6.61%／5.40%。本轮定向复算这两条，加原生40／200ns第440／480条、条件高斯长间隔第200条作对照。这5条按既有问题选择，不能当代表性恢复率样本。

保持纯α、固定95／351ns双τ、1024×2ns完整模拟电压和充分动态范围；均为已知K=2，Q2是按时刻排序的第二槽。均值、候选依赖协方差及logdet、Q域50–50000pe和观测K1触发[−96,+512]ns的局部时刻域均保持。没有用真实Q／t构造起点。

先在原网格重拟合跨界节点与最小值邻点：6起点、最多250迭代。随后仅将跨界及最小值相邻单元逐次几何二分，新点4起点／150迭代，直到端点所在lnQ括区不超过1%。再用8个观测起点、300迭代重新核验边界节点；每层级需区段数和域边界状态相同、配对端点变化≤1%、参考最小值改善≤0.01，最多12轮，未通过照实保留。所有搜索均正反扫描并保留较好的已计算点。

“括区1%”指相邻已计算Q点的比值范围，不是电荷误差1%，也不是连续剖面已全局收敛。重拟合和网格曲率的影响分别记录；波形模板本身的插值规则未改。
''')
    code(f'''from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Image
contour_root = Path.cwd().resolve()
while not (contour_root / "pileup/experiments/{NAME}").is_dir():
    if contour_root.parent == contour_root: raise RuntimeError("请从项目根目录或jupyter目录运行")
    contour_root = contour_root.parent
contour_output = contour_root / "jupyter/artifacts/{NAME}"
contour_record = contour_root / "jupyter/run_records/{NAME}"
contour_report = json.loads((contour_output / "report.json").read_text())
display(pd.DataFrame(contour_report["refinement"]))
''')
    fig('0004-108-1','两条问题波形的Q2剖面；(a,c)全Q域，(b,d)ΔJ≤9端点局部。比较旧网格、同网格重拟合及逐步细化，点为实际计算。')
    md('''<a id="contour-numeric-validation-20261027"></a>
## 109. 数值检查与恢复精度分别记录

图109-1显示5条的每步变化。新增检查要求边界所在网格本身足够细，避免仅因两次插值端点碰巧相近就认为已稳定。更强的局部搜索仍是有限搜索，其他未采样区间可能有遗漏的低代价解。

首轮5条算到verify阶段后，父进程因NumPy整数无法直接写JSON而中断，没有保存完整结果行。旧源码、清单和空输出已归档；修复诊断计数的保存类型并增加回归检查后，在相同配置上重新计算。物理模型、起点数、阈值未变，独立新噪声尚未生成时已修复。
''')
    md(boundary_explanation(report))
    fig('0004-109-1','(a)全部5条的最大边界括区，虚线1%；(b)参考最小代价变化。横轴0为同网格重拟合，其后依次细化与核验。')
    code('''display(pd.DataFrame(contour_report["refinement_scores"]))
print("这是已知K2的定向数值诊断；不能把5条的正确数当作总体盲恢复率。")
''')
    if 'coverage' in report:
        md('''<a id="independent-model-coverage-20261027"></a>
## 110. 冻结后新噪声：范围究竟多常包含真值？

5条数值诊断通过后，另冻结完整计算规则、源码、模型、参数设计和种子，再生成44条新波形。双脉冲固定Q=[18000,800]pe、第一时刻200.3ns，40／200ns间隔各16条；另有空窗、800pe孤立和18000pe孤立各4条。各条独立噪声，不同间隔也不复用噪声。Q是模型条件下的等效pe，未标定为沉积能量；受控间隔不等于源注入率。

抽样严格来自冻结JointNoise在FIR变换空间中的相关高斯模型，头三个残差样本条件固定为零，再逆变换回完整float64电压。沿用既有均值／协方差和冻结relative-MLP起点，没有新训练、原生光学母事件或ROOT。新噪声与此前开发抽样独立，但模型由旧资料标定，属于模型内部重复抽样检验。

对每条新双脉冲，先作33点对数网格加观测候选节点，4起点／120迭代、正反扫描；再将所有单元加几何中点，4起点／150迭代，最后执行第108节的自适应边界程序。这里检验的是这套完整冻结的有限算法，不能声称校准了此前的每条旧曲线。MLP仅提供观测起点，未加入学习型Q惩罚；真值文件在主推断进程中禁止读取。

主要“包含”按**实际输出的区段集合是否包含真Q2**计数。主结果哈希固定后，另只固定真实Q2、8观测起点／300迭代拟合其余参数作为评分探针；探针不修改区段或参考最小值，两种判定若不一致单独报告。全部32条保留在主要分母，数值不稳定样本不删除。

以下95%误差条是“包含比例”本身的精确二项抽样区间（SciPy binomtest的exact方法），不是每个脉冲Q2的95%置信区间。ΔJ=1／4／9也未宣称对应68%／95%／99%；本轮不据结果调阈值。每组只有16条，是先导检查。
''')
        md(coverage_table(report))
        md(result_explanation(report))
        fig('0004-110-1','新条件高斯样本：(a)区段包含真Q2的比例及精确二项95%区间；(b)跨度中位与16–84%分位；(c)原盲基线和已知K2的严格整对恢复数。')
        code('''display(pd.DataFrame(contour_report["coverage"]))
with pd.option_context("display.max_columns", None):
    display(pd.DataFrame(contour_report["conditional_metrics"]))
print("盲控制使用原baseline候选和原门槛；新剖面最优点只作已知K2评分，未接入盲门槛。")
print("Q偏差和sigma68来自时间匹配后的条件样本，sigma68不是能谱FWHM。")
''')
        md('严格整对恢复要求每个Q误差≤20%、每个时刻误差≤min(10ns,gap/4)，且没有额外输出。这里两档gap的时限均为10ns。表中missed／extra按10ns时间匹配定义；target−recovered还包括已匹配但电荷精度不合格的事件，不能把它全部叫作漏检。')
        md('''<a id="all-coverage-repetitions-20261027"></a>
## 111. 展示全部32条的相容区段

按预定编号展示所有新双脉冲，各层级的分离区段和边界截断全部保留。横轴除以800pe只为评分展示，真实电荷不参与区段构造；总跨度包含区段间不相容的空隙。范围很宽时容易包含真值，这不代表已经精确恢复弱脉冲。

两面板均用对数横轴，但显示范围不同；跨组比较宽度时使用第110节数值表。
''')
        fig('0004-111-1','(a)40ns、(b)200ns各16条的全部区段；三色为ΔJ≤1／4／9，竖线为真Q2。')
        code('''display(pd.DataFrame(contour_report["frames"]))
print({k:contour_report[k] for k in ["stable_frames", "probe_below_reference", "failed_truth_probes", "inclusion_disagreements"]})
print("数值筛查、实际真Q探针与区段插值分别记录；负Δ探针不能解释为更好的统计覆盖。")
''')
        md('''<a id="coverage-full-waveforms-20261027"></a>
## 112. 新独立噪声样本的完整波形与分量

例子在查看新结果前固定为各组前两条：0／1／16／17。左列为完整0–2048ns总电压（采样点0–2046ns），右列为后槽分量。黑线是模型在真参数处的条件均值；高斯残差由联合协方差生成，未声称可唯一拆成两个真实随机分量。全部拟合分解仍假定K=2。
''')
        fig('0004-112-1','四条预定新样本的全记录总波形及第二分量；(a–d)40ns，(e–h)200ns。')
        code('''contour_examples = json.loads((contour_output / "examples.json").read_text())
display(pd.DataFrame([dict(index=e["index"], true_Q=e["true_charge_pe"], true_t=e["true_times_ns"],
    fitted_Q=e["fit"]["charge_pe"], fitted_t=e["fit"]["times_ns"]) for e in contour_examples]))
''')
    md(f'''<a id="contour-coverage-conclusion-20261027"></a>
## 113. 结论、限制与下一步

{summary(report)}

区段数值稳定、包含率高、参数估计准确是三个不同问题。未通过数值检查、探针低于参考、边界截断和评分不一致均保留；不能靠删失败样本提高包含率。对于两个参数点之外的Q比、时距、free-τ、α/γ混合和高率K≥3，本轮没有验证。

下一步应先依据这里的数值失败／探针差异决定是否需继续处理搜索；数值问题处理后扩大固定参数点和独立重复次数。若要选择“95%”阈值，需另设校准样本并在新独立样本上检验，同时保留检测选择及范围宽度；原生模型迁移另做评价。

{LITERATURE}

本地完整[结果](artifacts/{NAME}/report.json)、[边界变化](artifacts/{NAME}/refinement-summary.csv)及运行记录保留。此前518个单元和图像保持，所有新增代码输出已执行。
''')
    for i,c in enumerate(cells):c.id=f'contour-20261027-{i:02d}'
    InteractiveShell.instance();space={};count=0;cwd=os.getcwd()
    try:
        os.chdir(PROJECT/'jupyter')
        for c in cells:
            if c.cell_type!='code':continue
            with capture_output(stdout=True,stderr=True,display=True) as cap:exec(compile(c.source,'<contour-notebook>','exec'),space)
            assert not cap.stderr,cap.stderr;c.outputs=[]
            if cap.stdout:c.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=cap.stdout))
            for output in cap.outputs:c.outputs.append(nbformat.v4.new_output('display_data',data=output.data,metadata=output.metadata))
            count+=1;c.execution_count=count
    finally:os.chdir(cwd)
    book.cells.extend(cells);assert book.cells[:518]==before.cells;nbformat.validate(book)
    encoded=(nbformat.writes(book)+'\n').encode()
    with tempfile.NamedTemporaryFile(dir=NOTEBOOK.parent,prefix=NOTEBOOK.stem+'.',suffix='.partial',delete=False) as stream:
        for i in range(0,len(encoded),256*1024):
            chunk=encoded[i:i+256*1024];assert stream.write(chunk)==len(chunk)
        stream.flush();os.fsync(stream.fileno());temp=type(NOTEBOOK)(stream.name)
    assert temp.read_bytes()==encoded and NOTEBOOK.read_bytes()==original
    temp.replace(NOTEBOOK);assert NOTEBOOK.read_bytes()==encoded
    save(RECORD/'notebook-delivery.json',dict(passed=True,old_cells=518,new_cells=len(cells),executed_cells=count,
        full_readback=True,history_preserved=True,notebook_sha256=digest(NOTEBOOK),source_sha256=digest(__file__)))
    print('CONTOUR_NOTEBOOK_DELIVERED',len(cells),count,flush=True)


def documents():
    assert load(RECORD/'notebook-delivery.json')['passed'];report=load(OUTPUT/'report.json');result=summary(report)
    root=f'<a id="{MARKER}"></a>\n2026-10-27完成[边界细化与新条件高斯重复抽样](pileup/experiments/{NAME}/README.md)：{result} 本地Notebook0004第108–113节保存固定图号、全部区段与完整波形；方法条件和文献采用见[L65](docs/literature-use.md#{MARKER})。原默认盲算法保持。'
    text=f'''<a id="{MARKER}"></a>
### L65：2026-10-27——边界括区细化与固定模型覆盖率先导

**依据与审核范围**：{LITERATURE}

**实际结果**：{result}

{boundary_explanation(report)}

{coverage_table(report)}

{result_explanation(report)}

**方法及区别**：两条上轮ΔJ≤9端点变化超过5%的高斯波形160／190，加原生440／480和高斯200为定向数值诊断。先同网格6起点／250迭代重拟合跨界与最小值邻点，再几何二分，新节点4起点／150迭代；括区≤1%后以8起点／300迭代核验，要求端点≤1%、参考代价改善≤0.01、区段数及边界状态不变，最多12轮。均正反扫描、保留旧更优点及阶段最优固定其余参数切片。有限搜索不保证全局最优，括区1%不代表电荷精度1%。波形模板插值规则、均值和协方差未改。

数值门通过后冻结完整新样本算法：33点基础logQ网格加观测候选，4起点／120迭代；所有单元中点4起点／150迭代，随后上述边界程序。Q=[18000,800]pe、t1=200.3ns、gap40／200ns各16次；另空／弱单／强单各4次。NumPy SeedSequence[202610271801,10000+index,994]独立噪声；固定模型FIR空间相关高斯、前三残差条件为零、逆变换成1024×2ns float64电压。纯α固定95／351ns、充分动态范围、Q域50–50000pe、观测K1局部时域和已知K2保持。冻结旧relative-MLP只给起点，没有学习型Q惩罚或新训练。无新ROOT／光学母事件；不能称原生物理独立确认。

主要包含率以实际输出区段集合包含真Q计数，所有32条含数值失败都留在分母。真实Q2探针只在主结果锁定哈希后评分：固定Q2、8观测起点／300迭代拟合未知Q1/t，不反馈主曲线。插值与探针不一致、探针低于参考、未收敛和域边界单列。表中精确二项95%区间只描述包含比例的抽样不确定性，不把ΔJ=1／4／9命名为68／95／99%电荷区间；每组16条只属先导，不用它反调阈值。盲空／单控制沿用原baseline及原门槛，未拿新剖面替换盲候选。完整指标同时记录Q偏差／sigma68、时误、漏检、额外、误拆／误并与严格恢复；Q非沉积能量，sigma68非FWHM，受控gap非源率。

**运行与复现**：实际源码、模型和协议在[pileup/experiments/{NAME}](../pileup/experiments/{NAME}/README.md)；主数值推断阻断truth／cases／ROOT。首轮仅因NumPy计数JSON序列化失败，旧源码/清单/空输出保留，修复类型并加回归后按同配置重算，在新抽样前完成。使用最多8核、BLAS1。Notebook路径jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb，第108–113节、图0004-108-1至112-1；旧518单元及结果保留。

**待验证**：处理残余数值问题后扩大独立重复次数和Q／gap参数点；另立校准与确认样本再选择覆盖阈值；原生迁移、free-τ、α/γ与高率K≥3分别检验。有限样本／非规则剖面统计的一手证据及Baker DOI仍需共享库另行补充，本项目不修改文献库。
'''
    for name,backup,addition in (('README.md','README.before.md',root),('docs/literature-use.md','literature.before.md',text)):
        p=PROJECT/name;old=(RECORD/backup).read_text();assert p.read_text()==old and MARKER not in old
        p.write_text(old.rstrip()+'\n\n'+addition.rstrip()+'\n')
    readme=f'''# Adaptive charge contours and fixed-model coverage pilot

{result}

{boundary_explanation(report)}

{coverage_table(report)}

{result_explanation(report)}

The numerical protocol targets the two previously changing level-9 Gaussian boundaries plus three prespecified controls. Actual algorithms, bounds, tolerances, budgets and selection are frozen in protocol.json; source and model hashes are in freeze.json. Initial serialization failure, original source in inference-source-v1, original freeze-v1.json and the empty interrupted output are retained. The fix only changes result serialization, not the physical model or numerical settings.

After all five diagnostics pass, coverage-protocol.json and coverage-freeze.json freeze a new, complete algorithm and independent-noise design before generation. Sixteen repeated pairs at each40/200ns gap have fixed Q=[18000,800]pe and t1=200.3ns. Twelve empty/single controls use the original blind baseline/gate. All32 pairs, including numerical failures, stay in the primary denominator. The main inclusion metric is membership in the computed, interpolated interval set; downstream fixed-trueQ2 probes are separate diagnostics and never alter the profiles. Exact binomial95% intervals quantify Monte Carlo uncertainty of inclusion proportions, not nominal charge intervals. Small per-condition samples cannot certify95% coverage. No threshold tuning, new optical parents, ROOTs, training or default change.

Run from the project root with the project Python environment, OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1. Up to8 workers; refinement5, pilot2, main7. Commands use frozen source/model snapshots rather than live code, with hashes checked; no Git checkout needed. Arrays, weights, Notebook and run records are local ignored inputs and require a separate reproducibility bundle.

```bash
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py tests
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/coverage.py tests
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/report.py analyze
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/report.py figures
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/deliver.py validate
```

The run/generate/pilot/infer/probe commands refuse existing outputs. Recompute in a separate project tree with all frozen inputs. coverage.py freeze requires the completed five-frame numerical gate. Source and seeds are fixed before any new-noise result, and pilot frames are reused exactly in the full sample. Full-waveform examples0/1/16/17 are prespecified. Inference workers cannot open truth/cases/ROOT; scoring-only trueQ2 probes run after primary hashes freeze. Full cost replay, sampler covariance/FIR checks, sixteenth-sample denominator checks, old-file hashes, Notebook preservation and embedded figures are validated. Notebook writes use256KiB chunks and full readback on shared storage.

Scope: pure alpha, fixed95/351ns, full1024x2ns analog voltage, ample dynamic range, knownK2 ordered later charge, observed-trigger-local time domain and50-50000pe charge box. This is frozen-model repeated sampling, not independent nativeSiPM validation, deposited-energy resolution, a universal minimum spacing or high-rateK>=3 recovery. Numerical stability, interval width/inclusion and individual-event accuracy are reported separately. Existing Baker/Cousins1984 and Liddick2012 notes are reused with detector/statistical conditions and AI-pending-review status; see docs/literature-use.md L65. Expand parameter points and independent samples after resolving numeric diagnostics; use separate calibration and confirmation samples for future threshold selection.

Local outputs: jupyter/artifacts/{NAME}/; records: jupyter/run_records/{NAME}/. Notebook0004 sections108–113 retain all518 preceding cells and add executed diagnostics and numbered figures. No Git commit or push in this study.
'''
    p=EXPERIMENT/'README.md';assert not p.exists();p.write_text(readme)
    save(RECORD/'document-delivery.json',dict(passed=True,files={str(x.relative_to(PROJECT)):digest(x) for x in (PROJECT/'README.md',PROJECT/'docs/literature-use.md',p)}))
    print('CONTOUR_DOCUMENTS_DELIVERED',flush=True)


def validate():
    check_freeze();from .alpha_pair_coverage_study import check_coverage
    check_coverage()
    for name in ('numerical-tests','refinement-validation','sampler-tests','coverage-generation','coverage-pilot-validation','coverage-inference',
                 'coverage-probe-validation','science-validation','figure-delivery','figure-visual-review','notebook-delivery','document-delivery'):
        assert load(RECORD/(name+'.json'))['passed'],name
    for manifest in ('freeze','coverage-freeze','report-freeze','delivery-freeze'):
        for name,h in load(EXPERIMENT/(manifest+'.json'))['files'].items():assert digest(EXPERIMENT/name)==h,name
    for name,h in load(EXPERIMENT/'freeze-v1.json')['files'].items():
        p=EXPERIMENT/name.replace('inference-source/','inference-source-v1/',1)
        assert digest(p)==h,name
    source={}
    for folder,names in (('inference-source/pileup/analysis',('alpha_pair_profile_contours.py','test_alpha_pair_profile_contours.py','alpha_pair_contour_study.py')),
                         ('coverage-source',('alpha_pair_coverage_study.py',)),('report-source',('alpha_pair_contour_report.py',)),('delivery-source',('alpha_pair_contour_delivery.py',))):
        for name in names:
            p=PROJECT/'pileup/analysis'/name;ast.parse(p.read_text());assert digest(p)==digest(EXPERIMENT/folder/name);source[str(p.relative_to(PROJECT))]=digest(p)
    provenance=load(RECORD/'input-provenance.json');prior=provenance['prior_validation']
    for name,h in provenance['inputs'].items():assert digest(PROJECT/name)==h,name
    for name,h in prior['source_sha256'].items():assert digest(PROJECT/name)==h,name
    for name,h in provenance['prior_freezes'].items():assert digest(PREVIOUS/name)==h,name
    assert digest(OLD_RECORD/'input-provenance.json')==provenance['inherited_sha256']
    for n,h in load(RECORD/'science-validation.json')['output_hashes'].items():assert digest(OUTPUT/n)==h,n
    run=load(RECORD/'coverage-inference.json');probe=load(RECORD/'coverage-probe-validation.json')
    assert digest(OUTPUT/'coverage/predictions.jsonl')==run['sha256']==probe['primary_sha256']
    assert digest(OUTPUT/'coverage/truth-probes.jsonl')==probe['sha256']
    assert max(run['maximum_cost_replay_error'],probe['maximum_cost_replay_error'],load(RECORD/'refinement-validation.json')['maximum_cost_replay_error'])<1e-9
    gen=load(RECORD/'coverage-generation.json');freeze=load(EXPERIMENT/'coverage-freeze.json')
    assert gen['source_freeze_sha256']==digest(EXPERIMENT/'coverage-freeze.json') and gen['created']>freeze['created']
    for n,h in gen['files'].items():assert digest(OUTPUT/'coverage'/n)==h,n
    pred=rows(OUTPUT/'coverage/predictions.jsonl');pilot=rows(OUTPUT/'coverage/pilot-predictions.jsonl')
    assert len(pred)==44 and len({r['index'] for r in pred})==44 and {r['index'] for r in pilot}=={0,16}
    lookup={r['index']:r for r in pred}
    for p in pilot:assert lookup[p['index']]==p
    report=load(OUTPUT/'report.json');assert report['all32_retained'] and all(c['frames']==16 for c in report['coverage'])
    assert len(load(OUTPUT/'coverage-intervals.json'))==96 and len(load(OUTPUT/'coverage-scores.json'))==76
    program='''import sys,json
sys.path.insert(0,sys.argv[1])
import pileup.analysis
pileup.analysis.__path__.insert(0,sys.argv[2])
from pileup.analysis.alpha_pair_coverage_study import initialize
initialize()
blocked=[]
for name in ('forbidden.truth.npz','forbidden.cases.json','forbidden.root','/tmp/truth/forbidden.npy','native-records.json'):
    try:open(name,'rb')
    except PermissionError:blocked.append(name)
    else:raise AssertionError(name)
print(json.dumps(blocked))
'''
    import json
    blocked=json.loads(subprocess.check_output([sys.executable,'-c',program,str(EXPERIMENT/'inference-source'),str(EXPERIMENT/'coverage-source')],cwd=PROJECT,text=True));assert len(blocked)==5
    docs=load(RECORD/'document-delivery.json')['files']
    for n,h in docs.items():assert digest(PROJECT/n)==h,n
    for n,b in (('README.md','README.before.md'),('docs/literature-use.md','literature.before.md')):
        assert digest(RECORD/b)==prior['documentation_sha256'][n]
        assert (PROJECT/n).read_text().startswith((RECORD/b).read_text().rstrip()+'\n\n')
    before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4);book=nbformat.read(NOTEBOOK,as_version=4);nbformat.validate(book)
    assert len(before.cells)==518 and book.cells[:518]==before.cells
    assert digest(RECORD/'notebook.before.ipynb')==prior['notebook_sha256']
    delivery=load(RECORD/'notebook-delivery.json');assert digest(NOTEBOOK)==delivery['notebook_sha256']
    new=book.cells[518:];codes=[c for c in new if c.cell_type=='code'];assert len(new)==delivery['new_cells'] and len(codes)==delivery['executed_cells']
    assert all(c.outputs and c.execution_count and not any(o.output_type=='error' for o in c.outputs) for c in codes)
    assert len({c.id for c in book.cells})==len(book.cells)
    anchors=re.findall(r'id="codex-figure-([^"]+)"','\n'.join(c.source for c in book.cells if c.cell_type=='markdown'));assert len(anchors)==len(set(anchors))
    for c in new:
        if c.cell_type=='markdown':assert not re.search(r'\]\([^)]*\.(?:pdf|ipynb|md)(?:[#?)][^)]*)?\)',c.source,re.I)
    figs=load(RECORD/'figure-delivery.json')['figures'];assert figs==load(RECORD/'figure-visual-review.json')['figures']
    for n,h in figs.items():assert digest(OUTPUT/'figures'/(n+'.png'))==h
    embedded=[hashlib.sha256(base64.b64decode(o.data['image/png'])).hexdigest() for c in codes for o in c.outputs if 'image/png' in o.get('data',{})]
    assert len(embedded)==5 and set(embedded)==set(figs.values())
    for p in (NOTEBOOK,OUTPUT/'report.json',EXPERIMENT/'coverage-relative-model.npz'):
        assert subprocess.run(['git','check-ignore','-q',str(p.relative_to(PROJECT))],cwd=PROJECT).returncode==0
    subprocess.run(['git','-c','core.fileMode=false','diff','--check','--','README.md','docs/literature-use.md'],cwd=PROJECT,check=True,timeout=120)
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip()==load(RECORD/'environment.json')['git_head']
    save(RECORD/'final-validation.json',dict(passed=True,numerical_tests=load(RECORD/'numerical-tests.json')['tests'],sampler_tests=5,
        refinement_frames=5,independent_noise_frames=44,knownK2_profiles=32,scoring_probes=32,contour_scores=96,recovery_scores=76,
        original518_cells_preserved=True,new_cells=len(new),executed_cells=len(codes),figures=figs,notebook_sha256=digest(NOTEBOOK),
        source_sha256=source,documentation_sha256=docs,truth_open_probes_blocked=blocked,old_results_preserved=True,coverage_thresholds_untuned=True,
        maximum_cost_replay_error=max(run['maximum_cost_replay_error'],probe['maximum_cost_replay_error']),
        new_ROOTs=0,new_training=0,no_default_change=True,no_git_commit_or_push=True,
        limitations=['16 repeats per gap','fixed-parameter conditional Gaussian draws','knownK2,fixedtau,localtime domain',
                     'finite search only','not a nominal-confidence calibration or native confirmation']))
    print('CONTOUR_FINAL_VALIDATION_PASSED',len(new),len(codes),flush=True)
