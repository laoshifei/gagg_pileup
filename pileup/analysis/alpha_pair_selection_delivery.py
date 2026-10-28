"""Local executed Notebook and audit records for the selection experiment."""
from __future__ import annotations

import ast
import base64
import hashlib
import os
import re
import subprocess
import tempfile

import nbformat
import numpy as np

from .alpha_pair_selection_study import (PROJECT, EXPERIMENT, OUTPUT, RECORD, OLD_OUT, NOTEBOOK, NAME,
    METHODS, check_freeze, inclusion_partition, design)
from .alpha_pair_selection_report import gates
from .alpha_pair_information_study import load, save, digest, rows
from .alpha_pair_joint_likelihood import calibrate_gate

MARKER='alpha-pair-selection-20261027'
LITERATURE='''本轮按paper-search只读检索共享库的主题／审核索引、数字脉冲处理／PSD／光产额主题页，复用Liddick、Darby、Grzywacz（2012）Algorithms for pulse shape analysis using silicon detectors完整既有精读笔记，未重读PDF或联网。来源papers/sources/2012_Liddick_NIMA_Algorithms-for-pulse-shape-analysis-using-silicon-detectors.md；DOI10.1016/j.nima.2011.12.034；据笔记PDF3–8、图3–13。对象为约65μm DSSD、40MSPS／12bit、25μs波形及8000余人工叠加双脉冲，采用“识别效率、已识别样本参数精度及假双背景分开报告”的方法思路；不移植100ns、幅比5或ADC精度结论。笔记仍为AI草稿待人工审核，门选记号与自由参数个数疑问保留，没有新的科学采用确认。

本轮分母分解、固定单脉冲校准点、秩门槛和样本数是项目的诊断设计，不是该论文提供的选择后置信区间算法。ΔJ=1／4／9仍不赋予名义置信水平。选择后区间校准和具有未知参数时的误报控制仍需专门的一手统计依据；共享库未修改。'''


def primary_table(report):
    lines=['| 间隔/ns | 门槛校准来源 | 输出双事件 | 严格整对恢复 | 时间漏检事件 | 额外输出事件 |',
           '| --- | --- | --- | --- | --- | --- |']
    for band in ('pair40','pair200'):
        for method,label in (('original_025','原门槛'),('matched_025','新条件高斯单脉冲')):
            r=next(r for r in report['summary'] if r['band']==band and r['method']==method)
            lines.append(f"| {band[4:]} | {label} | {r['K2_outputs']}/64 | {r['success']}/64 | {r['missed']}/128 | {r['extra']} |")
    return '\n'.join(lines)


def conclusions(report):
    parts=[]
    for band in ('pair40','pair200'):
        a=next(r for r in report['summary'] if r['band']==band and r['method']=='original_025')
        b=next(r for r in report['summary'] if r['band']==band and r['method']=='matched_025')
        cap=next(r for r in report['forced_candidate_limit'] if r['band']==band)
        parts.append(f"{band[4:]}ns：同为0.025门槛参数，改用新单脉冲校准后，输出两个事件由{a['K2_outputs']}/64变为{b['K2_outputs']}/64，严格整对恢复由{a['success']}/64变为{b['success']}/64。即使直接给这批真实双脉冲选K=2，当前候选中严格成功也只有{cap['success']}/64；这是同一有限候选的诊断上限，不是物理可辨识极限。")
        if band=='pair40':
            old=next(r for r in report['selection_bias'] if r['band']==band and r['method']=='original_025' and r['state']=='selected')
            new=next(r for r in report['selection_bias'] if r['band']==band and r['method']=='matched_025' and r['state']=='selected')
            parts.append(f"40ns被选中候选的有序Q2误差中位由{100*old['ordered_Q2_error_p16_p50_p84'][1]:+.1f}%变为{100*new['ordered_Q2_error_p16_p50_p84'][1]:+.1f}%，这是选择群体改变，不是同一条拟合参数变准；包含时间未匹配候选，不称能量分辨率。时间未匹配的额外输出由{a['extra']}增为{b['extra']}个，多检出没有转化成相应的准确恢复。额外输出表示10ns内找不到匹配真事件，不必然表示整窗原本无信号。")
    a=[r for r in report['summary'] if r['method']=='original_025' and r['kind']=='single']
    b=[r for r in report['summary'] if r['method']=='matched_025' and r['kind']=='single']
    parts.append(f"四档单脉冲共256条，原／新主门槛误报双事件分别为{sum(r['K2_outputs'] for r in a)}/{sum(r['K2_outputs'] for r in b)}条；含任意时间未匹配输出的帧分别为{sum(r['extra_frames'] for r in a)}/{sum(r['extra_frames'] for r in b)}条。各Q档结果另列，不能把这个等权混合比例当作任意输入谱的误报率。")
    return '\n\n'.join(parts)


def notebook():
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    assert load(RECORD/'science-validation.json')['passed'] and load(RECORD/'figure-visual-review.json')['passed']
    original=NOTEBOOK.read_bytes();book=nbformat.reads(original.decode(),as_version=4)
    before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4)
    assert len(book.cells)==543 and book.cells==before.cells
    report=load(OUTPUT/'report.json');cells=[]
    def md(s):cells.append(nbformat.v4.new_markdown_cell(s))
    def code(s):cells.append(nbformat.v4.new_code_cell(s))
    def fig(number,title):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{title}**')
        code(f'display(Image(filename=str(selection_output / "figures/{number}.png")))')
    md('''<a id="selection-denominators-20261027"></a>
## 114. 把“选出事件”和“区间包含真值”拆开

承接第108–113节，本轮首先复用原32条条件高斯双脉冲及其已算好的Q2区段，检查事件数选择。令S表示盲算法输出K=2，C表示已知K=2计算的后槽Q2范围包含真值。分别记录P(C)、P(S)、P(C|S)、P(C且S)、P(C|非S)。没有选中事件时P(C|S)未定义，程序保留空值。

原40ns组ΔJ≤4：全体14/16包含真Q，原门槛只选出6/16，其中包含真Q为4/6；“被选中且包含”是4/16。200ns组相应为16/16、16/16、16/16。这三个分母不能互换。原40ns被选中的6条，其原K2候选Q2全部偏高；这提示选择偏差，但这批旧样本已经看过，仅作探索诊断。新门槛作用于旧区间的重评分也不作为独立确认。

区间本身保持原样，仍是有序后槽Q2的条件K2范围；即使包含真Q，也不自动满足两时刻和两电荷的逐事件恢复要求。本轮没有计算新的区间或宣称校准了选择后覆盖率。
''')
    code(f'''from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Image
selection_root = Path.cwd().resolve()
while not (selection_root / "pileup/experiments/{NAME}").is_dir():
    if selection_root.parent == selection_root: raise RuntimeError("请从项目根目录或jupyter目录运行")
    selection_root = selection_root.parent
selection_output = selection_root / "jupyter/artifacts/{NAME}"
selection_report = json.loads((selection_output / "report.json").read_text())
selection_old = [r for r in selection_report["old_selection_summary"] if r["method"] == "original_025"]
display(pd.DataFrame([{{"gap_ns": r["band"], "level": r["level"],
    **{{k: str(r[k]["numerator"])+"/"+str(r[k]["denominator"]) for k in ["all","selected","unselected","selected_and_included"]}}}} for r in selection_old]))
''')
    fig('0004-114-1','(a,b)旧40／200ns全部候选的电荷与证据分数；红色为被选中。(c)ΔJ≤4的三个不同分母。绿色带为Q2误差±20%，真值仅用于评分。')
    md('''<a id="selection-calibration-20261027"></a>
## 115. 只改变选门：先校准，再独立确认

物理拟合器waveform_candidates保持不变：仅输入完整电压，生成K=0／1／2候选，仍使用固定95／351ns双τ、候选依赖协方差与logdet。没有新网络训练，也没有用新剖面最优点替换盲候选。

先冻结代码、模型、设计和种子，生成256条校准单脉冲：Q=800／3000／8000／18800pe各64条，起点200.3ns。按拟合K1电荷的1500／5000／12500pe边界分四档，保存每档单脉冲的D=max(J1−J2,0)；若K2任一Q<500pe则D=0。观测波形的秩分数为r=(1+校准中D_i≥D的个数)/(n+1)。只有K1存在条件、两个Q≥500pe、D>0且r≤门槛参数时才接受K2。相等值保守计入，不能把r误当电荷置信度。

预先固定六种流程：原校准／新条件高斯校准，门槛参数各0.025／0.05／0.10。主比较为两种0.025；其余四种只展示检出与误报取舍。新校准冻结之后才生成448条确认波形：空窗、上述四种单脉冲、40ns双脉冲、200ns双脉冲各64条。双脉冲Q=[18000,800]pe、t1=200.3ns。每个样本独立噪声，不同gap不复用抽样；没有看确认真值再选择门槛。

两种校准的每档样本数不同，秩分辨率也不同：n=64时最小非零秩为1/65，0.025只接受这一档；n=80时可接受至2/81。因此同名0.025不能解释成完全相同的实际误报预算。下表保存每档n、最小秩、门槛内最大可取秩与实际D切点；差异同时包括噪声来源、有限抽样和秩离散度，不能全归因于模型匹配。

抽样仍是冻结模型内的相关高斯：FIR空间生成，头三点残差固定为零，逆变换为1024×2ns float64电压。新种子SeedSequence[202610272001,split,index,995]，split为1／2。充分动态范围，无ADC削顶、原生光学母事件或新ROOT。四个固定单脉冲点及拟合Q分档不能保证所有未知Q/t条件下都控制误报；0.025是选门参数，不是已经证明的2.5%物理误报率。

运行中，6 worker校准在保存200条后约20分钟未继续写入，原因未确定。停止该批次，完整保留200条原始记录，另冻结续算脚本，以2 worker仅完成缺少的56条；后续确认也用2 worker。物理拟合、模型、波形、随机种子和选门规则不变，续算没有重算或删去已完成样本。中断记录及逐行一致性检查保留。
''')
    code('''display(pd.DataFrame(selection_report["thresholds"]))
print("strict_evidence_cutoff表示D必须严格大于该值；此外仍要求存在性、两个Q>=500和D>0。")
''')
    fig('0004-115-1','四个拟合Q档的校准单脉冲D经验分布；灰色原校准、蓝色新条件高斯校准。竖虚线为0.025对应的严格D门槛。')
    md('''<a id="selection-confirmation-20261027"></a>
## 116. 新确认样本：检出、严格恢复和假双分别比较

以下每档64条，全部保留。严格整对恢复要求两个Q误差各≤20%、两个时刻各≤min(10ns,gap/4)，且无额外时间未匹配输出；这里时限均为10ns。计数正确不等于严格恢复。missed是10ns内时间未匹配的真实事件；target−recovered还含电荷误差不合格，不能都称漏检。

单脉冲控制同时报告“错误输出K=2”与“有任意时间未匹配输出的帧”，两者定义不同。全指标另含Q偏差、sigma68、时间误差、漏检、额外、误拆、误并和事件恢复数；Q尚未标定成沉积能量，sigma68=(p84−p16)/2不是FWHM。

图中的精确二项95%区间只描述固定门槛下、这些固定参数点的确认样本计数不确定性，不包含有限校准样本的不确定性，也不是电荷的95%区间。所有方法在相同确认波形上比较，配对得失另外保存。
''')
    md(primary_table(report))
    fig('0004-116-1','(a)真实双脉冲中输出K=2的比例；(b)严格恢复两事件的比例。原／新校准各三档门槛参数，误差条为精确二项95%抽样区间。')
    fig('0004-116-2','新空／单控制：(a)错误双事件输出；(b)至少一个时间未匹配输出。每个格子分母都是64，空窗也保留。')
    code('''selection_metrics = pd.DataFrame(selection_report["summary"])
with pd.option_context("display.max_columns", None):
    display(selection_metrics[selection_metrics.method.isin(["original_025", "matched_025"])])
display(pd.DataFrame(selection_report["paired"]))
''')
    md('''<a id="selection-bias-20261027"></a>
## 117. 门槛偏向什么样的弱脉冲？

每条真实双脉冲的原K2候选在所有门槛中相同，只有是否接受改变。图117-1保留全部128条，包括未通过门槛、时间不准和Q<500的候选；纵轴是按拟合时刻排序的第二槽Q2/800−1，并非只在正确时间匹配后计算的能量分辨率。绿色带为±20%，纵轴在中心线性、远处对称对数，避免大误拆值被截掉。

下面同时列出全体、被选中和未选中候选的Q2及t2误差分位数。再给“真实双脉冲全部强制选K=2”的诊断：门槛只能放行已有候选，不能把电荷／时刻错误的候选变正确。这个上限仅针对当前候选和本批样本，不能解释成别的算法也无法恢复。
''')
    fig('0004-117-1','(a,b)40ns，(c,d)200ns；左列原0.025门槛、右列新0.025门槛。同一点的拟合参数保持不变，颜色表示选择结果。')
    code('''display(pd.DataFrame([r for r in selection_report["selection_bias"] if r["method"] in ["original_025","matched_025"]]))
display(pd.DataFrame(selection_report["forced_candidate_limit"]))
''')
    md('''<a id="selection-full-waveforms-20261027"></a>
## 118. 预先指定的完整波形

确认样本例子在生成前固定为40ns第320／321条、200ns第384／385条，不按效果挑选。左列显示0–2048ns全记录总波形；右列显示真实参数的弱分量条件均值，以及原／新主门槛实际输出的第二分量。若只输出一个事件，图例明确写no second output，不把隐藏K2候选画成已恢复事件。

黑线是无噪声条件均值。这里的联合高斯残差没有唯一的逐粒子分量，不声称画出了各粒子真实随机涨落。两个门槛均接受K2时，两条拟合曲线会完全重合，因为物理拟合参数没有变化。
''')
    fig('0004-118-1','四条预定确认波形的全记录及后槽输出；(a–d)40ns、(e–h)200ns。')
    md('''<a id="selection-conclusion-20261027"></a>
## 119. 本轮结论与下一步
''')
    md(conclusions(report))
    md('''本轮把旧区间的选择影响算清，并完成独立的新选门检验。原默认算法保持不变。下一步需要解决被接受候选中的电荷／时刻误差，并在选定的完整推断流程下另设区间校准及确认样本，把“未选中”和“选中但范围不覆盖”都保留；不能用本轮的检出改善代替覆盖验证。随后扩大Q／时距参数点并做原生SiPM迁移，free-τ、α/γ、高率K≥3和计数率扫描分别研究。

受控40／200ns间隔不是源注入率，704条新条件高斯波形也不是704个新光学母事件。本轮所有确认结果都局限于固定模型和参数点，有限多起点优化不保证全局最优。
''')
    md(LITERATURE)
    md(f'''完整[结果](artifacts/{NAME}/report.json)和[逐条件指标](artifacts/{NAME}/summary.csv)保存在本地。旧543单元及图像保持；新增代码均执行，图号114–118固定保留。''')
    for i,c in enumerate(cells):c.id=f'selection-20261027-{i:02d}'
    InteractiveShell.instance();space={};count=0;cwd=os.getcwd()
    try:
        os.chdir(PROJECT/'jupyter')
        for c in cells:
            if c.cell_type!='code':continue
            with capture_output(stdout=True,stderr=True,display=True) as cap:exec(compile(c.source,'<selection-notebook>','exec'),space)
            assert not cap.stderr,cap.stderr;c.outputs=[]
            if cap.stdout:c.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=cap.stdout))
            for output in cap.outputs:c.outputs.append(nbformat.v4.new_output('display_data',data=output.data,metadata=output.metadata))
            count+=1;c.execution_count=count
    finally:os.chdir(cwd)
    book.cells.extend(cells);nbformat.validate(book);assert book.cells[:543]==before.cells
    encoded=(nbformat.writes(book)+'\n').encode()
    with tempfile.NamedTemporaryFile(dir=NOTEBOOK.parent,prefix=NOTEBOOK.stem+'.',suffix='.partial',delete=False) as stream:
        for i in range(0,len(encoded),256*1024):
            chunk=encoded[i:i+256*1024];assert stream.write(chunk)==len(chunk)
        stream.flush();os.fsync(stream.fileno());temp=type(NOTEBOOK)(stream.name)
    assert temp.read_bytes()==encoded and NOTEBOOK.read_bytes()==original
    temp.replace(NOTEBOOK);assert NOTEBOOK.read_bytes()==encoded
    save(RECORD/'notebook-delivery.json',dict(passed=True,old_cells=543,new_cells=len(cells),executed_cells=count,
        notebook_sha256=digest(NOTEBOOK),history_preserved=True,full_readback=True))
    print('SELECTION_NOTEBOOK_DELIVERED',len(cells),count,flush=True)


def documents():
    assert load(RECORD/'notebook-delivery.json')['passed'];report=load(OUTPUT/'report.json');result=conclusions(report)
    root=f'''<a id="{MARKER}"></a>
2026-10-27完成[事件选择与条件高斯门槛确认](pileup/experiments/{NAME}/README.md)：固定拟合器，256条新单脉冲校准后另生成448条独立空／单／双确认波形；同时分解旧32条区间的选择分母。{result} 本轮未计算新的区间覆盖，默认算法保持；本地Notebook0004第114–119节及[L66](docs/literature-use.md#{MARKER})记录条件、全部指标和固定编号完整波形。'''
    literature=f'''<a id="{MARKER}"></a>
### L66：2026-10-27——事件数选择、弱电荷偏差与独立选门确认

**依据及审核**：{LITERATURE}

**已完成**：旧32条区间作探索性分母分解，40ns、ΔJ≤4的P(C)=14/16、P(C|S)=4/6、P(C且S)=4/16；不修改旧区间，不称新覆盖校准。主算法waveform_candidates、均值／协方差／logdet、固定双τ与全模拟电压保持，仅比较选门。冻结源码／模型／种子后生成256条校准单脉冲，800／3000／8000／18800pe各64，按观测K1电荷分档。新秩门槛冻结后才生成448条确认：空／四种单／40ns双／200ns双各64，双Q=18000／800pe、t1=200.3ns。

预先规定原／新校准各门槛参数0.025／0.05／0.10，主比较是0.025；没有根据确认真值反调参数或选择获胜门槛。分数D=max(J1−J2,0)，K2任一Q<500置零；存在性、两Q≥500、D>0及秩条件共同决定接受。有限固定Q点和观测分档不保证未知参数上统一控制误报；秩参数不是电荷置信度。校准样本数不同也改变秩离散度，例如n64的0.025仅接受1/65，n80可到2/81；因此比较同时包含噪声来源、有限抽样和秩网格的差异，不能全归因于模型匹配。逐档n和实际切点已保存。

{primary_table(report)}

{result}

**解释边界**：包含率、选出双事件和严格恢复三者分开。强制K2诊断只指同一有限候选的成功上限，不是信息论或物理极限。单脉冲假双与时间未匹配输出另报。Q为等效pe非沉积能量，sigma68非FWHM；误差分位含按时匹配版本和全体有序后槽版本，后者保留错误时刻／未选候选。确认误差条为固定校准后计数的精确二项95%区间，不含有限校准不确定性、不等于电荷区间。受控时距非源率。

**复现**：实验[pileup/experiments/{NAME}](../pileup/experiments/{NAME}/README.md)含实际源／模型哈希、设计和六门槛。新SeedSequence[202610272001,split,index,995]；同一冻结JointNoise高斯、FIR头三残差为零、1024×2ns float64完整电压，充分动态范围。无ROOT、原生光学母事件或训练，纯α固定95／351ns、局部K≤2。推断worker阻断truth／cases／ROOT；评分在预测锁定哈希后进行。6worker校准保存200条后约20分钟无新写入，原因未定；保留原200条及中断记录，2worker仅续算56条，确认也用2worker；计算参数和数据不变。原Notebook543单元及已有结果保留，新第114–119节，图0004-114-1、115-1、116-1／116-2、117-1、118-1。

**下一步／文献待办**：把电荷／时刻正确性作为候选改进目标，在冻结完整选择流程后另做选择后区间校准及独立确认；扩大Q／gap及原生SiPM迁移。选择后区间、未知参数下误报控制的一手统计依据需共享库另行整理。free-τ、α/γ、K≥3和源率扫描未验证。原默认不变，本轮无commit／push。
'''
    for name,backup,addition in (('README.md','README.before.md',root),('docs/literature-use.md','literature.before.md',literature)):
        p=PROJECT/name;old=(RECORD/backup).read_text();assert p.read_text()==old and MARKER not in old
        p.write_text(old.rstrip()+'\n\n'+addition.rstrip()+'\n')
    readme=f'''# Event selection and independent matched-Gaussian gate confirmation

{result}

{primary_table(report)}

The physical waveform fitter is unchanged. All six decision rules are prespecified: original calibration or fresh matched-Gaussian calibration, each with rank parameter0.025/0.05/0.10. The primary contrast is0.025. Use256 independent single pulses for calibration (64 at each800/3000/8000/18800pe); freeze gates before generating448 new confirmation waves (64 each empty,the four singles,and18000/800pe pairs at40/200ns). First time200.3ns. Actual source/model/design snapshots are recorded in freeze.json,with report/delivery manifests separate. Inference workers cannot open truth/cases/ROOT. Test truth is used only after prediction hashes freeze; no test-based gate tuning or default update.

Previous32 charge profiles are diagnostic data only. Retain all interval components and distinguish P(C),P(S),P(C|S),P(C and S),P(C|not S). No new charge interval profiles or nominal coverage claims. Counterfactual forcedK2 uses the same observed candidates and true-pair evaluation only; its strict-success count bounds gate-only improvement for these candidates,not physical resolvability.

Purealpha,fixed95/351ns,full1024x2ns float64 voltage,ample dynamic range,localK<=2. Exact frozen JointNoise Gaussian in FIR space,conditioned zero first-three residuals. SeedSequence[202610272001,split,index,995],split1calibration/2confirmation; all704 draws independent. No optical mothers,new ROOTs or training. Fitted-Q-bin calibration at four fixed true charges is not a uniform nuisance-parameter false-positive guarantee. Calibration sizes also differ: at parameter0.025,n64 accepts only rank1/65 whereas n80 can accept2/81. Report sample sizes,attainable ranks and actual evidence cutoffs; differences are not solely due to noise-model matching. Exact binomial confirmation bars condition on the frozen calibration and do not include calibration uncertainty. Q is equivalentpe,not deposited energy; controlled gaps are not injection rates.

Run from the project root with the project environment and BLAS1. Frozen scripts use snapshots rather than requiring Git checkout. Pilot2workers,initial calibration6,resumed calibration2,confirmation2; maximum8. Calibration stopped writing after200 rows for about20minutes,cause undetermined. The owned run was interrupted; all200 rows preserved,then only56 missing tasks recomputed using the unchanged fitter with2workers. The separate resume.py/source/manifest and interruption record retain this history. Source/config/run scripts are shareable; ignored local arrays,models and Notebook require a separate reproducibility bundle. Numeric phases refuse existing outputs; replay generation/inference in a separate tree with all frozen inputs. Pilot rows are reused byte-for-byte in calibration.

```bash
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py tests
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py generate-calibration
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py pilot --workers 2
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py infer-calibration --workers 2
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py gates
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py generate-test
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py infer-test --workers 2
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/report.py analyze
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/report.py figures
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/deliver.py validate
```

Local outputs jupyter/artifacts/{NAME}/; records jupyter/run_records/{NAME}/. Notebook0004 sections114–119 preserve543 prior cells and add executed tables and six numbered figures. Full-waveform examples320/321/384/385 are prespecified. Check relativeQ bias/sigma68,time errors,missed/extra/split/merge,strict pair recovery and per-control false outputs separately. Existing Liddick2012 notes support separating detection,parametric accuracy and false pairs; detector conditions and pending-review status are retained in docs/literature-use.md L66. Future work: improve accepted candidates,then independent selection-aware interval calibration/confirmation and native transfer. No commit or push.
'''
    p=EXPERIMENT/'README.md';assert not p.exists();p.write_text(readme)
    save(RECORD/'document-delivery.json',dict(passed=True,files={str(x.relative_to(PROJECT)):digest(x) for x in (PROJECT/'README.md',PROJECT/'docs/literature-use.md',p)}))
    print('SELECTION_DOCUMENTS_DELIVERED',flush=True)


def validate():
    check_freeze()
    for manifest in ('report-freeze.json','delivery-freeze.json','resume-freeze.json'):
        for name,h in load(EXPERIMENT/manifest)['files'].items():assert digest(EXPERIMENT/name)==h,name
    for name,h in load(RECORD/'input-provenance.json')['old_results'].items():assert digest(OLD_OUT/name)==h,name
    tests=load(RECORD/'tests.json');assert tests['passed'] and tests['tests']==11
    stages=[load(EXPERIMENT/'freeze.json')['created']]
    for split,n in (('calibration',256),('test',448)):
        generation=load(RECORD/(split+'-generation.json'));inf=load(RECORD/(split+'-inference.json'))
        assert generation['passed'] and inf['passed'] and generation['frames']==inf['frames']==n
        for name,h in generation['files'].items():assert digest(OUTPUT/split/name)==h
        assert digest(OUTPUT/split/'predictions.jsonl')==inf['sha256']
        assert inf['truth_inputs']==[] and inf['label_access_blocked'] and inf['max_cost_replay_error']<1e-9
        if split=='calibration':stages.extend([generation['created'],inf['created'],load(RECORD/'gate-validation.json')['created']])
        else:
            assert generation['gates_sha256']==digest(EXPERIMENT/'gates.json')
            stages.extend([generation['created'],inf['created']])
    assert stages==sorted(stages)
    a=rows(OUTPUT/'calibration/predictions.jsonl');pilot=rows(OUTPUT/'calibration/pilot.jsonl')
    assert a[:4]==pilot and load(RECORD/'calibration-inference.json')['inherited_pilot']==4
    assert digest(OUTPUT/'calibration/pilot.jsonl')==load(RECORD/'pilot-validation.json')['sha256']
    checkpoint=OUTPUT/'calibration/predictions.interrupted-01.jsonl'
    interruption=load(RECORD/'calibration-interruption.json')
    assert digest(checkpoint)==interruption['checkpoint_sha256']
    assert len(rows(checkpoint))==interruption['preserved_rows']==200
    assert a[:200]==rows(checkpoint)
    recovery=load(RECORD/'calibration-recovery-validation.json')
    assert recovery['passed'] and recovery['retained200_exact_bytes'] and recovery['pilot4_exact_bytes']
    assert digest(OUTPUT/'calibration/predictions.jsonl')==recovery['full_sha256']
    actual=calibrate_gate([r['candidates'] for r in a]);assert actual==gates()['matched_025']
    allseeds=[tuple(r['seed_words']) for split in ('calibration','test') for r in design(split)]
    assert len(set(allseeds))==704
    assert not set(allseeds)&{tuple(r['seed_words']) for r in load(OLD_OUT/'coverage/labels.cases.json')}
    assert load(RECORD/'truth-guard-preflight.json')['passed']
    cutoffs=load(RECORD/'rank-cutoff-validation.json')
    assert cutoffs['passed'] and cutoffs['checks']==96
    science=load(RECORD/'science-validation.json');assert science['passed'] and science['scored']==2688 and science['groups']==42
    assert digest(OUTPUT/'report.json')==science['report_sha256']
    scores=rows(OUTPUT/'scores.jsonl');report=load(OUTPUT/'report.json')
    for s in report['summary']:
        rr=[r for r in scores if r['method']==s['method'] and r['band']==s['band']]
        assert len(rr)==64 and sum(r['success'] for r in rr)==s['success']
        assert sum(r['count']==2 for r in rr)==s['K2_outputs']
    before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4);book=nbformat.read(NOTEBOOK,as_version=4);nbformat.validate(book)
    delivery=load(RECORD/'notebook-delivery.json');assert delivery['passed'] and digest(NOTEBOOK)==delivery['notebook_sha256']
    assert len(before.cells)==543 and book.cells[:543]==before.cells
    new=book.cells[543:];assert len(new)==delivery['new_cells']
    numbers=load(RECORD/'figure-delivery.json')['figures'];assert len(numbers)==6
    visual=load(RECORD/'figure-visual-review.json');assert visual['passed'] and visual['figures']==numbers
    alltext='\n'.join(c.source for c in book.cells)
    newtext='\n'.join(c.source for c in new if c.cell_type=='markdown')
    assert not re.search(r'\]\([^\n)]*\.(?:pdf|ipynb|md)(?:[)#])',newtext,re.I)
    images=[];code_count=0
    for c in new:
        if c.cell_type=='code':
            ast.parse(c.source);code_count+=1;assert c.execution_count and c.outputs
            for out in c.outputs:
                assert out.output_type!='error'
                if 'image/png' in out.get('data',{}):images.append(hashlib.sha256(base64.b64decode(out.data['image/png'])).hexdigest())
    for number,h in numbers.items():
        assert alltext.count('id="codex-figure-'+number+'"')==1
        assert digest(OUTPUT/'figures'/f'{number}.png')==h and images.count(h)==1
    assert code_count==delivery['executed_cells']
    for name,h in load(RECORD/'document-delivery.json')['files'].items():assert digest(PROJECT/name)==h,name
    source_names=('alpha_pair_selection_study.py','test_alpha_pair_selection.py','alpha_pair_selection_report.py','alpha_pair_selection_delivery.py','alpha_pair_selection_resume.py')
    sources={}
    for name in source_names:
        target=EXPERIMENT/('report-source' if name.endswith('_report.py') else 'delivery-source' if name.endswith('_delivery.py') else 'resume-source' if name.endswith('_resume.py') else 'source/pileup/analysis')/name
        assert digest(PROJECT/'pileup/analysis'/name)==digest(target)
        sources['pileup/analysis/'+name]=digest(target)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip();assert head==load(RECORD/'environment.json')['git_head']
    subprocess.run(['git','-c','core.fileMode=false','diff','--check','--','README.md','docs/literature-use.md'],cwd=PROJECT,check=True,timeout=60)
    for p in (NOTEBOOK,PROJECT/'AGENTS.md',EXPERIMENT/'joint-response.npz'):
        assert subprocess.run(['git','check-ignore','-q',str(p.relative_to(PROJECT))],cwd=PROJECT).returncode==0
    save(RECORD/'final-validation.json',dict(passed=True,tests=11,rank_cutoff_checks=cutoffs['checks'],calibration_frames=256,independent_confirmation_frames=448,
        scored=2688,forced_K2_diagnostics=128,old_interval_partitions=36,figures=numbers,new_cells=len(new),executed_cells=code_count,
        original543_cells_preserved=True,old_results_preserved=True,freeze_before_generation=True,gate_before_confirmation=True,
        interrupted200_rows_preserved_exactly=True,
        source_sha256=sources,documentation_sha256=load(RECORD/'document-delivery.json')['files'],notebook_sha256=digest(NOTEBOOK),
        max_cost_replay_error=max(load(RECORD/(s+'-inference.json'))['max_cost_replay_error'] for s in ('calibration','test')),
        no_new_interval_coverage=True,no_nominal_CI_claim=True,no_default_change=True,new_ROOTs=0,new_training=0,no_git_commit_or_push=True))
    print('SELECTION_FINAL_VALIDATION_PASSED',len(new),code_count,flush=True)


def main():
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('notebook','documents','validate'));args=p.parse_args()
    dict(notebook=notebook,documents=documents,validate=validate)[args.command]()


if __name__=='__main__':main()
