"""Preserve the existing Notebook and append executed independent diagnostics."""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import re
from pathlib import Path

import nbformat

from .alpha_pair_blind_uncertainty_study import (
    PROJECT, EXPERIMENT, OUTPUT, RECORD, NOTEBOOK, NAME,
    load, save, digest, write_new, check_freeze, now,
)
from .alpha_pair_blind_uncertainty_report import conclusions

MARKER='alpha-pair-blind-uncertainty-20261028'
LITERATURE='''本轮承接第126节及文献采用记录L68，复用Haselman等（2012，DOI10.1109/TNS.2012.2207403，笔记定位PDF3–7图2–11/表I）的能量、定时、计数分开评价；Žugec等（2016，DOI10.1016/j.nima.2015.12.054，PDF7–9式13–23/图9–10）的模板拟合质量记录；Usman、Patil（2018，DOI10.1016/j.net.2018.06.014，PDF1–3/5–9）的计数补偿与事件恢复区分。文献的LFS-3/MAPD-N合成流、n_TOF多类探测器和综述模型条件均不移植为GAGG数值保证。沿用既有AI待审状态；本轮没有重新检索或阅读PDF，没有修改共享库。

联合代价集合、经验上分位标定、八个容差外约束搜索及数值预算是项目试验设计，不是上述论文算法的复现。三篇均未提供本项目选择后、有限样本、边界/多解条件下的置信区间保证。相应一手统计证据仍列为共享库独立整理待办。'''


def freeze_extra():
    check_freeze(); folder=EXPERIMENT/'report-source'; folder.mkdir(exist_ok=False)
    for name in ('alpha_pair_blind_uncertainty_report.py','alpha_pair_blind_uncertainty_delivery.py'):
        write_new(folder/name,(PROJECT/'pileup/analysis'/name).read_bytes())
    files={str(p.relative_to(EXPERIMENT)):digest(p) for p in folder.glob('*.py')}
    for script,module in [('report.py','alpha_pair_blind_uncertainty_report'),('deliver.py','alpha_pair_blind_uncertainty_delivery')]:
        wrapper=f'''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
os.environ['GAGG_UNCERTAINTY_PROJECT']=str(folder.parents[2])
sys.path.insert(0,str(folder/'source'))
import pileup.analysis
pileup.analysis.__path__.insert(0,str(folder/'report-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_blind_uncertainty_delivery import check_extra
    check_extra()
    from pileup.analysis.{module} import main
    main()
'''
        write_new(EXPERIMENT/script,wrapper.encode()); files[script]=digest(EXPERIMENT/script)
    save(EXPERIMENT/'report-freeze.json',dict(created=now(),files=files))
    print('REPORT_SOURCE_FROZEN',len(files),flush=True)


def check_extra():
    for name,h in load(EXPERIMENT/'report-freeze.json')['files'].items(): assert digest(EXPERIMENT/name)==h,name


def replace_verified(path,payload,expected):
    path=Path(path); temporary=path.with_name(path.name+'.blind-uncertainty-writing')
    write_new(temporary,payload)
    assert digest(path)==expected, f'Concurrent change; prepared file retained: {path}'
    os.replace(temporary,path)
    assert path.read_bytes()==payload


def table(report):
    lines=['| 间隔/ns | 盲选双事件 | 严格整对恢复 | 选双后联合包含 | 全体中选双且包含 | 时间漏检/128 | 额外输出 |',
           '| --- | --- | --- | --- | --- | --- | --- |']
    for r in report['primary_containment']:
        s=next(x for x in report['summary'] if x['band']==r['band'])
        c=r['conditional_selected']; a=r['selected_and_contained']
        lines.append(f"| {r['band'][4:]} | {s['selected_K2']['numerator']}/{s['frames']} | {s['strict']['numerator']}/{s['frames']} | {c['numerator']}/{c['denominator']} | {a['numerator']}/{a['denominator']} | {s['missed']}/128 | {s['extra']} |")
    return '\n'.join(lines)


def notebook():
    check_freeze(); report=load(OUTPUT/'report.json')
    assert load(RECORD/'science-validation.json')['passed']
    figures=load(RECORD/'figures-validation.json'); assert figures['passed']
    old=NOTEBOOK.read_bytes(); before=(RECORD/'notebook.before.ipynb').read_bytes(); assert old==before
    book=nbformat.reads(old.decode(),as_version=4); assert len(book.cells)==600
    assert MARKER not in old.decode(); cells=[]
    def md(s): cells.append(nbformat.v4.new_markdown_cell(s))
    def code(s): cells.append(nbformat.v4.new_code_cell(s))
    def fig(number,caption):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{caption}**')
        code(f'display(Image(filename=str(uncertainty_output / "figures/{number}.png")))')
    md('''<a id="alpha-pair-blind-uncertainty-20261028"></a>
## 127. 独立盲重建验证：选出两个事件以后，参数有多可信？

承接第126节，本轮完成固定模型下的最小独立可靠性研究。原waveform_candidates及原生控制校准的0.025事件数门完全保留；新增诊断不改电荷、时刻和事件数输出。第一轮快照漏带analysis包，导入测试即失败，未生成数据；旧准备目录和备份保留，v2补全依赖后才运行。

校准144条：40/200ns、Q=[18000,800]pe、t1=200.3ns的先强后弱双脉冲各32条，另空窗及800/3000/8000/18800pe孤立控制各16条。确认288条：相同两档双脉冲各64条，五档控制各32条。两组噪声种子独立，确认波形在阈值与源码冻结后才生成；校准先导9条原样复用于校准总样本。

输入仍是固定95/351ns双衰减、1024×2ns完整float64电压、充分动态范围。均值与相关协方差完全沿用既有模型，在FIR空间抽样后反变换，前三个残差为零。这是条件高斯模型内独立验证，不是新光学母事件或原生SiPM实验；Q为等效pe，未标定为沉积能量，40/200ns不是源率。

推断仅接触观测和冻结模型，worker实际测试并阻断truth/cases/ROOT访问。真值只用于生成以及预测文件哈希冻结后的校准和评分；真实K不进入盲重建。失败、未收敛、漏检和额外输出保留在相应分母中。
''')
    code(f'''from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Image
uncertainty_root = Path.cwd().resolve()
while not (uncertainty_root / "pileup/experiments/{NAME}").is_dir():
    if uncertainty_root.parent == uncertainty_root: raise RuntimeError("请从项目根目录或jupyter目录运行")
    uncertainty_root = uncertainty_root.parent
uncertainty_output = uncertainty_root / "jupyter/artifacts/{NAME}"
uncertainty_report = json.loads((uncertainty_output / "report.json").read_text())
print("校准144条；独立确认288条；未知事件数；新ROOT和训练均为0。")
display(pd.DataFrame(uncertainty_report["summary"]))
''')
    md(table(report))
    fig('0004-127-1','(a)40ns、(b)200ns。输出双事件、两事件Q误差均≤20%且t误差均≤10ns的严格恢复、以及输出联合区域包含真参数，三项以同组全部64条真实双脉冲作分母。')
    md(r'''## 128. 联合代价区域与独立包含率

本轮不把单参数Q2曲线直接解释成四参数置信区间。为同时检查Q1、Q2、t1、t2，定义

$$R_c(y)=\{\theta\in D(y):J(y;\theta)-J(y;\widehat\theta_{K2})\le c\}.$$

J为包含logdet的完整高斯代价；参考点是冻结的原K2候选，不声称全局最优。D使用观测K1时刻前96ns至后512ns，并仅为包含原K2候选扩展、裁至0–2046ns；两时刻排序，每个Q限制50–50000pe。只有盲选K2后才实际输出这一条件模型的区域；其他事件数照原判决。

预测冻结后，在校准集中被选为K2的真实双脉冲上计算真参数相对参考代价。取第ceil(0.95n)个排序分数并与0取大，至少要求20条且阈值有限；真参数落在D外记为无穷大，不丢弃。这里的95%是经验标定目标，不是覆盖率定理，原参考并非全局最小值，不能套χ²阈值。校准混合权重由固定分层设计和选门共同形成，不是任意源谱的保证。

冻结c后才生成确认波形。报告：(1)在盲选K2子集中的联合包含；(2)全体真实双脉冲中选出K2且联合包含；(3)原点估计严格恢复。另列未选中K2的候选区域仅作已知真实双事件的评分诊断，不算实际输出效率。

图中的精确二项95%区间只适用于各固定参数组的重复抽样比例，条件于已经冻结的校准，不含校准不确定性和物理模型误差。跨固定分层的合并比例只给计数，不冒称同分布二项精确区间。当前联合区域的1/4/9对照与旧Q2剖面的参数维数、参考代价不同，不能直接比较数值或沿用旧置信解释。
''')
    code('''display(pd.DataFrame([uncertainty_report["calibration"]]))
display(pd.DataFrame(uncertainty_report["primary_containment"]))
display(pd.DataFrame(uncertainty_report["level_sensitivity"]))
''')
    fig('0004-128-1','(a)40ns、(b)200ns。蓝线为选双后的联合包含，橙线以全体真实双脉冲作分母；竖线为校准后冻结的c，横虚线只是标定目标。误差条描述固定模型重复抽样，不证明名义95%覆盖。')
    md(r'''## 129. 同样能解释波形的超差替代解

对盲选出的双事件，搜索八个达到或超出容差边界的参数区域：两个时刻各自低于/高于拟合值10ns，以及两个Q各自越过允许范围。由于严格恢复的相对误差分母是真实Q，对可能真实Q的允许范围应为

$$[\widehat Q/1.2,\ \widehat Q/0.8],$$

不能误写为拟合Q的±20%。八个搜索区均与D及t1≤t2相交；含等号的边界使诊断保守，边界本身仍可能恰好满足严格容差。

起点由原K2及观测K1的电荷/时刻生成，含8/40/200ns试探间隔、强弱两种顺序和五种份额；这些是固定搜索配置，不读取真实间隔。先每区2起点、最多100迭代，再每区4起点、最多180迭代，并复用第一阶段的解。保留所有已得更优点；最小外部代价变化≤0.05称本次搜索稳定，另记录各区最优点收敛、参数域限制和实际代价重放。它不是全局最优或完整拓扑证明。

若找到实际外部参数满足J−Jref≤c，记为“找到替代解”。若没找到，但搜索不稳定、存在未收敛最优点或容差区撞到参数域边界，记为“数值未定”。其余只记“有限搜索未找到替代解”，不能写成已证明可靠或唯一。辅助搜索发生数值异常时保留原盲拟合，标为数值未定。

本轮完成四维隐式联合区域的包含检验与有限容差外搜索；尚未计算各参数完整的投影区间、宽度和分离拓扑，不用八个外部最优点冒充完整区间。
''')
    fig('0004-129-1','(a)40ns、(b)200ns的盲选双事件。横轴≤0表示存在区域内的容差外替代点；纵轴为全部Q/t误差中最大的归一化值，≤1才满足参数容差。时间未匹配案例仍保留。')
    code('''display(pd.DataFrame([uncertainty_report["numerical"]]))
''')
    md('''## 130. 完整记录中的原波形、盲分解和替代分解

按生成前指定的确认序号0、1、64、65展示，不按误差挑成功案例。上排完整0–2046ns波形，下排各分量；真值虚线只用于展示与评分。原门未选K2时不强行补出第二个事件，也没有该帧的外部搜索。替代参数只用于诊断，不替换默认输出。
''')
    for j,index in enumerate((0,1,64,65),1):
        fig(f'0004-130-{j}',f'预定确认序号{index}的完整记录：(a)观测及总波形；(b)真实、盲拟合及可用时的替代分量。原门未选双时按实际输出展示。')
    md('''## 131. 质量标记的接受率、准确率与损失

把难事件标为有歧义，不会自动增加正确恢复事件。下面同时报告每组的接受数、严格正确数、控制样本中的假双，以及按全体真实双脉冲计算的正确产出。零接受对应0/0，准确率未定义，不能记为0%或100%。类型不同的固定分层合并结果只报计数；逐参数组的抽样区间另列。

匹配误差分位数只覆盖按时间匹配的事件，必须与漏检和额外输出一同解释。sigma68=(p84−p16)/2，不是能谱FWHM；误合并/误拆分沿用既有时间邻近定义，不被概括成所有失败的物理原因。
''')
    fig('0004-131-1','(a)所有盲选K2（含控制中的假双）的诊断组数量；(b)各组严格正确比例。混合分层不加同分布二项误差条；未找到替代解也不自动代表高准确率。')
    code('''display(pd.DataFrame(uncertainty_report["quality"]))
display(pd.DataFrame(uncertainty_report["quality_by_band"]))
display(pd.DataFrame([uncertainty_report["control_summary"]]))
''')
    md('## 132. 本轮结论、适用范围与复现\n\n'+conclusions(report)+'\n\n'+LITERATURE+f'''

实际源码、模型、原事件数门、协议及随机种子冻结在pileup/experiments/{NAME}/；观察、预测、评分、CSV和PNG在jupyter/artifacts/{NAME}/；测试、环境、SHA256及旧笔记备份在jupyter/run_records/{NAME}/。每个阶段保护已有输出，重跑用独立项目副本和所记录的本地模型输入，不依赖切换Git版本。

下一步应根据本轮“被选中但有歧义”的比例，决定先扩充参数点/窗口与相位，还是先补可实现的独立幅度或响应约束。各参数完整区间、原生SiPM母事件、α/γ、free-τ及连续K≥3尚未由本轮验证；不把40ns困难推广为普适物理极限。
''')
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    shell=InteractiveShell.instance(); shell.user_ns.clear(); shell.user_ns['__name__']='__main__'; executed=0
    for i,cell in enumerate(cells):
        cell.id=f'blind-uncertainty-{i:03d}'; cell.metadata['gagg_section_marker']=MARKER
        if cell.cell_type=='markdown':
            assert not re.search(r'\]\([^)]*\.(?:md|ipynb|pdf)(?:[#?][^)]*)?\)',cell.source,re.I)
        else:
            executed+=1
            with capture_output(display=True) as captured: result=shell.run_cell(cell.source,store_history=False)
            assert result.error_before_exec is None and result.error_in_exec is None
            outputs=[]
            if captured.stdout: outputs.append(nbformat.v4.new_output('stream',name='stdout',text=captured.stdout))
            if captured.stderr: outputs.append(nbformat.v4.new_output('stream',name='stderr',text=captured.stderr))
            for output in captured.outputs: outputs.append(nbformat.v4.new_output('display_data',data=output.data,metadata=output.metadata))
            cell.outputs=outputs; cell.execution_count=executed
    book.cells.extend(cells); nbformat.validate(book)
    payload=nbformat.writes(book).encode(); assert json.loads(payload)['cells'][:600]==json.loads(before)['cells']
    replace_verified(NOTEBOOK,payload,hashlib.sha256(old).hexdigest())
    after=nbformat.read(NOTEBOOK,as_version=4); assert after.cells[:600]==nbformat.reads(before.decode(),as_version=4).cells
    save(RECORD/'notebook-validation.json',dict(passed=True,old_cells=600,added_cells=len(cells),total_cells=len(after.cells),
        executed_cells=executed,figures=len(figures['figures']),old_cells_outputs_preserved=True,notebook_sha256=digest(NOTEBOOK)))
    print('NOTEBOOK_UPDATED',len(after.cells),executed,flush=True)


def documents():
    report=load(OUTPUT/'report.json'); text=conclusions(report)
    for path,backup in ((PROJECT/'README.md',RECORD/'README.before.md'),(PROJECT/'docs/literature-use.md',RECORD/'literature.before.md')):
        assert digest(path)==digest(backup)
    readme=(RECORD/'README.before.md').read_text()+f'\n\n### 2026-10-28：盲重建的联合参数可靠性独立验证\n\n{text}\n\n校准144条、独立确认288条；原NLL和事件数门不变，使用完整四参数代价区域及有限容差外搜索。模型条件、固定分层及数值限制详见 [实验记录](pileup/experiments/{NAME}/README.md)。本地0004笔记第127–132节保留原600单元，并给出八幅固定图号的新结果。没有新ROOT、训练或默认算法替换。\n'
    literature=(RECORD/'literature.before.md').read_text(); assert not re.search(r'^#{2,3}\s+L69\b',literature,re.M)
    literature+='\n\n## L69 — 2026-10-28：盲选后的联合参数代价区域与独立歧义验证\n\n'+LITERATURE+'\n\n'+text+f'''

**方法及边界**：参考代价为原K2候选，不声称全局最小值。校准144条，在盲选K2的真实双事件上取ceil(0.95n)排序分数并与0取大，域外真值保留为无穷大，至少20个有效选中对且阈值有限才推进；冻结后生成独立288条。分别报告选双后的四参数联合包含及全体真实双脉冲中的输出且包含。95%为经验目标，没有名义覆盖或未知参数统一保证；不把旧Q2剖面阈值等同当前四维区域。逐固定参数组给精确二项抽样区间，固定分层合并只报计数，区间不含校准及物理模型误差。

**歧义标记**：以Qhat/1.2、Qhat/0.8及时间±10ns定义八个闭合容差外约束，保持时刻排序和原参数域。观测起点2次/100迭代后4次/180迭代复核，保留旧更优点，最小外部代价变化≤0.05为有限数值稳定。实际找到区域内替代点说明存在歧义；没找到不证明全局唯一，收敛/边界不满足时数值未定。原输出不改变；辅助诊断失败仍保留原拟合。没有生成各参数完整投影区间、宽度和分离拓扑。所有失败、漏检、额外输出和质量筛选损失保留；Q是等效pe，sigma68不是能谱FWHM。

**复现与验证**：固定95/351ns、1024×2ns完整电压和充足动态范围；条件高斯FIR首三残差为零。SeedSequence[202610280701,split,index,997]，校准/确认分别split1/2。实际源码、模型、门、协议、运行器快照：pileup/experiments/{NAME}/；本地Notebook jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb第127–132节；本地运行记录jupyter/run_records/{NAME}/。首轮依赖快照导入失败发生在数据生成前，旧目录完整保留，v2补全absolute analysis依赖。未修改旧生产数据或共享库，未commit/push。

**待办**：各参数完整投影、扩大幅比/相位/窗口、原生独立母事件及模型失配、α/γ、free-τ、连续K≥3和计数率仍需分阶段验证。经验标定和有限搜索不解决已有一手统计依据缺口，也不能证明物理不可分辨。
'''
    for path,content,backup in ((PROJECT/'README.md',readme,RECORD/'README.before.md'),
                              (PROJECT/'docs/literature-use.md',literature,RECORD/'literature.before.md')):
        replace_verified(path,content.encode(),digest(backup))
    experiment_readme=f'''# Independent blind joint-parameter compatibility study

{text}

{table(report)}

This is a frozen-model diagnostic, not a new reconstruction algorithm or a theorem of nominal confidence coverage. The original blind candidate search and native count gate alpha0.025 are unchanged. Joint cost regions use the original K2 candidate as reference, not a claimed global optimum. Eight constrained exterior searches look for actual alternative Q/t sets at or beyond the20%/10ns accuracy boundary; absent finite-search witnesses do not prove containment or identifiability.

Calibration144 frames precedes independently generated confirmation288 frames. See protocol.json for the fixed stratified designs, SeedSequence words, complete numerical settings and stopping rules. The selected-pair empirical95% quantile is frozen before confirmation, with domains retained and no test tuning. Per-gap exact binomial sampling intervals condition on this frozen calibration; pooled fixed strata do not receive iid-binomial intervals. Truth is blocked in inference workers and used only for simulation and downstream calibration/scoring. Numerical failures, missed/extra/split/merge and quality-acceptance losses remain explicit.

Purealpha,fixed95/351ns,full1024x2ns float64 voltage,ample dynamic range,localK<=2. No new ROOT,optical mothers,training,free-tau,mixedalpha/gamma or rate-scan claims. Q is equivalentpe,not depositedenergy. No marginal profile interval widths or complete topology are computed. Source/model/protocol snapshots avoid a required Git checkout; local ignored arrays/models/Notebook must be supplied separately.

From the project root, using the project environment and BLAS1:

```bash
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py tests
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py generate-calibration
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py pilot --workers 2
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py infer-calibration --workers 4
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py calibrate
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py generate-test
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py infer-test --workers 4
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/report.py analyze
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/report.py figures
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/deliver.py validate
```

Existing output paths are protected. For clean replay use an isolated project copy and the recorded input assets; the pilot's nine rows are reused unchanged. Primary and report/delivery sources have separate frozen manifests. The sibling preparation without_v2 failed to import an omitted absolute analysis dependency before tests/generation and is preserved as an unsuccessful preparation, not a result. Local Notebook0004 sections127–132 append to600 preserved cells and include eight explicitly numbered figures. No commit or push.
'''
    write_new(EXPERIMENT/'README.md',experiment_readme.encode())
    save(RECORD/'documentation-validation.json',dict(passed=True,files={str(p.relative_to(PROJECT)):digest(p) for p in
         (PROJECT/'README.md',PROJECT/'docs/literature-use.md',EXPERIMENT/'README.md')}))
    print('DOCUMENTS_UPDATED',flush=True)


def validate():
    check_freeze(); check_extra()
    checks={name:load(RECORD/(name+'.json')) for name in ('tests','pilot-validation','calibration-generation',
        'calibration-inference','calibration-freeze','test-generation','test-inference','science-validation',
        'figures-validation','notebook-validation','documentation-validation')}
    assert all(r['passed'] for r in checks.values())
    assert checks['calibration-inference']['frames']==144 and checks['test-inference']['frames']==288
    for split in ('calibration','test'):
        assert digest(OUTPUT/split/'predictions.jsonl')==checks[split+'-inference']['sha256']
        for name,h in checks[split+'-generation']['files'].items(): assert digest(OUTPUT/split/name)==h
    assert digest(EXPERIMENT/'calibration.json')==checks['calibration-freeze']['sha256']==checks['test-generation']['calibration_sha256']
    assert digest(OUTPUT/'report.json')==checks['science-validation']['report_sha256']
    book=nbformat.read(NOTEBOOK,as_version=4); before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4)
    nbformat.validate(book); assert len(before.cells)==600 and book.cells[:600]==before.cells
    ids=[c.id for c in book.cells]; assert len(ids)==len(set(ids))
    current=book.cells[600:]; codes=[c for c in current if c.cell_type=='code']
    assert all(c.execution_count is not None and not any(o.output_type=='error' for o in c.outputs) for c in codes)
    pngs=sum('image/png' in output.get('data',{}) for cell in codes for output in cell.outputs)
    assert pngs==8
    for name,h in checks['figures-validation']['files'].items(): assert digest(OUTPUT/'figures'/name)==h
    for name,h in checks['documentation-validation']['files'].items(): assert digest(PROJECT/name)==h
    report=load(OUTPUT/'report.json'); assert report['baseline_unchanged']
    for path in (EXPERIMENT/'source').rglob('*.py'): ast.parse(path.read_text())
    for path in (EXPERIMENT/'report-source').glob('*.py'): ast.parse(path.read_text())
    for path in EXPERIMENT.glob('*.json'): load(path)
    for name,h in load(RECORD/'environment.json')['inputs'].items(): assert digest(PROJECT/name)==h
    result=dict(passed=True,created=now(),calibration_frames=144,confirmation_frames=288,confirmation_pairs=128,
        tests=checks['tests']['tests'],figures=pngs,previous_cells=600,total_cells=len(book.cells),
        executed_new_cells=len(codes),notebook_sha256=digest(NOTEBOOK),report_sha256=digest(OUTPUT/'report.json'),
        source_freeze_sha256=digest(EXPERIMENT/'freeze.json'),report_freeze_sha256=digest(EXPERIMENT/'report-freeze.json'),
        independent_splits=True,original_blind_outputs_unchanged=True,all_frames_retained=True,
        old_notebook_cells_preserved=True,new_ROOTs=0,new_training=0,new_native_mothers=0,
        no_default_change=True,no_commit_or_push=True)
    save(RECORD/'final-validation.json',result); print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['freeze-extra','notebook','documents','validate'])
    args=parser.parse_args()
    {'freeze-extra':freeze_extra,'notebook':notebook,'documents':documents,'validate':validate}[args.command]()


if __name__=='__main__': main()
