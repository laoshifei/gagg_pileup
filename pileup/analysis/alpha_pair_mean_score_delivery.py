"""Executed local Notebook appendix for the frozen mean-score diagnostics."""
from __future__ import annotations

import ast
import hashlib
import json
import os
import re
from pathlib import Path

import nbformat
import numpy as np

from .alpha_pair_mean_score_study import PROJECT,EXPERIMENT,OUTPUT,RECORD,INPUT,NOTEBOOK,NAME,check_freeze
from .alpha_pair_mean_score_report import METHODS,LABELS
from .alpha_pair_information_study import load,save,digest,rows

MARKER='alpha-pair-mean-score-20261028'
LITERATURE='''本轮只读复用共享库的主题／审核索引、数字脉冲处理／PSD／光产额主题页和Fixsen等（2004）Pulse estimation in nonlinear detectors with nonstationary noise既有精读笔记，来源papers/sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md，DOI10.1016/j.nima.2003.11.313。据笔记PDF2式(1)–(2)、PDF3图1–3及第4–5节，研究对象为非线性TES及非平稳噪声，0.5–5eV，每个能量点5000次模型噪声抽样；未给出可移植到本项目的采样率或独立真实探测器pile-up恢复验证。笔记仍为AI草稿待人工审核；符号及RMS/Variance术语疑问保留，不把2026-10-08用户提供正式PDF视为科学采用确认，也不沿用其他论文的审核状态。

采用的是条件均值与完整协方差建模的思路。原文插值均值与权重，并不等于本项目的“内层固定C、外层更新C”算法；本轮IRGLS、阻尼、搜索预算和推进标准均为项目试验设计。全高斯似然在当前模型下仍有正确的统计含义，不能把去掉协方差导数或logdet说成理论纠错，也不宣称普遍无偏／最优。未重读PDF、未联网、未修改共享库。'''


def summary_row(report,method,band):return next(r for r in report['summary'] if r['method']==method and r['band']==band)


def primary_table(report):
    text=['| 间隔/ns | 方法 | 严格整对恢复 | 两时刻匹配 | Q2误差中位 | 未数值收敛 |',
          '| --- | --- | --- | --- | --- | --- |']
    for band in ('pair40','pair200'):
        for m in METHODS:
            r=summary_row(report,m,band)
            text.append(f"| {band[4:]} | {LABELS[m]} | {r['success']}/16 | {r['both_times_strict']}/16 | {100*r['Q2_error_p16_p50_p84'][1]:+.1f}% | {r['not_converged']}/16 |")
    return '\n'.join(text)


def conclusions(report):
    parts=[]
    for band in ('pair40','pair200'):
        s=[summary_row(report,m,band) for m in METHODS]
        parts.append(f"{band[4:]}ns的严格整对恢复，原NLL／扩大搜索NLL／固定一次权重GLS／迭代GLS分别为"+'、'.join(str(r['success'])+'/16' for r in s)+'。所有失败与未收敛样本均保留在分母中。')
        oracle=[summary_row(report,m,band) for m in ('known_times','known_charges','truth_initialized_free')]
        parts.append('同一组在固定真实时刻／固定真实电荷／从全部真实参数启动后自由拟合时，严格成功分别为'+'、'.join(str(r['success'])+'/16' for r in oracle)+'；这三项使用额外真值信息，只是误差来源诊断。')
    search=report['search']
    parts.append(f"扩大搜索使{sum(r['baseline_minus_wide_cost']>.05 for r in search)}/32条的J降低超过0.05；从真值出发仍有{sum(r['wide_minus_truthstart_cost']>.05 for r in search)}/32条找到比扩大搜索更低超过0.05的J。有限起点和有限预算不能证明全局最优；J更低也不等于Q和t更接近真值。0.05只用于描述数值差异，不是统计显著性阈值。")
    counterexamples=[r for r in search if r['truthstart_success'] and not r['wide_success'] and r['wide_minus_truthstart_cost']<-.05]
    parts.append(f"共有{len(counterexamples)}/32条出现：真值启动后的自由拟合严格成功，扩大搜索候选失败，但成功解的J反而高出0.05以上。其中40ns为{sum(r['band']=='pair40' for r in counterexamples)}条、200ns为{sum(r['band']=='pair200' for r in counterexamples)}条。这些反例说明同一代价的排名与本次逐事件误差并不一致；增加搜索起点无法保证选中更准确的那个解，也不能据此判定物理上一定不可分辨。")
    if report['advancement']['chosen'] is None:
        parts.append('没有方法同时通过预先冻结的四项推进条件，因此本轮不增加112条确认样本、不改变默认重建流程。这是44条已暴露模型样本上的开发筛查结论，不是物理不可分辨极限，也没有证明其他优化器无效。')
    else:parts.append('有方法通过开发推进条件；需要另行完成冻结后的独立确认，不将开发成绩当作最终改进。')
    return '\n\n'.join(parts)


def atomic_bytes(path,payload,expected):
    """Readback-verified explicit-offset write followed by atomic replacement."""
    path=Path(path);temporary=path.with_name(path.name+'.mean-score-writing')
    fd=os.open(temporary,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o644);anomalies=[]
    try:
        for offset in range(0,len(payload),256*1024):
            block=payload[offset:offset+256*1024];n=os.pwrite(fd,block,offset)
            if n!=len(block):anomalies.append(dict(offset=offset,expected=len(block),returned=n))
        os.ftruncate(fd,len(payload));os.fsync(fd)
    finally:os.close(fd)
    assert temporary.read_bytes()==payload
    assert digest(path)==expected,'Concurrent edit detected; prepared file preserved'
    os.replace(temporary,path);assert path.read_bytes()==payload
    return anomalies


def notebook():
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    assert load(RECORD/'science-validation.json')['passed'] and load(RECORD/'figure-visual-review.json')['passed']
    report=load(OUTPUT/'report.json');assert report['advancement']['chosen'] is None,'Independent confirmation must be delivered first'
    original=NOTEBOOK.read_bytes();book=nbformat.reads(original.decode(),as_version=4)
    before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4)
    assert len(book.cells)==570 and book.cells==before.cells
    cells=[]
    def md(s):cells.append(nbformat.v4.new_markdown_cell(s))
    def code(s):cells.append(nbformat.v4.new_code_cell(s))
    def fig(number,title):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{title}**')
        code(f'display(Image(filename=str(mean_score_output / "figures/{number}.png")))')
    md('''<a id="alpha-pair-mean-score-20261028"></a>
## 120. 弱脉冲估不准：搜索不足，还是估计本身不稳定？

承接第114–119节，本轮先不放宽选门，而是复用已看过的条件高斯样本定位Q/t误差。从原测试序号320–383、384–447各每四条抽一条，得到40ns、200ns各16条；再取800pe／18800pe孤立脉冲各4条和空窗4条，共44条。双脉冲Q=[18000,800]pe，t1=200.3ns。它们是暴露过的开发数据，不是新独立验证。

本节明确给定真实脉冲数K=0／1／2，比较同一波形上的参数估计；结果不能当作盲计数效率。主要拟合进程禁止打开真值、cases或ROOT文件，仅从已保存的观测K1/K2候选生成起点。空窗K=0不拟合，四条全对是条件K=0的一致性检查，不是误报率测量。

沿用固定95／351ns双衰减、1024×2ns完整模拟电压及同一校准均值/协方差；动态范围充足，无ADC削顶。Q是等效pe，尚不是沉积能量；40／200ns是受控事件间隔，不是源注入率。没有新ROOT、网络训练或本征光学母事件。

预先选定320、384、64、0四条作运行检查，再原样复用到44条主计算。冻结源码、模型、配置、旧输入SHA256和环境后，用两个worker、每个BLAS线程为1运行。首次准备时共享盘备份写入长度异常，未进入拟合；该未完成备份完整保留，随后采用明确偏移分块写入及全文件校验。原Notebook没有因此改动。
''')
    code(f'''from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Image
mean_score_root = Path.cwd().resolve()
while not (mean_score_root / "pileup/experiments/{NAME}").is_dir():
    if mean_score_root.parent == mean_score_root: raise RuntimeError("请从项目根目录或jupyter目录运行")
    mean_score_root = mean_score_root.parent
mean_score_output = mean_score_root / "jupyter/artifacts/{NAME}"
mean_score_report = json.loads((mean_score_output / "report.json").read_text())
print("开发44条；明确已知K；其中真实双脉冲32条；新生成观测0条。")
''')
    md(primary_table(report))
    fig('0004-120-1','已知K的条件恢复：(a,b)40／200ns各16条严格整对恢复；(c)8条孤立控制。严格要求每个Q误差≤20%、每个t误差≤10ns且无额外输出。')
    md('''## 121. 四种估计流程与全样本误差

Baseline为上一轮原K候选。Wider NLL search扩大观测起点搜索，仍最小化完整高斯代价J=rᵀC⁻¹r+logdetC。One frozen-C GLS把原候选处的C固定，最小化加权残差平方。Iterated GLS从这一步出发，每次内层求解固定C，外层才用更新后的Q/t重建C；最多12轮，阻尼依次尝试0.5／0.25／0.125／0.0625，要求本轮固定C的残差代价不增加。

三个新流程使用同一时间参数化和观测种子：K1锚点两侧4／8／16／32／64／128／256／512ns及多种电荷份额，合并由原NLL和固定C各选出的3个分散起点，并包含原候选。时间域为观测K1起点前96ns至后512ns，必要时仅为包含原观测候选扩大，裁剪至0–2046ns；每个Q范围50–50000pe。每次L-BFGS-B最多150迭代，迭代GLS的每轮再使用当前／原候选／宽搜NLL三个起点。内层预算相同，但IRGLS的总预算更大，不能称等计算量比较。

下图保留全部32条，不仅画时间匹配成功或通过选门的样本。Q2/t2按候选时刻排序，可能落到错误事件，因此不是条件匹配后的能量分辨率。完整结果另列匹配Q偏差、sigma68=(p84−p16)/2、时间误差、漏检、额外输出、误拆／误并；sigma68不是FWHM，额外输出是时间找不到匹配真事件，不等于整窗没有信号。
''')
    fig('0004-121-1','四种方法的全部有序后槽Q2/t2误差；红圆40ns、蓝三角200ns，绿色为各自±20%／±10ns范围。双对数式对称刻度保留大误差点；进入绿色交集仍须检查第一个事件。')
    code('''mean_score_metrics = pd.DataFrame(mean_score_report["summary"])
with pd.option_context("display.max_columns", None):
    display(mean_score_metrics[~mean_score_metrics.truth_assisted])
display(pd.DataFrame(mean_score_report["paired"]))
''')
    md(r'''## 122. 降低拟合代价，是否就更接近真值？

对r=y−μ(θ)，完整高斯代价的一个参数方向导数为

$$\partial_\theta J=-2(\partial_\theta\mu)^T C^{-1}r+\mathrm{tr}(C^{-1}\partial_\theta C)-r^T C^{-1}(\partial_\theta C)C^{-1}r.$$

这是本项目对既有模型的推导。在模型正确且满足可微等条件时，后两项在真实参数处的期望抵消；完整NLL并不是错误目标。IRGLS尝试的是只令均值项接近零的估计方程，不是把同一个随θ变化的二次型直接删掉logdet；有限样本偏差和均方误差能否改善，只能靠验证。

解析小例子专门验证了完整NLL、固定权重均值估计和“变化权重但删logdet”三个最优值不同。另有坐标边界、白化二次型与密集矩阵等价、固定协方差一致性、阻尼不能伪造收敛等检查。下面把代价降低与全体Q2误差并列；代价下降本身不作为成功条件。
''')
    fig('0004-122-1','(a)扩大搜索降低J与弱电荷误差的关系；(b)从真实Q/t启动的自由NLL是否找到更低J。后者使用真值，仅定位搜索问题；零线不是置信阈值。')
    code('''display(pd.DataFrame(mean_score_report["search"]))
''')
    md('''## 123. 分别给出真实时刻或真实电荷的诊断

主要预测文件及其SHA256固定以后，另跑40条非空样本的真值辅助诊断。固定真实时刻时，仅从观测Q候选启动；固定真实电荷时，仅从观测时刻启动。每种保留最多6个观测候选起点、每次最多300迭代。第三项把真实Q和t作为起点后全部释放；第四项只计算真实参数的代价，不做拟合。结果写入独立oracle文件，不反馈给主要流程。

固定一个参数组之后仍失败，可以帮助定位剩余的不确定性；给出真值后成功，则表示额外信息有帮助，但不能证明盲算法也能获得这些信息。这些仍是有限优化器的诊断结果，不是给定信息条件下的理论上限。固定真实电荷的Q误差本来为零、固定真实时刻的时间误差本来为零，不能把它们当作学习或优化恢复的能力。真实参数无拟合行只作代价基准，没有画成100%成功的算法。
''')
    fig('0004-123-1','(a,b)40／200ns原NLL、扩大搜索，与三种真值辅助诊断的严格整对成功数。紫色柱使用额外真值信息。')
    code('''with pd.option_context("display.max_columns", None):
    display(mean_score_metrics[mean_score_metrics.truth_assisted & (mean_score_metrics.method != "true_parameters")])
''')
    md('''## 124. 全记录波形与弱分量

固定展示预先选作运行检查的序号320、384，分别属于40、200ns；不是按误差挑出的成功案例。上排给出完整观测、真实独立分量和四种拟合总波形，下排单独画弱分量。横轴均为原始记录0–2046ns，没有局部截取或更短拟合窗口。总波形看起来接近，仍可能把电荷分配给错误分量；需要结合上面逐事件Q/t指标判断。
''')
    fig('0004-124-1','(a,b)两条预先指定波形的完整观测、真实分量和拟合总和；(c,d)同一条的弱分量。真值虚线仅用于评分／展示，均明确已知K=2。')
    md('''## 125. 收敛、推进条件与本轮判断

迭代GLS要求内层优化器成功，且**未阻尼的候选更新**满足全部|ΔlnQ|≤0.002、全部|Δt|≤0.02ns，才记为固定点收敛。不能靠缩小阻尼让实际步长看起来很小而宣称收敛；达到12轮仍不满足时标为未收敛，并照常评分。

主计算前冻结的下一轮推进条件为：40ns严格成功比原候选至少增加2/16；200ns最多下降1/16；32条双脉冲至少80%数值收敛；8条孤立脉冲严格成功不减少。通过者按32条总严格成功数、全体|ln(Q2/真Q2)|中位及固定方法顺序择一。开发筛选不是显著性检验，且数值收敛不等于物理可辨识。若通过，才能冻结后在另112条独立噪声上验证；选门仍必须保持原样。
''')
    fig('0004-125-1','全部32条IRGLS的未阻尼电荷／时刻更新量；红40ns、蓝200ns，虚线为数值收敛条件。点重叠时包含多条样本，单轮也显示标记；仅为对数作图，把小于1e−8的数值放在1e−8。达到预算而未收敛的曲线照常保留；即使更新为零，优化器未成功仍不算收敛。')
    code('''display(pd.DataFrame(mean_score_report["advancement"]["candidates"]))
print("独立确认候选：", mean_score_report["advancement"]["chosen"])
''')
    md(conclusions(report)+'\n\n这批40ns数据中，固定真实时刻后仍只有4/16严格成功，支持继续检查弱电荷在强脉冲尾部噪声下的分配不确定性；它并未证明某个间隔的物理不可恢复性。下一步宜在冻结的完整盲选择流程下，用新独立样本校准Q/t歧义与区间覆盖，并把误差明确标记；若再改点估计，应验证其是否改善这些低代价错误分配。当前没有得到选择后覆盖率结论，原生光学迁移、自由τ、α/γ混合或K≥3也尚未验证。\n\n'+LITERATURE)
    shell=InteractiveShell.instance();shell.user_ns['__name__']='__main__';executed=0
    for cell in cells:
        if cell.cell_type!='code':continue
        executed+=1
        with capture_output(display=True) as captured:result=shell.run_cell(cell.source,store_history=False)
        assert result.error_before_exec is None and result.error_in_exec is None
        outputs=[]
        if captured.stdout:outputs.append(nbformat.v4.new_output('stream',name='stdout',text=captured.stdout))
        if captured.stderr:outputs.append(nbformat.v4.new_output('stream',name='stderr',text=captured.stderr))
        for output in captured.outputs:outputs.append(nbformat.v4.new_output('display_data',data=output.data,metadata=output.metadata))
        cell.outputs=outputs;cell.execution_count=executed
    for cell in cells:cell.metadata['gagg_section_marker']=MARKER
    book.cells.extend(cells);nbformat.validate(book)
    payload=nbformat.writes(book).encode();anomalies=atomic_bytes(NOTEBOOK,payload,hashlib.sha256(original).hexdigest())
    after=nbformat.read(NOTEBOOK,as_version=4);assert after.cells[:570]==before.cells
    save(RECORD/'notebook-validation.json',dict(passed=True,previous_cells=570,added_cells=len(cells),total_cells=len(after.cells),
        executed_cells=executed,old_cells_preserved=True,notebook_sha256=digest(NOTEBOOK),write_count_anomalies=anomalies))
    print('MEAN_SCORE_NOTEBOOK',len(after.cells),executed,flush=True)


def docs():
    report=load(OUTPUT/'report.json')
    paths=(PROJECT/'README.md',PROJECT/'docs/literature-use.md')
    for path,backup in zip(paths,('README.before.md','literature.before.md')):assert digest(path)==digest(RECORD/backup)
    main=paths[0].read_text();main+='\n\n### 2026-10-28：弱脉冲Q/t与均值估计方程诊断\n\n'+conclusions(report)+'\n\n结果为已知K、暴露过的44条条件高斯开发样本；未作为盲效率或物理分辨极限。笔记0004第120–125节含六幅固定编号图，旧570个单元及输出保留。源码及模型快照：`pileup/experiments/'+NAME+'/`；本地结果与运行记录：`jupyter/artifacts/'+NAME+'/`、`jupyter/run_records/'+NAME+'/`。本轮不改变默认流程，无新ROOT、训练、提交或推送。\n'
    literature=paths[1].read_text();assert '## L67' not in literature
    literature+='\n\n## L67 — 2026-10-28：条件协方差与弱脉冲均值估计方程\n\n'+LITERATURE+'\n\n'+conclusions(report)+'\n\n项目位置：`pileup/analysis/alpha_pair_mean_score.py`、`pileup/experiments/'+NAME+'/`；本地Notebook：`jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb`第120–125节。采用层级为方法启发，不是原论文算法复现或物理性能迁移。待办：若继续使用均值估计方程，需另找参数相关协方差下的有限样本偏差／效率、收敛及区间校准一手统计依据；新精读需求在共享库单独处理。本轮没有名义区间覆盖率、原生迁移或新独立确认结论。\n'
    for path,text,backup in zip(paths,(main,literature),('README.before.md','literature.before.md')):
        atomic_bytes(path,text.encode(),digest(RECORD/backup))
    readme=EXPERIMENT/'README.md'
    readme.write_text('# Weak-pulse mean-score diagnostics\n\nKnown-count development on 44 previously exposed frozen-model observations. See protocol.json, freeze.json, input provenance, separate primary/oracle predictions and advancement-decision.json. No default change.\n\n'+conclusions(report)+'\n\nReproduction: run.py tests; run.py pilot --workers 2; run.py infer --workers 2; run.py oracle --workers 2. Existing outputs are protected; use an isolated copy with the recorded relative input tree for a clean rerun. report.py regenerates scoring/figures from saved predictions; delivery.py appends to the preserved 570-cell baseline and checks it before writing. Models/source are archived; switching Git revisions is unnecessary.\n\nThe interrupted pre-inference Notebook backup remains under the sibling *_prepare_interrupted_005415 run record. The original Notebook was unchanged. No new ROOT or random draws were generated.\n')
    save(RECORD/'documentation-validation.json',dict(passed=True,files={str(p.relative_to(PROJECT)):digest(p) for p in (*paths,readme)}))


def validate():
    check_freeze();report=load(OUTPUT/'report.json');assert report['advancement']['chosen'] is None
    for manifest in ('report-freeze.json','delivery-freeze.json'):
        for n,h in load(EXPERIMENT/manifest)['files'].items():assert digest(EXPERIMENT/n)==h,n
    for n,h in load(RECORD/'input-provenance.json')['inputs'].items():assert digest(INPUT/n)==h
    pred=rows(OUTPUT/'predictions.jsonl');pilot=rows(OUTPUT/'pilot.jsonl');indexed={r['index']:r for r in pred}
    assert all(indexed[r['index']]==r for r in pilot)
    book=nbformat.read(NOTEBOOK,as_version=4);before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4)
    assert book.cells[:570]==before.cells and len(before.cells)==570;nbformat.validate(book)
    added=book.cells[570:];anchors=[]
    for cell in added:
        if cell.cell_type=='code':
            ast.parse(cell.source);assert cell.execution_count is not None and all(o.output_type!='error' for o in cell.outputs)
        else:
            assert not re.search(r'\]\([^)]*\.(?:pdf|ipynb|md)(?:[#?][^)]*)?\)',cell.source,re.I)
            anchors+=re.findall(r'id="(codex-figure-[^"]+)"',cell.source)
    assert len(anchors)==6 and len(set(anchors))==6
    images=[o for c in added if c.cell_type=='code' for o in c.outputs if o.output_type=='display_data' and 'image/png' in o.data]
    assert len(images)==6
    full='\n'.join(c.source for c in book.cells)
    assert all(full.count('id="'+a+'"')==1 for a in anchors)
    for n,h in load(RECORD/'figure-manifest.json')['figures'].items():assert digest(OUTPUT/'figures'/(n+'.png'))==h
    numerical=load(RECORD/'tests.json');assert numerical['passed'] and numerical['tests']==11
    for name in ('pilot','inference','oracle','science','notebook','documentation','figure-visual-review'):
        assert load(RECORD/(name+'-validation.json' if name!='figure-visual-review' else name+'.json'))['passed']
    source_names=('alpha_pair_mean_score.py','alpha_pair_mean_score_study.py','test_alpha_pair_mean_score.py',
                  'alpha_pair_mean_score_report.py','alpha_pair_mean_score_delivery.py')
    for n in source_names:assert digest(PROJECT/'pileup/analysis'/n)==digest(EXPERIMENT/'source/pileup/analysis'/n)
    save(RECORD/'final-validation.json',dict(passed=True,notebook_sha256=digest(NOTEBOOK),previous_cells=570,total_cells=len(book.cells),
        source_sha256={str((PROJECT/'pileup/analysis'/n).relative_to(PROJECT)):digest(PROJECT/'pileup/analysis'/n) for n in source_names},
        documentation_sha256=load(RECORD/'documentation-validation.json')['files'],primary_frames=44,oracle_frames=40,
        primary_scores=176,oracle_scores=160,tests=11,figures=6,pilot_reused=4,inputs_unchanged=True,
        exposed_development=True,known_count=True,new_ROOTs=0,new_observations=0,new_training=0,no_default_change=True,no_commit_or_push=True,
        report_sha256=digest(OUTPUT/'report.json'),advancement_sha256=digest(EXPERIMENT/'advancement-decision.json')))
    print('MEAN_SCORE_FINAL_VALIDATION_PASSED',flush=True)


if __name__=='__main__':notebook();docs();validate()
