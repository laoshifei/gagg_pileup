"""Append executed stability diagnostics to local Notebook 0004."""
from __future__ import annotations

import ast
import base64
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import nbformat
import numpy as np

from .alpha_pair_stability_study import PROJECT,EXPERIMENT,OUTPUT,RECORD,INPUT,OLD_RECORD,PREVIOUS,NOTEBOOK,NAME,check_freeze
from .alpha_pair_information_study import load,save,digest,rows

MARKER='alpha-pair-profile-stability-20261027'
GROUPS=(('native','40'),('native','200'),('continuous','middle'),('continuous','long'))
LITERATURE='''本轮通过paper-search只读复用共享库的主题、审核索引及两篇完整精读笔记，没有重新读取PDF或联网。Baker与Cousins（1984），Clarification of the use of chi-square and likelihood functions in fits to histograms，NIM221,437–442，已有笔记定位PDF1、PDF4§5–7、PDF5§8：点估计、区间覆盖及拟合优度应分开，有限样本需按具体模型研究抽样分布。其对象为Poisson／multinomial直方图，不是本项目的相关电压；不移用其计数似然或渐近覆盖结论。来源papers/sources/1984_Baker_NIM_Clarification-of-the-use-of-chi-square-and-likelihood-functions-in-fits-to-histograms.md，仍为AI草稿待人工审核，暂无科学采用确认；DOI未由已有笔记核实，式3概率记号疑问保留。

Liddick、Darby与Grzywacz（2012），Algorithms for pulse shape analysis using silicon detectors，NIMA669,70–78，DOI10.1016/j.nima.2011.12.034，已有笔记定位PDF3–8／图3–13，约65μm DSSD、40MSPS／12bit、25μs记录，8000余人工叠加样本。仅沿用检测、参数精度与假双控制分开报告的思路；不移植100ns、幅比5或ADC限制。来源papers/sources/2012_Liddick_NIMA_Algorithms-for-pulse-shape-analysis-using-silicon-detectors.md，仍为AI草稿待审，门选记号及自由参数数目缺口保留。增加起点、网格加密及本轮0.05／5%数值筛查门槛是项目预定诊断设计，非这两篇论文提供的定理或性能保证。'''


def table(report):
    lines=['| 条件 | 最后加密三层级均稳定 | ΔJ≤4宽度中位：旧→新／真实Q2 | 多区段：旧→新 | 已知K2严格正确：旧→新 |',
           '| --- | ---: | ---: | ---: | ---: |']
    for domain,band in GROUPS:
        a=next(r for r in report['profile_summary'] if (r['domain'],r['band'],r['stage'],r['level'])==(domain,band,'original',4.))
        b=next(r for r in report['profile_summary'] if (r['domain'],r['band'],r['stage'],r['level'])==(domain,band,'contour',4.))
        stable=sum(r['final_stable_all_levels'] for r in report['frames'] if (r['domain'],r['band'])==(domain,band))
        lines.append(f"| {domain}／{band} | {stable}/8 | {a['width_over_trueQ'][1]:.3f}→{b['width_over_trueQ'][1]:.3f} | {a['disconnected']}→{b['disconnected']}/8 | {a['conditional_best_success']}→{b['conditional_best_success']}/8 |")
    return '\n'.join(lines)


def interpretation(report):
    """Post-scoring consistency checks; never used to choose or refit profiles."""
    curves={(r['domain'],r['index']):r for r in rows(OUTPUT/'profiles.jsonl')}
    changes=load(OUTPUT/'changes.json');scores=load(OUTPUT/'profile-scores.json')
    old=load(INPUT/'report.json')['profiles'];context=[];consistency=[]
    for domain,band in GROUPS:
        prior=next(r for r in old if (r['domain'],r['band'],r['level'])==(domain,band,4.))
        selected=next(r for r in report['profile_summary'] if (r['domain'],r['band'],r['stage'],r['level'])==(domain,band,'original',4.))
        final=next(r for r in report['profile_summary'] if (r['domain'],r['band'],r['stage'],r['level'])==(domain,band,'contour',4.))
        context.append(dict(domain=domain,band=band,prior_all_frames=prior['frames'],selected_frames=selected['frames'],
            prior_all_width_median=prior['width_over_trueQ'][1],selected_original_width_median=selected['width_over_trueQ'][1],
            selected_final_width_median=final['width_over_trueQ'][1],prior_all_strict_success=prior['conditional_best_success'],
            selected_original_strict_success=selected['conditional_best_success'],selected_final_strict_success=final['conditional_best_success']))
    for s in scores:
        if s['stage']!='contour':continue
        c=curves[s['domain'],s['index']]['stages']['contour'];q=s['truth_Q2_pe']
        delta=float(np.interp(np.log(q),np.log(c['q_grid_pe']),[f['cost']-c['reference_cost'] for f in c['fits']]))
        grid_included=any(lo<=q<=hi for lo,hi in s['segments'])
        consistency.append(dict(domain=s['domain'],index=s['index'],level=s['level'],
            interpolated_delta=delta,probe_delta=s['truth_delta'],grid_minus_probe=delta-s['truth_delta'],
            grid_included=grid_included,probe_included=s['truth_included'],agrees=grid_included==s['truth_included']))
    failures=[dict(domain=c['domain'],index=c['index'],level=c['level'],
        endpoint_fraction_change=float(np.expm1(c['endpoint_log_change'])) if c['endpoint_log_change'] is not None else None,
        same_topology=c['same_topology'],same_boundaries=c['same_boundaries'],reference_gain=c['reference_gain'])
        for c in changes if c['comparison']=='dense_contour' and not c['stable']]
    trials=[c for r in curves.values() for stage,c in r['stages'].items() if stage!='original']
    files={'output/'+n:digest(OUTPUT/n) for n in ('profiles.jsonl','profile-scores.json','changes.json','report.json')}
    files['prior/report.json']=digest(INPUT/'report.json')
    return dict(passed=True,post_scoring_only=True,primary_changed=False,sample_context=context,
        final_stable_by_level={str(int(level)):sum(c['stable'] for c in changes if c['comparison']=='dense_contour' and c['level']==level) for level in (1.,4.,9.)},
        final_failures=failures,additional_optimizer_trials=sum(c['trials'] for c in trials),
        nonconverged_optimizer_trials=sum(c['nonconverged_trials'] for c in trials),
        new_nonconverged_truth_probes=sum(not f['new_probe_converged'] for f in report['frames']),
        selected_nonconverged_truth_probes=sum(not f['probe_converged'] for f in report['frames']),
        interpolation_probe_checks=consistency,inclusion_disagreements=sum(not r['agrees'] for r in consistency),
        maximum_absolute_grid_probe_delta_difference=max(abs(r['grid_minus_probe']) for r in consistency),input_sha256=files)


def interpretation_text(a):
    native=next(r for r in a['sample_context'] if (r['domain'],r['band'])==('native','40'))
    failed='、'.join(f"{r['domain']}第{r['index']}条（ΔJ≤{r['level']:g}，端点变化{100*r['endpoint_fraction_change']:.2f}%）" for r in a['final_failures'])
    return f'''末次加密在ΔJ≤1／4分别为{a['final_stable_by_level']['1']}/32和{a['final_stable_by_level']['4']}/32通过；未通过的是{failed}。这两条的区段数、边界标志及参考代价不变，只是端点变化超过预定5%。

**必须按同一批样本比较**：原生40ns组本轮8条的范围跨度／真实Q2中位数为{native['selected_original_width_median']:.2f}→{native['selected_final_width_median']:.2f}；上一轮全部40条的旧中位数为{native['prior_all_width_median']:.2f}。前两者才是加密前后比较，不能说从{native['prior_all_width_median']:.2f}恶化到{native['selected_final_width_median']:.2f}。系统抽取8条不保证统计代表性；跨度也包含分离区段之间的空隙。四组已知K2的严格恢复数仍分别为1／1／1／6（各8条），增加搜索没有改善这批样本的逐事件准确度。

本轮额外{a['additional_optimizer_trials']}次局部优化中，{a['nonconverged_optimizer_trials']}次保留未收敛标志；这不是失败波形条数。新真实Q2评分探针有{a['new_nonconverged_truth_probes']}条未收敛，取旧／新较低代价后实际选用的探针有{a['selected_nonconverged_truth_probes']}条未收敛。有限搜索仍可能遗漏更低解。

另做了**评分后的一致性检查**：图中区段来自网格的lnQ插值，真Q包含计数来自实际拟合探针，两者不能定义混用。本轮32条×3层级共96个包含判定有{a['inclusion_disagreements']}个不一致；插值ΔJ与探针ΔJ的最大绝对差为{a['maximum_absolute_grid_probe_delta_difference']:.3f}。包含判定相同不等于整条连续曲线已收敛。这项检查不调整主曲线，也不是新的独立覆盖率检验。'''


def conclusion(report):
    passed4=sum(r['stable'] for r in report['stability_summary'] if r['comparison']=='dense_contour' and r['level']==4.)
    text=f"32条预定波形中，{passed4}条在重点ΔJ≤4层级通过末次加密的数值筛查；{report['stable_all_levels']}条在ΔJ≤1／4／9三层级同时通过。增加搜索和加密网格没有改善这批样本的已知K2严格恢复数，原生40ns组3/8条的分离区段仍存在。"
    text+=f"{report['probe_improves_reference']}条真实Q2评分探针把最终参考代价降低超过0.01；该探针只用于评分，不反馈主剖面。"
    if report['stable_all_levels']<32 or report['probe_improves_reference']:
        text+='仍有数值变化或搜索遗漏，不能把当前范围当作已稳定的置信区间。本轮完成数值审核，独立覆盖率校准留待处理这些变化后开展。'
    else:
        text+='通过这个有限数值筛查仍不等于全局收敛，也不提供覆盖保证；下一阶段需冻结区间算法后独立重复抽样校准。'
    return text


def notebook():
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    assert load(RECORD/'science-validation.json')['passed'] and load(RECORD/'figure-visual-review.json')['passed']
    original=NOTEBOOK.read_bytes();book=nbformat.reads(original.decode(),as_version=4)
    before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4)
    assert len(book.cells)==len(before.cells)==497 and book.cells==before.cells
    report=load(OUTPUT/'report.json');audit=interpretation(report);save(RECORD/'result-interpretation.json',audit);cells=[]
    def md(text):cells.append(nbformat.v4.new_markdown_cell(text))
    def code(text):cells.append(nbformat.v4.new_code_cell(text))
    def fig(number,title):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{title}**')
        code(f'display(Image(filename=str(stability_output / "figures/{number}.png")))')
    md('''<a id="profile-stability-design-20261027"></a>
## 103. 上一轮电荷相容范围是否数值稳定？

上一轮160条有限剖面出现宽范围和分离区段。先验证增加搜索／加密网格会不会明显改变这些结论，再考虑独立覆盖率校准。这里不是增加机器学习训练，也不是另一次盲恢复效率试验。

从原生40／200ns和条件高斯24–64／128–256ns四个强→弱组中，预先取各组第0、5、10、15、20、25、30、35条，共32条。四组各第0条先导，原样并入主结果。样本全部是已公开开发资料，未按本轮好坏挑选；该子集不能代表未重算的另外128条。下文旧→新比较均使用同一32条，不能直接与第100节各40条中位数相减。

纯α、固定95／351ns双τ、1024×2ns完整模拟电压、充分动态范围保持。强制K=2，Q2为按时间排序的第二槽；既不传入真时刻／真Q，也不强制Q2<Q1。Q域50–50000pe、两时刻在观测K1触发附近[−96,+512]ns并裁至[0,2046]，与上一轮条件剖面一致。J仍为原相关高斯的残差二次型加logdet。Q为等效pe，未标定为沉积能量；这些受控双事件条件不对应某一随机注入率。

四阶段依次为：原始冻结曲线；同一网格增加搜索；所有相邻Q点插入几何中点；再在ΔJ≤1／4／9的跨界单元及最小值两侧插入几何中点。增加搜索使用6个起点、每次最多150迭代；新网格点用4个起点、150迭代，均正反扫描。起点来自原候选、此前剖面、邻点延续及观测触发偏移−8／0／+8ns和预定间隔；按代价排序后，在(t1/8ns,t2/8ns,lnQ1/0.25)坐标中优先选间距≥1的起点，不足再补。

每阶段保留旧点最小代价，也将本阶段最优点的固定其余参数切片纳入所有节点包络。因此网格阶段除了增加节点，也可能用新发现的切片降低旧点代价；这不是等计算量算法比较，更不是全局最优证明。所有未收敛标志保留。
''')
    code(f'''from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Image
stability_root = Path.cwd().resolve()
while not (stability_root / "pileup/experiments/{NAME}").is_dir():
    if stability_root.parent == stability_root:
        raise RuntimeError("请从项目根目录或jupyter目录运行")
    stability_root = stability_root.parent
stability_output = stability_root / "jupyter/artifacts/{NAME}"
stability_record = stability_root / "jupyter/run_records/{NAME}"
stability_report = json.loads((stability_output / "report.json").read_text())
print("32条预定开发波形，4阶段，已知K2；真值只用于下游评分。")
display(pd.DataFrame(stability_report["frames"]))
''')
    fig('0004-103-1','四组各8条的最小代价变化，显示全部样本；(a,b)原生40／200ns，(c,d)高斯中／长间隔。降低代价不代表参数更准确。')
    md('''<a id="profile-stability-curves-20261027"></a>
## 104. 曲线和每个分离区段分别比较

图104-1使用四条预定示例，所有阶段减去同一最终参考最小值，便于看搜索改进。图104-2则按每阶段自己的最小值定义ΔJ≤4，画出全部32条的每个相容区段。两图的参考不能混用：最小值降低时，相容范围也可能收缩，即使所有旧节点的绝对代价都没有变差。

图104-2横轴除以真实Q2仅用于展示；真值不参与主网格、起点或相容范围构造。区段间隙保留；范围总跨度=(最高端点−最低端点)，可能含不相容间隙。端点按lnQ线性插值，不是解析解。
''')
    fig('0004-104-1','四条预定示例的完整电荷域剖面；同一最终参考，横轴Q2为pe，竖点线为真实Q2，仅供展示。')
    fig('0004-104-2','四组全部32条ΔJ≤4分离区段；每阶段使用自己的最小值。真Q只作横轴归一，竖线为1。')
    md('''<a id="profile-stability-scores-20261027"></a>
## 105. 数值筛查与物理参数恢复分别评价

计算前约定：最后“全域中点→边缘加密”的参考最小值改善≤0.05、区段数相同、触及50／50000pe边界状态相同，且每个配对端点的乘性变化≤1.05，才称该层级通过本轮数值筛查。1／4／9分别检查；全部三层级通过才记入“均稳定”。这些容差是项目诊断约定，不是已证明的误差界；未穷尽起点、没有再次全域四分网格，也没有认证全局最优。

主结果锁定哈希后另运行真实Q2评分探针，只固定真实Q2，重新拟合Q1及两时刻，6个观测起点、150迭代。评分保留旧／新探针较低的实际代价，避免仅因探针搜索不同造成虚假包含率变化。真t／真Q1没有作为起点输入，探针不反馈主结果。如果探针比最终参考低超过0.01，明确标为数值遗漏；不能把负Δ当作更好的统计覆盖。

严格恢复仍为每Q±20%、每t±min(10ns,gap/4)、无额外输出。完整条件K2评分另报普通10ns匹配下的漏检、额外、误拆／误并、Q偏差与sigma68、时间误差以及未收敛；这些不替代原未知K盲效率。sigma68不是能谱FWHM。
''')
    md(table(report))
    md(interpretation_text(audit))
    fig('0004-105-1','(a)最后加密各层级通过数；(b)电荷范围总跨度中位；(c)有分离区段的条数；(d)已知K2最优点严格恢复数。各组分母均为8。')
    code('''display(pd.DataFrame(stability_report["stability_summary"]))
with pd.option_context("display.max_columns", None):
    display(pd.DataFrame(stability_report["conditional_metrics"]))
print("表中恢复、漏检、额外及误拆／误并均为已知K2条件评分，原盲输出未改变。")
''')
    md('''<a id="profile-stability-waveforms-20261027"></a>
## 106. 加密前后的完整波形分解

沿用四组各第0条预定波形，画完整0–2048ns（实际采样点0–2046ns），左列总电压，右列后槽分量。灰／绿线为原始／最终剖面的最优格点；黑线为真参数的条件均值，不等于原生随机雪崩分量。这里所有分解都假定K=2，不是最终盲算法的输出。
''')
    fig('0004-106-1','完整总电压与第二分量；(a,b)原生40ns，(c,d)原生200ns，(e,f)高斯中间隔，(g,h)高斯长间隔。')
    code('''stability_examples = json.loads((stability_output / "examples.json").read_text())
display(pd.DataFrame([dict(domain=e["domain"], index=e["index"], true_Q=e["truth_charge_pe"], true_t=e["truth_times_ns"],
    stage=s, Q=e[s]["charge_pe"], times=e[s]["times_ns"], cost=e[s]["cost"])
    for e in stability_examples for s in ("original", "contour")]))
''')
    md(f'''<a id="profile-stability-conclusion-20261027"></a>
## 107. 本轮结论与下一步

{conclusion(report)}

已有497个单元与历史图像原样保留。本轮不替换原160条剖面、不修改盲候选／门槛，也没有新训练、随机样本或ROOT。计算范围只覆盖32条系统抽取的开发样本、纯α固定双τ和局部K=2；不能据此推出普适最短可分间隔、最高计数率、free-τ或α/γ混合性能。

下一步首先针对仍变化的区段检查局部搜索与响应插值，再冻结数值规则后设计独立重复抽样，区分模型内覆盖与原生迁移。即使区间最终校准，也应同时报告检测选择、区间宽度与逐事件准确度，不以覆盖较高替代电荷分辨能力。

{LITERATURE}

完整[稳定性比较](artifacts/{NAME}/stability-summary.csv)、[条件范围汇总](artifacts/{NAME}/profile-summary.csv)、[条件K2全部指标](artifacts/{NAME}/conditional-metrics.csv)和[逐条结果](artifacts/{NAME}/report.json)保留。图号103-1、104-1／104-2、105-1、106-1固定不随执行顺序变化。
''')
    code('''stability_checks = json.loads((stability_record / "science-validation.json").read_text())
print({k:v for k,v in stability_checks.items() if k != "output_hashes"})
print("这些层级不是68%／95%／99%置信区间；没有新的独立覆盖率样本。")
''')
    for i,c in enumerate(cells):c.id=f'stability-20261027-{i:02d}'
    InteractiveShell.instance();namespace={};count=0;cwd=os.getcwd()
    try:
        os.chdir(PROJECT/'jupyter')
        for c in cells:
            if c.cell_type!='code':continue
            with capture_output(stdout=True,stderr=True,display=True) as capture:
                exec(compile(c.source,'<stability-notebook>','exec'),namespace)
            assert not capture.stderr,capture.stderr;c.outputs=[]
            if capture.stdout:c.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=capture.stdout))
            for item in capture.outputs:c.outputs.append(nbformat.v4.new_output('display_data',data=item.data,metadata=item.metadata))
            count+=1;c.execution_count=count
    finally:os.chdir(cwd)
    book.cells.extend(cells);assert book.cells[:497]==before.cells;nbformat.validate(book)
    encoded=(nbformat.writes(book)+'\n').encode()
    with tempfile.NamedTemporaryFile(dir=NOTEBOOK.parent,prefix=NOTEBOOK.stem+'.',suffix='.partial',delete=False) as stream:
        for i in range(0,len(encoded),256*1024):
            chunk=encoded[i:i+256*1024];assert stream.write(chunk)==len(chunk)
        stream.flush();os.fsync(stream.fileno());temp=Path(stream.name)
    assert temp.read_bytes()==encoded;nbformat.validate(nbformat.reads(encoded.decode(),as_version=4))
    assert NOTEBOOK.read_bytes()==original,'Preserve concurrent user edits'
    temp.replace(NOTEBOOK);assert NOTEBOOK.read_bytes()==encoded
    save(RECORD/'notebook-delivery.json',dict(passed=True,previous_cells=497,new_cells=len(cells),executed_cells=count,
        history_preserved=True,full_readback_verified=True,chunk_bytes=256*1024,notebook_sha256=digest(NOTEBOOK),source_sha256=digest(__file__)))
    print('STABILITY_NOTEBOOK_DELIVERED',len(cells),count,flush=True)


def documents():
    assert load(RECORD/'notebook-delivery.json')['passed']
    report=load(OUTPUT/'report.json');result=conclusion(report);numbers=table(report)
    interpretation_notes=interpretation_text(load(RECORD/'result-interpretation.json'))
    root=f'''<a id="{MARKER}"></a>
2026-10-27完成[条件电荷剖面的数值稳定性审核](pileup/experiments/{NAME}/README.md)：四组各系统抽取8条，共32条已公开开发波形，对比原曲线、增加观测起点、全域中点加密和相容边缘再加密。最后加密在1／4／9三层级同时通过预定端点／区段／参考代价筛查的有{report['stable_all_levels']}/32；{report['probe_improves_reference']}条真实Q2评分探针将最终参考再降低>0.01。结果为已知K2的有限数值诊断，不是盲效率改善或覆盖率校准；原160条剖面及默认流程保留。6项数值检查、32条四阶段曲线和独立标记的评分探针核验完成；本地Notebook0004第103–107节、五幅固定图及[L64](docs/literature-use.md#{MARKER})保留逐条变化和完整波形。'''
    literature=f'''<a id="{MARKER}"></a>
### L64：2026-10-27——条件电荷剖面的数值稳定性审核

**文献依据与审核范围**：{LITERATURE}

**项目实际计算**：{result}

{numbers}

{interpretation_notes}

每组原有40条，预先按第0／5／10／15／20／25／30／35条系统取8条，不按本轮拟合好坏选样。旧→新均为同一32条；不能直接与上一轮全部160条的中位数比较，也不认证未重算的128条。四个条件为原生先强后弱40／200ns及条件高斯24–64／128–256ns；两个域的Q／t分布不相同，不能把差异全归模型失配。

**数值方案**：保留原条件均值、候选依赖相关协方差、残差二次型加logdet，以及Q为50–50000pe、两时刻在观测单事件触发[−96,+512]ns裁至[0,2046]的参数域。先在原网格每方向用6个分离的观测起点、最多150迭代；再在每单元加几何中点，最后在1／4／9跨界单元及最小值两侧继续加几何中点，新节点每方向4起点／150迭代。起点来自观测缓存、旧剖面、邻点及预定触发偏移／间隔，按代价与数值坐标距离选择；没有传入真t／真Q。所有旧点的更优代价保留，阶段最优点的固定其余参数切片也纳入全节点包络，所以网格步骤不只是插值、也可能改善旧点。各阶段计算量不同，未证明全局最优。

末次加密的每层级仅在参考最小值改善≤0.05、相容区段数及边界标志不变、每个配对端点乘性变化≤1.05时记数值筛查通过。该规则先于计算冻结，是项目容差而非文献定理或严格数值误差界；没有穷举起点或再次全域四分网格。区段间隙保留，lnQ线性插值端点仍为近似。真值只用于评分和画图归一。

主结果哈希冻结后，另固定真实后槽Q2做评分探针，6观测起点／150迭代，未知Q1和两时刻；保留旧／新评分探针的较低实际代价。探针不反馈到主曲线、默认候选或门槛；低于最终参考>0.01单独记搜索遗漏，不能隐藏或解读成更好的覆盖。包含率不称置信水平。完整已知K2评分同时保留Q偏差／sigma68、时间误差、漏检、额外、误拆／误并和严格恢复；sigma68非能谱FWHM，Q未标定为沉积能量。

**复现与验证**：纯α固定95／351ns、完整1024×2ns电压、充分动态范围不变。本轮没有新ROOT、随机样本或训练；最多8核，先导4worker、主计算7worker且BLAS1。6项数值检查覆盖解析相关二次型、网格、分离区段与保留旧点；主曲线所有实际节点及评分探针重放代价，worker阻断truth／cases／ROOT访问。共享盘全仓库Git差异扫描的准备中断已记录；保留备份和协议，改为限定路径的版本差异记录后完成快照，不影响数值配置。

实验[pileup/experiments/{NAME}](../pileup/experiments/{NAME}/README.md)保存实际源码、模型／参数与协议哈希；本地Notebook为jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb第103–107节，图0004-103-1、104-1／104-2、105-1、106-1。原497单元、既有160条剖面和盲判决完整保留。

**待办**：先针对变化区段检查数值搜索及响应插值，冻结规则后再做独立重复抽样覆盖校准与原生迁移。小样本、条件K2、边界与选择不能忽略；尚未研究free-τ、α/γ、高率K≥3或普适分辨极限。正式有限样本／非规则剖面区间的一手证据及Baker DOI仍可由共享库另行整理，本项目不修改库。
'''
    experiment=f'''# Numerical stability of conditional later-charge profiles

{result}

{numbers}

{interpretation_notes}

## Fixed experiment

Pure alpha, fixed95/351ns,1024x2ns analog waveform, ample dynamic range. Known K=2, chronological later Q2, equivalent pe rather than calibrated deposited energy. All32 frames are previously exposed development data: systematic positions0,5,...,35 in each of the four40-frame native40/200ns and conditional-Gaussian24-64/128-256ns groups. The four first cases are pilot and plotting examples, retained in the final32. No new ROOT, training, seed or random draw. Inherited model/data/seed provenance is recorded locally.

Four stages: unchanged original finite profile; same grid with6 separated observed starts and150 iterations per direction; every geometric midpoint with4 starts/150 iterations; further midpoints at1/4/9 crossings and adjacent to the minimum,4 starts/150 iterations. Retain better old nodes and a fixed-nuisance slice from each stage minimum. Thus grid stages can also lower old-node costs through the slice. This is not an equal-compute algorithm comparison. The fixed response, parameter-dependent covariance and logdet remain unchanged, as do the charge50-50000pe and observed-trigger-local time domains. Details and numeric seed-selection rules are in protocol.json and the frozen source.

The last-step numerical screen requires reference improvement<=.05, unchanged interval count and boundary flags, and every matched interval endpoint within a factor1.05, separately at levels1/4/9. These are prespecified engineering tolerances, not certified error bounds or confidence levels. Finite multistart, incomplete convergence and interpolated contour endpoints remain limitations. No generalization to unexamined128 frames or physical independent tests is claimed.

After primary hashes are frozen, a separate true-Q2-only scoring probe refits nuisance from observed starts. Keep the lower old/new probe cost for scoring, never as a primary input. A cost below the primary reference by>.01 flags a missed point. Full known-K2 recovery, Q bias/sigma68, timing, missing/extra and split/merge are saved separately from interval inclusion. Original blind count candidates/gates and all160 earlier profiles remain untouched.

## Reproduction

Run from the project root using the project Python environment with OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1. The limit is8 workers; pilot uses4 and main7. Frozen source and model hashes are in freeze.json, report-freeze.json and delivery-freeze.json. run.py/report.py/deliver.py load these snapshots rather than current live code; no Git checkout is needed. Models/arrays/local records are ignored by Git and must be supplied alongside a reproducibility bundle.

```bash
bash jupyter/run_gagg_python.sh pileup/experiments/{NAME}/run.py tests
bash jupyter/run_gagg_python.sh pileup/experiments/{NAME}/report.py analyze
bash jupyter/run_gagg_python.sh pileup/experiments/{NAME}/report.py figures
bash jupyter/run_gagg_python.sh pileup/experiments/{NAME}/deliver.py validate
```

pilot/infer/probe are one-time commands refusing existing output. To recompute the complete run, use a separate project tree with the frozen inputs. Preparation interruption from an unscoped Git diff on shared storage was recorded, then resumed with identical protocol and verified backups. Scoped Git stat has a60s timeout; exact source and prior-document backups preserve provenance. Notebook writes use256KiB chunks, full readback and concurrent-edit checks.

Local outputs: jupyter/artifacts/{NAME}/; records: jupyter/run_records/{NAME}/. Notebook0004 sections103-107, five fixed figures including full0-2048ns comparisons. Original497 cells and outputs preserved; new cells executed. Validation checks6 analytic/numerical tests,32 four-stage profiles,32 true-Q2 probes,128 conditional-fit scores,384 contour scores,288 paired comparisons, actual cost replay, source/data hashes, truth-open guards, figures and Notebook preservation. No new default, Git commit or push.

## Evidence and next step

Only existing local scientific notes were reused: Baker/Cousins1984,NIM221437-442(PDF1/4/5; DOI unverified) distinguishes estimation, coverage and goodness of fit for Poisson/multinomial histograms, not this correlated voltage. Liddick2012,NIMA66970-78,DOI10.1016/j.nima.2011.12.034(PDF3-8;65umDSSD,40MSPS/12bit,25us,artificial doubles) motivates separate detection/accuracy/background reporting; its numerical limits are not transferred. Both notes remain AI drafts pending human scientific review. See docs/literature-use.md L64 for precise adoption and unresolved issues.

Address remaining grid/search/interpolation changes before independent repeated-sampling coverage calibration; distinguish model-internal calibration and native transfer. No free-tau,alpha/gamma,high-rateK>=3 or universal timing-limit result is claimed.
'''
    for name,backup,addition in (('README.md','README.before.md',root),('docs/literature-use.md','literature.before.md',literature)):
        path=PROJECT/name;old=(RECORD/backup).read_text();assert path.read_text()==old and MARKER not in old
        path.write_text(old.rstrip()+'\n\n'+addition.rstrip()+'\n')
    p=EXPERIMENT/'README.md';assert not p.exists();p.write_text(experiment)
    files={str(p.relative_to(PROJECT)):digest(p) for p in (PROJECT/'README.md',PROJECT/'docs/literature-use.md',p)}
    save(RECORD/'document-delivery.json',dict(passed=True,files=files,source_sha256=digest(__file__)))
    print('STABILITY_DOCUMENTS_DELIVERED',flush=True)


def validate():
    check_freeze()
    for name in ('numerical-tests','pilot-validation','inference','probe-validation','science-validation',
                 'figure-delivery','figure-visual-review','result-interpretation','notebook-delivery','document-delivery'):
        assert load(RECORD/(name+'.json'))['passed'],name
    manifests={key:load(EXPERIMENT/(key+'.json')) for key in ('freeze','report-freeze','delivery-freeze')}
    for manifest in manifests.values():
        for name,value in manifest['files'].items():assert digest(EXPERIMENT/name)==value,name
    source={}
    groups={'inference-source/pileup/analysis':('alpha_pair_profile_stability.py','alpha_pair_stability_study.py','test_alpha_pair_profile_stability.py'),
            'report-source':('alpha_pair_stability_report.py',),'delivery-source':('alpha_pair_stability_delivery.py',)}
    for folder,names in groups.items():
        for name in names:
            p=PROJECT/'pileup/analysis'/name;ast.parse(p.read_text());assert digest(p)==digest(EXPERIMENT/folder/name)
            source[str(p.relative_to(PROJECT))]=digest(p)
    provenance=load(RECORD/'input-provenance.json');prior=provenance['prior_validation']
    for name,value in provenance['inputs'].items():assert digest(PROJECT/name)==value,name
    for name,value in provenance['prior_freezes'].items():assert digest(PREVIOUS/name)==value,name
    assert digest(OLD_RECORD/'input-provenance.json')==provenance['inherited_provenance_sha256']
    for name,value in prior['source_sha256'].items():assert digest(PROJECT/name)==value,name
    for name,value in prior['documentation_sha256'].items():
        if name not in ('README.md','docs/literature-use.md'):assert digest(PROJECT/name)==value,name
    for old in ('alpha_pair_information_20261027','alpha_pair_charge_allocation_20261027',
                'alpha_pair_risk_training_v2_20261027','alpha_pair_physical_gate_20261027'):
        for name,value in load(PROJECT/'jupyter/run_records'/old/'final-validation.json')['source_sha256'].items():
            assert digest(PROJECT/name)==value,name
    science=load(RECORD/'science-validation.json')
    for name,value in science['output_hashes'].items():assert digest(OUTPUT/name)==value,name
    audit=load(RECORD/'result-interpretation.json');assert audit==interpretation(load(OUTPUT/'report.json'))
    run=load(RECORD/'inference.json');probe=load(RECORD/'probe-validation.json')
    assert digest(OUTPUT/'profiles.jsonl')==run['sha256']==probe['primary_profiles_sha256']
    assert digest(OUTPUT/'truth-probes.jsonl')==probe['sha256']
    expected={(t['domain'],t['index']) for t in load(EXPERIMENT/'tasks.json')}
    old_tasks=load(PREVIOUS/'tasks.json');assert load(EXPERIMENT/'tasks.json')==old_tasks[::5]
    curves=rows(OUTPUT/'profiles.jsonl');probes=rows(OUTPUT/'truth-probes.jsonl')
    assert len(curves)==len(probes)==32
    assert {(r['domain'],r['index']) for r in curves}=={(r['domain'],r['index']) for r in probes}==expected
    original={(r['domain'],r['index']):r for r in rows(INPUT/'profiles.jsonl')};assert len(original)==160
    for r in curves:
        old=dict(r['stages']['original']);assert old.pop('stage')=='original'
        assert old==original[r['domain'],r['index']]
        previous=None
        for stage in ('original','search','dense','contour'):
            c=r['stages'][stage];assert c['reference_cost']==min(f['cost'] for f in c['fits'])
            assert all(f['cost']<=s+1e-8 for f,s in zip(c['fits'],c['slice_cost']))
            current=dict(zip(c['q_grid_pe'],c['fits']))
            if previous:
                for q,f in previous.items():assert current[q]['cost']<=f['cost']+1e-9
            previous=current
    assert run['primary_truth_inputs']==[] and run['label_files_blocked_in_workers']
    assert run['maximum_cost_replay_error']<1e-9 and science['maximum_probe_cost_replay_error']<1e-9
    assert science['known_K2_scored_fits']==128 and science['level_scores']==384 and science['comparisons']==288
    program='''import sys,json
sys.path.insert(0,sys.argv[1])
from pileup.analysis.alpha_pair_stability_study import initialize
initialize()
blocked=[]
for name in ('forbidden.truth.npz','forbidden.cases.json','forbidden.root','/tmp/truth/forbidden.npy','native-records.json'):
    try:open(name,'rb')
    except PermissionError:blocked.append(name)
    else:raise AssertionError(name)
print(json.dumps(blocked))
'''
    blocked=json.loads(subprocess.check_output([sys.executable,'-c',program,str(EXPERIMENT/'inference-source')],cwd=PROJECT,text=True))
    assert len(blocked)==5
    documents=load(RECORD/'document-delivery.json')
    for name,value in documents['files'].items():assert digest(PROJECT/name)==value,name
    for name,backup in (('README.md','README.before.md'),('docs/literature-use.md','literature.before.md')):
        assert digest(RECORD/backup)==prior['documentation_sha256'][name]
        assert (PROJECT/name).read_text().startswith((RECORD/backup).read_text().rstrip()+'\n\n')
    before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4)
    book=nbformat.read(NOTEBOOK,as_version=4);nbformat.validate(book)
    assert len(before.cells)==497 and book.cells[:497]==before.cells
    assert digest(RECORD/'notebook.before.ipynb')==prior['notebook_sha256']==provenance['notebook_sha256']
    record=load(RECORD/'notebook-delivery.json');new=book.cells[497:];codes=[c for c in new if c.cell_type=='code']
    assert digest(NOTEBOOK)==record['notebook_sha256'] and len(new)==record['new_cells'] and len(codes)==record['executed_cells']
    assert all(c.id.startswith('stability-20261027-') for c in new) and len({c.id for c in book.cells})==len(book.cells)
    assert all(c.execution_count and c.outputs and not any(o.output_type=='error' for o in c.outputs) for c in codes)
    anchors=re.findall(r'id="codex-figure-([^"]+)"','\n'.join(c.source for c in book.cells if c.cell_type=='markdown'))
    assert len(anchors)==len(set(anchors))
    for c in new:
        if c.cell_type=='markdown':assert not re.search(r'\]\([^)]*\.(?:pdf|ipynb|md)(?:[#?)][^)]*)?\)',c.source,re.I)
    figures=load(RECORD/'figure-delivery.json')['figures'];assert figures==load(RECORD/'figure-visual-review.json')['figures']
    for number,value in figures.items():assert digest(OUTPUT/'figures'/(number+'.png'))==value
    embedded=[hashlib.sha256(base64.b64decode(o.data['image/png'])).hexdigest() for c in codes for o in c.outputs if 'image/png' in o.get('data',{})]
    assert len(embedded)==5 and set(embedded)==set(figures.values())
    for name in (str(NOTEBOOK.relative_to(PROJECT)),f'jupyter/artifacts/{NAME}/report.json',f'pileup/experiments/{NAME}/joint-response.npz'):
        assert subprocess.run(['git','check-ignore','-q',name],cwd=PROJECT).returncode==0,name
    subprocess.run(['git','-c','core.fileMode=false','diff','--check','--','README.md','docs/literature-use.md'],cwd=PROJECT,check=True,timeout=120)
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip()==load(RECORD/'environment.json')['git_head']
    save(RECORD/'final-validation.json',dict(passed=True,numerical_tests=6,profile_frames=32,stages=4,truth_Q2_scoring_probes=32,
        original160_profiles_preserved=True,known_K2_scored_fits=128,contour_scores=384,comparisons=288,
        maximum_cost_replay_error=run['maximum_cost_replay_error'],maximum_probe_cost_replay_error=science['maximum_probe_cost_replay_error'],
        frozen_files={k:len(v['files']) for k,v in manifests.items()},source_sha256=source,truth_open_probes_blocked=blocked,
        original497_cells_preserved=True,new_cells=len(new),executed_new_cells=len(codes),figures=figures,
        notebook_sha256=digest(NOTEBOOK),documentation_sha256=documents['files'],prior_research_preserved=True,
        result_interpretation_sha256=digest(RECORD/'result-interpretation.json'),
        final_stable_by_level=audit['final_stable_by_level'],inclusion_disagreements=audit['inclusion_disagreements'],
        new_ROOTs=0,new_training=0,new_random_draws=0,no_default_change=True,no_git_commit_or_push=True,
        limitations=['32 systematic exposed development frames','conditional knownK2 profiles','finite search and nested grid screen only',
                     'not a confidence-coverage calibration','purealpha fixedtau local trigger domain']))
    print('STABILITY_FINAL_VALIDATION_PASSED',len(new),len(codes),flush=True)
