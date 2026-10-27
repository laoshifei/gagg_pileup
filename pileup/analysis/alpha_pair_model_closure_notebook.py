"""Append executed model-closure diagnostics while preserving Notebook0004 history."""
from __future__ import annotations

import base64
import hashlib
import json
import os
import re
import shutil

import nbformat

from .alpha_pair_model_closure_study import PROJECT, OUTPUT, RECORD
from .alpha_pair_resolution_data import load, save, digest
from .alpha_pair_conditional_notebook import atomic_notebook


def deliver():
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output

    assert load(RECORD / 'science-validation.json')['passed']
    assert load(RECORD / 'figure-delivery.json')['passed']
    report = load(OUTPUT / 'report.json')
    path = PROJECT / 'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    previous = path.read_bytes()
    notebook = nbformat.reads(previous.decode(), as_version=4)
    before = nbformat.read(RECORD / 'notebook.before.ipynb', as_version=4)
    assert len(before.cells) == 364
    assert notebook.cells[:len(before.cells)] == before.cells
    if len(notebook.cells) > len(before.cells):
        assert digest(path) == load(RECORD / 'notebook-delivery.json')['notebook_sha256']
        assert all(c.get('id', '').startswith('model-closure-20261026-') for c in notebook.cells[len(before.cells):])
        notebook.cells = notebook.cells[:len(before.cells)]
    shutil.copyfile(path, RECORD / 'notebook.before-delivery.ipynb')
    original = json.loads(json.dumps(notebook.cells))
    cells = []

    def md(source):
        cells.append(nbformat.v4.new_markdown_cell(source))

    def code(source):
        cells.append(nbformat.v4.new_code_cell(source))

    def figure(number, title):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{title}**')
        code(f'display(Image(filename=str(closure_output / "figures/{number}.png")))')

    md('''<a id="model-closure-design-20261026"></a>
## 71. 模型内对照：生成过程与拟合假设一致时，还会失败吗？

第70节提出的问题在本轮实际检验：保持原生事件的数值Q/t和整个盲拟合流程，另由拟合所用的条件均值、协方差生成高斯波形。比较这两类输入，检查现有困难有多少能在这个统计模型内重现。这里没有训练新网络，也没有修改物理拟合算法。

**原生数据（native）：**已公开latest批中的760帧，包括mid_mid、high_weak、weak_high三种组合×8/40/200 ns×40对，以及80空窗和四档各80个单事件。全部沿用原生光子／共同微胞模拟，不叠加独立电压。它们已是开发和诊断资料，本轮没有生成新光学母ROOT。

**匹配模型数据（matched Gaussian）：**对上述每帧的相同Q/t生成两个独立高斯噪声实现，共1520帧；另在280个独立门槛母事件的Q/t上各生成一帧，只用于标定门槛。因此实际新盲拟合1800帧，而评价帧数为1520。两次噪声抽样不等于增加了光学母事件数量。原生评价320个母事件与训练280、门槛280个母事件互不交叉，条件和间隔之间仍复用母事件。

原生Q是该次模拟实现的雪崩等效电荷pe；高斯模型中的Q是数值相同的**条件参数**，不是新实现的一批雪崩的积分电荷，也不是沉积能量。t沿用有效SiPM面首光子时刻。固定95/351 ns闪烁双τ、2 ns×1024点、0–2048 ns完整模拟电压和充足动态范围，未知局部K≤2。

模型只给出固定FIR变换后1021维数据的分布，协方差为C=B−VVᵀ。本轮用带状Cholesky和秩≤2的矩阵平方根精确采样该协方差，再逆FIR得到电压；前三个电压残差样点显式设零，没有假设它们的平稳联合分布。16个预先指定的FIR零空间边界扰动另行检查。这里的读出2/20/2 ns常数与95/351 ns发光衰减时间不同。

推断worker只读电压和冻结模型，禁止读取truth、case标签及ROOT。真Q/t只用于生成、评分和明确标注的辅助诊断。没有按评价结果调门槛或选择有利噪声种子；所有未收敛输出都保留。

文献采用沿用Fixsen等2004的条件均值／非平稳协方差思想（DOI：10.1016/j.nima.2003.11.313，PDF2式1–2及PDF3；TES模拟、每能量5000次，既有AI精读稿待人工审核）。原文是协方差加权最小二乘；本项目的logdet、高斯精确采样及此项对照不是该文的GAGG结论。本轮复用已有精读，没有重新核查原PDF。
''')
    code('''from pathlib import Path
import json
import numpy as np
import pandas as pd
from IPython.display import display, Image
closure_root = Path.cwd().resolve()
while not (closure_root / "pileup/experiments/alpha_pair_model_closure_20261026").is_dir():
    if closure_root.parent == closure_root:
        raise RuntimeError("请从项目根目录或jupyter目录运行")
    closure_root = closure_root.parent
closure_output = closure_root / "jupyter/artifacts/alpha_pair_model_closure_20261026"
closure_record = closure_root / "jupyter/run_records/alpha_pair_model_closure_20261026"
closure_report = json.loads((closure_output / "report.json").read_text())
closure_summary = pd.DataFrame(closure_report["summary"])
closure_preflight = json.loads((closure_record / "sampler-preflight.json").read_text())
print("原生评价760帧；匹配模型评价1520帧，另280帧只标定门槛。")
print("采样检查通过：", closure_preflight["passed"], "；新ROOT：", closure_report["new_roots"])
''')
    md('''<a id="model-closure-blind-20261026"></a>
## 72. 同一盲算法、两类输入，以及门槛的交叉检查

主比较对两类输入都使用原280个原生单事件标定的增峰门槛。第二组比较用280个匹配模型单事件重新标定门槛，再分别施用于两类评价输入。除门槛来源外，条件均值、候选依赖协方差、logdet、全部候选搜索和拟合Q≥500保护均不变。

增峰统计为max(NLL₁−NLL₂,0)，按拟合单事件Q的1500/5000/12500 pe分档，加一并包含并列的经验上尾秩p≤0.025。高斯门槛没有用评价帧的真值调优；交叉使用只是一项诊断，本轮不据此替换默认方案。

严格整对成功要求两真实Q均≥500 pe、各Q误差≤20%、各t误差≤min(10 ns,间隔/4)，且没有额外输出。表中同时保留漏检、误拆／误并标志、普通10 ns匹配和匹配后Q/t误差。未匹配到真事件的输出也计作额外，不能一概解释为完全无信号时凭空检出。控制组每类原生80帧、高斯160帧；双事件每点原生40对、高斯80对。
''')
    figure('0004-72-1', '三种双事件组合的严格恢复率；(a)中／中，(b)先强后弱，(c)先弱后强。四条曲线区分输入与门槛来源；两种门槛在同类输入上的曲线完全重合。')
    code('''display(closure_summary[closure_summary.kind == "pair"][[
    "method", "condition", "gap_ns", "success", "eligible_frames", "extra_frames",
    "q_error_p16_p50_p84", "time_error_p16_p50_p84_ns"]])
closure_paired = json.loads((closure_output / "paired-comparison.json").read_text())
display(pd.DataFrame([dict(方法=r["method"], 参照=r["reference"], 组合=r["condition"], 间隔_ns=r["gap_ns"],
    恢复率差_pp=100*r["difference"], 配对95区间_pp=[100*v for v in r["paired_block95"]])
    for r in closure_paired["comparisons"] if r["method"] == "matched/native_gate"]))
print("2000次配对重抽保留母事件块的跨条件复用；区间条件于这些参数、已拟合模型和两次噪声抽样，不是独立探测器性能区间。")
print("全成／全败点的bootstrap差值区间退化为[0,0]，不代表总体差异已被精确确定；逐帧Wilson区间也仅作补充描述。")
''')
    figure('0004-72-2', '空窗和单事件控制；(a)额外输出帧比例，(b)严格成功率，(c)匹配后Q偏差中位，空窗无Q偏差。')
    code('''display(closure_summary[closure_summary.kind != "pair"][[
    "method", "condition", "success", "eligible_frames", "extra_frames", "q_error_p16_p50_p84"]])
''')
    md('''<a id="model-closure-truth-diagnostics-20261026"></a>
## 73. 区分增峰选择、有限搜索与参数精度

在两种强弱方向、三个间隔上，对原生240对沿用上一轮诊断，对高斯480对执行相同诊断。四层依次为：原生门槛的盲结果；评分端指定K=2并取已有候选；真实Q/t初始化后联合优化；固定真实时刻仅拟合Q。

真实初始化允许各t在真值±48 ns内变化、最多200迭代；固定真实时刻使用观测总Q的七种预定分配起点、每个最多100迭代。后两项有更多真值信息和不同搜索预算，只能定位问题。即使某项结果很好，也不能称作可实现的盲恢复率或理论最优界限。

代价比较使用相同rᵀC⁻¹r＋logdet C。真值初始化找到更低代价可提示原搜索没有到达同一解；它找到更准确但代价更高的参数则说明单纯扩大同目标搜索未必改善恢复。两者都不单独证明全局最优或不可分辨。
''')
    figure('0004-73-1', '四层诊断；(a,b)原生，(c,d)匹配高斯；左列先强后弱，右列先弱后强。除盲结果外均用到真实K或更多真值。')
    code('''closure_diagnostic = json.loads((closure_output / "diagnostic-summary.json").read_text())
display(pd.DataFrame(closure_diagnostic["summary"]))
closure_costs = pd.DataFrame(closure_diagnostic["cost_comparisons"])
display(closure_costs.groupby(["domain", "condition", "gap_ns"])[[
    "true_start_lower_cost", "new_correct_but_higher_cost"]].sum())
print("new_correct_but_higher_cost：真实初始化使严格参数由错变对，但同一代价没有比原K=2候选更低。")
''')
    md('''<a id="model-closure-residuals-20261026"></a>
## 74. 在真实参数处检查噪声模型

先减去真实Q/t对应的条件均值，再用同一条件协方差白化。若输入确为这个高斯模型，1021个白化坐标应为独立标准正态，平方和／1021应在1附近。这是在**真实参数**处的检查，没有从1021减去拟合参数数目，也不把它当作拟合残差的约化χ²。

原生数据可能同时包含条件均值误差、协方差误差和非高斯涨落。平方范数、残差偏度／峰度、相邻相关与强事件对齐后的平均残差用于描述差异，不能仅据一项统计量归因于某个物理过程。双事件使用共同微胞，弱分量的统计也不能直接从孤立单事件照搬。

图中直方图按全部样点数归一化，显示区间外的比例另列，未把重尾样点删除后重新归一化。汇总和图中的帧来自复用母事件，不把所有时间点或两次高斯噪声实现当成独立实物曝光。
''')
    figure('0004-74-1', '真参数白化检查；(a)平方范数／1021的帧间16–50–84百分位，(b)白化残差分布，(c)滞后相关，(d)先强后弱40 ns按强事件对齐的平均残差。')
    code('''closure_residual = json.loads((closure_output / "residual-diagnostics.json").read_text())
display(pd.DataFrame(closure_residual["summary"]))
closure_boundary = json.loads((closure_output / "boundary-comparison.json").read_text())
print("16个预定FIR零空间边界扰动中，输出K改变：", closure_boundary["changed_selected_count"],
      "；严格成功判定改变：", closure_boundary["changed_strict_success"])
print("这是有限搜索的边界敏感性检查；不补充真实硬件的窗前历史模型。")
''')
    md('''<a id="model-closure-full-waveforms-20261026"></a>
## 75. 完整波形、弱分量与计算记录

沿用预先指定的high_weak第0对，分别展示40、200 ns；高斯取第0次噪声抽样，不按成功与否挑图。所有面板展示0–2048 ns完整窗口，全部拟合曲线使用盲结果和原生门槛。

原生行的弱分量来自原联合微胞模拟保存的雪崩归属；高斯模型没有唯一的逐雪崩弱分量，所以高斯行只能显示条件弱均值。这两种参考曲线含义不同。若盲流程只选一个输出，拟合弱分量显示为零；若选两个，画较小Q输出，它也未必正确对应真实弱事件。真时刻和分量只在展示／评价端读取。
''')
    figure('0004-75-1', '全窗对照，四行依次为原生40 ns、高斯40 ns、原生200 ns、高斯200 ns；三列依次为总波形、弱分量、总残差，子图(a)至(l)。')
    code('''closure_examples = json.loads((closure_output / "examples.json").read_text())
display(pd.DataFrame([dict(输入=e["domain"], 间隔_ns=e["gap_ns"], 严格成功=e["score"]["success"],
    真t_ns=e["score"]["true_t"], 拟合t_ns=e["score"]["fit"]["times_ns"],
    真Q_pe=e["score"]["true_q"], 拟合Q_pe=e["score"]["fit"]["charge_pe"])
    for e in closure_examples]))
display(pd.DataFrame([dict(输入与门槛=k, 最终选中候选未收敛=v)
    for k, v in closure_report["selected_nonconverged"].items()]))
print("高斯盲拟合耗时16/50/84/100百分位（秒）：", closure_report["runtime_seconds_quantiles"])
print("包括校准帧和小样本阶段，并发不同，仅作运行记录，不与缓存原生结果做速度排名。")
''')
    counts = []
    for method in ('native/native_gate', 'matched/native_gate', 'matched/matched_gate'):
        values = [next(r for r in report['summary'] if r['method'] == method
                       and r['condition'] == 'high_weak' and r['gap_ns'] == gap) for gap in (8., 40., 200.)]
        counts.append(f"{method}：" + '、'.join(f"{int(r['gap_ns'])} ns为{r['success']}/{r['eligible_frames']}" for r in values))
    md('''<a id="model-closure-conclusion-20261026"></a>
## 76. 本轮结果的解释范围

先强后弱的实际盲结果为：

''' + '\n\n'.join(counts) + '''。

**主要发现：在这个完全匹配的高斯统计模型内，先强后弱8、40 ns的盲失败仍然存在。** 共用原生门槛时，200 ns为原生50%、高斯52.5%；改用高斯门槛后，三点严格恢复数均没有增加。因此，仅修正原生响应／噪声近似或重新标定门槛，不能被当作当前短间隔困难的已验证解法。

门槛迁移还有代价：高斯门槛用于原生high单事件时，额外输出由0/80升至7/80；该门槛在高斯high单事件中只有2/160。这个差异提示门槛尾部统计不能直接跨生成域照搬，但它没有解释掉匹配模型自身的双事件失败。本轮不升级默认方法。

40 ns高斯先强后弱中，已有K=2候选仅2/80严格正确，真实初始化为24/80、固定真实时刻拟合Q为29/80。后两项虽改善，仍远未达到可靠逐事件恢复，也都用了不可在盲推断中取得的信息。原生对应为3/40、9/40、12/40，不能把盲结果全归因于增峰门槛。

两域的总白化平方范数接近，仍不能说噪声模型已完全准确：原生先强后弱40 ns的逐帧白化残差偏度中位约0.251，高斯约−0.0049；原生其他条件也有正偏和较重正尾。现有高斯模型仍是近似，只是这项近似不是短间隔失败的充分解释。

下一步应把“检出第二事件”和“分别估准两个Q/t”拆开，先量化弱事件电荷与时刻的联合不确定性，再检验有针对性的估计或搜索改进；不把当前点估计失败直接命名为硬性分辨极限。这一后续分析尚未执行。

原生与高斯的同Q/t对照改变了整个噪声分布及条件均值误差，并非单一器件效应消融。真值辅助诊断与白化检查仍需同时阅读，不能从总波形拟合得好直接推断两个事件的Q/t已经分别恢复。

即使真实初始化或固定时刻仍失败，也不能断言“任何算法都不可能分开”；即使匹配模型恢复很好，也不能证明实物探测器可达到相同性能。这里没有全局优化证明、信息论下界、连续Q／位置失配研究或独立实物验证。结论仅限固定双τ、当前几何／噪声、纯α、局部K≤2，不回答free-τ最短窗、α/γ、随机高率K≥3或普适最小间隔。

完整[逐点结果](artifacts/alpha_pair_model_closure_20261026/summary.csv)、[对照报告](artifacts/alpha_pair_model_closure_20261026/report.json)、[配对区间](artifacts/alpha_pair_model_closure_20261026/paired-comparison.json)和[真参数残差统计](artifacts/alpha_pair_model_closure_20261026/residual-diagnostics.json)保留。
''')
    code('''closure_checks = json.loads((closure_record / "science-validation.json").read_text())
print("通过：", closure_checks["passed"], "；盲评分：", closure_checks["scores"],
      "；真值辅助评分：", closure_checks["diagnostic_scores"], "；真参数残差帧：", closure_checks["residual_frames"])
print("原生结果与上一轮逐组一致：", closure_checks["native_scores_match_previous"],
      "；原输入及冻结源码保持：", closure_checks["original_inputs_unchanged"] and closure_checks["freeze_unchanged"])
print("盲候选／诊断候选代价重放最大相对差：", closure_checks["maximum_candidate_relative_cost_error"],
      closure_checks["maximum_diagnostic_relative_cost_error"])
''')
    for i, cell in enumerate(cells):
        cell['id'] = f'model-closure-20261026-{i:02d}'
    InteractiveShell.instance()
    namespace = {}
    executed = 0
    cwd = os.getcwd()
    try:
        os.chdir(PROJECT / 'jupyter')
        for cell in cells:
            if cell.cell_type != 'code':
                continue
            with capture_output(stdout=True, stderr=True, display=True) as capture:
                exec(compile(cell.source, '<model-closure-notebook>', 'exec'), namespace)
            assert not capture.stderr, capture.stderr
            cell.outputs = []
            if capture.stdout:
                cell.outputs.append(nbformat.v4.new_output('stream', name='stdout', text=capture.stdout))
            for output in capture.outputs:
                cell.outputs.append(nbformat.v4.new_output('display_data', data=output.data, metadata=output.metadata))
            executed += 1
            cell.execution_count = executed
    finally:
        os.chdir(cwd)
    notebook.cells.extend(cells)
    assert notebook.cells[:len(original)] == original
    anchors = re.findall(r'id="codex-figure-([^\"]+)"', '\n'.join(c.source for c in notebook.cells if c.cell_type == 'markdown'))
    assert len(anchors) == len(set(anchors))
    for cell in cells:
        if cell.cell_type == 'markdown':
            assert not re.search(r'\]\([^)]*\.(?:pdf|ipynb|md)(?:[#?)][^)]*)?\)', cell.source, re.I)
    images = [output['data']['image/png'] for cell in cells if cell.cell_type == 'code'
              for output in cell.outputs if 'image/png' in output.get('data', {})]
    figures = load(RECORD / 'figure-delivery.json')['figures']
    assert len(images) == len(figures) == 5
    assert {hashlib.sha256(base64.b64decode(image)).hexdigest() for image in images} == set(figures.values())
    nbformat.validate(notebook)
    atomic_notebook(path, notebook, previous)
    save(RECORD / 'notebook-delivery.json', dict(passed=True, previous_cells=len(original), new_cells=len(cells),
        executed_cells=executed, history_preserved=True, figures=figures, notebook_sha256=digest(path),
        source_sha256=digest(__file__)))
    print('MODEL_CLOSURE_NOTEBOOK_DELIVERED', len(cells), executed, flush=True)


if __name__ == '__main__':
    deliver()
