"""Append executed profile/search results without rewriting Notebook0004 history."""
from __future__ import annotations

import base64
import hashlib
import json
import os
import re
import shutil

import nbformat

from .alpha_pair_weak_search_study import PROJECT, EXPERIMENT, OUTPUT, RECORD
from .alpha_pair_resolution_data import load, save, digest
from .alpha_pair_conditional_notebook import atomic_notebook


def deliver():
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    assert load(RECORD / 'science-validation.json')['passed']
    assert load(RECORD / 'figure-delivery.json')['passed']
    report = load(OUTPUT / 'report.json')
    path = PROJECT / 'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    previous = path.read_bytes(); notebook = nbformat.reads(previous.decode(), as_version=4)
    before = nbformat.read(RECORD / 'notebook.before.ipynb', as_version=4)
    assert notebook.cells[:len(before.cells)] == before.cells
    if len(notebook.cells) > len(before.cells):
        assert digest(path) == load(RECORD / 'notebook-delivery.json')['notebook_sha256']
        assert all(c.get('id', '').startswith('weak-search-20261026-') for c in notebook.cells[len(before.cells):])
        notebook.cells = notebook.cells[:len(before.cells)]
    shutil.copyfile(path, RECORD / 'notebook.before-delivery.ipynb')
    original = json.loads(json.dumps(notebook.cells)); cells = []

    def md(source): cells.append(nbformat.v4.new_markdown_cell(source))
    def code(source): cells.append(nbformat.v4.new_code_cell(source))
    def figure(number, title):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{title}**')
        code(f'display(Image(filename=str(search_output / "figures/{number}.png")))')

    md('''<a id="weak-search-profiles-20261026"></a>
## 66. 弱事件时间—电荷剖面：检查现有搜索遗漏了什么

上一轮先强后弱40 ns盲恢复0/40，即使评分端直接取已有K=2候选也只有3/40。本轮保持面积守恒均值、联合协方差与logdet、95/351 ns双τ、2 ns采样和完整2048 ns模拟电压，先检查局部目标函数，再比较两种扩大搜索的方法。动态范围仍假定充足，不按ADC削顶拒绝窗口。

**数据身份改变须明确：**第60节1600帧和第63节公开的1600帧，在本轮都属于已公开开发资料，分别记为older与latest。第63节冻结后的独立结果保持原样；用该批研究新搜索后，不能把本轮结果再次称为独立性能。另用原来的280个门槛母样本；它们与均值／噪声训练280个母样本及两批评价各320个母样本互不交叉。不同条件与间隔仍复用母事件。

weak、low、mid、high为四档单α源能量对应的电荷类别；A_B表示先A后B，全部是α。Q为实际雪崩等效pe，t为SiPM有效面首光子；源能量、沉积能量及Q不混用。

这里的剖面把弱事件时间和电荷逐个固定，每一点重新优化强事件Q/t，并重新计算候选依赖协方差及logdet。预先选high_weak每间隔的第0、1对，共六例。弱时间是真值±48 ns内13点，弱Q为真值0.25–3倍的9点；强时间限于真值±48 ns，使用真实强事件和原单事件拟合两个初值。

**这是有真值辅助的局部诊断。** 已知K=2、强弱身份、扫描范围及一个强事件初值均来自真值，不进入盲方法；局部最低值不等于全局最优或物理分辨极限。图中的ΔNLL只是代价差，没有按Wilks近似解释为置信区域。噪声数据的拟合代价低于精确真值点本身正常，不能仅凭这一点断言模型有偏。

文献层面沿用Fixsen等2004的条件均值／协方差思路（DOI：10.1016/j.nima.2003.11.313，PDF2式1–2、PDF3，TES模型每能量5000次；既有AI笔记待人工审核）。本轮只复用此前精读，未重新核PDF；时间—电荷剖面和搜索方案是项目诊断，不是该文报告的GAGG堆积结果。
''')
    code('''from pathlib import Path
import json
import numpy as np
import pandas as pd
from IPython.display import display, Image
search_root=Path.cwd().resolve()
while not (search_root / "pileup/experiments/alpha_pair_weak_search_20261026").is_dir():
    if search_root.parent==search_root:
        raise RuntimeError("请从项目根目录或jupyter目录运行")
    search_root=search_root.parent
search_output=search_root / "jupyter/artifacts/alpha_pair_weak_search_20261026"
search_record=search_root / "jupyter/run_records/alpha_pair_weak_search_20261026"
search_report=json.loads((search_output / "report.json").read_text())
search_diagnostics=json.loads((search_output / "diagnostic-summary.json").read_text())
print("两批均为开发资料；本轮选中：",search_report["chosen"])
''')
    figure('0004-66-1', '六个预先指定的局部剖面；星号为真值，圆圈／叉号为盲搜索产生的K=2候选中较小Q输出。超出显示范围的拟合已标注。')
    figure('0004-66-2', '40 ns第0对的网格密度核验；(a)13×9点，(b)25×17点，共用较细网格最低值作为色标零点。')
    code('''display(pd.DataFrame(search_diagnostics["profiles"]))
print("共用网格点代价最大差：",search_diagnostics["dense_grid_common_max_cost_error"])
print("加密网格最低值下降：",search_diagnostics["dense_grid_gain"])
print("boundary_minimum=True表示最低值位于本次扫描边界；不能据此宣称已找到全局最低值。")
''')
    md('''<a id="weak-search-oracles-20261026"></a>
## 67. 全部强弱组合的真值辅助对照

对latest批high_weak和weak_high各三个间隔、每点40对，共240对，保留三种不同层级：评分端取已有双事件候选（只知K）、从真实Q/t初始化局部联合拟合、固定真实时刻仅拟合Q。最后一种的Q初值仅来自观测单事件拟合总电荷的预设份额，没有真实Q初值。

真实初始化的时间搜索半径48 ns、最多200迭代；固定真实时刻时七个电荷分配初值、每个最多100迭代。它们的预算和真值信息均与盲流程不同，只用于定位问题，不能当作可实现的盲恢复率或理论上限。

所有表仍用完整严格标准：两真实Q≥500 pe、每个Q相对误差≤20%、每个t误差≤min(10 ns,间隔/4)、无额外输出。普通10 ns时间匹配另列，不能把匹配到两个时刻等同于完整恢复。
''')
    figure('0004-67-1', '全部240对的诊断结果；(a)先强后弱，(b)先弱后强。所有曲线均至少使用真实K，两个橙／紫色方法还使用真实时刻或初值。')
    code('''display(pd.DataFrame(search_diagnostics["summary"]))
search_costs=pd.DataFrame(search_diagnostics["cost_comparisons"])
display(search_costs.groupby(["condition","gap_ns"])[[
    "original_minus_truth_start","weak_search_minus_truth_start","exact_truth_minus_truth_start"]].agg(["median","min","max"]))
print("搜索代价差用于比较同一似然下的候选；真值点代价较高本身不证明模型失配。")
search_rank_audit=json.loads((search_output / "exploratory-rank-audit.json").read_text())
display(pd.DataFrame(search_rank_audit["summary"]))
print("补充的事后排序诊断：把已有真值辅助解加入候选，再按同一NLL取最小值；没有新增拟合或改变盲实验。它不是盲恢复上限。")
''')
    md('''<a id="weak-search-blind-20261026"></a>
## 68. 优化预算相同的两种盲搜索与重新标定的门槛

两种方法都保留原K=0/1候选和原最佳K=2候选，均只使用观测电压及同一电压产生的旧候选。`legacy_extra`继续使用15%/50%/85%电荷分配；`weak_extra`另加入1%/3%/5%/10%及90%/95%/97%/99%。两种方向都考虑，间隔仍为4/8/12/20/40/80/160/320 ns。

每个间隔／方向先按初始NLL取最好的份额，再从每个方向取四个起点；最多八次20迭代筛选、两次70迭代精修。两者额外非线性优化预算相同，弱份额方案有更多初始代价计算；实际迭代次数和耗时另报。原最佳候选始终保留，所以最低代价不会变差，但严格参数恢复仍可能变差。

扩大搜索会改变单事件的增峰证据分布，因此每个方案各自用相同280个独立原生单事件重新标定门槛：max(NLL₁−NLL₂,0)、加一／并列上尾秩p≤0.025，拟合Q分档1500/5000/12500 pe；第二个事件的拟合Q须≥500 pe。门槛不使用评价批，也没有按真实强弱顺序逐帧择优。

每个新方案共3480帧，其中280帧只标定门槛，两批各1600帧评价。进入新母样本确认的条件提前固定：latest批high_weak、40 ns至少增加10个百分点（40对时为4对）；older批40 ns不退步，两批200 ns不退步；两批weak_high、8 ns最多下降5个百分点，弱组合平均不下降、较强组合平均最多下降2个百分点。每组空／单额外输出≤5%，各单成功≥95%，匹配Q偏差中位绝对值≤5%。这是开发筛选标准，不是总体误报保证。
''')
    code('''search_overview=[]
for method,datasets in search_report["results"].items():
    for dataset,r in datasets.items():
        search_overview.append(dict(方法=method,开发批=dataset,弱组合平均_percent=100*r["weak_mean"],
            较强组合平均_percent=100*r["strong_mean"],控制通过=r["controls_pass"],
            强后弱_40ns_percent=100*r["high_weak"]["40"],强后弱_200ns_percent=100*r["high_weak"]["200"]))
display(pd.DataFrame(search_overview))
display(pd.DataFrame(search_report["eligibility"]))
search_selection_audit=json.loads((search_output / "candidate-selection-audit.json").read_text())
display(pd.DataFrame([r for r in search_selection_audit["summary"] if r["condition"]=="high_weak"]))
print("事后候选审计：known_K_success为评分端已知K的诊断；blind_success是原固定门槛选出的实际结果。")
print("真Q≥500定义评价目标，拟合Q≥500是另外一个既有筛选条件；本轮没有修改它们。")
''')
    figure('0004-68-1', '两批全部十种组合；上两排older，下两排latest，均为开发资料，采用各自重新标定的门槛。')
    figure('0004-68-2', '两批空窗／原生单事件控制；各组80例，依次为额外输出帧比例、严格成功率、匹配后Q偏差中位。')
    code('''search_table=pd.DataFrame(search_report["summary"])
display(search_table[search_table.kind!="pair"][["dataset_role","method","condition","success","eligible_frames","extra_frames","q_error_p16_p50_p84"]])
search_paired=json.loads((search_output / "paired-development.json").read_text())
display(pd.DataFrame([dict(开发批=r["dataset"],方法=r["method"],基线=r["baseline"],组合=r["group"],
    平均差_pp=100*r["difference"],母事件配对95区间_pp=[100*v for v in r["paired_block95"]]) for r in search_paired["aggregate"]]))
print("配对区间按母事件块分层重抽，不能消除已公开开发资料上的方法选择偏差。")
''')
    md('''<a id="weak-search-full-waveforms-20261026"></a>
## 69. 完整波形与额外计算量

以下仍为预先指定的latest批high_weak第0对、40和200 ns；横轴展示整个0–2048 ns。上排观测和总拟合，中排真实弱分量及选中K=2输出中较小Q的拟合分量，下排总残差。若未选中两个输出，中排拟合弱分量画为零；较小Q输出也未必正确对应真实弱事件。真实分量和虚线时刻仅用于展示。

初次完整运行的两个进程收到SIGTERM，原因未确定；逐帧检查点均保留，随后减少并发继续预定未完成任务，没有按结果选择重跑。所有最终未收敛拟合也计入指标，不选择性删除。
''')
    figure('0004-69-1', '完整波形对照；(a,b)总电压，(c,d)弱／较小Q输出分量，(e,f)总残差，全部拟合曲线来自盲选择。')
    code('''search_examples=json.loads((search_output / "examples.json").read_text())
display(pd.DataFrame([dict(间隔_ns=e["case"]["gap_ns"],方法=method,严格成功=r["success"],
    真t_ns=r["true_t"],拟合t_ns=r["fit"]["times_ns"],真Q_pe=r["true_q"],拟合Q_pe=r["fit"]["charge_pe"])
    for e in search_examples for method,r in e["methods"].items()]))
search_runtime=json.loads((search_output / "runtime-and-convergence.json").read_text())
display(pd.DataFrame(search_runtime["runtimes"]).T)
display(pd.DataFrame([dict(方法与开发批=k,选中候选未收敛帧=v,该批帧数=1600) for k,v in search_runtime["selected_not_converged"].items()]))
print("这里的运行时间只包含原候选缓存之后的新增计算，不能当作完整流程耗时。")
print("耗时含一致性复核，并来自不同并发阶段；不据此给出严格速度排名。")
''')
    chosen = report['chosen']
    conclusion = ('两个搜索候选均未同时满足预定的开发约束，本轮不升级默认方法，也不启动新ROOT独立确认。所有局部收益和退步都保留。'
                  if chosen is None else
                  f'预定规则选中{chosen}进入下一阶段。这里仍只有开发证据；须冻结候选和门槛，再在新光学母事件上独立确认。')
    latest = report['results']
    detail = '；'.join(f"{name}的latest先强后弱8/40/200 ns恢复率为" + '/'.join(f"{100*latest[name]['latest']['high_weak'][str(gap)]:.1f}%" for gap in (8,40,200)) for name in ('baseline','legacy_extra','weak_extra'))
    md('''<a id="weak-search-conclusion-20261026"></a>
## 70. 本轮采用结论与边界

'''+conclusion+'\n\n'+detail+'''。

诊断结果应分开回答三个问题：已有搜索是否漏掉代价更低的候选；降低代价是否同时改善逐事件Q/t；重新标定增峰门槛后是否仍保留收益并控制误拆。最低代价更低、优化器成功停止、总波形残差更小，均不能单独替代严格恢复率。

**40 ns失败的具体拆分：**最新批high_weak的原K=2候选仅3/40满足严格参数精度，另外37/40在已知K诊断中已经失败。这3个正确候选的两个拟合Q均≥500 pe，被拒绝来自增峰证据：上尾秩p分别约0.383、0.321和0.222，高于固定0.025要求。因此这3例不能靠单独降低电荷阈值补回；把证据门槛放宽到这些值也须重新评价误报，不能直接作为改进。补充弱起点后仍只有3/40个参数正确的K=2候选，最终盲恢复仍0/40。

**下一步待验证：**用当前拟合模型自己的均值和协方差生成匹配模型的模拟数据，与现有原生光子／共同微胞数据在相同Q、间隔和评价标准下比较。前者用于检查该统计模型内的信息与搜索限制，二者差异用于检查响应或噪声近似；两类结果都不能直接宣称实物探测器的绝对分辨极限。该对照尚未运行，本轮结论不预支其结果。

当前剖面只是有限真值邻域内、有限起点和预算的数值诊断。扫描边界、离散网格、有限四档源能量、固定τ／几何／噪声，以及高斯条件似然近似仍限制结论。这里不能给出普适最小时间间隔、free-τ最短窗、α/γ混合或随机高率K≥3性能。

全部[逐点结果](artifacts/alpha_pair_weak_search_20261026/summary.csv)、[开发报告](artifacts/alpha_pair_weak_search_20261026/report.json)、[诊断统计](artifacts/alpha_pair_weak_search_20261026/diagnostic-summary.json)和[配对分析](artifacts/alpha_pair_weak_search_20261026/paired-development.json)保留。
''')
    code('''search_checks=json.loads((search_record / "science-validation.json").read_text())
print("评分重放：",search_checks["scores_replayed"],"；汇总：",search_checks["summary_rows"])
print("诊断记录：",search_checks["diagnostic_records"],"；检查通过：",search_checks["passed"])
print("原输入保持：",search_checks["original_inputs_unchanged"],"；各搜索重新标定门槛：",search_checks["gate_recalibrated_for_each_search"])
''')
    for i, cell in enumerate(cells): cell['id'] = f'weak-search-20261026-{i:02d}'
    InteractiveShell.instance(); namespace = {}; executed = 0; cwd = os.getcwd()
    try:
        os.chdir(PROJECT / 'jupyter')
        for cell in cells:
            if cell.cell_type != 'code': continue
            with capture_output(stdout=True, stderr=True, display=True) as capture:
                exec(compile(cell.source, '<weak-search-notebook>', 'exec'), namespace)
            assert not capture.stderr, capture.stderr
            cell.outputs = []
            if capture.stdout: cell.outputs.append(nbformat.v4.new_output('stream', name='stdout', text=capture.stdout))
            for output in capture.outputs:
                cell.outputs.append(nbformat.v4.new_output('display_data', data=output.data, metadata=output.metadata))
            executed += 1; cell.execution_count = executed
    finally: os.chdir(cwd)
    notebook.cells.extend(cells); assert notebook.cells[:len(original)] == original
    anchors = re.findall(r'id="codex-figure-([^"]+)"', '\n'.join(c.source for c in notebook.cells if c.cell_type == 'markdown'))
    assert len(anchors) == len(set(anchors))
    for cell in cells:
        if cell.cell_type == 'markdown':
            assert not re.search(r'\]\([^)]*\.(?:pdf|ipynb|md)(?:[#?)][^)]*)?\)', cell.source, re.I)
    images = [output['data']['image/png'] for cell in cells if cell.cell_type == 'code'
              for output in cell.outputs if 'image/png' in output.get('data', {})]
    figures = load(RECORD / 'figure-delivery.json')['figures']
    assert len(images) == len(figures) == 6
    assert {hashlib.sha256(base64.b64decode(image)).hexdigest() for image in images} == set(figures.values())
    atomic_notebook(path, notebook, previous)
    save(RECORD / 'notebook-delivery.json', dict(passed=True, previous_cells=len(original), new_cells=len(cells),
        executed_cells=executed, history_preserved=True, figures=figures, notebook_sha256=digest(path)))
    print('WEAK_SEARCH_NOTEBOOK_DELIVERED', len(cells), executed, flush=True)


if __name__ == '__main__': deliver()
