"""Executed Notebook 0004 delivery for blind charge-allocation comparisons."""
from __future__ import annotations

import base64
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile

import nbformat

from .alpha_pair_allocation_study import PROJECT,OUTPUT,RECORD
from .alpha_pair_information_study import load,save,digest


def deliver():
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    assert load(RECORD/'science-validation.json')['passed']
    assert load(RECORD/'figure-delivery.json')['passed']
    report=load(OUTPUT/'report.json')
    path=PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    previous=path.read_bytes();notebook=nbformat.reads(previous.decode(),as_version=4)
    before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4)
    assert len(before.cells)==410 and notebook.cells[:410]==before.cells
    if len(notebook.cells)>410:
        assert digest(path)==load(RECORD/'notebook-delivery.json')['notebook_sha256']
        assert all(c.get('id','').startswith('allocation-20261027-') for c in notebook.cells[410:])
        notebook.cells=notebook.cells[:410]
    shutil.copyfile(path,RECORD/'notebook.before-delivery.ipynb')
    original=json.loads(json.dumps(notebook.cells));cells=[]
    def md(text):cells.append(nbformat.v4.new_markdown_cell(text))
    def code(text):cells.append(nbformat.v4.new_code_cell(text))
    def figure(number,title):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{title}**')
        code(f'display(Image(filename=str(allocation_output / "figures/{number}.png")))')

    md(r'''<a id="allocation-design-20261027"></a>
## 82. 实际实现总电荷与分配比例的盲搜索

第81节提出的优化在本轮执行。问题是：原来的$\ln Q_1,\ln Q_2$坐标是否让有限搜索难以沿强弱电荷反相关方向移动？此外，对每个候选时间／份额先优化总Q，是否能找到更好的起点？这两项分开比较，不把增加计算量的收益都归给换坐标。

**四种对照与原基线：**

| 代码名 | 电荷坐标 | 候选与优化预算 |
| --- | --- | --- |
| baseline | 原$\ln Q_1,\ln Q_2$ | 已有冻结结果 |
| component_local | 原分量坐标 | 对原最佳K=2候选再优化至多70次迭代 |
| allocation_local | 总量／分配坐标 | 同一候选、同样至多70次迭代 |
| component_profile | 原分量坐标 | 加入共同的总Q剖面起点，再筛选与精修 |
| allocation_profile | 总量／分配坐标 | 与上行相同起点及迭代上限 |

总量坐标为$S=Q_1+Q_2$和$r$：$Q_1=L(S)+[U(S)-L(S)]r$，$Q_2=S-Q_1$，其中$L=\max(50,S-50000)$、$U=\min(50000,S-50)$、$0\le r\le1$。因此完整保留每个Q在50–50000 pe内的原范围。内部用$100\ln S,10r$作数值缩放；时刻仍为$t_1/10,t_2/10$，每次局部时刻范围为起点±24 ns。Q的单位为pe，实际对数采用相应数值。

这里的Q范围是拟合参数域，与ADC位数或削顶不同；采集动态范围足够、完整电压可用的前提保持。

这是对同一物理参数空间换坐标和缩放；均值、协方差及包含logdet的目标完全相同，最大似然目标不添加变量变换Jacobian。有限差分在物理Q中的实际步长会随坐标改变，故比较的是有限数值优化实现，不是宣称纯代数变换增加了信息。两个坐标使用同一L-BFGS-B设置，实际函数调用数和耗时仍可能不同，另行报告。

**共同剖面起点：**只以观测单脉冲拟合的时刻为锚，在前／后4/8/12/20/40/80/160/320 ns放第二候选。每个时间对扫描第一事件份额0.01/0.03/0.05/0.10/0.15/0.50/0.85/0.90/0.95/0.97/0.99；在完整可行电荷范围数值优化lnS（最多24次、容差10⁻⁴），并保留原未优化值。各间隔／方向取代价最低份额，每方向前4个，共最多8个起点。各做20次迭代筛选，最好的两个再做70次精修，始终保留原最佳K=2候选和local结果。这里的“剖面”是有限数值最小化，不保证全局最优。

280个独立原生单事件只标定门槛；760帧已公开原生数据与其同Q/t的**高斯draw0共760帧**用于开发比较。合计1800帧，每帧产生四种变体候选，共7200组；高斯draw1本轮未运行，不能将40对的表直接当上一节80对的重复结果。没有新母ROOT、网络训练或电荷先验，推断worker禁止读取truth、case和ROOT。相同物理输入与已公开母样本复用均保留，不能称新独立确认。

固定95/351 ns发光双τ、2 ns×1024点完整模拟电压及充分动态范围、局部未知K≤2。原生Q为实现雪崩等效pe，高斯Q为相同数值的条件参数；t为有效面首光子，不作沉积能量。真值只在评分、已知K诊断和图像对照中使用。

主表为每种新搜索在相同280个原生单事件上重新计算增峰参考分布；分档仍是拟合单事件Q的1500/5000/12500 pe，上尾加一并保留并列的秩p≤0.025。单事件原模型与Q≥500保护保持。另保存直接沿用原数值门槛的敏感性结果，避免扩大搜索后误把更容易增峰当成可靠改善。

文献依据沿用Fixsen等2004的条件均值／协方差思想（DOI：10.1016/j.nima.2003.11.313，PDF2式1–2、PDF3；TES模型0.5–5 eV，每能量5000次）及Liddick等2012的检测与参数精度分开评价（DOI：10.1016/j.nima.2011.12.034，PDF3–8、图3–13；DSSD、40 MSPS、人工合成双脉冲）。据既有精读笔记，两篇科学内容保持AI草稿待人工审核；本轮没有重查PDF。当前总Q参数化、数值缩放和剖面搜索为项目方案，不冒充两篇原算法或已有GAGG性能结论。
''')
    code('''from pathlib import Path
import json
import numpy as np
import pandas as pd
from IPython.display import display, Image
allocation_root = Path.cwd().resolve()
while not (allocation_root / "pileup/experiments/alpha_pair_charge_allocation_20261027").is_dir():
    if allocation_root.parent == allocation_root:
        raise RuntimeError("请从项目根目录或jupyter目录运行")
    allocation_root = allocation_root.parent
allocation_output = allocation_root / "jupyter/artifacts/alpha_pair_charge_allocation_20261027"
allocation_record = allocation_root / "jupyter/run_records/alpha_pair_charge_allocation_20261027"
allocation_report = json.loads((allocation_output / "report.json").read_text())
allocation_summary = pd.DataFrame(allocation_report["summary"])
print("原生760帧、高斯draw0 760帧；另280个原生单事件标定门槛。")
print("真值不进入推断；新母ROOT：", allocation_report["new_roots"])
display(allocation_summary[allocation_summary.kind == "pair"][["domain", "method", "condition", "gap_ns", "frames", "success", "K2_outputs", "both_times10", "both_times_strict", "extra_frames"]])
''')
    figure('0004-82-1','五种流程的严格整对恢复率；(a–c)原生中／中、先强后弱、先弱后强，(d–f)高斯draw0对应组合；各点40对，同规则分别校准增峰门槛。')

    md('''<a id="allocation-controls-20261027"></a>
## 83. 空／单控制、门槛敏感性与计算量

严格整对成功要求两真实Q≥500 pe，各Q误差≤20%，各时刻误差≤min(10 ns,真实间隔/4)，且没有额外输出。单／空和双事件都保留所有失败与未收敛结果。误拆／误并为原时间匹配图的操作性标志；额外输出也可包含时刻未匹配的候选，不一定是完全无信号位置的假脉冲。

各空／单组、每种输入均80帧。预定开发控制为额外输出≤5%，各单事件严格成功≥95%、匹配后Q偏差中位绝对值≤5%。Q分位及sigma68=(p84−p16)/2是匹配条件下的电荷误差宽度，没有沉积能标，不能写成能谱FWHM。

计算量为已有缓存基线以外的增量：两个profile变体各自独立使用时都需加共同剖面耗时，本次实际只计算一次共同剖面。local结果被各自profile复用，不重复记作新拟合。四进程条件下逐帧wall time不能当作实时吞吐量；“相同迭代上限”也不等于相同函数调用数。
''')
    figure('0004-83-1','空／单控制；(a–c)原生，(d–f)高斯draw0；三列为额外输出、严格成功、匹配后Q偏差中位。虚线为预定开发控制要求。')
    code('''display(allocation_summary[allocation_summary.kind != "pair"][["domain", "method", "condition", "frames", "success", "extra_frames", "q_error_p16_p50_p84", "q_error_sigma68", "not_converged"]])
allocation_original_gate = pd.DataFrame(allocation_report["original_gate_summary"])
display(allocation_original_gate[(allocation_original_gate.condition == "high_weak") | (allocation_original_gate.kind != "pair")][["domain", "method", "gate", "condition", "gap_ns", "frames", "success", "extra_frames"]])
print("上面第二表直接沿用原门槛数值，仅作敏感性对照；主表采用同规则重新校准后的门槛。")
''')
    figure('0004-83-2','(a,d)先强后弱的K=2候选代价下降比例；(b,e)强制K=2后的严格成功；(c,f)各方法平均增量耗时，上排原生、下排高斯。C为分量坐标，S为总量／分配。')
    code('''display(pd.DataFrame(allocation_report["timing"]))
allocation_diagnostic = pd.DataFrame(allocation_report["diagnostics"])
display(allocation_diagnostic[allocation_diagnostic.condition == "high_weak"][["domain", "method", "gap_ns", "frames", "positive_cost_gain", "known_K_success", "below500", "gap_under2ns", "not_converged", "cost_gain_p16_p50_p84", "gained_correct", "lost_correct"]])
print("强制K=2属于评分诊断；below500表示最低代价双候选含阈下Q，不能据其候选成功数冒充盲输出。gap_under2ns只是描述，不是新增拒绝规则。")
''')

    md('''<a id="allocation-errors-20261027"></a>
## 84. 电荷分配误差与成对比较

下图画已有K=2候选中弱事件与强事件的绝对Q误差。按拟合时刻先后与真事件对应，包含时刻错误、Q低于500、增峰门槛拒绝及未收敛的候选。它用于诊断电荷是否沿近似恒定总Q方向转移，不能作为盲恢复率。虚线ΔQ弱+ΔQ强=0表示总电荷误差抵消；横纵轴保留所有点，不裁掉极端失败。

成对差值另按母事件块重抽2000次，0–19与20–39两层分别重抽20块，同一组权重用于全部方法、条件、间隔和两类输入。中／中40对涉及40块，强弱各40对涉及前20块、每块2对；因此不能把所有帧独立重抽。区间描述已公开开发资料下的差异，未包含完整物理／模型不确定性；全成／全败点的[0,0]区间只是样本退化。
''')
    figure('0004-84-1','40 ns先强后弱的强制K=2电荷误差；(a–c)原生的基线、分量profile、总量profile，(d–f)高斯对应方案。颜色为该K=2候选严格正确／失败，不等于被盲门槛接受。')
    code('''display(allocation_diagnostic[(allocation_diagnostic.condition == "high_weak") & (allocation_diagnostic.gap_ns == 40)][["domain", "method", "known_K_success", "weak_q_error_p16_p50_p84", "total_q_error_p16_p50_p84", "weak_time_error_p16_p50_p84_ns", "fitted_gap_p16_p50_p84_ns"]])
allocation_paired = pd.DataFrame(allocation_report["paired"])
display(allocation_paired[(allocation_paired.condition == "high_weak") & allocation_paired.method.isin(["allocation_local", "allocation_profile"])])
print("difference与paired_block95以比例表示，乘100为百分点；gained/lost保留逐对变化。")
''')

    md('''<a id="allocation-waveforms-20261027"></a>
## 85. 完整波形和各分量

在计算前指定原生先强后弱40/200 ns的第0对作示例，不按改善程度挑选。三列分别为完整含噪总电压、较大Q输出、较小Q输出，全部显示0–2048 ns。真实原生雪崩分量与真参数条件均值分别标注，两者不是同一条波形。

曲线采用各方法经各自校准门槛最终选择的输出。若只选K=1，弱分量图中的该方法画零线，不把评分端强制K=2的候选偷偷当作解堆积结果。分量的“较大／较小”仅按盲输出Q排序，真实强弱对应是否正确还需看评分表。
''')
    figure('0004-85-1','预定原生第0对、全长0–2048 ns；(a–c)40 ns，(d–f)200 ns。总电压、强分量和弱分量分别显示，比较原基线与allocation_profile实际盲输出。')
    code('''allocation_examples = json.loads((allocation_output / "examples.json").read_text())
display(pd.DataFrame([dict(间隔_ns=e["gap_ns"], 方法=m, K=len(v["fit"]["times_ns"]), 拟合Q_pe=v["fit"]["charge_pe"], 拟合时刻_ns=v["fit"]["times_ns"], 严格成功=v["score"]["success"], 额外输出=v["score"]["extra"]) for e in allocation_examples for m,v in e["methods"].items()]))
''')

    lookup={(r['domain'],r['method'],r['condition'],r['gap_ns']):r for r in report['summary']}
    base40=lookup['native','baseline','high_weak',40.]['success']
    total40=lookup['native','allocation_profile','high_weak',40.]['success']
    base200=lookup['native','baseline','high_weak',200.]['success']
    total200=lookup['native','allocation_profile','high_weak',200.]['success']
    chosen=report['chosen']
    decision=(f'预定开发要求下，暂选{chosen}作为后续候选；这不是独立确认或默认算法升级。' if chosen else '四种变体均未同时通过预定开发要求，本轮不选择新默认方案。')
    md(f'''<a id="allocation-conclusions-20261027"></a>
## 86. 本轮取舍与后续方向

实际完成四种盲搜索对照。原生先强后弱40 ns严格恢复：原基线{base40}/40，allocation_profile为{total40}/40；200 ns为{base200}/40→{total200}/40。其他变体、反向强弱、中／中、控制和高斯draw0结果均在前表完整保留。

{decision}

预定选择要求同时检查原生空／单控制，先强后弱40 ns至少提高10个百分点、200 ns不下降，先弱后强8 ns降幅不超过5个百分点，中／中各点≥95%，以及六个强弱条件平均不下降。数据已公开且只有三种组合和三档间隔；即使通过，也需连续Q、幅比外推、模型失配和冻结后的新母样本验证。

本轮严格区分似然代价和恢复准确度：更低代价可以对应另一组强弱电荷、错误时刻或阈下分量；换坐标不会创造新的波形信息。参数域、响应、协方差和门槛规则固定，所有未收敛输出计入。这项对照用于判断有限搜索的作用，不是任意算法的理论上界。

后续若收益不足，应转向检验估计目标、响应／噪声模型及训练先验，而不是仅增加同一似然的起点。若引入先验或学习型电荷估计，先验只能从训练母样本获得，须用连续Q检验，不能把当前四档真实能量当作测试答案。若已有收益满足开发要求，再扩大条件并冻结独立确认；两种情况下均不直接推广到α/γ、free-τ或随机高率K≥3。

完整[主表与门槛敏感性](artifacts/alpha_pair_charge_allocation_20261027/summary.csv)、[候选诊断](artifacts/alpha_pair_charge_allocation_20261027/candidate-summary.csv)、[成对区间](artifacts/alpha_pair_charge_allocation_20261027/paired-comparison.json)与[报告](artifacts/alpha_pair_charge_allocation_20261027/report.json)保留。旧410个单元和图像不变，新图号固定。
''')
    code('''display(pd.DataFrame(allocation_report["decisions"]))
allocation_checks = json.loads((allocation_record / "science-validation.json").read_text())
print("通过：", allocation_checks["passed"], "；输入帧：", allocation_checks["frames"], "；变体候选组：", allocation_checks["variant_candidate_sets"])
print("主评分：", allocation_checks["primary_scores"], "；含原门槛敏感性：", allocation_checks["all_scores"], "；强制K=2诊断：", allocation_checks["forced_K2_scores"])
print("基线与既有缓存一致：", allocation_checks["baseline_matches_previous"], "；最大代价重放相对差：", allocation_checks["maximum_cost_replay_error"])
print("标量剖面点：", allocation_checks["scalar_profile_points"], "；未获优化器成功状态：", allocation_checks["scalar_profile_not_converged"], "，全部保留。")
''')
    for i,c in enumerate(cells):c['id']=f'allocation-20261027-{i:02d}'
    InteractiveShell.instance();namespace={};executed=0;cwd=os.getcwd()
    try:
        os.chdir(PROJECT/'jupyter')
        for cell in cells:
            if cell.cell_type!='code':continue
            with capture_output(stdout=True,stderr=True,display=True) as capture:
                exec(compile(cell.source,'<allocation-notebook>','exec'),namespace)
            assert not capture.stderr,capture.stderr
            cell.outputs=[]
            if capture.stdout:cell.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=capture.stdout))
            for out in capture.outputs:cell.outputs.append(nbformat.v4.new_output('display_data',data=out.data,metadata=out.metadata))
            executed+=1;cell.execution_count=executed
    finally:os.chdir(cwd)
    notebook.cells.extend(cells);assert notebook.cells[:410]==original
    anchors=re.findall(r'id="codex-figure-([^\"]+)"','\n'.join(c.source for c in notebook.cells if c.cell_type=='markdown'))
    assert len(anchors)==len(set(anchors))
    for c in cells:
        if c.cell_type=='markdown':assert not re.search(r'\]\([^)]*\.(?:pdf|ipynb|md)(?:[#?)][^)]*)?\)',c.source,re.I)
    images=[o['data']['image/png'] for c in cells if c.cell_type=='code' for o in c.outputs if 'image/png' in o.get('data',{})]
    expected=load(RECORD/'figure-delivery.json')['figures'];assert len(images)==len(expected)==5
    assert {hashlib.sha256(base64.b64decode(s)).hexdigest() for s in images}==set(expected.values())
    nbformat.validate(notebook);encoded=(nbformat.writes(notebook)+'\n').encode()
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=path.stem+'.',suffix='.partial',delete=False) as stream:
        stream.write(encoded);stream.flush();os.fsync(stream.fileno());temporary=Path(stream.name)
    assert temporary.read_bytes()==encoded
    nbformat.validate(nbformat.reads(encoded.decode(),as_version=4))
    assert path.read_bytes()==previous,'Notebook changed during delivery; retain user edits'
    temporary.replace(path);assert path.read_bytes()==encoded
    save(RECORD/'notebook-delivery.json',dict(passed=True,previous_cells=410,new_cells=len(cells),executed_cells=executed,
        history_preserved=True,figures=expected,notebook_sha256=digest(path),source_sha256=digest(__file__)))
    print('ALLOCATION_NOTEBOOK_DELIVERED',len(cells),executed,flush=True)


if __name__=='__main__':deliver()
