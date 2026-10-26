"""Deliver weak-pulse diagnoses and covariance development to Notebook0004."""
from __future__ import annotations

import base64
import hashlib
import json
import os
import re
import shutil

import nbformat

from . import alpha_pair_weak_diagnostic_study as diagnostic
from . import alpha_pair_empirical_study as empirical
from .alpha_pair_conditional_notebook import atomic_notebook
from .alpha_pair_resolution_data import load, save, digest


def deliver():
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    assert load(empirical.RECORD/'science-validation.json')['passed']
    selection=load(empirical.EXPERIMENT/'development-selection.json')
    path=empirical.PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    previous=path.read_bytes();nb=nbformat.reads(previous.decode(),as_version=4)
    assert not any(c.get('id','').startswith('weak-noise-20261026-') for c in nb.cells)
    shutil.copyfile(path,empirical.RECORD/'notebook.before-delivery.ipynb')
    original=json.loads(json.dumps(nb.cells));cells=[]
    def md(text):cells.append(nbformat.v4.new_markdown_cell(text))
    def code(text):cells.append(nbformat.v4.new_code_cell(text))
    def figure(number,title,variable):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{title}**')
        code(f'display(Image(filename=str({variable} / "figures/{number}.png")))')
    md('''<a id="weak-pulse-diagnosis-20261026"></a>
## 53. 强脉冲后的弱事件为什么失败：先分开搜索和模型问题

**研究问题**：第52节的high_weak、40 ns组，盲恢复0/40，已知真实事件数后也仅4/40。这一轮先检查：找到正确起点是否足够；已知真实时刻后能否恢复电荷；强脉冲的实际形状涨落是否影响弱事件。

本轮没有新增G4仿真或神经网络训练。上一轮320个独立母ROOT的结果已经公开，此后只称为**开发数据**。预先选定十种电荷组合的8/40/200 ns、每点40对，以及全部80空窗和320原生单事件，共1600帧。跨条件／间隔复用母样本，1600帧不等于1600个独立光学母样本。

仍为纯α、固定95/351 ns双τ、2 ns采样、完整2.048 μs局部模拟电压和充分动态范围。Q为实际雪崩等效pe；时间是SiPM有效面首光子时刻。整对成功要求两个真实Q均≥500 pe，各Q误差≤20%，各t误差≤min(10 ns,Δ/4)，且无额外输出。时间间隔不是脉冲长度，也不是最高注入率的倒数。

依次增加真值信息，所有辅助流程都与盲结果分开：

| 流程 | 使用的真实信息 | 目的 |
| --- | --- | --- |
| 旧盲流程 | 无 | 保留实际可用基线 |
| 已知K | 真实阈上事件数 | 排除事件数选择因素；仍用原候选搜索 |
| 固定真实时刻 | K、t；Q初值仅来自观测总量的五种比例 | 查看时刻—电荷耦合 |
| 从真实参数启动 | K、Q/t初值；随后联合优化，时刻局部±24 ns | 查看目标函数是否把真值拉向错误解 |
| 减去真实另一分量 | K、真实强分量及其Q/t、弱事件真实初值 | 移除强事件的实际涨落；仅作oracle诊断 |

固定真实时刻与真值启动均另做OLS对照。Q依赖的条件模板是非线性的，固定t后仍通过非线性优化求Q。oracle弱分量仍保留真实邻居的条件响应及观测中的暗噪声／电子噪声；返回的另一事件直接使用真值，所以40/40不能当作盲恢复率或理论分辨上限。
''')
    code('''from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Image
weak_root = Path.cwd().resolve()
while not (weak_root / "pileup/experiments/alpha_pair_weak_diagnostics_20261026").is_dir():
    if weak_root.parent == weak_root:
        raise RuntimeError("请从项目根目录或jupyter目录运行")
    weak_root = weak_root.parent
weak_output = weak_root / "jupyter/artifacts/alpha_pair_weak_diagnostics_20261026"
noise_output = weak_root / "jupyter/artifacts/alpha_pair_empirical_covariance_20261026"
noise_experiment = weak_root / "pileup/experiments/alpha_pair_empirical_covariance_20261026"
weak_report = json.loads((weak_output / "report.json").read_text())
noise_report = json.loads((noise_output / "development-report.json").read_text())
weak_table = pd.DataFrame(weak_report["summary"])
rows = []
for method, group in weak_table[weak_table.condition=="high_weak"].groupby("method",sort=False):
    row = {"方法":method}
    for gap in (8.,40.,200.):
        r = group[group.gap_ns==gap].iloc[0]
        row[f"{gap:g} ns成功/对数"] = f"{r.success}/{r.eligible_frames}"
    rows.append(row)
display(pd.DataFrame(rows))
print("1600帧均已诊断：", weak_report["complete"], "；数据用途：已公开开发数据")''')
    figure('0004-53-1','强→弱组的诊断；(a)8 ns，(b)40 ns，(c)200 ns。仅第一根柱为盲结果，其余用真值；误差线为逐点Wilson95区间。','weak_output')
    md('''40 ns时，固定真实时刻得到13/40，从真实参数启动联合GLS得到8/40，仍有大量失败；减去真实强分量后，8、40、200 ns均为40/40。这支持继续研究强分量涨落、条件响应和时间—电荷耦合，不能把失败全部归于候选没有找到，也不能据此断言8 ns盲解一定可行。

40/40的逐点Wilson95下限约91.2%，0/40的上限约8.8%。它们不包含跨条件复用母样本、模型失配及选择方法带来的不确定性。
''')
    md('''<a id="weak-pulse-cost-20261026"></a>
## 54. 总波形拟合得像，弱事件参数仍可能错误

继续展示预先指定的high_weak第0对、40及200 ns，全横轴0–2048 ns。中间一行单独放大弱分量纵轴；原生真分量由共同微胞数字化的雪崩祖先归属构造，不能在实测中直接取得。

GLS诊断重放原盲流程每帧固定精度矩阵，全部候选代价的最大相对差约9.73×10⁻¹⁶。只有共享矩阵、使用完整原始波形的代价可以直接比较；OLS及oracle减去另一分量后的代价不与它混比。
''')
    figure('0004-54-1','预先指定的完整波形；(a,b)总信号，(c,d)弱分量，(e,f)完整残差。图中拟合曲线均有真值辅助。','weak_output')
    code('''weak_examples = json.loads((weak_output / "examples.json").read_text())
display(pd.DataFrame([dict(间隔_ns=e["case"]["gap_ns"],真Q_pe=e["true_q"],真t_ns=e["true_t"],
    已知K恢复Q_pe=e["diagnostics"]["known_count_previous"]["charge_pe"],
    已知K恢复t_ns=e["diagnostics"]["known_count_previous"]["times_ns"],
    真值处GLS代价=e["diagnostics"]["exact_truth_cost"],
    真值启动后GLS代价=e["diagnostics"]["truth_start_joint"]["cost"])
    for e in weak_examples]))
display(pd.DataFrame(json.loads((weak_output / "cost-diagnosis.json").read_text())["counts"]))''')
    figure('0004-54-2','真值启动后的弱事件误差；(a)8 ns使用对称对数纵轴，(b)40 ns，(c)200 ns。绿色区域仅表示弱事件Q/t容限，整对成功还要求强事件正确。颜色为拟合代价／真值处代价。','weak_output')
    md('''散点按拟合后的时间排序，把第二个输出与后来的真实弱事件比较；短间隔发生分量交换时仍保留，不能把图中的第二个输出预先认定为正确弱事件。整对评价使用一对一匹配。

40 ns组真值启动后失败32/40，这32例的代价都低于真值处代价；200 ns也有14例。优化随机数据后残差下降本身是正常现象，不能单凭这一点证明模型有偏。但与严重的Q/t误差、固定时刻和精确分量扣除的对照一起看，当前代价较小并不足以判定逐事件恢复正确。更多起点可以改善搜索，无法保证消除这类失败。
''')
    md('''<a id="empirical-pulse-covariance-20261026"></a>
## 55. 实现新候选：用孤立脉冲涨落估计相关噪声

前述诊断后，新增一个实际盲候选：保持条件均值响应，改变残差的协方差权重。这里没有释放双τ，也没有利用当前待测波形的真实分量。

从此前两批较强原生单事件和一批弱原生单事件取得560个母样本，按原生replicate奇偶分为280个协方差训练母样本和互不重叠的280个增峰门槛标定母样本。协方差训练的weak/low/mid/high数量为40/80/80/80；另一个320母样本、1600帧的已公开集合用于开发验证。三组光学母样本种子无交叉。

训练时，按真实Q和首光子时刻对齐“实际单脉冲减条件均值模板”的残差。读取原始20 μs记录，在−64到2112 ns、2 ns网格上作线性插值；整个相对网格都有实际观测，不把未观测尾部填零。每档残差中心化后提取经验相关模式，残差均值只保存作诊断，不加回模板。

令U为选定的残差模式、v为逐样点经验方差，用C=diag(max(v−rowSum(U²),0.01 mV²))+UUᵀ保证正定。试四种预先登记的方案：仅经验对角方差，16模式×0.5协方差份额，16模式×0.9份额，32模式×0.9份额。模式幅度乘份额平方根；“份额”不是相关系数ρ。

推断先只用观测电压做一个OLS单事件初估，按初估Q选档、按初估t平移噪声模式，并把该帧C固定给所有K/Q/t候选。原残差候选搜索和优化预算保持。小矩阵白化与稠密逆矩阵二次型已作数值对照。每个推断worker禁止读取ROOT、truth或case标签；真值仅在前面的训练、明确标注的诊断及事后评价中使用。

各方案再用独立的280个单事件标定同一(C1−C2)/C1经验上尾秩，门槛p≤0.025、Q分档边界1500/5000/12500 pe。另把旧AR模型权重用同样280例重标定，作为相同门槛标定样本量的对照。原先560例标定的基线也保留，不能把样本量变化的作用全归于新协方差。

**选用规则在计算前固定**：各空／单控制额外输出≤5%，各单严格成功≥95%，各单Q偏差中位绝对值≤5%；包含弱事件的15点等权平均必须超过两种旧盲基线，较强组合15点平均相对原观测权重最多下降2个百分点。所有帧为开发用途，控制阈值指样本点估计。

文献只借用Fixsen等2004的条件均值／协方差思想（DOI：[10.1016/j.nima.2003.11.313](https://doi.org/10.1016/j.nima.2003.11.313)，共享库AI精读稿待审；TES模拟，每能量5000次，原PDF2式1–2及PDF3）。本文不含这里的低秩截断、正则化或GAGG双事件验证；这些是项目实现。本轮复用既有笔记，未重新精读PDF。
''')
    figure('0004-55-1','经验协方差训练；(a)四档残差标准差，(b)高Q档16模式、0.9份额的相关矩阵。图中的训练对齐使用真值，盲推断对齐使用观测初估。','noise_output')
    md('''<a id="empirical-covariance-results-20261026"></a>
## 56. 四种新权重的开发对照：恢复与误拆一起看

每种新权重实际拟合1880帧：280个门槛标定单事件加1600个开发验证帧；四种合计7520次盲帧分析。下面全部成功率来自后1600帧，不把门槛标定数据混进验证分母。每个双事件点生成40对，阈下数量在完整结果中保留。

两类平均各覆盖五种组合×三个间隔，共15点，与第52节七间隔的35点平均覆盖不同；只能在本节的相同点集上作方法比较。
''')
    code('''noise_rows = []
for name, result in noise_report["results"].items():
    noise_rows.append(dict(方法=name,弱组合平均_percent=100*result["weak_mean"],
        较强组合平均_percent=100*result["strong_mean"],控制通过=result["controls_pass"],
        满足选用条件=result["eligible"]))
display(pd.DataFrame(noise_rows))
print("同点集旧基线：",noise_report["baselines"])
print("选中方案：",noise_report["chosen"],"；独立性能声明：",noise_report["independent_performance"])
noise_points = []
for method in noise_report["results"]:
    row = {"方法":method}
    for gap in (8.,40.,200.):
        r = next(r for r in noise_report["summary"] if r["method"]==method and r["condition"]=="high_weak" and r["gap_ns"]==gap)
        row[f"strong→weak {gap:g} ns"] = f'{r["success"]}/{r["eligible_frames"]}'
    noise_points.append(row)
display(pd.DataFrame(noise_points))''')
    figure('0004-56-1','十种组合的完整开发恢复曲线；旧两种盲流程、旧AR同样本量门槛对照及四种新权重均保留。','noise_output')
    figure('0004-56-2','空／单事件控制，每组80例；(a)额外输出帧比例，(b)严格成功率，(c)匹配Q误差中位。虚线为预定选用条件，完整报告另存逐点Wilson区间。','noise_output')
    code('''noise_controls = []
for method,result in noise_report["results"].items():
    for condition,r in result["controls"].items():
        noise_controls.append(dict(方法=method,条件=condition,
            成功=f'{r["success"]}/{r["frames"]}',额外输出=f'{r["extra"]}/{r["frames"]}',
            Q误差中位_percent=None if r["median_q_bias"] is None else 100*r["median_q_bias"]))
display(pd.DataFrame(noise_controls))''')
    if selection['chosen'] is None:
        conclusion='四种新权重没有同时达到预先约定的选用条件，本轮不把它们升级为默认方法，也未启动新的独立测试。结果作为完整的开发否定性证据保留；局部点的改善不能代替整体恢复与误拆约束。'
    else:
        conclusion=f"开发选择了{selection['chosen']}。这是已公开数据上的选择结果；只有在冻结该方案后用新母样本检验，才能报告独立改善。当前章节不把开发收益称为独立性能。"
    md('''<a id="weak-noise-conclusion-20261026"></a>
## 57. 本轮结论、限制与下一步

'''+conclusion+'''

目前证据支持：强→弱失败包含时间—电荷耦合，以及强分量涨落／条件模型的影响；单纯继续增加候选起点不够。真值扣除的高成功率只定位问题，不能转化为可实现的盲算法性能。

经验协方差仍受四档离散Q、每档40或80个母样本、真值对齐后的线性插值以及仅用单脉冲初估的限制。两个事件共同推进SiPM微胞时，联合涨落不一定等于孤立强事件涨落；四种份额／截断也没有覆盖所有正则化选择。

下一步应先在开发数据上检验**双事件条件均值及联合协方差**：以多个独立母样本的已知组合建立开发基准，对比真实强弱分量的均值、方差和互相关；再检验仅由观测候选Q/t推定这些量的盲实现。如果C随候选变化，高斯近似的目标还须包含log det C，不能只比较加权残差。可以继续使用模拟真分量训练／标定，但推断不能直接读取它们。候选通过后再冻结并生成新母样本，随后才回到随机高率、K≥3及α/γ比较。

完整[诊断结果](artifacts/alpha_pair_weak_diagnostics_20261026/report.json)、[逐点开发指标](artifacts/alpha_pair_empirical_covariance_20261026/development-summary.csv)、[开发选择](artifacts/alpha_pair_empirical_covariance_20261026/development-report.json)及[同门槛样本量配对变化](artifacts/alpha_pair_empirical_covariance_20261026/paired-development.json)均保留。本轮没有新ROOT，Run_Catalog沿用原有登记；原实验、笔记旧计算、图像和编号不变。
''')
    for i,cell in enumerate(cells):cell['id']=f'weak-noise-20261026-{i:02d}'
    InteractiveShell.instance();namespace={};executed=0;old_cwd=os.getcwd()
    try:
        os.chdir(empirical.PROJECT/'jupyter')
        for cell in cells:
            if cell.cell_type!='code':continue
            with capture_output(stdout=True,stderr=True,display=True) as captured:
                exec(compile(cell.source,'<weak-noise-notebook>','exec'),namespace)
            assert not captured.stderr,captured.stderr
            cell.outputs=[]
            if captured.stdout:cell.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=captured.stdout))
            for out in captured.outputs:cell.outputs.append(nbformat.v4.new_output('display_data',data=out.data,metadata=out.metadata))
            executed+=1;cell.execution_count=executed
    finally:os.chdir(old_cwd)
    nb.cells.extend(cells);assert nb.cells[:len(original)]==original
    anchors=re.findall(r'id="codex-figure-([^"]+)"','\n'.join(c.source for c in nb.cells if c.cell_type=='markdown'))
    assert len(anchors)==len(set(anchors))
    for c in cells:
        if c.cell_type=='markdown':assert not re.search(r'\]\([^)]*\.(?:pdf|ipynb|md)(?:[#?)][^)]*)?\)',c.source,re.I)
    atomic_notebook(path,nb,previous)
    images=[out['data']['image/png'] for c in cells if c.cell_type=='code' for out in c.outputs
        if out.get('output_type')=='display_data' and 'image/png' in out.get('data',{})]
    figures={str(p.relative_to(empirical.PROJECT)):digest(p) for folder in (diagnostic.OUTPUT,empirical.OUTPUT)
        for p in (folder/'figures').glob('*.png')}
    assert len(images)==len(figures)==6
    assert {hashlib.sha256(base64.b64decode(v)).hexdigest() for v in images}==set(figures.values())
    save(empirical.RECORD/'notebook-delivery.json',dict(passed=True,previous_cells=len(original),new_cells=len(cells),
        executed_cells=executed,history_preserved=True,figures=figures,notebook_sha256=digest(path)))
    print('WEAK_NOISE_NOTEBOOK_DELIVERED',len(cells),'cells',executed,'executed',flush=True)


if __name__=='__main__':deliver()
