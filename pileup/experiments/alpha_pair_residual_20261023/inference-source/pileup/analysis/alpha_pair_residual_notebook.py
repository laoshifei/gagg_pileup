"""Append the residual-search study and execute only its new display cells."""
from __future__ import annotations

import base64
import hashlib
import json
import os
import re

import nbformat

from . import alpha_pair_residual_study as study
from .alpha_pair_conditional_notebook import atomic_notebook


def interpretation():
    """Describe all predeclared outcomes without selecting a new procedure."""
    report=study.load(study.OUTPUT/'test-report.json');summary=report['summary']
    def row(method,condition,gap=None):
        return next(r for r in summary if r['method']==method and r['condition']==condition and r['gap_ns']==gap)
    old_low=row('selective_1','low');new_low=row('residual_regularized','low')
    old_mid=row('selective_1','mid_mid',8.);new_mid=row('residual_regularized','mid_mid',8.)
    weak=[]
    for gap in (40.,200.):
        old=row('selective_1','high_low',gap);new=row('residual_regularized','high_low',gap)
        weak.append(f"{gap:g} ns：{old['success']}/{old['frames']} → {new['success']}/{new['frames']}")
    controls=[row('residual_regularized',name) for name in ('all','low','mid','high')]
    false='、'.join(f"{name}{r['extra_frames']}/{r['frames']}" for name,r in zip(('空窗','低Q','中Q','高Q'),controls))
    budget='全部通过' if all(r['budget_pass'] for r in controls) else '未全部通过'
    return f'''**本批独立结果**：低Q原生单事件的匹配后电荷偏差中位数，原GLS为{100*old_low['q_error_p16_p50_p84'][1]:.2f}%，完整新流程为{100*new_low['q_error_p16_p50_p84'][1]:.2f}%；严格整帧成功为{old_low['success']}/{old_low['frames']} → {new_low['success']}/{new_low['frames']}。中／中组8 ns整对成功为{old_mid['success']}/{old_mid['frames']} → {new_mid['success']}/{new_mid['frames']}。先强后弱的整对成功为{'；'.join(weak)}。这里的箭头均比较同一批新测试，不拼接上一批测试数值。

**误报代价**：完整新流程各控制的额外输出为{false}，5%点估计预算{budget}。成功率与误报共同保留；即使样本点估计通过，也不意味着总体误报率已被严格认证低于5%。下面列出弱第二脉冲的事件数分布和全部关键对照，不能仅根据某一例波形外观判断恢复正确。'''


def deliver():
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    assert study.load(study.RECORD/'science-validation.json')['passed']
    path=study.PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    before=path.read_bytes();nb=nbformat.reads(before.decode(),as_version=4)
    assert not any(c.get('id','').startswith('residual-20261023-') for c in nb.cells)
    original=json.loads(json.dumps(nb.cells));cells=[]
    def md(text):cells.append(nbformat.v4.new_markdown_cell(text))
    def code(text):cells.append(nbformat.v4.new_code_cell(text))
    def figure(number,title):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{title}**')
        code(f'display(Image(filename=str(residual_output / "figures/{number}.png")))')
    md('''<a id="residual-search-20261023"></a>
## 45. 相关权重正则化与残差新增候选

第44节留下两类问题：低电荷单事件被经验GLS系统性低估；网络若只输出一个事件，保持K的重拟合无法补出第二个事件。本轮分别检验权重与候选搜索，沿用纯α、固定95/351 ns双τ、2 ns采样、2.048 μs局部电压和充分动态范围。Q仍指实际雪崩等效pe，不直接作为源能量或沉积能量。

均值响应和两个网络权重均沿用上一轮；没有重新训练网络。先在原响应训练母样本上比较GLS／OLS精度矩阵混合及ρ缩减，再用上一轮已公开的独立样本作为本轮开发选择数据：240个母ROOT、每条件／间隔前20对和全部空／单控制，共1320帧。它们不再称为本轮独立测试。响应标定的训练母样本与这些选择母样本分离。

**权重调整**：从原经验AR(1)的ρ≈0.9783出发，比较保留原值、缩减至0.9及0.8。观测波形给出的平滑先验在每帧固定，所有K/t/Q候选共用相同精度矩阵；将零模型代价统一缩放至约10⁶只用于保持优化器的数值尺度，选择K时再除以零模型代价。这是经验GLS正则化，没有重新测量器件的相关系数，也不是精确微观似然。

**新增候选**：拟合单事件后，用固定线性字典扫描残差；先投影掉现有事件的Q和时刻局部切向，减少把已有脉冲的小幅修正误作新事件。残差局部峰提供最多4个起点，另在离散间隔与电荷分配上筛出4个分裂起点。它们只供初始化，随后都用同一非线性条件响应联合拟合Q/t。新起点先作最多20次函数评估，再精修最好的两个；数值Jacobian调用另计。局部搜索半径24 ns，最多仍为K=2，尚非全局搜索或任意多事件重建。

未知K由归一残差代价＋每事件惩罚决定。选择时要求各空／单控制额外输出点估计≤5%，各单事件成功率≥95%，各档单事件Q相对误差中位数绝对值≤5%；在满足条件的方案中提高五种组合等权的整对成功率。网络概率不作为强制事件数，真K、Q、t、源能量或注入率不进入推断。

上述约束用于选择主方案。原权重等消融即使未达到电荷偏差目标也保留，并按开发误报预算选择其惩罚；各版本的惩罚可不同。它们比较经过同样开发程序的完整流程，不将最终差异全部归因于某一个ρ值。

继承的波形检测器若只提供零模型，本轮跳过新增候选并返回空输出；这是观测波形的预筛，不是真值判断。同K网络对照如没有可用的同K精修则保留原网络输出。尚未检验阈附近物理弱事件因此不能忽略这项预筛的漏检可能。

文献层面仅复用已有精读笔记：Fixsen等（TES模拟，DOI 10.1016/j.nima.2003.11.313）支持分别检查条件均值与非平稳权重；Zhang等（LaBr₃／PMT、500 MHz，DOI 10.1016/j.nima.2024.169273）采用幅度相关模板和残差增峰，并提醒模型失配可产生幽灵脉冲。两篇科学结论均保持AI待审，本轮没有重读PDF；当前ρ缩减、切向投影和门槛均为本项目实现，不搬用论文的ns或计数率数字。
''')
    code('''from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Image
residual_root = Path.cwd().resolve()
while not (residual_root / "pileup/experiments/alpha_pair_residual_20261023").is_dir():
    if residual_root.parent == residual_root:
        raise RuntimeError("请从项目根目录或其jupyter目录运行")
    residual_root = residual_root.parent
residual_output = residual_root / "jupyter/artifacts/alpha_pair_residual_20261023"
residual_report = json.loads((residual_output / "test-report.json").read_text())
residual_summary = pd.DataFrame(residual_report["summary"])
residual_selection = residual_report["residual_selection"]
display(pd.DataFrame([{k:v for k,v in residual_selection["chosen"].items() if k!="controls"}]))
display(pd.DataFrame(residual_selection["chosen"]["controls"]).T)''')
    figure('0004-45-1','开发集的权重与K选择；(a)低电荷单事件偏差，(b)双事件恢复与额外输出取舍。')
    md('''<a id="residual-independent-20261023"></a>
## 46. 冻结后重新生成母样本的检验

算法、阈值、源码、权重及配置冻结后，新生成2.5／5.486／8 MeV中央入射单α母ROOT各80个。双事件由新光子共同推进SiPM微胞，未缩放光子或叠加独立电压；另保留原生单事件与空窗。五种组合、4/8/12/20/30/40/60/80/120/200 ns，每点40对，另有80空窗及240单事件，总2320帧。各点复用母样本，不能把跨间隔帧数当作独立曝光数。

严格整对成功要求每个Q误差≤20%，每个时刻误差≤min(10 ns, Δ/4)，一对一恢复且无额外输出；真Q门槛500 pe。空／单控制另列。比较14种处理流程，两个网络仍为同架构的已有权重种子，不能称为14个新机器学习模型。

t定义为进入SiPM有效面的首光子时刻，Δ是两个这样的到达时刻之差，不是源发射时刻或波形峰顶的间隔。母光子按该时刻对齐后平移以构造受控双事件；真实事件的发射时刻分辨率尚未由本实验评价。

`regularized_1/2`保持各网络K；`regularized_count`重新比较既有候选的K；`residual_regularized`同时用波形与网络起点并新增候选；`residual_waveform`完全不使用网络起点；`residual_original`保留旧相关权重作搜索消融。仅波形流程仍是有限起点的局部搜索，不是纯拟合能力上限。旧GLS局部半径16 ns，新流程为24 ns，跨轮比较包含搜索设置变化；本轮权重消融共享相同搜索设置。
''')
    code('''residual_points = []
for condition, gap in [("low_low",30.),("mid_mid",8.),("high_high",8.),("low_high",8.),("high_low",40.),("high_low",200.)]:
    row = dict(条件=condition,间隔_ns=gap)
    for method in ["old_hybrid_1","selective_1","conditional_fit","regularized_1","residual_regularized","residual_waveform"]:
        r = residual_summary[(residual_summary.condition==condition)&(residual_summary.gap_ns==gap)&(residual_summary.method==method)].iloc[0]
        row[method] = f"{r.success}/{r.frames}"
    row["真实Q1中位_pe"] = r.q1_p16_p50_p84_pe[1]
    row["真实Q2中位_pe"] = r.q2_p16_p50_p84_pe[1]
    residual_points.append(row)
display(pd.DataFrame(residual_points))''')
    figure('0004-46-1','新独立双事件恢复曲线；(a–e)五种电荷组合，40对/点，右下图例。')
    figure('0004-46-2','原生空／单事件控制；(a)匹配后电荷误差，(b)额外输出帧，(c)整帧正确且无额外输出。')
    code('''residual_controls = []
for r in residual_summary.itertuples():
    if r.kind=="pair" or r.method not in ["blind_fit","old_hybrid_1","selective_1","regularized_1","residual_regularized","residual_waveform"]:
        continue
    residual_controls.append(dict(方法=r.method,条件=r.kind+"_"+r.condition,
        严格成功=f"{r.success}/{r.frames}",额外输出=f"{r.extra_frames}/{r.frames}",
        额外输出Wilson95上限=r.extra_wilson95[1],点估计预算通过=r.budget_pass,
        Q偏差中位_percent=None if r.q_error_p16_p50_p84 is None else 100*r.q_error_p16_p50_p84[1]))
display(pd.DataFrame(residual_controls))''')
    md('图46-2(a)的p16–p84是匹配条件下的误差分布，不是中位数置信区间，须与成功率一起看。0/80的额外输出Wilson95上限约4.58%；40/40整对成功的逐点下界约91.2%，均非所有条件的同时保证。图中连线只连接已测点。全部漏配、额外输出、误拆／误并标志、Q误差标准差及分位数保存在完整结果表。')
    md('''<a id="residual-examples-20261023"></a>
## 47. 完整波形：强脉冲之后的弱脉冲

预先固定high_low第0对、间隔40和200 ns，完整展示0–2048 ns记录。未按成功／失败挑图。真实分量来自雪崩祖先归属，只用于评价；恢复分量由最终Q/t经条件响应绘制。总波形残差变小不自动表示各事件电荷正确。
''')
    figure('0004-47-1','完整波形与新增事件；(a,b)总波形，(c,d)真实及恢复分量，(e,f)残差。')
    code('''residual_examples = json.loads((residual_output / "examples.json").read_text())
display(pd.DataFrame([dict(间隔_ns=e["case"]["gap_ns"],真Q_pe=e["true_q"],恢复Q_pe=e["fit"]["charge_pe"],
    真t_ns=e["true_t"],恢复t_ns=e["fit"]["times_ns"],原流程事件数=len(e["previous"]["times_ns"]),
    新流程事件数=len(e["fit"]["times_ns"]),严格整对成功=e["score"]["success"]) for e in residual_examples]))''')
    md('''<a id="residual-conclusions-20261023"></a>
## 48. 结果、限制与下一步

下表直接对比同一批新测试的原GLS、单独权重调整与完整新增候选流程；所有方法及失败条件均保留。是否改善由逐事件评分、单事件电荷偏差和误报预算共同判断，未在测试后换门槛或按条件拼接方法。
''')
    md(interpretation())
    code('''residual_diagnostics = json.loads((residual_output / "test-diagnostics.json").read_text())
display(pd.DataFrame(residual_diagnostics["high_low_200ns_counts"]))
residual_paired = pd.DataFrame(residual_report["residual_comparisons"])
display(residual_paired[(residual_paired.right=="residual_regularized") &
    (residual_paired.condition=="high_low") & (residual_paired.gap_ns.isin([40.,200.]))])''')
    code('''residual_conclusions = []
for condition in ["low_low","mid_mid","high_high","low_high","high_low"]:
    for gap in [8.,20.,40.,200.]:
        row = dict(条件=condition,间隔_ns=gap)
        for method in ["selective_1","regularized_1","regularized_count","residual_original","residual_regularized","residual_waveform"]:
            r=residual_summary[(residual_summary.condition==condition)&(residual_summary.gap_ns==gap)&(residual_summary.method==method)].iloc[0]
            row[method]=f"{r.success}/{r.frames}"
        residual_conclusions.append(row)
display(pd.DataFrame(residual_conclusions))''')
    md('''这些结果仍局限于当前中央入射、三档源能量、纯α和固定双τ响应。真实弱事件Q主要在数千pe，本轮没有新增500–1000 pe物理事件；没有free-τ、混合α/γ、随机高率K≥3或实物验证，也不能把可恢复间隔的倒数解释为最高计数率。下一步根据本表保留的失败类型，检查不同位置／连续电荷、未扫描的间隔与三事件复杂度；不得把本测试再当作未见数据进行门槛优化。

完整[逐点评分](artifacts/alpha_pair_residual_20261023/test-summary.csv)及[指标与冻结选择](artifacts/alpha_pair_residual_20261023/test-report.json)保留全部流程；实验源码／宏／模型标识在pileup/experiments/alpha_pair_residual_20261023，ROOT及种子已登记Run_Catalog。旧笔记单元、图像、旧实验和本轮早期开发输出完整保留。''')
    for i,cell in enumerate(cells):cell['id']=f'residual-20261023-{i:02d}'
    InteractiveShell.instance();namespace={};executed=0;old_cwd=os.getcwd()
    try:
        os.chdir(study.PROJECT/'jupyter')
        for cell in cells:
            if cell.cell_type!='code':continue
            with capture_output(stdout=True,stderr=True,display=True) as captured:
                exec(compile(cell.source,'<residual-notebook>','exec'),namespace)
            assert not captured.stderr,captured.stderr
            cell.outputs=[]
            if captured.stdout:cell.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=captured.stdout))
            for out in captured.outputs:cell.outputs.append(nbformat.v4.new_output('display_data',data=out.data,metadata=out.metadata))
            executed+=1;cell.execution_count=executed
    finally:os.chdir(old_cwd)
    anchors=[];images=[]
    for cell in cells:
        if cell.cell_type=='markdown':
            anchors.extend(re.findall(r'id="codex-figure-([^"]+)"',cell.source))
            assert not re.search(r'\]\([^)]*\.(?:pdf|ipynb|md)(?:[?#][^)]*)?\)',cell.source,re.I)
        else:
            for out in cell.outputs:
                assert out.output_type!='error'
                if 'image/png' in out.get('data',{}):images.append(hashlib.sha256(base64.b64decode(out.data['image/png'])).hexdigest())
    assert len(anchors)==len(set(anchors))==4 and len(images)==4
    assert set(images)=={study.digest(p) for p in (study.OUTPUT/'figures').glob('*.png')}
    nb.cells.extend(cells);assert json.loads(json.dumps(nb.cells[:len(original)]))==original
    assert len({c.id for c in nb.cells})==len(nb.cells)
    atomic_notebook(path,nb,before)
    study.save(study.RECORD/'notebook-validation.json',dict(passed=True,preserved_cells=len(original),
        new_cells=len(cells),new_code_cells_executed=executed,new_figure_ids=anchors,historical_cells_unchanged=True,
        new_outputs_executed=True,notebook_sha256=study.digest(path),atomic_write=True))


def catalog():
    path=study.PROJECT/'jupyter/Run_Catalog.ipynb';before=path.read_bytes();nb=nbformat.reads(before.decode(),as_version=4)
    cell_id='residual-catalog-20261023';assert not any(c.get('id')==cell_id for c in nb.cells)
    runs=study.load(study.EXPERIMENT/'manifest.json')['runs'];assert len(runs)==240
    text='''## 2026-10-23：相关权重与残差新增候选的独立检验

240个新原生单α母ROOT，源2.5／5.486／8 MeV各80，中央入射，20 μs单WaveformEvent，G4内部串行、最多4个任务。实际宏、源／SiPM／readout种子已核对RunInfo。双事件为新光子的共同微胞离线数字化，未称作原生双事件ROOT。固定95/351 ns双τ，充分动态范围；分析2.048 μs、局部K≤2，非注入率扫描。

运行清单：[manifest.json](../pileup/experiments/alpha_pair_residual_20261023/manifest.json)。过程在codex_0004第45–48节，结果位于artifacts/alpha_pair_residual_20261023。参数在生成前冻结；同一母样本跨间隔复用。

| ROOT | 源能量 MeV | 源seed | SiPM seed | readout seed | 状态 |
| --- | ---: | ---: | ---: | ---: | --- |
'''
    for r in sorted(runs,key=lambda r:r['seed']):
        text+=f"| {r['root']} | {r['energy_MeV']} | {r['seed']} | {r['sipm_seed']} | {r['readout_seed']} | {r['status']} |\n"
    cell=nbformat.v4.new_markdown_cell(text);cell.id=cell_id;nb.cells.append(cell)
    atomic_notebook(path,nb,before)
    study.save(study.RECORD/'catalog-validation.json',dict(passed=True,registered_roots=240,sha256=study.digest(path)))


if __name__=='__main__':deliver();catalog()
