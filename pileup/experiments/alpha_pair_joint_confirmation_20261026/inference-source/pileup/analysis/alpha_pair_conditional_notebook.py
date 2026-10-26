"""Append checked results without renumbering or replacing historical outputs."""
from __future__ import annotations

import base64
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import tempfile

import nbformat
import pandas as pd

from . import alpha_pair_conditional_study as study
from .alpha_pair_resolution_data import load,save,digest


def atomic_notebook(path,notebook,previous):
    nbformat.validate(notebook)
    encoded=(nbformat.writes(notebook)+'\n').encode()
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=path.stem+'.',suffix='.partial',delete=False) as stream:
        stream.write(encoded);stream.flush();os.fsync(stream.fileno());temp=Path(stream.name)
    assert temp.read_bytes()==encoded
    nbformat.validate(nbformat.reads(encoded.decode(),as_version=4))
    assert path.read_bytes()==previous,'Notebook changed during delivery; preserve the user edit'
    temp.replace(path)
    assert path.read_bytes()==encoded


def deliver():
    result=load(study.OUTPUT/'test-report.json');summary=result['summary']
    selection=result['selection'];cal=load(study.EXPERIMENT/'calibration.json')
    figures=study.OUTPUT/'figures'
    path=study.PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    previous=path.read_bytes();nb=nbformat.reads(previous.decode(),as_version=4)
    assert not any(c.get('id','').startswith('conditional-20261022-') for c in nb.cells),'Already delivered'
    before_cells=json.loads(json.dumps(nb.cells));cells=[]
    def md(text):cells.append(nbformat.v4.new_markdown_cell(text))
    def code(text,outputs):cells.append(nbformat.v4.new_code_cell(text,execution_count=1,outputs=outputs))
    def table(frame,code_text):
        code(code_text,[nbformat.v4.new_output('display_data',data={'text/html':frame.to_html(index=False),
            'text/plain':frame.to_string(index=False)},metadata={})])
    def figure(number,title):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{title}**')
        image=figures/(number+'.png')
        source=f'from IPython.display import display, Image\ndisplay(Image(filename=str(conditional_output / "figures/{number}.png")))'
        code(source,[nbformat.v4.new_output('display_data',data={'image/png':base64.b64encode(image.read_bytes()).decode(),
            'text/plain':f'<图{number}：{title}>'},metadata={})])
    md('''<a id="conditional-response-20261022"></a>
## 41. 电荷／邻近事件相关模板：先修正响应，再决定是否重拟合

第40节已经发现：同一个网络输出经固定模板拟合后，短间隔电荷分配可能变差。这一轮直接实现响应校准和选择性重拟合，沿用纯α、固定95/351 ns双τ、2 ns采样、完整2.048 μs局部模拟电压和充分动态范围。Q仍为实际雪崩等效pe，不直接当作源能量或沉积能量。

**开发边界**：上一轮240个母ROOT现在全部归为开发数据；按母样本块分成互不重叠的120个训练母样本和120个验证母样本。此前看过的数据不再称为独立测试。训练用真实祖先分量、Q及首光子时刻标定响应；验证用于选K惩罚和重拟合接受规则。最终测试重新生成240个母ROOT，参数在生成前冻结。

**均值修正**：每个分量写成Q乘以经验响应；响应在原模板上加入9项正则化修正，依赖自身Q、估计的邻居Q及相对时刻。邻居状态是从候选Q/t计算的代理量，不是读取真实微胞状态。后来事件的修正不提前生效；起始修正采用0–8 ns平滑因果包络；未观测尾部不按零信号学习。没有逐事件释放发光寿命，也没有修改G4的物理参数。

**相关加权**：从训练残差标定异方差AR(1)近似，再用本帧观测的固定平滑先验建立协方差。同帧各K/Q/t候选使用相同C，故比较时logdet为常数。此处为经验GLS，不宣称精确微观似然或全局最优。相关涨落的建模和均值修正分别保留消融结果。

方法思想参考已读Fixsen等的条件均值／非平稳噪声估计工作（DOI: 10.1016/j.nima.2003.11.313；TES模拟，原笔记AI待审），这里的9项特征、AR(1)近似和选择规则是本项目方案，不是该论文在GAGG上的既有结论。

新均值与新GLS均在盲起点附近±16 ns联合优化时刻和log Q，max_nfev=60（数值Jacobian调用另计），Q边界50–50000 pe；原固定模板为±12 ns变投影。新均值／新GLS共享优化设置，可作加权消融；旧／新比较属于整个流程的比较，不能把全部变化归于单一因素。

**实际比较**：原固定模板盲拟合；原始网络两个种子；原网络＋固定模板；仅波形起点的条件GLS；加入网络起点的条件GLS；网络＋新均值／新GLS；选择性重拟合。两个网络沿用既有同架构权重，本轮没有新增神经网络架构。未知K仍由观测决定，推断不读真实K、Q、t或注入率。
''')
    init='''from pathlib import Path
import json
import pandas as pd
conditional_root = Path.cwd().resolve()
while not (conditional_root / "pileup/experiments/alpha_pair_conditional_20261022").is_dir():
    if conditional_root.parent == conditional_root:
        raise RuntimeError("请从项目根目录或其 jupyter 目录运行")
    conditional_root = conditional_root.parent
conditional_output = conditional_root / "jupyter/artifacts/alpha_pair_conditional_20261022"
conditional_report = json.loads((conditional_output / "test-report.json").read_text())
conditional_summary = pd.DataFrame(conditional_report["summary"])
conditional_selection = conditional_report["selection"]
display(pd.DataFrame([dict(seed=k, **v["chosen"]) for k,v in conditional_selection["gates"].items()]))'''
    table(pd.DataFrame([dict(seed=k,**v['chosen']) for k,v in selection['gates'].items()]),init)
    md(f'''上表只用开发验证集选择。`gap_threshold_ns`是**网络估计**的间隔下限，`max_log_charge_change`是允许的最大绝对log电荷改变量；优化不收敛时保留原始输出。保持网络预测的事件数，既不读取真实间隔，也不按测试真值逐例选优。空窗和三档单事件各自额外输出≤5%才允许入选。

**实际选择的含义**：两个种子均选新GLS、间隔下限0、最大log电荷改变量10。在本实验范围，这几乎就是接受所有已收敛重拟合，只对不收敛结果保留原始输出；验证并没有支持一个更严格的间隔／改变量门槛。因此不能把这轮收益归功于选择性门控本身。

训练残差得到AR(1)系数ρ={cal['residual_calibration']['rho']:.4f}、shot系数={cal['residual_calibration']['shot_coefficient_mV']:.4f} mV。它们是本模拟域的残差参数，不是另测出的器件参数。首次开发尝试发现对称平滑产生非因果前沿，已保留旧尝试并在生成独立测试前修正；图和表均使用修正后重跑结果。''')
    figure('0004-41-1','响应标定与开发验证；(a)电荷相关单脉冲均值，(b)邻近事件状态代理，(c)开发验证中能组。')
    md('''<a id="conditional-test-20261022"></a>
## 42. 新母ROOT的独立检验

新生成源能量2.5、5.486、8 MeV各80个单α光学母ROOT，中央入射。双事件用新光子离线共同推进微胞，未缩放光子，未相加两个独立数字化电压。单事件控制保留原生电压。Δ=4/8/12/20/30/40/60/80/120/200 ns，五组条件每点40对；另有80空窗和240单事件，共2320帧。跨间隔及部分条件复用母样本，帧数不是独立物理曝光数。

严格双事件成功要求每个Q误差≤20%、每个时刻误差≤min(10 ns,Δ/4)，且一对一恢复、无额外输出；真Q门槛500 pe。下表“成功/n”是**整对成功**，不是仅匹配率。完整结果同时保存匹配后Q偏差／标准差、时刻误差、漏配、额外输出及误拆／误并标志。条件误差分布必须和效率一起解读。
''')
    chosen=['blind_fit','ml_1','old_hybrid_1','conditional_fit_only','conditional_fit','selective_1','ml_2','selective_2']
    subset=[r for r in summary if r['condition']=='mid_mid' and r['gap_ns'] in (4.,8.,12.,20.,40.,60.) and r['method'] in chosen]
    table_frame=pd.DataFrame([dict(方法=r['method'],间隔_ns=r['gap_ns'],成功=f"{r['success']}/{r['frames']}",
        Wilson95=f"{r['wilson95'][0]:.3f}–{r['wilson95'][1]:.3f}",额外输出=r['extra'],漏配=r['missed']) for r in subset])
    table(table_frame,'''conditional_mid = conditional_summary[
    (conditional_summary.condition == "mid_mid") & conditional_summary.gap_ns.isin([4,8,12,20,40,60]) &
    conditional_summary.method.isin(["blind_fit","ml_1","old_hybrid_1","conditional_fit_only","conditional_fit","selective_1","ml_2","selective_2"])]
display(pd.DataFrame([dict(方法=r.method, 间隔_ns=r.gap_ns, 成功=f"{r.success}/{r.frames}",
    Wilson95=f"{r.wilson95[0]:.3f}–{r.wilson95[1]:.3f}",额外输出=r.extra,漏配=r.missed)
    for r in conditional_mid.itertuples()]))''')
    figure('0004-42-1','独立双事件恢复曲线；(a–e)五种电荷组合，右下图例。40对/点，条件GLS曲线包含网络起点。')
    md('`conditional_fit_only`只使用波形候选；`conditional_fit`还比较两个网络给出的起点，均以观测建立候选并选择K。不能把后者的全部收益单独归因于均值模板。`selective_1/2`分别对应两套冻结网络及其开发选择规则；未在测试后逐点选取较好种子。图中线条连接已测点，不表示中间间隔已经验证。')
    figure('0004-42-2','独立空窗／原生单事件；(a)额外输出帧比例及Wilson95区间，(b)恢复正确且无额外输出的整帧成功率。')
    controls=[r for r in summary if r['kind']!='pair' and r['method'] in chosen]
    frame=pd.DataFrame([dict(方法=r['method'],条件=r['kind']+'_'+r['condition'],
        正确且无额外=f"{r['success']}/{r['frames']}",额外输出帧=f"{r['extra_frames']}/{r['frames']}",
        点估计预算通过=r['budget_pass'],额外输出区间上限=r['extra_wilson95'][1]) for r in controls])
    table(frame,'''display(pd.DataFrame([dict(方法=r.method,条件=r.kind+"_"+r.condition,
    正确且无额外=f"{r.success}/{r.frames}",额外输出帧=f"{r.extra_frames}/{r.frames}",
    点估计预算通过=r.budget_pass,额外输出区间上限=r.extra_wilson95[1])
    for r in conditional_summary.itertuples() if r.kind!="pair" and r.method in
    ["blind_fit","ml_1","old_hybrid_1","conditional_fit_only","conditional_fit","selective_1","ml_2","selective_2"]]))''')
    md('“点估计预算通过”仅表示该组实测额外输出≤5%；区间上限单列。0/80的Wilson95上限约4.58%，不是绝对零误报保证。逐点区间及同一条件／间隔的配对bootstrap差值也不是所有间隔的同时保证。原来的真事件被找回但另报事件，不能算整帧成功。')
    md('''<a id="conditional-examples-20261022"></a>
## 43. 完整堆积波形和分量

固定取新测试mid_mid第0对、Δ=8及40 ns，未按成功率挑例。每幅完整展示0–2048 ns记录。真实分量来自雪崩祖先归属，只用于评价和绘图；恢复分量用最终Q/t通过校准响应呈现。若接受规则保留原始网络Q/t，这条渲染曲线也不是该网络直接输出的真实波形。残差较小仍不等于逐事件Q/t都正确。
''')
    figure('0004-43-1','完整局部波形与解堆积分量；(a,b)总波形，(c,d)真实及重建分量，(e,f)残差。')
    examples=load(study.OUTPUT/'examples.json')
    frame=pd.DataFrame([dict(间隔_ns=e['case']['gap_ns'],真Q_pe=e['true_q'],恢复Q_pe=e['fit']['charge_pe'],
        真t_ns=e['true_t'],恢复t_ns=e['fit']['times_ns'],接受重拟合=e['fit']['accepted_refit'],严格成功=e['score']['success']) for e in examples])
    table(frame,'''conditional_examples = json.loads((conditional_output / "examples.json").read_text())
display(pd.DataFrame([dict(间隔_ns=e["case"]["gap_ns"],真Q_pe=e["true_q"],恢复Q_pe=e["fit"]["charge_pe"],
    真t_ns=e["true_t"],恢复t_ns=e["fit"]["times_ns"],接受重拟合=e["fit"]["accepted_refit"],严格成功=e["score"]["success"])
    for e in conditional_examples]))''')
    md('''<a id="conditional-conclusion-20261022"></a>
## 44. 本轮结果、限制和下一步

以下将原始网络、无条件固定模板重拟合和冻结的选择性流程直接并列。全部为这批新母样本的实际结果；未事后替换接受阈值。是否改进须同时看双事件恢复和各空／单控制，不用总波形残差排名代替逐事件评分。
''')
    points=[]
    for condition in ('mid_mid','high_low','low_high'):
        for gap in (8.,20.,40.,80.,200.):
            row=dict(条件=condition,间隔_ns=gap)
            for method in ('ml_1','old_hybrid_1','conditional_fit_only','selective_1','ml_2','selective_2'):
                r=next(r for r in summary if (r['condition'],r['gap_ns'],r['method'])==(condition,gap,method))
                row[method]=f"{r['success']}/{r['frames']}"
            points.append(row)
    frame=pd.DataFrame(points)
    code_text='''conditional_points = []
for condition in ["mid_mid","high_low","low_high"]:
    for gap in [8.,20.,40.,80.,200.]:
        row = dict(条件=condition,间隔_ns=gap)
        for method in ["ml_1","old_hybrid_1","conditional_fit_only","selective_1","ml_2","selective_2"]:
            r = conditional_summary[(conditional_summary.condition==condition)&(conditional_summary.gap_ns==gap)&(conditional_summary.method==method)].iloc[0]
            row[method] = f"{r.success}/{r.frames}"
        conditional_points.append(row)
display(pd.DataFrame(conditional_points))'''
    table(frame,code_text)
    md('''单次独立验证并不保证所有波形都能正确分离。当前经验响应仅覆盖三档源能量、中央入射和这个器件模型；实际输入不读取源能量。相同母光子跨间隔复用，不能把跨点汇总的样本数用于独立曝光置信区间。Q误差标准差是匹配条件下的电荷指标，不直接声称沉积能量分辨率。

本轮只扫描4–200 ns，不能替代第36–40节0–800 ns的完整条件扫描；没有得出正数“绝对不可分间隔”，也不能取恢复阈值的倒数当作最高计数率。500–1000 pe物理弱事件、不同位置、free-τ、α/γ混合及随机高率K≥3仍待验证。

下一步应检查新的源能量和光学位置，检验响应插值及模型失配，再扩展到α/γ和多事件流。选择性重拟合出现的收益不意味着应对所有网络输出强制重拟合。完整[测试评分表](artifacts/alpha_pair_conditional_20261022/test-summary.csv)和[配对比较及指标](artifacts/alpha_pair_conditional_20261022/test-report.json)保留全部方法及条件。

冻结源码、配置、二进制和权重标识在pileup/experiments/alpha_pair_conditional_20261022；原生ROOT在g4/data，实际宏／种子／检查记录已登记Run_Catalog。观测与真值分开存储；推断工作进程安装标签文件读取拦截。旧开发尝试、原笔记计算与图像完整保留。
''')
    for i,cell in enumerate(cells):cell['id']=f'conditional-20261022-{i:02d}'
    nb.cells.extend(cells)
    assert json.loads(json.dumps(nb.cells[:len(before_cells)]))==before_cells
    (study.RECORD/'notebook.before-final-delivery.ipynb').write_bytes(previous)
    atomic_notebook(path,nb,previous)
    save(study.RECORD/'notebook-validation.json',dict(passed=True,preserved_cells=len(before_cells),
        appended_cells=len(cells),new_figure_ids=['0004-41-1','0004-42-1','0004-42-2','0004-43-1'],
        original_sha256=hashlib.sha256(previous).hexdigest(),notebook_sha256=digest(path),
        notebook_format_valid=True,atomic_write=True,generated_from='Frozen test-report.json and saved images'))


def catalog():
    path=study.PROJECT/'jupyter/Run_Catalog.ipynb';previous=path.read_bytes()
    nb=nbformat.reads(previous.decode(),as_version=4)
    assert not any(c.get('id')=='conditional-catalog-20261022' for c in nb.cells)
    manifest=load(study.EXPERIMENT/'manifest.json');runs=manifest['runs']
    assert len(runs)==240 and all(r['status']=='completed' for r in runs)
    text='''## 2026-10-22：条件响应／选择性重拟合的独立母ROOT

240个新单α母ROOT，源能量2.5/5.486/8 MeV各80、中央入射，固定双τ95/351 ns，每文件一个20 μs WaveformEvent。G4每进程串行、最多4任务、启动至少相隔1秒。实际配置、宏、源种子和SiPM/readout种子已逐文件核对RunInfo；ROOT保留于g4/data。双事件分析来自这些新光子离线共同数字化，不是把它们称为原生双事件ROOT，也不是源注入率扫描。

完整清单：[运行清单](../pileup/experiments/alpha_pair_conditional_20261022/manifest.json)。波形／真值与冻结分析结果保存于artifacts/alpha_pair_conditional_20261022；母样本跨间隔复用。codex_0004第41–44节记录方法和结果。

| ROOT | 源能量 MeV | G4 seed | SiPM seed | readout seed | 状态 |
| --- | ---: | ---: | ---: | ---: | --- |
'''
    for r in sorted(runs,key=lambda r:r['seed']):
        text+=f"| {r['root']} | {r['energy_MeV']} | {r['seed']} | {r['sipm_seed']} | {r['readout_seed']} | {r['status']} |\n"
    cell=nbformat.v4.new_markdown_cell(text);cell['id']='conditional-catalog-20261022';nb.cells.append(cell)
    atomic_notebook(path,nb,previous)
    save(study.RECORD/'catalog-validation.json',dict(passed=True,registered_roots=240,sha256=digest(path)))


if __name__=='__main__':deliver();catalog()
