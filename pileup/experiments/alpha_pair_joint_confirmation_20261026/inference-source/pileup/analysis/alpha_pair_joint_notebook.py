"""Append executed joint-likelihood results to the existing Notebook0004."""
from __future__ import annotations

import base64
import hashlib
import json
import os
import re
import shutil

import nbformat

from .alpha_pair_joint_data import PROJECT,EXPERIMENT,OUTPUT,RECORD
from .alpha_pair_conditional_notebook import atomic_notebook
from .alpha_pair_resolution_data import load,save,digest


def deliver():
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    assert load(RECORD/'science-validation.json')['passed']
    report=load(OUTPUT/'development-report.json')
    path=PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    previous=path.read_bytes();nb=nbformat.reads(previous.decode(),as_version=4)
    assert not any(c.get('id','').startswith('joint-nll-20261026-') for c in nb.cells)
    shutil.copyfile(path,RECORD/'notebook.before-delivery.ipynb')
    original=json.loads(json.dumps(nb.cells));cells=[]
    def md(text):cells.append(nbformat.v4.new_markdown_cell(text))
    def code(text):cells.append(nbformat.v4.new_code_cell(text))
    def figure(number,title):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{title}**')
        code(f'display(Image(filename=str(joint_output / "figures/{number}.png")))')
    md('''<a id="joint-component-diagnosis-20261026"></a>
## 58. 联合分量诊断：互相关是否是主要缺口

第57节提出检查双事件条件均值、联合协方差和候选依赖的似然。本轮完成这三项，继续使用**纯α、固定95/351 ns双τ、2 ns采样、充分动态范围、完整2.048 μs局部模拟电压、未知K≤2**。没有重训网络或新增G4 ROOT；研究对象仍是实际雪崩电荷Q（pe），并非源能量或沉积能量。t仍取SiPM有效面首光子时刻。

标定资料由此前三批开发记录组成：5400个共同微胞双事件、560个原生单事件和160个空窗，共6120帧。按光学母事件replicate整除2再模2分折，使同档成对母事件留在同一折；训练和留出各280个母事件，各3060帧。跨间隔、条件重用母样本，6120不是独立母事件数。留出折只有280个原生单事件用于增峰门槛；另一个此前已公开的320母事件集合提供1600帧开发验证。三组母ROOT无交叉，后者不能再叫独立测试。

**输入核验**：旧打包单事件的components是零占位，而原生truth中的signal_mV包含暗计数。初次汇总文件保留，但未用于标定；正式输入按avalanche_owner重放粒子分量，暗链留在背景。560例“粒子分量＋暗分量＝原无噪声信号”的最大差约1.48×10⁻¹² mV；逐例核对实际导出metadata中的读出配置，未用当前宏猜历史条件。

联合双事件分量来自两个粒子的光子共同推进SiPM微胞后的真实雪崩归属。统计在原采样网格上进行，按实验条件中心化；图中均值／协方差诊断使用真Q/t，不作为盲推断输入。

读出为2 ns上升、20 ns下降和另一个2 ns低通，故用固定离散滤波器A(z)=(1−e^(−2/20)z⁻¹)(1−e^(−2/2)z⁻¹)²消去电路极点。这三个时间尺度是**电子学参数，不是释放或改动闪烁双τ**。固定变换后舍弃前3个瞬态样点，拟合1021维边缘高斯近似；当前事件前有96–256 ns基线，不把此实现外推到事件紧贴记录起点的情况。
''')
    code('''from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Image
joint_root=Path.cwd().resolve()
while not (joint_root / "pileup/experiments/alpha_pair_joint_likelihood_20261026").is_dir():
    if joint_root.parent==joint_root:
        raise RuntimeError("请从项目根目录或jupyter目录运行")
    joint_root=joint_root.parent
joint_output=joint_root / "jupyter/artifacts/alpha_pair_joint_likelihood_20261026"
joint_record=joint_root / "jupyter/run_records/alpha_pair_joint_likelihood_20261026"
joint_report=json.loads((joint_output / "development-report.json").read_text())
joint_diagnosis=json.loads((joint_output / "joint-diagnosis.json").read_text())
joint_noise=json.loads((joint_record / "joint-noise-calibration.json").read_text())
display(pd.DataFrame([dict(间隔_ns=r["gap_ns"],留出帧=r["frames"],
    互协方差占总方差_percent=100*r["cross_variance_fraction"],
    相邻点相关=r["innovation_lag_correlation"][1])
    for r in joint_diagnosis["summary"] if r["model"]=="joint" and r["fold"]==1 and
    r["condition"]=="high_weak" and r["gap_ns"] in (8.,40.,200.)]))''')
    figure('0004-58-1','(a)留出强→弱组的互协方差贡献；(b)消去电子学极点后的残差相关。仅为使用真实分量的诊断。')
    md('''在这三个high_weak点，互协方差对总方差的贡献绝对值均不足1%，不能继续把强弱互相关当作已证实的主要原因。变换后主要剩相邻样点相关，适合检验带状协方差；短程近似仍不能代表完整微观SiPM分布。
''')
    md('''<a id="joint-normalized-likelihood-20261026"></a>
## 59. 实际改动：面积守恒、联合条件噪声与归一化似然

重新标定9特征条件均值时，发现一个物理约束缺口：旧经验形状修正会改变模板积分。这里Q定义为实际雪崩总电荷，线性读出完整积分应为Q×13.618501389 mV·ns/pe；微胞竞争可以改变形状和已实现Q，却不能在给定Q后任意改变这个面积。因此新增方案对每个候选的完整条件模板归一化，观测窗外尾部仍由模板预测，不把截断窗口的实测积分当完整Q。

留出分量的均方根误差：旧均值约4.486 mV；只重新标定约4.470 mV；再加入面积约束约4.327 mV。只是均值诊断，不能替代未知K逐事件成功率。积分残差的波动明显减小，支持在协方差中加入有限强度的总电荷条件约束。

协方差模型以**候选自身**的Q、t和滤波后均值确定局部噪声强度，用训练真实分量估计电荷依赖、邻居依赖、相邻相关及两个分量的交叉相关。电子噪声使用实际0.1 mV白噪声经过A后的带状矩阵，暗噪声与平均背景来自训练折。候选分量积分约束通过秩不超过2的修正实现，完整方案固定条件化份额0.95；训练估计约0.995，取0.95以保留余量。这是正定的高斯近似，仍未建模所有长程波动或离散雪崩统计。

对残差u=A(y−b−Σμ)，比较F=uᵀC⁻¹u＋log det C。**log det C已实际加入**，避免候选仅通过扩大预测方差降低代价。代价可能为负，增峰证据因此用max(F₁−F₂,0)，不再除以F₁。带状Cholesky和小矩阵更新已与独立构造的稠密联合协方差、逆矩阵和行列式核对。

| 新流程 | 均值 | 协方差及积分约束 |
| --- | --- | --- |
| old_mean_nll | 原条件模板 | 新候选依赖协方差＋log det；不加积分条件化 |
| joint_conserved_nll | 训练折重新标定＋完整面积守恒 | 新候选依赖协方差＋log det，0.95积分条件化 |

两者都是本轮新增盲流程。它们同时改变目标函数、候选排序和门槛；完整流程还改均值及积分约束，不能把最终变化全归于某一个因素。

文献沿用Fixsen等2004的条件均值／非平稳协方差思想，DOI：[10.1016/j.nima.2003.11.313](https://doi.org/10.1016/j.nima.2003.11.313)。其条件为TES模型、每能量5000次模拟，共享库AI稿仍待人工审核。原文加权最小二乘不能直接等同这里含log det的高斯近似；电路滤波、面积约束、带状加低秩实现、未知K搜索和GAGG性能均属本项目。此次复用既有笔记，没有重新精读或修改文献库。
''')
    figure('0004-59-1','(a)完整模板积分与实际电荷的关系；(b)训练分量积分残差的方差／未条件化短程模型方差。横轴Q扫描只检查模型性质，不是连续能量独立验证。')
    code('''display(pd.DataFrame([dict(方法=k,相邻相关=v["lag1"],交叉系数=v["cross"],
    积分条件化份额=v["kappa"],背景_mV=v["background_mV"],
    留出分量RMSE_mV=v["heldout_component_rmse_mV"])
    for k,v in joint_noise["models"].items()]))''')
    md('''<a id="joint-likelihood-blind-comparison-20261026"></a>
## 60. 未知K盲比较：恢复、额外输出与电荷偏差一起评价

每种新流程实际分析1880帧，其中280原生单事件只标定增峰门槛，另1600帧为开发评价：十种电荷组合×8/40/200 ns×40对，以及80空窗和320原生单事件。门槛按拟合单事件Q分档（1500/5000/12500 pe边界），采用含并列值和加一修正的经验上尾秩p≤0.025。Q≥500 pe且单事件NLL低于空窗NLL才输出首个事件；不读取真实K、Q、t、源率或粒子能量。

候选来自电压起点、残差扫描和固定间隔／电荷比例网格；NLL筛选后有限次数局部优化。此搜索不保证全局最优。推断worker禁止打开ROOT、truth及case文件。小样本首先发现K=0空分量轴的软件错误，修正后用run-v2.py执行；原冻结源码、空结果和错误原因保留，物理、标定参数和搜索预算没有改变。

整对成功要求两个真实Q均≥500 pe，各Q误差≤20%，各t误差≤min(10 ns,Δ/4)，且无额外输出。常规10 ns匹配同时用于漏检、额外输出、误拆／误并及Q/t误差；未匹配的时刻偏差输出也算额外，不能全部解释为凭空造峰。匹配Q误差是条件电荷精度，不是沉积能量分辨率。

预定选用条件不变：每组空／单额外输出≤5%；各单严格成功≥95%、Q误差中位绝对值≤5%；弱组合15点等权平均须超过两种原盲基线，较强组合15点平均相对原观测权重最多下降2个百分点。这些是开发样本点估计要求，并非总体误报保证。
''')
    code('''joint_rows=[]
for k,v in joint_report["baselines"].items():
    joint_rows.append(dict(方法=k,弱组合平均_percent=100*v["weak_mean"],较强组合平均_percent=100*v["strong_mean"],控制通过="见旧结果",满足选用条件="基线"))
for k,v in joint_report["results"].items():
    joint_rows.append(dict(方法=k,弱组合平均_percent=100*v["weak_mean"],较强组合平均_percent=100*v["strong_mean"],控制通过=v["controls_pass"],满足选用条件=v["eligible"]))
display(pd.DataFrame(joint_rows))
print("开发选中：",joint_report["chosen"],"；独立性能证据：",joint_report["independent_performance"])
joint_table=pd.DataFrame(joint_report["summary"])
display(joint_table[(joint_table.condition=="high_weak") | (joint_table.condition=="weak_high")][
    ["method","condition","gap_ns","success","eligible_frames","extra_frames","missed","split","merge"]]
    if all(k in joint_table for k in ["missed","split","merge"]) else
    joint_table[joint_table.condition.isin(["high_weak","weak_high"])][["method","condition","gap_ns","success","eligible_frames","extra_frames"]])''')
    figure('0004-60-1','全部十种组合的严格整对恢复率；每点40对，跨点复用母样本。所有方法使用同一批开发帧。')
    figure('0004-60-2','每组80个空／单控制；(a)额外输出帧，(b)严格成功，(c)匹配电荷偏差中位。虚线为预定要求。')
    code('''display(pd.DataFrame([dict(方法=k,条件=c,成功=f'{r["success"]}/{r["frames"]}',
    额外输出=f'{r["extra"]}/{r["frames"]}',Q误差中位_percent=None if r["median_q_bias"] is None else 100*r["median_q_bias"])
    for k,v in joint_report["results"].items() for c,r in v["controls"].items()]))
print("逐帧重放及收敛统计：")
display(pd.DataFrame(json.loads((joint_record / "science-validation.json").read_text())["diagnostics"]).T)''')
    md('''<a id="joint-likelihood-full-waveforms-20261026"></a>
## 61. 完整波形：拟合外观与逐事件恢复是否一致

仍展示先前已指定的high_weak第0对、40和200 ns，横轴覆盖0–2048 ns。上排总波形、中排弱分量、下排残差。拟合曲线均由盲流程产生；真实分量只用于画图和评价。中排在K=2时画时间排序后的第二个输出，它未必正确对应弱事件；K=1时不伪造第二条分量曲线。
''')
    figure('0004-61-1','预先指定的全窗例；(a,b)总信号，(c,d)真实弱分量与较晚拟合输出，(e,f)总残差。')
    code('''joint_examples=json.loads((joint_output / "examples.json").read_text())
display(pd.DataFrame([dict(间隔_ns=e["case"]["gap_ns"],方法=k,严格成功=v["success"],
    真Q_pe=v["true_q"],拟合Q_pe=v["fit"]["charge_pe"],真t_ns=v["true_t"],拟合t_ns=v["fit"]["times_ns"])
    for e in joint_examples for k,v in e["methods"].items()]))
print("真K辅助诊断：直接选用同一次盲搜索已有的双事件候选，跳过增峰门槛；不重新拟合。")
display(pd.DataFrame(json.loads((joint_output / "known-count-diagnosis.json").read_text())["summary"]))''')
    if report['chosen'] is None:
        conclusion='本轮两种新流程都未同时满足预定选用条件，chosen=null。保留完整实现、局部改善与失败，不升级默认流程，也未生成新的独立测试样本。'
    else:
        conclusion=f'开发阶段选中{report["chosen"]}。这只是已公开资料上的选择；本节不能据此宣称独立性能改善，冻结后的新母事件检验须另行记录。'
    md('''<a id="joint-likelihood-conclusion-20261026"></a>
## 62. 本轮结论与可继续推进的方向

'''+conclusion+'''

这一轮已经改变了拟合模型和目标函数：条件模板面积守恒、候选依赖联合协方差、log det归一化和重新标定的增峰门槛均已实际运行。对high_weak的诊断暂不支持把小的分量互相关作为主要瓶颈；完整积分的物理约束则应保留，但它本身不能保证把总电荷正确分给相近事件。

仍需面对强尾部涨落、时间—电荷耦合、短间隔的弱第二分量、有限候选搜索和高斯近似失配。四种离散源能量不是连续Q/位置覆盖；标定、门槛及开发评价都有母样本量限制。当前只检查局部K≤2、纯α、固定双τ，不给出自由τ最短窗、普适不可分辨间隔或最高注入率结论。

下一步优先把“总电荷正确、但分量划分错误”与“没有接受第二事件”分开，检查联合NLL在真实参数附近的形状和单／双模型证据；需要保留面积约束，再依据失败类型选择改进，避免只追求更小总残差。任何开发选择通过后，仍须冻结并用全新母事件确认，才能推进到随机高率、K≥3和α/γ。

完整[逐点指标](artifacts/alpha_pair_joint_likelihood_20261026/development-summary.csv)、[开发结果](artifacts/alpha_pair_joint_likelihood_20261026/development-report.json)、[同门槛标定样本的配对变化](artifacts/alpha_pair_joint_likelihood_20261026/paired-development.json)及[全窗例数据](artifacts/alpha_pair_joint_likelihood_20261026/examples.json)保留；原笔记输出、图号及历史实验不变。
''')
    for i,c in enumerate(cells):c['id']=f'joint-nll-20261026-{i:02d}'
    InteractiveShell.instance();namespace={};executed=0;old_cwd=os.getcwd()
    try:
        os.chdir(PROJECT/'jupyter')
        for c in cells:
            if c.cell_type!='code':continue
            with capture_output(stdout=True,stderr=True,display=True) as captured:
                exec(compile(c.source,'<joint-likelihood-notebook>','exec'),namespace)
            assert not captured.stderr,captured.stderr
            c.outputs=[]
            if captured.stdout:c.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=captured.stdout))
            for out in captured.outputs:c.outputs.append(nbformat.v4.new_output('display_data',data=out.data,metadata=out.metadata))
            executed+=1;c.execution_count=executed
    finally:os.chdir(old_cwd)
    nb.cells.extend(cells);assert nb.cells[:len(original)]==original
    anchors=re.findall(r'id="codex-figure-([^"]+)"','\n'.join(c.source for c in nb.cells if c.cell_type=='markdown'))
    assert len(anchors)==len(set(anchors))
    for c in cells:
        if c.cell_type=='markdown':assert not re.search(r'\]\([^)]*\.(?:pdf|ipynb|md)(?:[#?)][^)]*)?\)',c.source,re.I)
    atomic_notebook(path,nb,previous)
    images=[out['data']['image/png'] for c in cells if c.cell_type=='code' for out in c.outputs if out.get('output_type')=='display_data' and 'image/png' in out.get('data',{})]
    figures={str(p.relative_to(PROJECT)):digest(p) for p in (OUTPUT/'figures').glob('*.png')}
    assert len(images)==len(figures)==5
    assert {hashlib.sha256(base64.b64decode(v)).hexdigest() for v in images}==set(figures.values())
    save(RECORD/'notebook-delivery.json',dict(passed=True,previous_cells=len(original),new_cells=len(cells),executed_cells=executed,
        history_preserved=True,figures=figures,notebook_sha256=digest(path)))
    print('JOINT_NOTEBOOK_DELIVERED',len(cells),executed,flush=True)


if __name__=='__main__':deliver()
