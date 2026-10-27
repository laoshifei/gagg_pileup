"""Deliver independent joint-likelihood confirmation and its ROOT catalog."""
from __future__ import annotations

import base64
import hashlib
import json
import os
import re
import shutil

import nbformat

from .alpha_pair_joint_confirmation import PROJECT,EXPERIMENT,OUTPUT,RECORD
from .alpha_pair_resolution_data import load,save,digest
from .alpha_pair_conditional_notebook import atomic_notebook


def deliver():
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    assert load(RECORD/'science-validation.json')['passed']
    report=load(OUTPUT/'report.json')
    path=PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb';previous=path.read_bytes()
    nb=nbformat.reads(previous.decode(),as_version=4)
    assert not any(c.get('id','').startswith('joint-confirm-20261026-') for c in nb.cells)
    shutil.copyfile(path,RECORD/'notebook.before-delivery.ipynb');original=json.loads(json.dumps(nb.cells));cells=[]
    def md(text):cells.append(nbformat.v4.new_markdown_cell(text))
    def code(text):cells.append(nbformat.v4.new_code_cell(text))
    def figure(number,title):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{title}**')
        code(f'display(Image(filename=str(confirm_output / "figures/{number}.png")))')
    md('''<a id="joint-likelihood-independent-20261026"></a>
## 63. 冻结之后：320个新母事件的独立确认

第60节完整方案在开发数据上达到预定要求，故在查看本批拟合结果之前冻结均值、协方差、搜索预算和增峰门槛，随后生成320个新光学母ROOT：源能量1.625、2.5、5.486、8 MeV各80例。源能量不是沉积或可见能量；Q仍为实现雪崩电荷，首光子t定义不变。每个ROOT为一个20 μs单α原生窗口，实际物理源码和二进制与前批一致；三种G4／雪崩／读出随机流均事前核对不与旧试验重用。

weak、low、mid、high依次对应上述四档单α源能量；A_B表示先A后B。这些名称区分强弱，不表示α／γ类别，本批全部是α。

独立样本生成时发生一次进程SIGTERM：2个ROOT已完整核验，4次中断尝试保留原文件；随后使用单进程及原预定种子重跑这些中断项并完成其余样本，未按物理结果选择重跑。实际宏、种子、状态和原文件位置均登记Run_Catalog，未改写旧实验。

用这些新光子共同推进SiPM微胞，构造十种组合×8/40/200 ns×40对，加80空窗和320原生单事件，共1600帧。首对联合数字化与原标量实现逐数组完全一致。**没有直接叠加独立电压，也没有缩放光子数。** 不同间隔和组合仍复用母事件；1600帧不是1600个独立光学样本。

只确认事先选中的joint_conserved_nll，与原观测GLS、原模型GLS及使用共同280门槛母事件的原模型GLS对照。新样本不重新标定任何门槛、不选择新参数。推断仍为未知K≤2，只读电压；真实信息仅用于评分及明确标注的真K诊断。

整对恢复与控制规则沿用第60节。逐点报告Wilson95区间；两类平均各为15个条件点的等权平均。为处理母事件重复，另按replicate整除2的40个母事件块做2000次配对重采样：前20块与后20块分层，各层重抽20块，所有方法、能量档及间隔共享权重，同时重算阈上分母。该方式保留共同母事件相关性，不按1600帧独立重抽；仍只描述当前四档能量和固定物理模型。
''')
    code('''from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Image
confirm_root=Path.cwd().resolve()
while not (confirm_root / "pileup/experiments/alpha_pair_joint_confirmation_20261026").is_dir():
    if confirm_root.parent==confirm_root:
        raise RuntimeError("请从项目根目录或jupyter目录运行")
    confirm_root=confirm_root.parent
confirm_output=confirm_root / "jupyter/artifacts/alpha_pair_joint_confirmation_20261026"
confirm_record=confirm_root / "jupyter/run_records/alpha_pair_joint_confirmation_20261026"
confirm_report=json.loads((confirm_output / "report.json").read_text())
display(pd.DataFrame([dict(方法=k,弱组合平均_percent=100*v["weak_mean"],
    较强组合平均_percent=100*v["strong_mean"],控制通过=v["controls_pass"])
    for k,v in confirm_report["results"].items()]))
print("冻结前选中：",confirm_report["chosen_before_test"],"；本批是否满足确认要求：",confirm_report["test_requirements_passed"])
confirm_intervals=json.loads((confirm_output / "paired-confirmation.json").read_text())
display(pd.DataFrame([dict(基线=r["baseline"],组合=r["group"],平均差_pp=100*r["difference"],
    母事件分块95区间_pp=[100*x for x in r["mother_block_bootstrap95"]])
    for r in confirm_intervals["aggregate_intervals"]]))''')
    figure('0004-63-1','独立新母事件的全部十种组合；新方法误差线为逐点Wilson95区间，跨点相关性另由配对分块分析处理。')
    figure('0004-63-2','独立空／单控制，每组80例；(a)额外输出帧，(b)严格成功，(c)匹配Q偏差中位。')
    code('''confirm_table=pd.DataFrame(confirm_report["summary"])
display(confirm_table[(confirm_table.kind!="pair")][["method","condition","success","eligible_frames","extra_frames","q_error_p16_p50_p84"]])
print("单个输出未在容限内匹配也会计入额外，不能一概理解为无中生有的信号。")
display(confirm_table[confirm_table.condition.isin(["high_weak","weak_high"])][
    ["method","condition","gap_ns","success","eligible_frames","missed","extras","split","merge"]])''')
    md('''<a id="joint-confirmation-waveforms-20261026"></a>
## 64. 新样本完整波形及已知K诊断

以下仍是事前指定的high_weak第0对、40和200 ns，而不是按成功与否挑选的例子。横轴覆盖整个0–2048 ns。上排总波形、中排真实弱分量及较晚拟合输出、下排总残差。所有拟合曲线来自冻结盲流程；中排真实分量只作评价参照，第二个输出未必正确对应弱事件。
''')
    figure('0004-64-1','新独立样本的完整波形；(a,b)总信号，(c,d)弱分量，(e,f)残差。')
    code('''confirm_examples=json.loads((confirm_output / "examples.json").read_text())
display(pd.DataFrame([dict(间隔_ns=e["case"]["gap_ns"],方法=k,严格成功=v["success"],
    真Q_pe=v["true_q"],拟合Q_pe=v["fit"]["charge_pe"],真t_ns=v["true_t"],拟合t_ns=v["fit"]["times_ns"])
    for e in confirm_examples for k,v in e["methods"].items()]))
print("真K诊断：评分端强制取已有的双事件候选，跳过门槛；没有真实Q/t初值或新搜索。")
display(pd.DataFrame(json.loads((confirm_output / "known-count-diagnosis.json").read_text())["summary"]))''')
    conclusion=('本批满足事先的恢复与控制点估计要求。独立结果支持在当前纯α、固定双τ、局部K≤2条件下继续使用该研究候选；仍须结合上面的配对区间判断改善幅度，不能据此替换未检验的随机高率全窗流程。'
                if report['test_requirements_passed'] else
                '本批未同时达到事先的恢复与控制要求。开发选择仍完整保留，但本轮不将该方案升级为默认流程，也不根据新测试结果重新调门槛。应先依据失效类型开展下一轮开发，再另取新样本检验。')
    def success(method,condition,gap):
        r=next(r for r in report['summary'] if r['method']==method and r['condition']==condition and r['gap_ns']==gap)
        return f"{r['success']}/{r['eligible_frames']}"
    tradeoff=(f"独立样本中，先弱后强8 ns由原观测GLS的{success('observed_penalty','weak_high',8.)}提高到{success('joint_conserved_nll','weak_high',8.)}，是主要收益来源。"
        f"但先强后弱200 ns由{success('observed_penalty','high_weak',200.)}降为{success('joint_conserved_nll','high_weak',200.)}，原模型GLS为{success('model_relative_residual','high_weak',200.)}；"
        f"先强后弱8／40 ns仍分别为{success('joint_conserved_nll','high_weak',8.)}、{success('joint_conserved_nll','high_weak',40.)}。这是一种依赖条件的取舍，不能把平均改善写成所有场景改善。")
    md('''<a id="joint-confirmation-conclusion-20261026"></a>
## 65. 独立结果的适用范围与下一步

'''+conclusion+'\n\n'+tradeoff+'''

面积守恒约束的是总电荷；正确分配到两个事件还依赖时刻、相对强度和实际形状涨落。需要同时阅读强→弱、弱→强、空／单误拆和参数误差，不能仅凭平均恢复率或总波形残差下结论。

优化器未成功停止的候选也保留并计入评价，下表给出数量；有限搜索不等于全局最优。匹配Q误差仍是条件电荷精度，不是沉积能量分辨率；离散源能量、固定几何／噪声／τ、2 ns采样及局部K≤2限制全部保留。

接下来的研究应围绕弱事件的时间—电荷可辨识性及更准确的涨落模型，保留本轮面积守恒和严格控制作为基线；随后才扩展连续Q、位置与器件失配、随机高率K≥3和α/γ。当前没有得到“低于某个间隔必定无法解堆积”的普适定理，也没有将间隔直接换算为最高注入率。

不能根据真实强弱顺序在新旧方法之间逐帧挑选更好的结果；若尝试条件式融合，规则只能依赖观测，并须在开发数据上制定后另行冻结和检验。

完整[逐点结果](artifacts/alpha_pair_joint_confirmation_20261026/summary.csv)、[独立报告](artifacts/alpha_pair_joint_confirmation_20261026/report.json)、[配对分块分析](artifacts/alpha_pair_joint_confirmation_20261026/paired-confirmation.json)及[全窗例](artifacts/alpha_pair_joint_confirmation_20261026/examples.json)均保留。
''')
    code('''confirm_checks=json.loads((confirm_record / "science-validation.json").read_text())
print("ROOT核验数：",confirm_checks["completed_roots"],"；评分重放数：",confirm_checks["scores_replayed"])
print("完整性检查：",confirm_checks["passed"],"；门槛保持：",confirm_checks["gate_unchanged"])
display(pd.DataFrame([dict(方法=k,最终候选未收敛帧=v,总帧=1600)
    for k,v in confirm_checks["selected_not_converged"].items()]))''')
    for i,c in enumerate(cells):c['id']=f'joint-confirm-20261026-{i:02d}'
    InteractiveShell.instance();namespace={};executed=0;cwd=os.getcwd()
    try:
        os.chdir(PROJECT/'jupyter')
        for c in cells:
            if c.cell_type!='code':continue
            with capture_output(stdout=True,stderr=True,display=True) as captured:
                exec(compile(c.source,'<joint-confirmation-notebook>','exec'),namespace)
            assert not captured.stderr,captured.stderr;c.outputs=[]
            if captured.stdout:c.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=captured.stdout))
            for out in captured.outputs:c.outputs.append(nbformat.v4.new_output('display_data',data=out.data,metadata=out.metadata))
            executed+=1;c.execution_count=executed
    finally:os.chdir(cwd)
    nb.cells.extend(cells);assert nb.cells[:len(original)]==original
    anchors=re.findall(r'id="codex-figure-([^"]+)"','\n'.join(c.source for c in nb.cells if c.cell_type=='markdown'))
    assert len(anchors)==len(set(anchors))
    for c in cells:
        if c.cell_type=='markdown':assert not re.search(r'\]\([^)]*\.(?:pdf|ipynb|md)(?:[#?)][^)]*)?\)',c.source,re.I)
    atomic_notebook(path,nb,previous)
    images=[out['data']['image/png'] for c in cells if c.cell_type=='code' for out in c.outputs if out.get('output_type')=='display_data' and 'image/png' in out.get('data',{})]
    figures={str(p.relative_to(PROJECT)):digest(p) for p in (OUTPUT/'figures').glob('*.png')}
    assert len(images)==len(figures)==3
    assert {hashlib.sha256(base64.b64decode(v)).hexdigest() for v in images}==set(figures.values())
    save(RECORD/'notebook-delivery.json',dict(passed=True,previous_cells=len(original),new_cells=len(cells),executed_cells=executed,
        history_preserved=True,figures=figures,notebook_sha256=digest(path)))
    print('CONFIRMATION_NOTEBOOK_DELIVERED',len(cells),executed,flush=True)


def catalog():
    path=PROJECT/'jupyter/Run_Catalog.ipynb';previous=path.read_bytes();nb=nbformat.reads(previous.decode(),as_version=4)
    identifier='joint-confirm-runs-20261026';assert not any(c.get('id')==identifier for c in nb.cells)
    shutil.copyfile(path,RECORD/'catalog.before-delivery.ipynb');original=json.loads(json.dumps(nb.cells));manifest=load(EXPERIMENT/'manifest.json')
    assert len(manifest['runs'])==320
    lines=['## 2026-10-26：联合条件似然的独立母事件确认','',
        '冻结后新增320个完成ROOT，1.625/2.5/5.486/8 MeV源能量各80例，每个ROOT一个20 μs单α原生窗。双事件由保留光子共同推进微胞，局部分析窗2.048 μs；源能量不是沉积／可见能量。',
        '', '最初并行生产收到SIGTERM后保留2个已核验ROOT和4次中断尝试，后续单进程按相同预定种子恢复。每个完成ROOT的RunInfo、实际宏、源／雪崩／读出种子和哈希已核验。',
        '', '| 源能量 / MeV | ROOT | 源种子 | 雪崩种子 | 读出种子 | 实际宏 | 状态 |',
        '| ---: | --- | ---: | ---: | ---: | --- | --- |']
    for r in sorted(manifest['runs'],key=lambda r:r['seed']):
        status='中断后同种子重跑并核验' if r.get('retry_of') else '完成并核验'
        lines.append(f"| {r['energy_MeV']:g} | {r['root']} | {r['seed']} | {r['sipm_seed']} | {r['readout_seed']} | {r['macro']} | {status} |")
    lines+=['','### 保留的中断尝试','','| 源种子 | 原计划ROOT | 原文件状态 | 实际宏 | 日志 |','| ---: | --- | --- | --- | --- |']
    for r in manifest.get('failed_attempts',[]):
        state=f"保留未完成文件，{r['root_bytes']} bytes" if r['root_exists'] else '未创建'
        lines.append(f"| {r['plan']['seed']} | {r['root']} | {state} | {r['macro']} | {r['log']} |")
    cell=nbformat.v4.new_markdown_cell('\n'.join(lines));cell['id']=identifier;nb.cells.append(cell)
    assert nb.cells[:len(original)]==original;atomic_notebook(path,nb,previous)
    save(RECORD/'catalog-delivery.json',dict(passed=True,completed_roots=320,interrupted_attempts=len(manifest.get('failed_attempts',[])),
        previous_cells=len(original),history_preserved=True,notebook_sha256=digest(path)))
    print('CONFIRMATION_CATALOG_DELIVERED',flush=True)


if __name__=='__main__':deliver();catalog()
