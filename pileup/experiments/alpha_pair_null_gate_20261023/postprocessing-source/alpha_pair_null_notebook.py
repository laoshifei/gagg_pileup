"""Append the calibrated-birth experiment, preserving existing notebook cells."""
from __future__ import annotations

import base64
import hashlib
import json
import os
import re
import shutil

import nbformat

from . import alpha_pair_null_study as study
from .alpha_pair_conditional_notebook import atomic_notebook


def interpretation():
    report = study.load(study.OUTPUT/'test/report.json')
    primary = report['primary']; summary = report['summary']
    def row(method, condition, gap=None):
        return next(r for r in summary if r['method']==method and r['condition']==condition and r['gap_ns']==gap)
    weak_old = row('observed_penalty','weak'); weak_new = row(primary,'weak')
    low_old = row('observed_penalty','low'); low_new = row(primary,'low')
    status = '本批独立空／单事件控制全部通过预定要求' if report['all_control_requirements_passed'] else '本批独立空／单事件控制没有全部通过预定要求'
    text = f'''**同批独立结果**：约700 pe的原生弱单事件，匹配后Q相对误差中位数由{100*weak_old['q_error_p16_p50_p84'][1]:.2f}%变为{100*weak_new['q_error_p16_p50_p84'][1]:.2f}%；严格整帧成功为{weak_old['success']}/{weak_old['frames']} → {weak_new['success']}/{weak_new['frames']}。低Q单事件的额外输出为{low_old['extra_frames']}/{low_old['frames']} → {low_new['extra_frames']}/{low_new['frames']}。{status}。这些数字均比较本批相同观测，不拼接前一批结果。

**真实弱双事件**：以下分母只包含两个真实Q都≥500 pe的对；完整表同时保留所有生成帧及阈下数。'''
    for condition, gap in [('low_weak',40.),('high_weak',40.),('high_weak',200.),('weak_high',8.),('weak_weak',40.)]:
        old=row('observed_penalty',condition,gap); new=row(primary,condition,gap); known=row('known_count_model',condition,gap)
        text += f"\n\n{condition}、{gap:g} ns：原观测权重{old['success']}/{old['eligible_frames']}，新流程{new['success']}/{new['eligible_frames']}；使用真实K的模型权重诊断为{known['success']}/{known['eligible_frames']}。"
    text += '''

**结论边界**：降低假增峰、修正单事件电荷偏差与正确恢复每一个弱事件，是三个不同要求。残差较小或输出两个事件，不足以证明两者的Q/t都正确。已知K诊断仍使用同一有限候选搜索，也受模型失配和近邻参数相关性影响；其曲线不是理论上限，不能据此断言某个间隔“一定不可解”。

本轮不直接替换高率全窗默认流程。后续应针对新表保留的失败区间分析时间／电荷耦合、邻居条件响应和模型失配，再扩展不同位置、连续电荷、α/γ及K≥3；不能把本批测试继续当未见数据选择门槛，也不能将受控间隔的倒数当作最高计数率。'''
    return text


def deliver():
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    assert study.load(study.RECORD/'science-validation.json')['passed']
    path=study.PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    previous=path.read_bytes(); nb=nbformat.reads(previous.decode(),as_version=4)
    assert not any(c.get('id','').startswith('null-gate-20261023-') for c in nb.cells)
    shutil.copyfile(path,study.RECORD/'notebook.before-delivery.ipynb')
    original=json.loads(json.dumps(nb.cells)); cells=[]
    def md(text): cells.append(nbformat.v4.new_markdown_cell(text))
    def code(text): cells.append(nbformat.v4.new_code_cell(text))
    def figure(number,title):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{title}**')
        code(f'display(Image(filename=str(null_output / "figures/{number}.png")))')
    md('''<a id="null-gate-development-20261023"></a>
## 49. 单事件增峰标定与真实弱α对照

第48节的8个误拆案例把低Q单脉冲分出了约522–754 pe、相隔4–8 ns的额外事件。本轮同时检验两件事：单事件的残差改善能否给出有效的增峰筛选；这种筛选会损失多少真实500–1000 pe弱事件。原来的独立测试已经公开，现只作开发数据。

仍采用纯α、固定95/351 ns双τ、2 ns采样、2.048 μs局部电压与充分动态范围；未知K限于0/1/2。没有重训网络，也没有向推断提供真实K、Q、t或源能量。前一轮纯波形搜索与增加网络起点在整对成败上完全一致，本轮沿用纯波形搜索，保留其有限起点和局部优化限制。

**真实弱事件**：先做源能量1.5／1.75／2.0 MeV各4个ROOT，实际Q中位约423／1034／1583 pe；再补1.625 MeV的4个ROOT，中位约713 pe，固定该源能量。Q是实际雪崩等效pe，源能量不等于入射、沉积或可见能量。之后生成80个新光学母ROOT，与已公开的较强母光子一起进行共同微胞数字化，构造五种弱事件组合、七个间隔、40对/点，共1400双事件及80原生单事件。未缩放光子、未相加独立数字化电压。

**开发中发现的偏差**：只在旧搜索后收紧增峰判据，低Q单事件误拆可由12/160降到3/160；但约700 pe单事件的Q偏差中位仍为−10.10%，没有通过全部条件。20个开发弱单事件的对照中，普通最小二乘约+1%，由观测构造权重的GLS约−10%，以初步拟合的平滑模型构造权重后约−1%。这些是开发诊断，支持优先检查权重；没有单凭此对照证明唯一的微观偏差机制。

**模型权重**：先对观测电压作单事件OLS拟合，以其平滑均值响应形成经验GLS权重，ρ仍为0.9；随后所有K/Q/t候选共享该帧固定矩阵。这里的单事件是用于估计权重的假设，最终仍盲比较K=0/1/2。均值响应未重新标定，器件物理参数未改变；该矩阵仍是经验近似，不是精确的SiPM微观似然。

**增峰证据**：将完整搜索后的单／双事件代价记为C1、C2，比较S=(C1−C2)/C1及S=(C1−C2)/C0两个正改善评分；不合电荷门槛的双事件不给增峰证据。按单事件拟合Q在1500、5000、12500 pe处分四档，利用560个开发单事件的评分分布计算p=(1+分数≥当前S的标定数量)/(N+1)，并列保守处理，固定p≤0.025。这个经验上尾秩不是“第二个事件为假的概率”，也不是χ²显著性或跨模型误报保证。

最终模型权重方案用经验秩替代K1→K2的固定惩罚，允许恢复改善小于原固定门槛、但超过单事件涨落范围的弱脉冲；零／单事件检测保留旧判据。没有统一拒绝短间隔，也没有提高500 pe输出门槛。

模型权重开发共3120帧：两批已公开资料各取每条件／间隔前10对及全部控制，再加新1480帧弱事件资料。要求各控制额外输出≤5%、各单事件成功≥95%、各单事件Q偏差中位绝对值≤5%，然后按强／弱组合等权的平均整对成功选择评分。控制样本同时用于经验秩标定，这里的通过只是开发筛选；后续新独立测试才检验实际控制效果。

文献依据仅复用已有精读笔记：Zhang等2024（LaBr₃／R6233 PMT、500 MHz，DOI 10.1016/j.nima.2024.169273）提醒模板失配会产生幽灵脉冲，具体增峰阈值未充分报告。其AI笔记仍待人工审核，本轮未重读PDF。当前模型权重和经验秩是本项目实现，不是论文给出的GAGG处方。
''')
    code('''from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Image
null_root = Path.cwd().resolve()
while not (null_root / "pileup/experiments/alpha_pair_null_gate_20261023").is_dir():
    if null_root.parent == null_root:
        raise RuntimeError("请从项目根目录或jupyter目录运行")
    null_root = null_root.parent
null_output = null_root / "jupyter/artifacts/alpha_pair_null_gate_20261023"
null_experiment = null_root / "pileup/experiments/alpha_pair_null_gate_20261023"
null_report = json.loads((null_output / "test/report.json").read_text())
null_summary = pd.DataFrame(null_report["summary"])
null_primary = null_report["primary"]
null_selection = null_report["selection"]
display(pd.DataFrame(json.loads((null_output / "pilot-report.json").read_text())["summary"]))
display(pd.DataFrame(null_selection["results"][null_selection["chosen"]]["controls"]).T)
print("开发选择：", null_primary, "；标定单事件数：", null_selection["calibration_rows"])''')
    figure('0004-49-1','单事件经验上尾分布；(a,b)两种证据评分，颜色为拟合Q分档，竖虚线为各档对应门槛。')
    md('''<a id="null-gate-independent-20261023"></a>
## 50. 冻结后独立测试：电荷、误拆与真实弱事件

冻结源码、响应和筛选规则后，另生成1.625／2.5／5.486／8 MeV各80个母ROOT，共320个。十种电荷组合、4/8/12/20/40/80/200 ns、每点40对，再加80空窗与320原生单事件，总3200帧。双事件由新光子共同推进微胞；原生单事件按名义触发加独立抖动截取，未按真实首光子居中。

同点40对的母样本互不重叠，但跨间隔及部分条件复用母样本，不能把全部帧数当独立曝光数。t与Δ均指进入SiPM有效面的首光子时刻及其间隔，不是源发射时刻、峰值时刻或单脉冲持续长度。

整对指标要求两个真实Q都≥500 pe，各Q误差≤20%，各t误差≤min(10 ns,Δ/4)，一对一恢复且无额外输出。阈下对不混入两事件恢复分母，但保留全部生成数量及阈下数量。通常10 ns匹配下的逐事件效率、漏配、未正确恢复、额外输出及误拆／误并标志另报。Q误差分布条件于匹配，不等于能量分辨率。

比较五种流程：原观测权重＋固定惩罚、模型权重＋固定惩罚、模型权重＋两种经验秩，以及真实K辅助诊断。最后一种由评分端指定K，只用于判断计数选择与参数拟合各自的影响，不是盲重建结果，也不是真正理论上限。其余流程只接收观测波形和冻结参数。
''')
    code('''null_controls = []
for r in null_report["summary"]:
    if r["kind"]=="pair" or r["method"] not in ["observed_penalty","model_penalty",null_primary]:
        continue
    q = r["q_error_p16_p50_p84"]
    null_controls.append(dict(方法=r["method"],条件=r["condition"],
        严格成功=f'{r["success"]}/{r["frames"]}',额外输出=f'{r["extra_frames"]}/{r["frames"]}',
        额外输出Wilson95上限=r["extra_wilson95"][1],Q误差中位_percent=None if q is None else 100*q[1]))
display(pd.DataFrame(null_controls))
print("本批独立控制全部通过：", null_report["all_control_requirements_passed"])''')
    figure('0004-50-1','原三档电荷组合的独立恢复曲线；灰色虚线使用真实K，只作诊断。')
    figure('0004-50-2','包含真实弱α的五种组合；每点生成40对，效率分母要求两个真实Q≥500 pe。')
    figure('0004-50-3','80例/组的独立控制；(a)额外输出及Wilson95区间，(b)严格整帧成功，(c)匹配Q误差中位及p16–p84。')
    code('''null_points = []
for condition, gap in [("high_low",40.),("high_low",200.),("low_weak",40.),("high_weak",40.),("high_weak",200.),("weak_high",8.),("weak_weak",40.)]:
    point = dict(条件=condition,间隔_ns=gap)
    for method in ["observed_penalty","model_penalty",null_primary,"known_count_model"]:
        r = null_summary[(null_summary.condition==condition)&(null_summary.gap_ns==gap)&(null_summary.method==method)].iloc[0]
        point[method] = f"{r.success}/{r.eligible_frames}"
    point["生成对数"] = r.frames
    point["阈下对数"] = r.below_threshold_pairs
    point["真Q1中位_pe"] = r.q1_p16_p50_p84_pe[1]
    point["真Q2中位_pe"] = r.q2_p16_p50_p84_pe[1]
    null_points.append(point)
display(pd.DataFrame(null_points))''')
    md('逐点Wilson95区间描述有限样本：0/80额外输出的上限约4.58%，40/40成功的下限约91.2%，不是全曲线或总体同时保证。p16–p84是匹配误差的分布宽度，不是中位数置信区间。配对bootstrap仅在同一条件／间隔内比较，不跨复用母样本累计独立曝光。')
    md('''<a id="null-gate-waveforms-20261023"></a>
## 51. 完整波形与弱分量

预先固定high_weak第0对、间隔40和200 ns，不按成功／失败挑图。每列从上至下为总波形、第一个分量、第二个分量及残差；各分量独立设置纵轴，全部横轴均为完整0–2048 ns。真分量由共同微胞雪崩的祖先归属取得，只用于评价；恢复分量由盲拟合参数计算。
''')
    figure('0004-51-1','完整记录；(a,b)总波形，(c,d)第一分量，(e,f)真实弱分量与恢复，(g,h)残差。')
    code('''null_examples = json.loads((null_output / "test/examples.json").read_text())
display(pd.DataFrame([dict(间隔_ns=e["case"]["gap_ns"],真Q_pe=e["true_q"],恢复Q_pe=e["new"]["charge_pe"],
    真t_ns=e["true_t"],恢复t_ns=e["new"]["times_ns"],原事件数=len(e["old"]["times_ns"]),
    新事件数=len(e["new"]["times_ns"]),严格成功=e["new_score"]["success"],经验上尾秩=e["new"].get("null_rank_p"))
    for e in null_examples]))''')
    md('<a id="null-gate-conclusion-20261023"></a>\n## 52. 本轮结论与保留问题\n\n'+interpretation())
    code('''null_comparisons = pd.DataFrame(null_report["paired_comparisons"])
display(null_comparisons[(null_comparisons.baseline=="observed_penalty") &
    null_comparisons.condition.isin(["low_weak","high_weak","weak_high","weak_weak"]) &
    null_comparisons.gap_ns.isin([8.,40.,200.])])''')
    md('''完整[逐点结果](artifacts/alpha_pair_null_gate_20261023/test/summary.csv)、[完整指标与选择](artifacts/alpha_pair_null_gate_20261023/test/report.json)及[示例数值](artifacts/alpha_pair_null_gate_20261023/test/examples.json)保留全部方法与失败条件。实验位于pileup/experiments/alpha_pair_null_gate_20261023，实际ROOT／宏／种子已登记Run_Catalog。开发阶段重启的完整结果经核验后保留，未完成尝试保留原文件并以相同预定种子重跑；未根据物理结果选择重跑对象。

本轮限于当前中央入射、固定双τ、纯α和局部K≤2；未验证free-τ、混合α/γ、随机高率K≥3或实物器件。旧笔记计算、图像和图号完整保留。''')
    for i,cell in enumerate(cells): cell['id']=f'null-gate-20261023-{i:02d}'
    InteractiveShell.instance(); namespace={}; executed=0; old_cwd=os.getcwd()
    try:
        os.chdir(study.PROJECT/'jupyter')
        for cell in cells:
            if cell.cell_type!='code': continue
            with capture_output(stdout=True,stderr=True,display=True) as captured:
                exec(compile(cell.source,'<null-gate-notebook>','exec'),namespace)
            assert not captured.stderr,captured.stderr
            cell.outputs=[]
            if captured.stdout: cell.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=captured.stdout))
            for out in captured.outputs: cell.outputs.append(nbformat.v4.new_output('display_data',data=out.data,metadata=out.metadata))
            executed+=1; cell.execution_count=executed
    finally: os.chdir(old_cwd)
    nb.cells.extend(cells)
    assert nb.cells[:len(original)]==original
    anchors=re.findall(r'id="codex-figure-([^"]+)"','\n'.join(c.source for c in nb.cells if c.cell_type=='markdown'))
    assert len(anchors)==len(set(anchors))
    atomic_notebook(path,nb,previous)
    figures={p.name:study.digest(p) for p in (study.OUTPUT/'figures').glob('*.png')}
    images=[out['data']['image/png'] for c in cells if c.cell_type=='code' for out in c.outputs
            if out.get('output_type')=='display_data' and 'image/png' in out.get('data',{})]
    assert {hashlib.sha256(base64.b64decode(value)).hexdigest() for value in images}==set(figures.values())
    study.save(study.RECORD/'notebook-delivery.json',dict(passed=True,previous_cells=len(original),new_cells=len(cells),
        executed_cells=executed,figures=figures,notebook_sha256=study.digest(path),history_preserved=True))
    print('NOTEBOOK_DELIVERED',len(cells),'cells',executed,'executed',flush=True)


def catalog():
    path=study.PROJECT/'jupyter/Run_Catalog.ipynb'
    previous=path.read_bytes(); nb=nbformat.reads(previous.decode(),as_version=4)
    identifier='null-gate-runs-20261023'
    assert not any(c.get('id')==identifier for c in nb.cells)
    shutil.copyfile(path,study.RECORD/'catalog.before-delivery.ipynb')
    original=json.loads(json.dumps(nb.cells))
    lines=['## 2026-10-23：真实弱α、经验增峰筛选和独立测试',
        '', '本轮416个完成ROOT：先导12、补测4、弱事件开发80、冻结后独立测试320。均为一个20 μs单α原生窗口；双事件使用保留光子共同推进微胞，分析窗2.048 μs。源能量不是沉积／可见能量。',
        '', '各行实际宏、源／雪崩／读出种子与ROOT的对应关系如下；RunInfo已核对。恢复登记行原进程退出码不可得，保留该限制。',
        '', '| 阶段 | 源能量 / MeV | ROOT | 源种子 | 雪崩种子 | 读出种子 | 实际宏 | 完成方式 |',
        '| --- | ---: | --- | ---: | ---: | ---: | --- | --- |']
    count=0
    for stage in ('pilot','pilot-refinement','development','test'):
        manifest=study.load(study.EXPERIMENT/stage/'manifest.json')
        assert len(manifest['runs'])==len(manifest['planned_runs'])
        for r in sorted(manifest['runs'],key=lambda r:r['seed']):
            note='重启后核验完整ROOT，原退出码未知' if r.get('recovered_after_reboot') else '原生完成并核验'
            lines.append(f"| {stage} | {r['energy_MeV']:g} | {r['root']} | {r['seed']} | {r['sipm_seed']} | {r['readout_seed']} | {r['macro']} | {note} |")
            count+=1
    recovery=study.load(study.PROJECT/study.load(study.RECORD/'latest-recovery.json')['directory']/'audit.json')
    lines += ['', '### 重启中断保留记录', '', '原文件未删除；完整13个已登记ROOT哈希不变，另1个完整ROOT恢复登记。以下4次尝试按相同预定种子重跑，未按物理结果筛选。',
              '', '| 源种子 | 原计划ROOT | 原文件状态 | 保留的实际宏 |', '| ---: | --- | --- | --- |']
    for r in recovery['interrupted']:
        state=f"保留，{r['root_bytes']} bytes，未完成" if r['root_exists'] else '未创建'
        lines.append(f"| {r['plan']['seed']} | {r['root']} | {state} | {r['.mac']} |")
    cell=nbformat.v4.new_markdown_cell('\n'.join(lines));cell['id']=identifier;nb.cells.append(cell)
    assert nb.cells[:len(original)]==original and count==416
    atomic_notebook(path,nb,previous)
    study.save(study.RECORD/'catalog-delivery.json',dict(passed=True,completed_roots=count,
        interrupted_attempts=len(recovery['interrupted']),previous_cells=len(original),history_preserved=True,
        notebook_sha256=study.digest(path)))
    print('CATALOG_DELIVERED',count,flush=True)


if __name__=='__main__':
    deliver();catalog()
