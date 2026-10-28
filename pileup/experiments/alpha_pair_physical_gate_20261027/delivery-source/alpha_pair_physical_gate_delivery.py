"""Executed research record and final checks for physical count-gate ablation."""
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

from .alpha_pair_physical_gate_study import PROJECT,EXPERIMENT,OUTPUT,RECORD,PREVIOUS,NAME,check_freeze
from .alpha_pair_physical_gate import METHODS,PRIMARY
from .alpha_pair_information_study import load,save,digest,rows

MARKER='alpha-pair-physical-gate-20261027'
NOTEBOOK=PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
LITERATURE='''只读复用共享库已有精读及主题／审核记录，本轮未联网、未重查PDF、未修改共享库。Liddick、Darby、Grzywacz，2012，*Algorithms for pulse shape analysis using silicon detectors*，DOI 10.1016/j.nima.2011.12.034，来源为papers/sources/2012_Liddick_NIMA_Algorithms-for-pulse-shape-analysis-using-silicon-detectors.md，PDF3–8、图3–13：约65 μm DSSD、40 MSPS／12 bit、25 μs记录，以选出的实测单脉冲缩放／移位构造8000余双脉冲。仅借鉴“单／双判别与参数精度分别评价”和慢上升／噪声引起的假双控制；等幅100 ns与幅比5的250 ns不是同一条件，不移用为GAGG性能。该科学笔记仍为AI稿待人工审核，参数数目3／4和门选记号疑问保留，不移用其6.8／8.5数值门槛。

Morad等，2024，*Model-Based Deep Learning Algorithm for Detection and Classification at High Event Rates*，DOI 10.1109/TNS.2024.3371573，来源为papers/sources/2024_Morad_TNS_Model-Based-Deep-Learning-Algorithm-for-Detection-and-Classification-at-High-Event-Rates.md，PDF3–5方法、PDF10表II：25 mm CLYC、PMT800 V、1 GS/s、17 mV阈值，Poisson线性叠加＋AWGN、K≤3。仅采用检测、条件分类和联合成功须区分的原则，不把其分类率当逐事件电荷恢复率。人工审核第553项只是限定确认指标区分和适用条件，未全面确认方法；NoA posterior／likelihood、Poisson项和误报分母的疑问保持。当前NLL差、秩校准与候选消融是本项目实现，不是两篇论文的算法复现。'''


def result_table(report):
    lines=['| 流程 | 原生360对 | 已公开连续Q高斯600对 |','| --- | ---: | ---: |']
    for method in METHODS:
        values=[next(r for r in report['aggregates'] if (r['domain'],r['method'],r['part'])==(d,method,'all')) for d in ('native','continuous')]
        label=method+('（门槛敏感性）' if method not in PRIMARY else '')
        lines.append('| '+label+' | '+' | '.join(f"{r['success']}/{r['frames']}（{100*r['efficiency']:.2f}%）" for r in values)+' |')
    return '\n'.join(lines)


def notebook():
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    assert load(RECORD/'science-validation.json')['passed']
    assert load(RECORD/'figure-visual-review.json')['passed']
    previous=NOTEBOOK.read_bytes();book=nbformat.reads(previous.decode(),as_version=4)
    before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4)
    assert len(book.cells)==len(before.cells)==455 and book.cells==before.cells
    report=load(OUTPUT/'report.json');cells=[]
    def md(text):cells.append(nbformat.v4.new_markdown_cell(text))
    def code(text):cells.append(nbformat.v4.new_code_cell(text))
    def figure(number,title):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{title}**')
        code(f'display(Image(filename=str(gate_output / "figures/{number}.png")))')

    md('''<a id="physical-gate-design-20261027"></a>
## 93. 网络只提供起点，事件数交给物理似然

上一轮混合流程同时改变了K=2候选与事件数判决。本轮将二者分开，回答：正确候选是否被学习型分数拒绝？增加网络起点是否比只追加一次局部物理精修更有用？

范围仍为纯α、固定95／351 ns双τ、2 ns×1024完整模拟电压、充分动态范围、局部未知K≤2。Q为等效pe：原生数据使用实现雪崩电荷，高斯数据使用生成均值／协方差的条件Q；均未连接沉积能标。t为有效面首光子时刻。本轮没有新训练、新网络推断或新优化，只重放并选择上一轮缓存候选；所有物理NLL先转float64再扣基线，2040帧代价重放最大相对差为0。

| 流程 | K=2候选 | K=2判决 |
| --- | --- | --- |
| baseline | 原物理拟合 | 原物理NLL门槛 |
| learned_gate | 上轮未加Q正则的MLP＋NLL候选 | 上轮学习型logit门槛，精确重放对照 |
| local_physical | 原K2起步的一次70步局部精修，并保留原候选 | 物理NLL差，重新校准 |
| mlp_physical | 原起点与网络起点两次70步精修，并保留两个初值，按原NLL取最优 | 物理NLL差，重新校准 |
| mlp_fixed_gate | 与mlp_physical完全相同 | 沿用原数值门槛，仅作敏感性对照 |

保留K0／K1与有信号判据：单候选Q≥500 pe且NLL1<NLL0。物理增峰证据为max(NLL1−NLL2,0)，要求两个候选Q均≥500、证据>0；在同一280个原生单事件上，按观测初拟合单Q的1500／5000／12500 pe分档，用保留并列的加一上尾秩p≤0.025选择K2。搜索方式改变后重新校准同一种规则；“规则相同”不等于“数值门槛相同”，末行专门检查数值门槛影响。物理判决不读取网络logit，也不加入网络Q惩罚或额外Poisson先验。

若最低NLL候选含Q<500，仍按原规则拒绝，不暗中换成更差但过门槛的另一候选。局部对照只有一次额外起点，网络混合有两次；这是检验追加网络候选的实际价值，不是等计算预算的结构优劣试验。

输入为280门槛帧＋760已公开原生评价帧＋1000上一轮已公开的连续Q高斯帧。本轮两种评价域都属于开发复用，没有新的独立检验。原生360对覆盖三组合×8／40／200 ns；连续Q600对含三组合×三间隔各40对、另240对均匀Q分布改变。两域各400空／单控制。合计人为条件配比不能换算成源注入率。

严格整对成功要求每Q误差≤20%、每t误差≤min(10 ns,真实间隔/4)、无额外输出；所有目标真实Q≥500 pe。以下同时保留漏检、额外输出、误拆／误并、偏差及误差宽度。''')
    code(f'''from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Image
gate_root = Path.cwd().resolve()
while not (gate_root / "pileup/experiments/{NAME}").is_dir():
    if gate_root.parent == gate_root:
        raise RuntimeError("请从项目根目录或jupyter目录运行")
    gate_root = gate_root.parent
gate_output = gate_root / "jupyter/artifacts/{NAME}"
gate_record = gate_root / "jupyter/run_records/{NAME}"
gate_report = json.loads((gate_output / "report.json").read_text())
gate_summary = pd.DataFrame(gate_report["summary"])
print("门槛280帧；评价1760帧。所有评价数据本轮均已公开；无新训练／优化／光学母事件。")
display(pd.DataFrame(gate_report["aggregates"]).query("part == 'all'"))
''')
    figure('0004-93-1','原生开发复用的严格恢复率；(a)中／中，(b)先强后弱，(c)先弱后强，各点40对。斜线为沿用原数值门槛的敏感性结果。')
    figure('0004-93-2','已公开连续Q高斯样本；(a)近等强，(b)先强后弱，(c)先弱后强各点40对，(d)均匀Q分布改变240对。')

    md('''<a id="physical-gate-stages-20261027"></a>
## 94. 正确候选在哪一步丢失？

对每个K2候选，先仅在评分阶段用真值判断Q和时刻是否正确，再统计最终判决去向：接受、Q阈下、无单事件、无双候选、非正证据或秩门槛拒绝。强制K=2诊断不是盲效率，也不是所有算法的性能上界。

先强后弱200 ns：原候选、局部精修及网络混合都各有20/40个严格正确K2。学习型门槛只接受其中7个，另13个由秩门槛拒绝；相同网络混合候选换回物理判决后20个全部接受。因此7→20来自判决改变，不能归因于新拟合。

先强后弱40 ns：三种候选库均只有3/40严格正确K2，且全部被各自门槛拒绝；最终仍0/40。即使知道该保留这三个候选，也无法解释其余37个失败。当前候选存在定位／电荷分配误差，门槛也造成损失，两者应分开处理。原生8 ns三个库均0/40正确K2；这些只是本数据／有限搜索结果，不能断言物理上不可分。

物理门槛相对学习型门槛的原生200 ns严格成功差为+32.5个百分点，母事件分块配对95%区间为+17.5至+47.5个百分点。每条件2000次bootstrap：原生按实际母事件块（强弱20块、每块两对），高斯按独立抽样帧。它是已公开数据的描述区间，退化区间不代表总体确定；未将跨条件复用的原生帧当作独立样本。''')
    figure('0004-94-1','(a,c)依次为最终选K2、两时刻均正确、Q和时刻全正确；(b,d)仅统计严格正确K2候选的去向。上排40 ns，下排200 ns。')
    code('''gate_stages = pd.DataFrame(gate_report["stages"])
display(gate_stages.query("domain == 'native' and condition == 'high_weak'"))
gate_forced = pd.DataFrame(gate_report["forced_summary"])
display(gate_forced.query("domain == 'native' and condition == 'high_weak'"))
print("forced_correct/success为强制K2诊断；strict_success才是经最终判决的严格盲成功数。")
''')

    md('''<a id="physical-gate-controls-20261027"></a>
## 95. 控制误拆，以及网络是否贡献了额外正确事件

预定“修复门槛”要求：原生200 ns先强后弱≥20/40，原生总成功≥241/360，高斯总成功≥500/600；每个空／单组额外输出帧≤5%、单事件严格成功≥95%、匹配Q偏差中位绝对值≤5%。局部＋物理、网络＋物理及固定数值门槛敏感性均通过这些开发条件。

另外单列“网络有实际增益”的要求：网络＋物理相对局部＋物理，在两个域都至少多5个百分点且通过上述控制。本轮不满足：原生都是241/360，高斯500/600→501/600，仅多1对。新增正确对位于先强后弱128–256 ns组，组内增加1/40（+2.5个百分点，配对区间0至+7.5），未见原生正确对增加。不能把恢复到基线说成超越基线。

网络起点使原生有K2库的686帧中92帧、连续Q920帧中136帧的NLL低于局部对照超过10⁻⁶；两域下降量中位均为0。数值代价更低不保证真实Q/t更准。原生物理新流程的空／weak／low／mid／high额外输出为0／0／0／1／0（每组80）；高斯局部为0／0／0／2／1，网络＋物理为0／0／0／2／0。相对原高斯基线mid的1/80，新流程mid为2/80；虽然通过5%要求，不能称所有条件的误拆都下降。

匹配Q偏差、sigma68=(p84−p16)/2及时间误差在下表完整保留，sigma68是时间匹配条件下的相对Q误差宽度，非沉积能谱FWHM。普通漏检／额外输出使用10 ns时间匹配；严格整对与both_times使用min(10 ns,gap/4)。未收敛结果未删除。

本轮2040帧缓存组装与代价重放约11.48 s，单worker／BLAS1；这是重放开销，不包括原物理拟合、特征提取或网络训练，不能作为端到端吞吐。''')
    figure('0004-95-1','空／单控制，各组80帧；(a–c)原生，(d–f)高斯。三列分别为额外输出帧率、严格成功率和匹配Q偏差中位。')
    code('''with pd.option_context("display.max_rows", None, "display.max_columns", None):
    display(gate_summary[["domain", "method", "kind", "condition", "band", "frames", "success", "missed", "extra", "extra_frames", "split", "merge", "q_error_p16_p50_p84", "q_error_sigma68", "time_error_p16_p50_p84_ns", "not_converged"]])
gate_paired = pd.DataFrame(gate_report["paired"])
display(gate_paired.query("method == 'mlp_physical' and reference == 'local_physical'"))
print("difference／paired95乘100为百分点；比较条件内配对，不能视为独立物理确认。")
display(pd.DataFrame(gate_report["search"]).T)
''')

    md('''<a id="physical-gate-waveforms-20261027"></a>
## 96. 完整堆积波形和最终选出的分量

在评分前指定原生先强后弱40／200 ns的第0对，显示完整0–2048 ns记录；实际采样点为0–2046 ns、步长2 ns。左列总电压，中列较大Q输出，右列较小Q输出。真实原生分量与真参数条件均值只用于画图，未输入盲判决；较小Q输出缺失时明确画零，不能拿被拒绝的K2冒充恢复。所选输出按Q大小归列，不保证与真强／弱配对正确。

200 ns示例展示学习型判决选K1而物理判决选K2，但选K2本身仍不保证两个Q/t合格；以紧随图后的实际评分为准。全部记录与失败例保留，没有按成功与否挑图。''')
    figure('0004-96-1','全长波形；(a–c)40 ns，(d–f)200 ns；原基线、学习型判决、局部＋物理、网络＋物理。缺失第二输出画零线。')
    code('''gate_examples = json.loads((gate_output / "examples.json").read_text())
display(pd.DataFrame([dict(间隔_ns=e["gap_ns"], 方法=m, K=len(v["fit"]["times_ns"]), Q_pe=v["fit"]["charge_pe"], 时刻_ns=v["fit"]["times_ns"], 严格成功=v["score"]["success"], 额外输出=v["score"]["extra"]) for e in gate_examples for m,v in e["methods"].items()]))
''')

    md(f'''<a id="physical-gate-conclusions-20261027"></a>
## 97. 结论、依据与下一步

{result_table(report)}

本轮完成了学习型判决损失的定位与修复对照，未证明当前MLP能提高短间隔逐事件恢复。网络仅作初始化时基本回到原物理基线，保持原默认流程。原生40 ns的正确候选本就少，下一步优先量化这些帧的双事件证据与Q分配区间：在空／单误拆控制下报告完整的恢复率—额外输出取舍，检查原生弱事件偏差和区间覆盖；再判断是否需要多解输出、不同估计目标或响应模型。不能用已知真值挑候选／挑门槛，也不能把当前0/40当普适间隔极限。

仍未验证新光学母事件、训练随机性、位置／噪声变化、混合α/γ、free-τ及完整随机高率K≥3。连续Q高斯数据在上一轮虽于权重冻结后生成，本轮已经公开复用，不能再次标作新测试。没有新ROOT、随机抽样、模型参数或默认算法修改。

{LITERATURE}

已执行6项数值测试、32帧先导、全2040帧候选重放；7040条四主流程评分，另1760条数值门槛敏感性评分，2880条强制K2诊断。旧基线和学习型流程逐条评分精确复现，所有正确候选拒绝去向加和核对。之前455单元与输出保留，新增五幅固定图。

完整[逐组结果](artifacts/{NAME}/summary.csv)、[候选去向](artifacts/{NAME}/stage-diagnostics.csv)、[已知K2诊断](artifacts/{NAME}/forced-K2-summary.csv)、[配对区间](artifacts/{NAME}/paired-comparison.json)与[结果记录](artifacts/{NAME}/report.json)可直接复核。''')
    code('''display(pd.DataFrame(gate_report["decisions"]))
print("网络额外价值要求：", gate_report["network_practical_value"])
gate_checks = json.loads((gate_record / "science-validation.json").read_text())
print("核验通过：", gate_checks["passed"], "；主评分：", gate_checks["primary_scores"], "；含敏感性：", gate_checks["all_scores"], "；已知K2诊断：", gate_checks["forced_K2_scores"])
print("原基线／学习型流程精确重放：", gate_checks["old_baseline_and_learned_scores_exact"], "；NLL重放最大相对差：", gate_checks["input_cost_replay_relative_error"])
''')
    for i,cell in enumerate(cells):cell.id=f'physical-gate-20261027-{i:02d}'
    InteractiveShell.instance();namespace={};executed=0;cwd=os.getcwd()
    try:
        os.chdir(PROJECT/'jupyter')
        for cell in cells:
            if cell.cell_type!='code':continue
            with capture_output(stdout=True,stderr=True,display=True) as capture:
                exec(compile(cell.source,'<physical-gate-notebook>','exec'),namespace)
            assert not capture.stderr,capture.stderr
            cell.outputs=[]
            if capture.stdout:cell.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=capture.stdout))
            for out in capture.outputs:cell.outputs.append(nbformat.v4.new_output('display_data',data=out.data,metadata=out.metadata))
            executed+=1;cell.execution_count=executed
    finally:os.chdir(cwd)
    book.cells.extend(cells);assert book.cells[:455]==before.cells
    anchors=re.findall(r'id="codex-figure-([^"]+)"','\n'.join(c.source for c in book.cells if c.cell_type=='markdown'))
    assert len(anchors)==len(set(anchors))
    for c in cells:
        if c.cell_type=='markdown':assert not re.search(r'\]\([^)]*\.(?:pdf|ipynb|md)(?:[#?)][^)]*)?\)',c.source,re.I)
    images=[hashlib.sha256(base64.b64decode(o.data['image/png'])).hexdigest() for c in cells if c.cell_type=='code' for o in c.outputs if 'image/png' in o.get('data',{})]
    figures=load(RECORD/'figure-delivery.json')['figures']
    assert len(images)==len(figures)==5 and set(images)==set(figures.values())
    nbformat.validate(book);encoded=(nbformat.writes(book)+'\n').encode()
    with tempfile.NamedTemporaryFile(dir=NOTEBOOK.parent,prefix=NOTEBOOK.stem+'.',suffix='.partial',delete=False) as stream:
        stream.write(encoded);stream.flush();os.fsync(stream.fileno());temporary=Path(stream.name)
    assert temporary.read_bytes()==encoded
    nbformat.validate(nbformat.reads(encoded.decode(),as_version=4))
    assert NOTEBOOK.read_bytes()==previous,'Concurrent Notebook change; preserve user edits'
    temporary.replace(NOTEBOOK);assert NOTEBOOK.read_bytes()==encoded
    save(RECORD/'notebook-delivery.json',dict(passed=True,previous_cells=455,new_cells=len(cells),executed_cells=executed,
        figures=figures,history_preserved=True,notebook_sha256=digest(NOTEBOOK),source_sha256=digest(__file__)))
    print('PHYSICAL_GATE_NOTEBOOK_DELIVERED',len(cells),executed,flush=True)


def documents():
    report=load(OUTPUT/'report.json');assert load(RECORD/'notebook-delivery.json')['passed']
    table=result_table(report)
    summary='''在相同缓存拟合候选上分离事件数判决：原生先强后弱200 ns的学习型门槛拒绝13个严格正确K2，改用物理似然判决后严格恢复7/40→20/40；40 ns仍0/40。追加局部精修／网络起点＋物理判决，原生均241/360，高斯分别500/600、501/600。门槛修复通过开发控制，但网络相对局部精修的增益未达两域各5个百分点要求；保持原默认流程。'''
    root=f'''<a id="{MARKER}"></a>
2026-10-27完成[物理事件数判决与网络起点对照](pileup/experiments/{NAME}/README.md)：{summary} 全部数据本轮均为已公开开发复用，没有新训练／拟合／ROOT。6项数值检查、2040帧代价重放、7040主评分／8800含敏感性评分和2880已知K2诊断核验；Notebook0004第93–97节保存五幅固定图、完整波形与候选拒绝去向，依据及边界见[L62](docs/literature-use.md#{MARKER})。'''
    literature=f'''<a id="{MARKER}"></a>
### L62：2026-10-27——物理判决、学习型分数和网络候选的分离对照

**依据与审核边界**：{LITERATURE}

**项目实施**：{summary}

{table}

纯α、固定95／351 ns、2 ns×1024、充分动态范围、局部未知K≤2；Q为等效pe、未作沉积能标。复用relative网络的未正则化物理NLL候选；local只取原K2起点的70步精修，mlp额外含网络起点的70步精修，两者都保留初值。没有新优化或训练，也不是等计算预算的网络结构比较。K0／K1及单事件存在规则不变；物理增峰分数max(NLL1−NLL2,0)，每候选Q≥500、正证据、同280原生单事件按拟合单Q分档1500／5000／12500 pe作加一上尾秩≤.025。搜索改变后重新校准同一种规则；沿用原数值门槛单列敏感性，不混进四主流程。物理判决不使用网络logit或学习型Q惩罚；最低代价候选阈下时不偷换高代价过门候选。

**故障分解**：原生200 ns三个库各20/40正确K2，学习型门槛拒绝13个，物理门槛全部接受；配对恢复差+32.5个百分点，母事件块95%区间+17.5至+47.5。40 ns三个库均只有3/40正确K2且都被秩门槛拒绝，最终0/40，显示候选参数与判决两方面问题，不能据此推导普适分辨极限。网络相对局部对照的新增成功仅高斯先强后弱128–256 ns的1/40，组内差+2.5个百分点、区间0至+7.5。每条件2000次bootstrap，种子270327062；原生按母事件块，高斯按独立抽样。区间为已公开开发样本描述，跨条件复用及模型误差不作独立帧推断。

**控制与代价**：原生空／weak／low／mid／high额外输出帧，local和mlp均0／0／0／1／0；高斯local为0／0／0／2／1、mlp为0／0／0／2／0，每组80。高斯基线mid为1/80，故不可称所有误拆都降低。每组额外≤5%、单事件严格成功≥95%、匹配Q偏差中位绝对值≤5%的开发要求均通过。每Q±20%、每t±min(10 ns,gap/4)、无额外为严格整对成功；另保留10 ns匹配下的漏检／额外／误拆／误并及Q偏差、sigma68和时间误差。sigma68非能谱FWHM。网络降低NLL>10⁻⁶的帧为原生92/686、高斯136/920，两域下降量中位为0，不能据数值改善宣称真实恢复。

**记录与验证**：280校准＋760原生＋1000高斯共2040帧全部代价重放，最大相对差0；6项数值测试、32帧先导、7040主评分＋1760门槛敏感性、2880强制K2诊断，旧基线和学习型流程精确重现。采用上一轮缓存、全部未收敛结果保留，推断禁读truth／case／ROOT。单worker、BLAS1的11.48 s只为缓存重放，不作端到端吞吐。实验为[alpha_pair_physical_gate_20261027](../pileup/experiments/{NAME}/README.md)，本地jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb第93–97节，图0004-93-1／93-2／94-1／95-1／96-1。旧455单元、输入、模型与此前源码保留。

**下一步及待核查**：优先量化弱Q分配歧义及区间覆盖，在空／单控制下展示恢复率—额外输出曲线，再评估多解或不同估计目标；不再把单纯更多起点当作已证实方向。学习型分数需另行校准而非当作物理似然。新增统计／区间校准的一手文献整理、Liddick门选异文及Morad统计量疑问仍交共享库处理。当前没有新母事件、噪声／位置变化、free-τ、混合α/γ或随机高率K≥3验证，全部评价已公开，不更新默认流程。'''
    readme=f'''# Physical count decision and ML candidate ablation

{summary}

{table}

## 范围与数据边界

纯α，固定双τ 95／351 ns，完整1024×2 ns模拟电压，充分动态范围，未知局部K≤2。原生Q为实现雪崩等效pe，高斯Q是均值／协方差条件参数；没有沉积能标。保持物理响应、噪声模型、K0／K1与源数据不变。

输入来自上一实验alpha_pair_risk_training_v2_20261027：9600训练／960验证所得relative网络及已生成的候选；本轮未重新训练、推断网络、优化或生成样本。280原生门槛帧、760原生评价、1000连续Q高斯评价在本轮均为开发复用。训练／验证／高斯生成种子基数270310001／270310002／270310003；权重选择与母事件分组、输入哈希完整继承于本地input-provenance.json。没有新的G4宏或ROOT。当前代码HEAD及依赖记录于environment.json，未提交的先前研究与本轮源码分别冻结保留。

## 对照与判决

baseline为原候选与原物理门槛；learned_gate为与上轮mlp_nll逐条相同的混合候选及学习型门槛；local_physical从缓存trial0取原K2起点的70步局部精修；mlp_physical额外包含网络起点的70步精修，均保留原初值。两种新物理流程以同280校准帧重新校准同一种NLL差规则；mlp_fixed_gate沿用原数值门槛，仅作敏感性。一个与两个起点不是等预算结构比较。

单事件存在需要Q≥500且NLL1<NLL0。K2需两个Q均≥500、max(NLL1−NLL2,0)>0；观测初拟合单Q按1500／5000／12500分档，加一上尾秩≤.025。并列保留，不读取网络logit、不加Q正则或重复源率先验；阈下最低代价候选不被悄悄换掉。详见protocol.json与gate-calibrations.json。

严格整对正确：每Q±20%、每t±min(10 ns,gap/4)、无额外。普通漏检／额外等使用10 ns时间匹配；时间匹配后的相对Q偏差／sigma68和时间误差、误拆／误并、未收敛均保留。强制K2只诊断候选。预定门槛修复和网络增益要求分开，后者未通过。不改默认算法。

原生先强后弱200 ns：三个库各20/40正确K2，学习型门槛只保留7，物理保留20；40 ns三个库均3/40正确候选且被拒，最终0/40。原生control额外0/0/0/1/0（每组80），高斯local为0/0/0/2/1、mlp为0/0/0/2/0；高斯mid比基线多一帧，未声称全面降低误拆。网络相对local只在高斯长间隔强→弱新增1个成功，不能证明短间隔改善。

配对bootstrap2000次、种子270327062；原生条件内按实际母事件块，高斯独立抽样；已公开开发结果不作物理确认。只有局部K≤2／固定τ／纯α，不提供最高计数率或普适最小间隔。

## 复现与输出

从项目根目录使用项目Python环境运行；设置OPENBLAS_NUM_THREADS=1、OMP_NUM_THREADS=1。run.py使用inference-source，report.py使用report-source，deliver.py使用delivery-source；每层校验冻结哈希，不需要切Git版本。

```bash
bash jupyter/run_gagg_python.sh pileup/experiments/{NAME}/report.py analyze
bash jupyter/run_gagg_python.sh pileup/experiments/{NAME}/report.py figures
bash jupyter/run_gagg_python.sh pileup/experiments/{NAME}/deliver.py validate
```

原始prepare／pilot／infer为一次性运行器，拒绝覆盖原目录／缓存。已有缓存可重复analyze／figures；重建缓存时应复制实验及依赖输入到独立项目输出树，保持本记录。Notebook与文档交付亦拒绝覆盖用户后续修改。

本地结果jupyter/artifacts/{NAME}/，记录jupyter/run_records/{NAME}/；其report／summary／scores／forced-K2／stage-diagnostics／paired-comparison等文件保留全部结果。Notebook0004第93–97节、固定图0004-93-1／93-2／94-1／95-1／96-1，40／200 ns第0对完整0–2048 ns波形，缺失第二输出为零。

6项数值检查、32帧先导、2040帧代价重放（相对差0）、7040主评分／8800含敏感性评分、2880已知K2诊断；旧基线与学习型评分精确复现。五图目视检查、Notebook新单元实际执行，旧455单元原样保留；最终完整校验记录见final-validation.json。输入和各层源码哈希追溯，未收敛候选不删除，推断阻止读取标签。缓存重放11.48 s不是完整拟合吞吐。

## 文献采用与后续

已有Liddick2012 DOI 10.1016/j.nima.2011.12.034（硅DSSD、40 MSPS）仅支持将检测与参数精度分开、控制假双；AI精读待人工审核，门选记号有未决项。Morad2024 DOI 10.1109/TNS.2024.3371573（CLYC、1 GS/s）仅借鉴联合指标，人工审核553为限定确认。当前实验不是其复现，未重查PDF或修改文献库；具体位置、条件、疑问与采用边界见项目docs/literature-use.md的L62。

下一步分析弱Q歧义、证据与区间覆盖，公开误拆—恢复取舍；尚需新母事件、噪声／位置、α/γ、free-τ及高率多事件验证。当前网络不是已证实更好的默认解法。
'''
    for name,backup,addition in (('README.md','README.before.md',root),('docs/literature-use.md','literature.before.md',literature)):
        path=PROJECT/name;original=(RECORD/backup).read_text()
        assert path.read_text()==original and MARKER not in original
        path.write_text(original.rstrip()+'\n\n'+addition.rstrip()+'\n')
    path=EXPERIMENT/'README.md';assert not path.exists();path.write_text(readme)
    files={str(p.relative_to(PROJECT)):digest(p) for p in (PROJECT/'README.md',PROJECT/'docs/literature-use.md',path)}
    save(RECORD/'document-delivery.json',dict(passed=True,files=files,source_sha256=digest(__file__)))
    print('PHYSICAL_GATE_DOCUMENTS_DELIVERED',flush=True)


def validate():
    check_freeze()
    for name in ('numerical-tests','pilot-validation','inference','science-validation','figure-delivery','figure-visual-review','notebook-delivery','document-delivery'):
        assert load(RECORD/(name+'.json'))['passed'],name
    manifests={name:load(EXPERIMENT/(name+'.json')) for name in ('freeze','report-freeze','delivery-freeze')}
    for manifest in manifests.values():
        for name,expected in manifest['files'].items():assert digest(EXPERIMENT/name)==expected,name
    groups={'inference-source/pileup/analysis':('alpha_pair_physical_gate.py','alpha_pair_physical_gate_study.py','test_alpha_pair_physical_gate.py'),
        'report-source':('alpha_pair_physical_gate_report.py',),'delivery-source':('alpha_pair_physical_gate_delivery.py',)}
    source={}
    for folder,names in groups.items():
        for name in names:
            p=PROJECT/'pileup/analysis'/name;ast.parse(p.read_text());assert digest(p)==digest(EXPERIMENT/folder/name)
            source[str(p.relative_to(PROJECT))]=digest(p)
    provenance=load(RECORD/'input-provenance.json');prior=provenance['prior_validation']
    for name,expected in provenance['inputs'].items():assert digest(PROJECT/name)==expected,name
    for name,expected in provenance['previous_source_freezes'].items():assert digest(PREVIOUS/name)==expected,name
    for name,expected in prior['source_sha256'].items():assert digest(PROJECT/name)==expected,name
    for name,expected in prior['documentation_sha256'].items():
        if name not in ('README.md','docs/literature-use.md'):assert digest(PROJECT/name)==expected,name
    for old in ('alpha_pair_information_20261027','alpha_pair_charge_allocation_20261027'):
        for name,expected in load(PROJECT/'jupyter/run_records'/old/'final-validation.json')['source_sha256'].items():assert digest(PROJECT/name)==expected,name
    inference=load(RECORD/'inference.json');science=load(RECORD/'science-validation.json')
    for name,expected in inference['outputs'].items():assert digest(OUTPUT/name)==expected,name
    for name,expected in science['output_hashes'].items():assert digest(OUTPUT/name)==expected,name
    assert len(rows(OUTPUT/'scores.jsonl'))==8800 and len(rows(OUTPUT/'forced-K2.jsonl'))==2880
    assert science['primary_scores']==7040 and science['input_cost_replay_relative_error']==0.
    assert inference['truth_inputs']==[] and inference['label_access_blocked'] and inference['no_new_optimization'] and inference['no_new_network_inference']
    program='''import sys,json
sys.path.insert(0,sys.argv[1])
from pileup.analysis.alpha_pair_charge_allocation import forbid_labels
sys.addaudithook(forbid_labels)
blocked=[]
for name in ('forbidden.truth.npz','forbidden.cases.json','forbidden.root','/tmp/truth/forbidden.npy','native-records.json'):
    try:open(name,'rb')
    except PermissionError:blocked.append(name)
    else:raise AssertionError(name)
print(json.dumps(blocked))
'''
    blocked=json.loads(subprocess.check_output([sys.executable,'-c',program,str(EXPERIMENT/'inference-source')],cwd=PROJECT,text=True));assert len(blocked)==5
    for domain,size in (('native',1040),('continuous',1000)):
        cached=rows(OUTPUT/(domain+'-candidates.jsonl'));assert len(cached)==size
        for r in cached:
            for bank in r['sets'].values():
                assert bank[:2]==r['sets']['baseline'][:2]
                assert all('network_score' not in c for c in bank)
    doc=load(RECORD/'document-delivery.json')
    for name,expected in doc['files'].items():assert digest(PROJECT/name)==expected,name
    for name,backup in (('README.md','README.before.md'),('docs/literature-use.md','literature.before.md')):
        assert digest(RECORD/backup)==prior['documentation_sha256'][name]
        assert (PROJECT/name).read_text().startswith((RECORD/backup).read_text().rstrip()+'\n\n')
    book=nbformat.read(NOTEBOOK,as_version=4);nbformat.validate(book)
    before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4)
    assert digest(RECORD/'notebook.before.ipynb')==provenance['notebook_before_sha256']==prior['notebook_sha256']
    assert len(before.cells)==455 and book.cells[:455]==before.cells
    delivery=load(RECORD/'notebook-delivery.json');new=book.cells[455:]
    assert len(new)==delivery['new_cells'] and digest(NOTEBOOK)==delivery['notebook_sha256']
    assert all(c.id.startswith('physical-gate-20261027-') for c in new)
    assert len({c.id for c in book.cells})==len(book.cells)
    anchors=re.findall(r'id="codex-figure-([^"]+)"','\n'.join(c.source for c in book.cells if c.cell_type=='markdown'))
    assert len(anchors)==len(set(anchors))
    codes=[c for c in new if c.cell_type=='code'];assert len(codes)==delivery['executed_cells']
    assert all(c.execution_count and c.outputs and not any(o.output_type=='error' for o in c.outputs) for c in codes)
    figures=load(RECORD/'figure-delivery.json')['figures'];assert figures==load(RECORD/'figure-visual-review.json')['figures']
    for number,expected in figures.items():assert digest(OUTPUT/'figures'/(number+'.png'))==expected
    embedded=[hashlib.sha256(base64.b64decode(o.data['image/png'])).hexdigest() for c in codes for o in c.outputs if 'image/png' in o.get('data',{})]
    assert len(embedded)==5 and set(embedded)==set(figures.values())
    for c in new:
        if c.cell_type=='markdown':assert not re.search(r'\]\([^)]*\.(?:pdf|ipynb|md)(?:[#?)][^)]*)?\)',c.source,re.I)
    for name in (str(NOTEBOOK.relative_to(PROJECT)),f'jupyter/artifacts/{NAME}/report.json',f'pileup/experiments/{NAME}/joint-response.npz'):
        assert subprocess.run(['git','check-ignore','-q',name],cwd=PROJECT).returncode==0,name
    subprocess.run(['git','-c','core.fileMode=false','diff','--check','--','README.md','docs/literature-use.md','pileup/analysis'],cwd=PROJECT,check=True)
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip()==load(RECORD/'environment.json')['git_head']
    save(RECORD/'final-validation.json',dict(passed=True,numerical_tests=6,all_frames=2040,calibration_frames=280,evaluation_frames=1760,
        primary_scores=7040,all_scores=8800,forced_K2_scores=2880,maximum_cost_replay_relative_error=0.,
        frozen_files={k:len(v['files']) for k,v in manifests.items()},source_sha256=source,truth_open_probes_blocked=blocked,
        original_455_cells_preserved=True,new_cells=len(new),executed_new_cells=len(codes),figures=figures,notebook_sha256=digest(NOTEBOOK),
        documentation_sha256=doc['files'],prior_information_allocation_risk_preserved=True,physical_gate_independent_of_network_score=True,
        new_ROOTs=0,new_training=0,new_optimization=0,new_independent_samples=0,no_default_change=True,no_git_commit_or_push=True,
        limitations=['exposed development reuse in both domains','one calibrated conditional Gaussian model','one versus two cached refinement starts, not equal compute','local K<=2; fixed tau; pure alpha']))
    print('PHYSICAL_GATE_FINAL_VALIDATION_PASSED',len(new),len(codes),flush=True)
