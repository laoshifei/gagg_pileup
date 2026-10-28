"""Notebook 0004, literature-use record and reproducibility checks."""
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

from .alpha_pair_evidence_study import PROJECT,EXPERIMENT,OUTPUT,RECORD,PREVIOUS,OLD_RECORD,NAME,check_freeze
from .alpha_pair_information_study import load,save,digest,rows

NOTEBOOK=PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
MARKER='alpha-pair-evidence-profile-20261027'
LITERATURE='''本轮通过paper-search只读复用本地文献库的主题、审核及完整来源笔记，没有联网、重查PDF或修改共享库。Liddick、Darby、Grzywacz，2012，Algorithms for pulse shape analysis using silicon detectors，NIMA669,70–78，DOI 10.1016/j.nima.2011.12.034，PDF3–8／图3–13；约65 μm DSSD、40 MSPS／12 bit、25 μs记录，选出的实测单脉冲缩放移位构造8000余人工双脉冲。借鉴检测、参数精度及假双控制分别报告；不移用其数值门槛、100 ns等幅结果或250 ns幅比5结果。来源papers/sources/2012_Liddick_NIMA_Algorithms-for-pulse-shape-analysis-using-silicon-detectors.md，仍为AI精读稿待人工审核，门选记号和参数数目疑问保留。

Baker与Cousins，1984，Clarification of the use of chi-square and likelihood functions in fits to histograms，NIM221,437–442，来源papers/sources/1984_Baker_NIM_Clarification-of-the-use-of-chi-square-and-likelihood-functions-in-fits-to-histograms.md；PDF1、PDF4§5–7、PDF5§8区分点估计、区间覆盖与拟合优度，提醒有限样本需按具体模型研究抽样分布。其Poisson／multinomial直方图模型不是本项目相关电压的统计模型，不移植计数似然或声称本文保证GAGG区间覆盖。该笔记暂无人工科学确认，DOI未由已有PDF核实，不补造。式3模型概率记号等原文疑问保留。本轮的秩扫描、触发局部参数域、剖面与1／4／9层级均为项目诊断设计，非上述文献复现。'''


def tables(report):
    threshold=['| 方法 | 开发控制允许的最大α | 原生40 ns：默认→放宽 | 原生200 ns：默认→放宽 | 原生总成功 | 高斯总成功 |',
        '| --- | ---: | ---: | ---: | ---: | ---: |']
    for r in report['sweep']['envelopes']:
        a,b=r['default_native'],r['native'];c=r['continuous']
        threshold.append(f"| {r['method']} | {r['alpha']:.5f} | {a['target_success']}→{b['target_success']}/40 | {a['long_success']}→{b['long_success']}/40 | {b['strict_success']}/360 | {c['strict_success']}/600 |")
    profile=['| 输入／间隔 | ΔJ≤4包含真实Q2 | 范围宽度／真实Q2中位 | 触及50 pe下界 | 分离区段 | 真Q探针降低参考最小值 |',
        '| --- | ---: | ---: | ---: | ---: | ---: |']
    for domain,band in (('native','40'),('native','200'),('continuous','middle'),('continuous','long')):
        r=next(r for r in report['profiles'] if (r['domain'],r['band'],r['level'])==(domain,band,4.))
        profile.append(f"| {r['domain']}／{r['band']} | {r['truth_included']}/40 | {r['width_over_trueQ'][1]:.3f} | {r['lower_boundary']}/40 | {r['disconnected']}/40 | {r['probe_lower_than_reference']}/40 |")
    return '\n'.join(threshold),'\n'.join(profile)


def conclusion(report):
    r=next(r for r in report['sweep']['envelopes'] if r['method']=='mlp_physical')
    short=r['native']['target_success'];permissive=r['all_native']['target_success']
    text=f"网络起点＋物理判决在两个域所有空／单组额外输出帧≤5%的开发约束下，原生40 ns严格成功为{short}/40；即使允许所有合资格候选通过也只有{permissive}/40。"
    text+='放宽判决只能接纳已有正确候选，不能修正其余候选的电荷／时刻错误。'
    rows4=[r for r in report['profiles'] if r['level']==4.]
    lower=sum(r['lower_boundary'] for r in rows4);disconnected=sum(r['disconnected'] for r in rows4);probes=sum(r['probe_lower_than_reference'] for r in rows4)
    text+=f"ΔJ≤4的160条条件剖面中，{lower}条触及50 pe下界、{disconnected}条存在多个分离区段；{probes}条真实Q2评分探针把参考最小代价再降低超过0.01。这些数值及有限搜索误差需与包含率一起解释，不将该层级称为95%置信区间。"
    return text


def notebook():
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    assert load(RECORD/'science-validation.json')['passed'] and load(RECORD/'figure-visual-review.json')['passed']
    previous=NOTEBOOK.read_bytes();book=nbformat.reads(previous.decode(),as_version=4)
    before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4)
    assert len(book.cells)==len(before.cells)==475 and book.cells==before.cells
    report=load(OUTPUT/'report.json');threshold,profile_table=tables(report);cells=[]
    def md(text):cells.append(nbformat.v4.new_markdown_cell(text))
    def code(text):cells.append(nbformat.v4.new_code_cell(text))
    def figure(number,title):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{title}**')
        code(f'display(Image(filename=str(evidence_output / "figures/{number}.png")))')

    md('''<a id="evidence-sweep-20261027"></a>
## 98. 放宽事件数判决，需要付出多少误拆？

本轮沿上一轮的候选库，固定四种流程的拟合参数、存在判据、Q≥500条件和原生单事件标定，仅扫描允许的上尾秩α。物理判决仍需正NLL差；学习型分数仍按自己原来的秩校准。没有用双事件真值调出一套新默认门槛。

门槛样本各拟合Q档为40／80／80／80帧，加一秩因此离散；扫描两个分母41／81的全部可取值，另含0、0.025、0.05、0.1、0.25、0.5、1，共127点。α=0.025精确复现上一轮四方法的7040条逐帧评分。α=1也不能越过Q<500、无单事件或物理非正证据的限制，因此与不加任何资格约束的已知K2不同。

两评价域均为开发复用：原生360对＋400空／单，高斯600对＋400空／单，另280门槛帧。图中横轴为该域五个空／单组中最大的额外输出帧率；每组80，额外输出并非只发生在空白区域。下表的“最大允许α”回看两个域共十个空／单组，要求每组≤5%；这是开发数据上的取舍边界，不是新的独立验证或未来误拆保证。所有漏检、额外、误拆／误并、匹配Q偏差／sigma68、时间误差与未收敛数保存在逐组CSV中。

严格双事件恢复仍需每Q±20%、每t±min(10 ns,间隔/4)、无额外输出；普通漏检／额外输出采用10 ns匹配。条件合计不代表某个注入率或能谱。''')
    code(f'''from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Image
evidence_root = Path.cwd().resolve()
while not (evidence_root / "pileup/experiments/{NAME}").is_dir():
    if evidence_root.parent == evidence_root:
        raise RuntimeError("请从项目根目录或jupyter目录运行")
    evidence_root = evidence_root.parent
evidence_output = evidence_root / "jupyter/artifacts/{NAME}"
evidence_record = evidence_root / "jupyter/run_records/{NAME}"
evidence_report = json.loads((evidence_output / "report.json").read_text())
print("127个秩门槛；1760评价帧；4方法；候选参数保持。全部评价是已公开开发复用。")
display(pd.DataFrame([dict(method=r["method"], alpha=r["alpha"], native40=r["native"]["target_success"], native200=r["native"]["long_success"], native_controls=r["native"]["control_extra"], Gaussian_controls=r["continuous"]["control_extra"], native_single_min=r["native"]["min_single_success"], Gaussian_single_min=r["continuous"]["min_single_success"]) for r in evidence_report["sweep"]["envelopes"]]))
''')
    md(threshold)
    figure('0004-98-1','恢复率与最差空／单组额外输出的取舍；(a,b)原生40／200 ns，(c,d)高斯24–64／128–256 ns。圆点为原α=0.025，竖虚线为5%。')
    figure('0004-98-2','原生先强后弱；(a–c)40 ns、(d–f)200 ns。三列依次为选K2、两时刻正确、两Q及两时刻全正确。横轴α放宽不改变拟合参数。')
    code('''display(pd.DataFrame(evidence_report["sweep"]["short_correct_candidates"]))
print("表中needed_alpha为让已有正确候选通过所需的秩门槛，控制数每组分母80；只是回看开发结果。")
evidence_thresholds = pd.read_csv(evidence_output / "threshold-summary.csv")
display(evidence_thresholds[evidence_thresholds.alpha == .025])
''')

    md(r'''<a id="evidence-profile-20261027"></a>
## 99. 扫描第二个事件电荷，重新拟合其余参数

单个最优Q无法表达波形允许的其他电荷分配。对预先选定的原生40／200 ns强→弱各40条、高斯24–64／128–256 ns各40条，共160条，强制K=2，扫描**按时间排序的第二槽Q2**。在这四个评价组中真实第二事件较弱，但算法没有读取这个真值关系，也没有强制Q2<Q1。

使用原均值与候选依赖协方差：$J=r^T C^{-1}r+\ln|C|$。在固定线性电子学变换后的条件高斯模型下，J等于$-2\ln L$加参数无关常数；不是只比较残差，也不是把旧约化GLS当似然。

对每个Q2，重新拟合Q1、t1、t2。Q域50–50000 pe，基础33个对数格点，加500、原候选Q2和网络混合Q2的0.5／0.8／1／1.2／2倍，裁到同一域。起点包含原／混合候选、总Q重新分配、相邻格点延续，以及0／4／8／16／32／64／128／256／512 ns间隔；按起始代价选两处，各最多50次L-BFGS-B迭代，正向、反向各扫一次，保留初值和未收敛结果。

两时刻均限定在观测单事件初拟合时刻附近[max(0,t_pilot−96),min(2046,t_pilot+512)] ns，并保证t1≤t2。这个**预定触发局部参数域**避免把记录末端的截断事件混进本问题；不同于旧候选在整个记录上的搜索范围。旧候选落在域外、真实时刻是否落在域内分别统计。它是明确有K=2条件的诊断，不能当成完整盲恢复效率，也不据此悄悄修改第98节候选。

再固定剖面所得一组拟合时刻与Q1，只改变Q2，得到“固定其他参数的切片”。切片也加入有限搜索包络，保证数值剖面不会比同Q切片更差；两条曲线的宽度差用于展示其他参数可调整时的歧义，固定切片并非正确误差条。

展示ΔJ≤1／4／9的相容范围；这些是预定代价层级，**不命名为68%／95%／99%区间**。Q网格在lnQ中线性插值求交点，保留所有分离区段、50／50000 pe边界标志及区间宽度。范围是有限起点、有限格点的数值近似，不保证全局最优。图99-1用真实Q2归一横轴只发生在展示阶段。
''')
    figure('0004-99-1','(a,b)原生40／200 ns，(c,d)高斯中／长间隔。绿色为剖面中位及跨40帧的16–84%带，橙色为固定其他参数的切片中位；不是总体置信带。')

    md('''<a id="evidence-inclusion-20261027"></a>
## 100. 用真值评分包含率，同时检查搜索和选择影响

主剖面保存并锁定哈希后，单独运行评分探针：只固定真实Q2，使用观测候选及主剖面邻点作起点，重新拟合Q1与两时刻。真时刻、真Q1没有传入探针。此步骤明确使用已知K=2与真实Q2，仅用于检查该Q是否落在指定代价层级；不把结果反馈到主剖面、门槛或盲拟合。

因此表中“包含真实Q2”由真实Q处的实际代价求得，不依赖网格插值。但它仍是该有限条件／有限搜索上的经验包含率：若探针代价比主剖面最小值低超过0.01，说明主剖面的格点或搜索遗漏了更低点，单独报告；探针未收敛同样保留。范围端点和宽度仍为格点插值近似。不能将当前包含比例当作经过独立校准的通用覆盖保证。

既报告全部40条，也报告原α=0.025最终选K2的子集，避免只看较容易被检出的波形。报告宽度=(上界−下界)/真实Q2时，真值仅作评价尺度；分离区段之间的空隙仍列出，不能误作连续支持。下表边界数说明截断范围，不能把50 pe边界当作数据给出的精确下限。
''')
    md(profile_table)
    figure('0004-100-1','(a)真实Q评分探针包含率；(b)相容范围总跨度／真实Q2的中位；(c)触及50 pe下界的比例；(d)全部波形与最终选K2子集的ΔJ≤4包含率。')
    code('''with pd.option_context("display.max_rows", None, "display.max_columns", None):
    display(pd.DataFrame(evidence_report["profiles"]))
print("conditional_best_success为强制K2剖面最优点的诊断，不是盲效率。")
print("true_times_outside、cached_outside、failed_trueQ_probes和probe_lower_than_reference必须与范围宽度一起阅读。")
''')

    md('''<a id="evidence-full-waveforms-20261027"></a>
## 101. 同一条完整波形允许哪些不同解释？

在计算前指定四组各第0条波形，画完整0–2048 ns记录；实际点为0–2046 ns、2 ns步长。左列为总电压，右列为第二槽分量。实线为主剖面最优格点，另画ΔJ≤4内格点中最小和最大Q2对应的完整模型，各自使用自己拟合的Q1/t1/t2。端点取实际计算的格点而非插值交点，数值以随后的表为准。

这些线是**已知K=2的相容模型示例**，与前面最终盲选输出图不同；它们不代表真的恢复了两个事件。黑色真参数条件均值仅用于比较，也不等于原生事件的随机雪崩分量。总波形看起来接近时，第二分量仍可能差很多，应结合时间、电荷和代价评分判断。
''')
    figure('0004-101-1','(a,b)原生40 ns，(c,d)原生200 ns，(e,f)高斯中间隔，(g,h)高斯长间隔。完整总电压与第二槽的条件K2解释，非盲恢复曲线。')
    code('''evidence_examples = json.loads((evidence_output / "examples.json").read_text())
display(pd.DataFrame([dict(domain=e["domain"], index=e["index"], true_gap_ns=e["true_gap_ns"], true_Q2_pe=e["true_Q2_pe"], model=m, fitted_Q_pe=f["charge_pe"], fitted_times_ns=f["times_ns"], cost=f["cost"]) for e in evidence_examples for m,f in e["models"].items()]))
''')

    md(f'''<a id="evidence-conclusions-20261027"></a>
## 102. 结果、文献边界与下一步

{conclusion(report)}

保持原默认算法。这轮只完成固定候选门槛取舍和已知K2下电荷相容性的诊断，没有新的训练、随机抽样或光学母ROOT。纯α、固定95／351 ns双τ、充分动态范围不变，未研究free-τ、混合α/γ、位置／噪声外推和高率K≥3。Q仍为等效pe，不是已标定沉积能量；相对Q误差sigma68非能谱FWHM。

接下来应先通过格点加密／更多观测起点检查相容范围的数值稳定性，再用冻结后的独立模型重复抽样校准需要报告的覆盖率，并检查原生迁移。若逐事件电荷确有多解，应同时给出不确定性或不可靠标志；不能在不知道真值时挑选恰好正确的一支。只有在误拆控制与独立样本都通过后，才评估更换实际重建流程。

{LITERATURE}

完整[门槛扫描](artifacts/{NAME}/threshold-summary.csv)、[条件剖面汇总](artifacts/{NAME}/profile-summary.csv)、[逐条剖面评分](artifacts/{NAME}/profile-scores.json)和[结果记录](artifacts/{NAME}/report.json)保留。原475单元与图像不变，新增五幅固定图。
''')
    code('''evidence_sweep_checks = json.loads((evidence_record / "sweep-validation.json").read_text())
evidence_profile_checks = json.loads((evidence_record / "science-validation.json").read_text())
print("门槛核验：", evidence_sweep_checks)
print("剖面核验：", {k:v for k,v in evidence_profile_checks.items() if k != "output_hashes"})
''')
    for i,c in enumerate(cells):c.id=f'evidence-20261027-{i:02d}'
    InteractiveShell.instance();namespace={};executed=0;cwd=os.getcwd()
    try:
        os.chdir(PROJECT/'jupyter')
        for c in cells:
            if c.cell_type!='code':continue
            with capture_output(stdout=True,stderr=True,display=True) as cap:exec(compile(c.source,'<evidence-notebook>','exec'),namespace)
            assert not cap.stderr,cap.stderr;c.outputs=[]
            if cap.stdout:c.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=cap.stdout))
            for out in cap.outputs:c.outputs.append(nbformat.v4.new_output('display_data',data=out.data,metadata=out.metadata))
            executed+=1;c.execution_count=executed
    finally:os.chdir(cwd)
    book.cells.extend(cells);assert book.cells[:475]==before.cells;nbformat.validate(book)
    encoded=(nbformat.writes(book)+'\n').encode()
    with tempfile.NamedTemporaryFile(dir=NOTEBOOK.parent,prefix=NOTEBOOK.stem+'.',suffix='.partial',delete=False) as stream:
        for start in range(0,len(encoded),256*1024):
            chunk=encoded[start:start+256*1024];assert stream.write(chunk)==len(chunk)
        stream.flush();os.fsync(stream.fileno());temporary=Path(stream.name)
    assert temporary.read_bytes()==encoded;nbformat.validate(nbformat.reads(encoded.decode(),as_version=4))
    assert NOTEBOOK.read_bytes()==previous,'Notebook changed; preserve user edits'
    temporary.replace(NOTEBOOK);assert NOTEBOOK.read_bytes()==encoded
    save(RECORD/'notebook-delivery.json',dict(passed=True,previous_cells=475,new_cells=len(cells),executed_cells=executed,
        history_preserved=True,full_readback_verified=True,chunk_bytes=256*1024,notebook_sha256=digest(NOTEBOOK),source_sha256=digest(__file__)))
    print('EVIDENCE_NOTEBOOK_DELIVERED',len(cells),executed,flush=True)


def documents():
    assert load(RECORD/'notebook-delivery.json')['passed']
    report=load(OUTPUT/'report.json');threshold,profile_table=tables(report);result=conclusion(report)
    root=f'''<a id="{MARKER}"></a>
2026-10-27完成[事件数门槛与条件电荷剖面](pileup/experiments/{NAME}/README.md)：固定四流程候选，扫描127个可取秩门槛；两域十个空／单控制组均≤5%额外输出时，四流程原生40 ns先强后弱仍0/40。另计算160条已知K2的第二槽Q2剖面及单独真实Q2评分探针，报告范围宽度、分离区段、边界与数值搜索遗漏，不把ΔJ层级称为置信水平。6项数值测试、7040条默认评分精确重放及全量剖面代价检查通过。全部数据为开发复用，不改默认算法；Notebook0004第98–102节、五幅固定图及[L63](docs/literature-use.md#{MARKER})保存完整结果。'''
    literature=f'''<a id="{MARKER}"></a>
### L63：2026-10-27——检测门槛取舍与条件电荷相容范围

**证据及审核状态**：{LITERATURE}

**项目实际结果**：{result}

{threshold}

表中α为回看两个域十个空／单控制组、每组80帧均≤5%额外输出时可取最大值，是开发约束边界，不作为新冻结默认。127点为同280标定帧的上尾秩全部可取值及预定0／.025／.05／.1／.25／.5／1，分档样本40／80／80／80，加一分母41／81。除α外候选、Q≥500、存在规则和物理正证据均不变；α=1不是无条件强制K2。物理及学习型门槛分别保留原校准。严格恢复每Q±20%、每t±min(10 ns,gap/4)、无额外；普通10 ns匹配的漏检／额外／误拆／误并及Q偏差／sigma68、时间误差、未收敛完整保存。sigma68不是沉积能谱FWHM。

**条件剖面**：原生先强后弱40／200 ns及高斯24–64／128–256 ns各40条，共160条，全部预定、不按成功挑选。强制K=2，Q2为按时间排序的后槽；每固定Q2重新拟合Q1和两时刻，使用原相关高斯J=rᵀC⁻¹r+ln|C|。两时刻均在观测单拟合触发附近[−96,+512] ns、裁至记录内，明确不同于旧候选全记录域；旧候选域外与真时刻域外分别计数。基础33个对数Q点50–50000 pe，加500、原Q2及混合Q2的五个比例点，正反扫描、每点两个观测起点各50迭代；保留初值、未收敛及固定其他参数切片，未证明全局最优。分离相容区段保留，lnQ插值端点仅为近似。

{profile_table}

ΔJ≤1／4／9为预定代价层级，不对应本项目已验证的68／95／99%覆盖。主剖面冻结后，另以真实Q2固定的评分探针重新拟合其余参数；不提供真t或真Q1，不反馈到主剖面、门槛或盲输出。真实Q处实际代价用于经验包含率；如果低于主参考最小值>0.01，标明格点／搜索遗漏，不能隐藏。区间宽度仍需加密和独立重复抽样检验。报告全部40条与原物理门槛选K2子集，避免忽略选择影响；Q真值只用于评分／画图归一。

**条件与验证**：纯α、固定95／351 ns、2 ns×1024完整模拟电压、充分动态范围、Q为等效pe而非沉积能量。门槛沿用280原生单事件；原生760、高斯1000评价均已公开。本轮无新训练／随机抽样／母ROOT。6项数值检查、4条先导、160剖面与160评分探针，主曲线代价重放、原α=.025四流程7040条逐帧评分精确复现；127×7040只是同一端点缓存的894080次选择，不是新增独立样本。所有未收敛保留，主剖面worker禁读truth／case／ROOT；真实Q2仅在另行标记的评分阶段传入。最多8核，实际主剖面7worker／BLAS1。

实验为[alpha_pair_evidence_profile_20261027](../pileup/experiments/{NAME}/README.md)；本地jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb第98–102节、图0004-98-1／98-2／99-1／100-1／101-1。五幅图含四条预定完整0–2048 ns波形；图101-1为已知K2相容模型，不能当作最终盲恢复。旧475单元与所有既有输入／结果保留。

**后续与文献待办**：先检查格点／搜索稳定性，再在冻结后的独立模拟重复抽样中校准所需覆盖并检验原生迁移；必要时给出多解与不可靠标志，禁止真值挑解。正式区间的有限样本／非规则模型选择一手统计证据、Baker文献DOI及相关边界理论可由共享库另行整理，本项目不修改库。新光学母事件、训练随机性、位置／噪声失配、α/γ、free-τ和随机高率K≥3尚未验证。'''
    experiment=f'''# Detection tradeoffs and conditional later-charge profiles

{result}

{threshold}

{profile_table}

## Experiment definition

Pure alpha, fixed95/351 ns, full1024×2 ns analog waveform, ample dynamic range. Original blind count fits allow K≤2; profiles explicitly assume K=2. Charge is equivalent pe, not calibrated deposited energy. All samples are previously exposed development reuse:280 native gate singles,760 native evaluation frames and1000 conditional-Gaussian frames. No new ROOT, training or random draws.

The threshold sweep reuses all four prior candidate/gate workflows.127 alpha points cover the attainable plus-one upper ranks with40/80/80/80 calibration samples. Presence,Q≥500 and positive physical evidence guards stay fixed. At alpha=.025 all7040 scored decisions exactly match the previous experiment. The observed ≤5% worst-control envelope is retrospective; no default is changed. Full errors and controls appear in threshold-summary.csv.

Profiles use160 prespecified strong-then-weak frames: native40/200 ns and Gaussian24–64/128–256 ns,40 each. For fixed chronological Q2, optimize Q1 and ordered t1/t2. Both times lie in the observed K1 trigger region[−96,+512]ns clipped to[0,2046]. Flag cached candidates outside this region.33 logarithmic50–50000 pe points plus observed-candidate points; two best observation-derived starts,50 iterations, forward/reverse continuation. The original candidate-dependent covariance and logdet are retained. Frozen-nuisance slices are included in the finite envelope. Disconnected ranges and boundary truncation are reported; global optimality is not claimed.

Levels1/4/9 are compatibility contours, not confidence levels. After primary profiles are frozen, a separate scoring probe fixes true later Q2 and refits nuisance using observed starts only. It never updates profiles, gates or default fits. Negative probe delta below−.01 flags a better point missed by the finite grid/search. Inclusion and interpolated width are empirical diagnostics, not calibrated coverage. All failed convergence flags are retained.

## Reproduction and provenance

Run from the project root with the project Python environment and OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1. Actual inference uses7 workers; the maximum is8. Input observations, candidates and truth come from the risk-training-v2/physical-gate studies; hashes, prior model/parent provenance and Git HEAD are in local input-provenance.json/environment.json. No new macro or stochastic seed is needed for this deterministic diagnostic. Earlier uncommitted work is retained.

run.py uses inference-source; report.py uses report-source; deliver.py uses delivery-source. Source, models, protocol and task indices are frozen with hashes, so no Git checkout is needed. Local arrays/models are ignored by Git and must accompany a reproducibility bundle.

```bash
bash jupyter/run_gagg_python.sh pileup/experiments/{NAME}/report.py sweep
bash jupyter/run_gagg_python.sh pileup/experiments/{NAME}/report.py analyze
bash jupyter/run_gagg_python.sh pileup/experiments/{NAME}/report.py figures
bash jupyter/run_gagg_python.sh pileup/experiments/{NAME}/deliver.py validate
```

prepare/pilot/infer/probe are one-time commands that reject existing outputs. To recompute inference, copy the frozen experiment and its inputs into a separate project tree. Notebook/document delivery refuses to replace subsequent user edits. Notebook writes use256 KiB chunks and complete readback before replacement for shared-folder compatibility.

Local artifacts: jupyter/artifacts/{NAME}/. Local records: jupyter/run_records/{NAME}/. Notebook0004 sections98–102, figures0004-98-1/98-2/99-1/100-1/101-1. All original475 cells retained; new code cells executed. Full0–2048 ns examples are prespecified; plotted profile components are conditional K2 explanations, not final blind reconstructions.

6 numerical tests,4 pilot profiles,160 primary profiles,160 true-Q2 scoring probes,7040 default-score replays and127-point sweep verified. final-validation.json records hashes, label-access probes, executed Notebook and image checks. Gate sweep involves894080 cached choices, not that many independent samples or fresh fits.

## Evidence and limitations

Liddick2012 (NIMA66970–78; DOI10.1016/j.nima.2011.12.034, siliconDSSD40MSPS) motivates separate detection/accuracy/false-double evaluation; its numerical limits do not transfer. Baker/Cousins1984 (NIM221437–442; DOI unverified in existing library note), PDF1/4/5, distinguishes estimation, coverage and goodness of fit under Poisson/multinomial histograms. No histogram likelihood is transplanted to voltage. Both scientific notes remain AI drafts pending review; source conditions and unresolved details are preserved in docs/literature-use.md L63. Only local existing notes were consulted; no PDF reread or library edits.

Before confidence claims, check grid/search stability and independent repeated-sampling calibration. No new optical parents, detector measurements, free-tau, alpha/gamma or full high-rate K≥3 validation; keep the current default reconstruction.
'''
    for name,backup,addition in (('README.md','README.before.md',root),('docs/literature-use.md','literature.before.md',literature)):
        path=PROJECT/name;old=(RECORD/backup).read_text();assert path.read_text()==old and MARKER not in old
        path.write_text(old.rstrip()+'\n\n'+addition.rstrip()+'\n')
    path=EXPERIMENT/'README.md';assert not path.exists();path.write_text(experiment)
    files={str(p.relative_to(PROJECT)):digest(p) for p in (PROJECT/'README.md',PROJECT/'docs/literature-use.md',path)}
    save(RECORD/'document-delivery.json',dict(passed=True,files=files,source_sha256=digest(__file__)))
    print('EVIDENCE_DOCUMENTS_DELIVERED',flush=True)


def validate():
    check_freeze()
    for name in ('numerical-tests','pilot-validation','inference','probe-validation','sweep-validation','science-validation','figure-delivery','figure-visual-review','notebook-delivery','document-delivery'):
        assert load(RECORD/(name+'.json'))['passed'],name
    manifests={k:load(EXPERIMENT/(k+'.json')) for k in ('freeze','report-freeze','delivery-freeze')}
    for manifest in manifests.values():
        for name,value in manifest['files'].items():assert digest(EXPERIMENT/name)==value,name
    source={}
    groups={'inference-source/pileup/analysis':('alpha_pair_evidence_profile.py','alpha_pair_evidence_study.py','test_alpha_pair_evidence_profile.py'),
        'report-source':('alpha_pair_evidence_report.py',),'delivery-source':('alpha_pair_evidence_delivery.py',)}
    for folder,names in groups.items():
        for name in names:
            path=PROJECT/'pileup/analysis'/name;ast.parse(path.read_text());assert digest(path)==digest(EXPERIMENT/folder/name)
            source[str(path.relative_to(PROJECT))]=digest(path)
    provenance=load(RECORD/'input-provenance.json');prior=provenance['prior_validation']
    for name,value in provenance['inputs'].items():assert digest(PROJECT/name)==value,name
    for name,value in provenance['prior_freezes'].items():assert digest(PREVIOUS/name)==value,name
    assert digest(OLD_RECORD/'input-provenance.json')==provenance['inherited_provenance_sha256']
    for name,value in prior['source_sha256'].items():assert digest(PROJECT/name)==value,name
    for name,value in prior['documentation_sha256'].items():
        if name not in ('README.md','docs/literature-use.md'):assert digest(PROJECT/name)==value,name
    for old in ('alpha_pair_information_20261027','alpha_pair_charge_allocation_20261027','alpha_pair_risk_training_v2_20261027'):
        for name,value in load(PROJECT/'jupyter/run_records'/old/'final-validation.json')['source_sha256'].items():assert digest(PROJECT/name)==value,name
    for record in ('science-validation','sweep-validation'):
        for name,value in load(RECORD/(record+'.json'))['output_hashes'].items():assert digest(OUTPUT/name)==value,name
    inference=load(RECORD/'inference.json');probe=load(RECORD/'probe-validation.json')
    assert digest(OUTPUT/'profiles.jsonl')==inference['sha256']==probe['primary_profiles_sha256']
    assert digest(OUTPUT/'truth-probes.jsonl')==probe['sha256']
    profiles=rows(OUTPUT/'profiles.jsonl');probes=rows(OUTPUT/'truth-probes.jsonl')
    assert len(profiles)==len(probes)==160
    plan={(r['domain'],r['index']) for r in load(EXPERIMENT/'tasks.json')}
    assert {(r['domain'],r['index']) for r in profiles}=={(r['domain'],r['index']) for r in probes}==plan
    assert inference['primary_truth_inputs']==[] and inference['label_files_blocked_in_workers']
    assert inference['maximum_cost_replay_error']<1e-9
    for curve in profiles:
        assert curve['reference_cost']==min(r['cost'] for r in curve['fits'])
        assert all(f['cost']<=s+1e-8 for f,s in zip(curve['fits'],curve['slice_cost']))
        assert len(curve['fits'])==len(curve['q_grid_pe'])
    assert len(rows(OUTPUT/'decision-records.jsonl'))==7040
    sweep=load(RECORD/'sweep-validation.json');assert sweep['alpha_count']==127 and sweep['equivalent_selected_scores']==894080
    program='''import sys,json
sys.path.insert(0,sys.argv[1])
from pileup.analysis.alpha_pair_evidence_study import initialize
initialize()
blocked=[]
for name in ('forbidden.truth.npz','forbidden.cases.json','forbidden.root','/tmp/truth/forbidden.npy','native-records.json'):
    try:open(name,'rb')
    except PermissionError:blocked.append(name)
    else:raise AssertionError(name)
print(json.dumps(blocked))
'''
    blocked=json.loads(subprocess.check_output([sys.executable,'-c',program,str(EXPERIMENT/'inference-source')],cwd=PROJECT,text=True));assert len(blocked)==5
    document=load(RECORD/'document-delivery.json')
    for name,value in document['files'].items():assert digest(PROJECT/name)==value,name
    for name,backup in (('README.md','README.before.md'),('docs/literature-use.md','literature.before.md')):
        assert digest(RECORD/backup)==prior['documentation_sha256'][name]
        assert (PROJECT/name).read_text().startswith((RECORD/backup).read_text().rstrip()+'\n\n')
    book=nbformat.read(NOTEBOOK,as_version=4);nbformat.validate(book)
    before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4);assert len(before.cells)==475 and book.cells[:475]==before.cells
    assert digest(RECORD/'notebook.before.ipynb')==prior['notebook_sha256']==provenance['notebook_sha256']
    delivery=load(RECORD/'notebook-delivery.json');new=book.cells[475:];codes=[c for c in new if c.cell_type=='code']
    assert digest(NOTEBOOK)==delivery['notebook_sha256'] and len(new)==delivery['new_cells'] and len(codes)==delivery['executed_cells']
    assert all(c.id.startswith('evidence-20261027-') for c in new) and len({c.id for c in book.cells})==len(book.cells)
    assert all(c.execution_count and c.outputs and not any(o.output_type=='error' for o in c.outputs) for c in codes)
    anchors=re.findall(r'id="codex-figure-([^"]+)"','\n'.join(c.source for c in book.cells if c.cell_type=='markdown'))
    assert len(anchors)==len(set(anchors))
    for c in new:
        if c.cell_type=='markdown':assert not re.search(r'\]\([^)]*\.(?:pdf|ipynb|md)(?:[#?)][^)]*)?\)',c.source,re.I)
    figures=load(RECORD/'figure-delivery.json')['figures'];assert figures==load(RECORD/'figure-visual-review.json')['figures']
    for number,value in figures.items():assert digest(OUTPUT/'figures'/(number+'.png'))==value
    embedded=[hashlib.sha256(base64.b64decode(o.data['image/png'])).hexdigest() for c in codes for o in c.outputs if 'image/png' in o.get('data',{})]
    assert len(embedded)==5 and set(embedded)==set(figures.values())
    for name in (str(NOTEBOOK.relative_to(PROJECT)),f'jupyter/artifacts/{NAME}/report.json',f'pileup/experiments/{NAME}/joint-response.npz'):
        assert subprocess.run(['git','check-ignore','-q',name],cwd=PROJECT).returncode==0,name
    subprocess.run(['git','-c','core.fileMode=false','diff','--check','--','README.md','docs/literature-use.md','pileup/analysis'],cwd=PROJECT,check=True)
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=PROJECT,text=True).strip()==load(RECORD/'environment.json')['git_head']
    save(RECORD/'final-validation.json',dict(passed=True,numerical_tests=6,profile_frames=160,truth_Q2_scoring_probes=160,
        original_default_scores_replayed=7040,rank_thresholds=127,cached_threshold_choices=894080,
        maximum_cost_replay_error=inference['maximum_cost_replay_error'],frozen_files={k:len(v['files']) for k,v in manifests.items()},
        source_sha256=source,truth_open_probes_blocked=blocked,original_475_cells_preserved=True,new_cells=len(new),executed_new_cells=len(codes),
        figures=figures,notebook_sha256=digest(NOTEBOOK),documentation_sha256=document['files'],prior_research_preserved=True,
        new_ROOTs=0,new_training=0,new_random_draws=0,no_default_change=True,no_git_commit_or_push=True,
        limitations=['exposed development reuse','forced K2 trigger-local profiles','finite multistart and interpolated contour endpoints','reference levels are not calibrated confidence intervals','pure alpha, fixed tau, local K<=2']))
    print('EVIDENCE_FINAL_VALIDATION_PASSED',len(new),len(codes),flush=True)
