"""Executed Notebook 0004 record of continuous-charge risk learning."""
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

from .alpha_pair_risk_study import PROJECT,EXPERIMENT,OUTPUT,RECORD,NAME
from .alpha_pair_information_study import load,save,digest


def deliver():
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    assert load(RECORD/'science-validation.json')['passed'] and load(RECORD/'figure-delivery.json')['passed']
    report=load(OUTPUT/'report.json');chosen=report['selected_training_risk']
    path=PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    previous=path.read_bytes();notebook=nbformat.reads(previous.decode(),as_version=4)
    before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4)
    assert len(before.cells)==431 and notebook.cells[:431]==before.cells
    if len(notebook.cells)>431:
        assert digest(path)==load(RECORD/'notebook-delivery.json')['notebook_sha256']
        assert all(c.get('id','').startswith('risk-20261027-') for c in notebook.cells[431:])
        notebook.cells=notebook.cells[:431]
    shutil.copyfile(path,RECORD/'notebook.before-delivery.ipynb')
    original=json.loads(json.dumps(notebook.cells));cells=[]
    def md(text):cells.append(nbformat.v4.new_markdown_cell(text))
    def code(text):cells.append(nbformat.v4.new_code_cell(text))
    def figure(number,title):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{title}**')
        code(f'display(Image(filename=str(risk_output / "figures/{number}.png")))')

    md(r'''<a id="risk-design-20261027"></a>
## 87. 改变估计目标：连续Q监督训练

前几轮增加同一似然的搜索，仍难恢复先强后弱的电荷。本轮实际检验：使用连续Q训练网络，并改变电荷损失函数，能否改善逐事件结果？再比较网络直接输出、网络起点＋原物理NLL、网络电荷约束＋物理NLL。

**本轮范围：**纯α、固定95/351 ns双τ、2 ns×1024点完整模拟电压、充分动态范围、局部未知K≤2。未生成新ROOT，未研究随机高率K≥3或混合α/γ。原生Q为实现雪崩等效pe；新高斯样本Q为生成均值／协方差的条件参数，不能等同新测量的雪崩电荷或沉积能量。时刻定义为有效面首光子。

训练9600帧、验证960帧，K=0/1/2抽样概率为0.1/0.3/0.6。单事件Q在500–25000 pe对数均匀；双事件一半独立对数均匀，四分之一弱500–3000＋强6000–25000 pe随机先后，四分之一近等强；单Q均在500–25000、总Q≤48000 pe。间隔一半在4–320 ns均匀，一半对数均匀，首时刻96–288 ns均匀。此分布是明确声明的模拟设计先验，不是实测源谱，也没有拿旧测试四档能量作网络答案。

训练目标是模拟器的两事件Q/t及事件数；模拟器可给出每个条件分量，但本轮训练参数头，并未训练逐点原生雪崩波形重建。输入只用观测波形推得的单事件初拟合残差、Q和时刻：1021个完整FIR创新、385个相对初拟合时刻的创新、lnQ和时刻，共1408维。创新做asinh(z/3)，仅用训练集均值／标准差归一化并截在±8。该数值预处理不是ADC削顶；原始完整电压仍保留。

一个NumPy MLP，隐藏层128/96，SiLU；三种损失共用初始化、训练数据和小批次顺序，各训练3000步、batch64、Adam学习率0.0007。没有另装深度学习框架。模型输出两个按时间排序的Q/t候选及K=2分数；时刻用相对初拟合−384至384 ns、步长2 ns网格的峰附近五点局部加权，避免把两个时间峰取成中间值。

总Q头为$S=Q_{pilot}\exp[\tanh(z)]$，第一事件份额为sigmoid，$Q_1=Sf,Q_2=S(1-f)$。初版总Q修正限制为exp(±0.1)；先导及训练输入检查发现5757个训练双事件中740个所需对数修正超过0.1，范围−0.759至0.425。因此在正式训练及测试生成前，将上限改为±1，另建v2实验目录；初版源码／数据完整保留。训练、验证数组与生成函数逐项核对，复用的是相同样本，没有重抽到更容易的例子。

三个目标均为BCE(K=2)+3×双事件Q损失+0.4×双事件时间交叉熵，时间软标签宽度2 ns。Q损失分别为：

| 名称 | 两事件平均的Q损失 | 含义 |
| --- | --- | --- |
| absolute | $[(\hat Q-Q)/5000]^2$ | 固定尺度的绝对平方误差 |
| logarithmic | $[\ln(\hat Q/Q)]^2$ | 对数误差，对大小比例对称 |
| relative | $(\hat Q/Q-1)^2$ | 相对平方误差，加强小Q的权重 |

三种目标对应不同估计偏好；都不保证在歧义波形中恢复真实弱事件。只训练一个初始化种子，尚未量化训练随机性。

文献仅复用已有精读：Morad等2024，DOI 10.1109/TNS.2024.3371573，PDF3–5式4–9／算法1–2和PDF10表II，支持检测／估计／分类分开评价；CLYC、1 GS/s、Poisson叠加＋AWGN的检出结果不是本项目电荷恢复保证。人工审核第553项仅限定确认联合指标区分及条件，未全面认可方法。Kim等2023 NIMA，DOI 10.1016/j.nima.2023.168492，PDF2–6／表1–2，TexAT、40 ns采样、100点对齐CNN波形MSE，提示弱成分和合成／物理迁移应分开；其第二信号常作为噪声目标，不能等同本项目两个α都需恢复。Kim科学内容仍为AI精读稿待人工审核。本轮未重查原PDF；当前MLP、三损失、总量头与学习型惩罚为项目实现，非原论文复现。
''')
    code(f'''from pathlib import Path
import json
import numpy as np
import pandas as pd
from IPython.display import display, Image
risk_root = Path.cwd().resolve()
while not (risk_root / "pileup/experiments/{NAME}").is_dir():
    if risk_root.parent == risk_root:
        raise RuntimeError("请从项目根目录或jupyter目录运行")
    risk_root = risk_root.parent
risk_output = risk_root / "jupyter/artifacts/{NAME}"
risk_record = risk_root / "jupyter/run_records/{NAME}"
risk_experiment = risk_root / "pileup/experiments/{NAME}"
risk_report = json.loads((risk_output / "report.json").read_text())
risk_summary = pd.DataFrame(risk_report["summary"])
risk_protocol = json.loads((risk_experiment / "protocol.json").read_text())
print("训练9600，验证960；冻结后新高斯评价1000；已有原生评价760，另门槛280。")
print("新母ROOT：", risk_report["new_ROOTs"], "；训练初始化种子数：", risk_report["training_seeds"])
print("训练／验证／评价种子：", risk_protocol["training"]["seeds"])
''')
    figure('0004-87-1','(a)训练、(b)验证、(c)冻结后评价的连续Q分布；散点仅作参数分布展示，不是新原生雪崩电荷。')

    md(f'''<a id="risk-validation-20261027"></a>
## 88. 验证集选择与物理拟合

每500步及第0步保存公共验证指标。选择顺序为：强制K=2严格整对成功率最高，其次Q绝对对数误差中位最低，再次事件数BCE最低。三种损失各保存自己的最佳权重；同一规则选出`{chosen}`供两个混合方案使用。这里明确已知K=2，仅作训练选择，不写成盲算法效率；测试真值没有用于选权重或宽度。

六种流程：原物理NLL基线、三种MLP直接候选、所选MLP＋原NLL、所选MLP＋带学习型Q惩罚的NLL。两种混合拟合均从原最佳K=2及网络K=2两处起步，各至多70次L-BFGS-B迭代，时间限于各起点±24 ns，并保留初值。惩罚为Σ[(lnQ−lnQ网络)/σ]²；σ按验证集全部双事件（包含错时刻）的68%绝对lnQ误差取值，下限0.2。网络和拟合使用同一波形，这是数据依赖正则项，不能当作独立贝叶斯先验再重复计证据。

K=2网络分数不作为似然，不再乘第二份Poisson先验。各流程在相同280个原生单事件、同一拟合单Q分档1500/5000/12500 pe上分别标定；保留并列的加一上尾秩p≤0.025且两个候选Q均≥500。K0/K1沿用原基线。因此流程差异同时含候选估计和检测分数，不能把所有增益都归给电荷拟合；后面另列强制K=2诊断。未收敛及失败均保留。
''')
    code('''risk_training = []
for risk_name in ("absolute", "logarithmic", "relative"):
    selection = json.loads((risk_experiment / "models" / risk_name / "selection.json").read_text())
    training = json.loads((risk_record / ("training-" + risk_name + ".json")).read_text())
    risk_training.append(dict(risk=risk_name, step=selection["step"], forced_pair_success=selection["forced_pair_success"], median_abs_logQ=selection["median_abs_logQ"], count_BCE=selection["count_BCE"], seconds=training["seconds"]))
display(pd.DataFrame(risk_training))
print("混合拟合所选损失：", risk_report["selected_training_risk"], "；两个时间槽的lnQ正则宽度：", risk_report["regularization_sigma"])
''')
    figure('0004-88-1','(a)已知K=2验证严格成功，(b)全部验证双事件Q绝对对数误差中位，(c)计数BCE；同样3000步，各自保存最佳检查点。')

    md('''<a id="risk-continuous-20261027"></a>
## 89. 冻结后连续Q评价

权重、损失选择与正则宽度冻结后才生成1000个新高斯样本。360对为近等强／先强后弱／先弱后强，各含6–12、24–64、128–256 ns三档、每格40对；弱Q为600–1500、强Q为12000–24000 pe。另240对采用两个Q各自在500–23000 pe均匀、间隔4–320 ns均匀，检查训练分布改变后的表现。还含80空帧、四档连续Q单事件各80帧。

严格整对成功要求两真实Q≥500 pe，各Q误差≤20%，各时刻误差≤min(10 ns,真实间隔/4)，且无额外输出。Q偏差／sigma68只统计时间匹配后误差；漏检、误拆／误并、额外输出及效率一起报告。sigma68=(p84−p16)/2不是沉积能谱FWHM。额外输出包括没有与真时刻匹配的候选，也可能发生在确有信号的区域。

这些是新参数／噪声种子的模型内评价，仍条件于同一已标定高斯均值和协方差；不是新G4光学母事件或探测器实测确认。600对合计反映本轮人为条件配比，不代表某注入率的自然分布。
''')
    figure('0004-89-1','新高斯盲评价；(a)近等强，(b)先强后弱，(c)先弱后强，(d)均匀Q分布改变；前三格各点40对，最后240对。')
    code('''display(pd.DataFrame(risk_report["aggregates"]))
display(risk_summary[(risk_summary.domain == "continuous") & (risk_summary.kind == "pair")][["method", "condition", "band", "frames", "success", "K2_outputs", "both_times_strict", "missed", "extra", "split", "merge", "q_error_p16_p50_p84", "q_error_sigma68", "not_converged"]])
''')

    md('''<a id="risk-native-20261027"></a>
## 90. 原生迁移、空／单控制与失败环节

原生760帧来自此前已公开的确认批：中／中、先强后弱、先弱后强×8/40/200 ns×40对，加空和四档单事件各80帧。已有母事件在本轮是公开迁移诊断，不算未见物理确认。相同280个原生单事件用于门槛，母事件与模型标定／评价分组不重合；跨组合和间隔有母事件复用。

预定开发要求为连续Q双事件总成功率至少比基线高5个百分点；两类输入每个空／单组额外输出≤5%，各单组成功≥95%、匹配Q偏差中位绝对值≤5%；原生双事件合计不下降。这只是开发筛选标准，达到也不自动替换默认方法。

表中强制K=2候选包括被门槛拒绝、Q<500和错时刻的输出，用于区分定位／电荷估计与事件数选择。对两槽按时间先后与真事件比较，弱／强Q误差分别保留。它们不是盲恢复结果。

配对区间按条件计算2000次bootstrap：新高斯逐对重抽，原生按实际母事件块重抽（中／中40块，强弱20块、每块两对）。区间描述当前模型／公开样本，不含完整物理不确定性；全成／全败的退化区间不能视为总体确定性。本轮未对原生合并数给独立逐帧置信区间。
''')
    figure('0004-90-1','已公开原生迁移诊断；(a)中／中，(b)先强后弱，(c)先弱后强，各点40对。')
    figure('0004-90-2','空／单控制；(a–c)新高斯，(d–f)原生，三列依次为额外输出帧、严格成功、匹配Q偏差中位，各组80帧。')
    code('''display(risk_summary[risk_summary.kind != "pair"][["domain", "method", "condition", "frames", "success", "extra_frames", "q_error_p16_p50_p84", "q_error_sigma68", "not_converged"]])
risk_diagnostic = pd.DataFrame(risk_report["diagnostics"])
display(risk_diagnostic[risk_diagnostic.condition == "high_weak"])
display(pd.DataFrame(risk_report["paired"]).query("condition == 'high_weak'"))
print("上述诊断success为强制K=2候选成功数；主表另经盲门槛选择。配对difference和paired95乘100为百分点。")
display(pd.DataFrame(risk_report["timing"]))
print("耗时为已有基线之外的增量；网络耗时不含观测特征初拟合，不是端到端实时吞吐。")
''')

    md('''<a id="risk-waveforms-20261027"></a>
## 91. 完整堆积波形与所选分量

计算前指定原生先强后弱、40/200 ns第0对，显示完整0–2048 ns记录。三列依次为总电压、较大Q输出、较小Q输出。真实原生雪崩分量和真参数条件均值分开画，仅用于展示。所有算法曲线均经过最终盲门槛；若K=1，弱输出为零，不能用被拒绝的双候选冒充恢复结果。盲输出按Q大小归列不保证与真实强／弱对应正确。
''')
    figure('0004-91-1','预定完整波形；(a–c)40 ns，(d–f)200 ns；基线、验证选中的MLP及两种物理混合，缺失的第二输出画零线。')
    code('''risk_examples = json.loads((risk_output / "examples.json").read_text())
display(pd.DataFrame([dict(间隔_ns=e["gap_ns"], 方法=m, K=len(v["fit"]["times_ns"]), 拟合Q_pe=v["fit"]["charge_pe"], 拟合时刻_ns=v["fit"]["times_ns"], 严格成功=v["score"]["success"], 额外输出=v["score"]["extra"]) for e in risk_examples for m,v in e["methods"].items()]))
''')

    lines=['| 方法 | 新高斯600对 | 已公开原生360对 |','| --- | ---: | ---: |']
    for method in ('baseline','nn_absolute','nn_logarithmic','nn_relative','mlp_nll','mlp_regularized'):
        values=[next(r for r in report['aggregates'] if (r['domain'],r['method'],r['part'])==(d,method,'all')) for d in ('continuous','native')]
        lines.append('| '+method+' | '+' | '.join(f"{r['success']}/{r['frames']}（{100*r['efficiency']:.1f}%）" for r in values)+' |')
    promising=[r['method'] for r in report['decisions'] if r['passed']]
    decision='预定开发要求通过：'+('、'.join(promising) if promising else '无')+'。默认流程保持；是否继续确认应结合分条件收益与失败成本。'
    md(f'''<a id="risk-conclusions-20261027"></a>
## 92. 实际结果、限制与下一步

{chr(10).join(lines)}

{decision}

这里比较完整输出流程。若混合方法优于直接MLP，可能来自时间精修、Q精修及候选资格改变；若MLP优于某个似然结果，也不能推出它增加了波形信息。估计损失、训练设计先验和事件数判决都会改变偏差／方差取舍。更低训练损失、K=2候选代价或整波形残差均不能替代上述整对正确标准。

下一步应先按失败环节决定：时间候选不准时改善局部时间监督及候选覆盖；Q分配仍有歧义时检查弱事件偏差与不确定区间，而不是把旧四档Q塞进先验。若新流程满足开发要求，再冻结并使用连续能量、位置变化和新的光学母事件确认。当前只有一个训练初始化、一个高斯生成模型及有限条件；不能推出普适最小间隔、最高计数率、α/γ或free-τ结论。

完整[逐组评分](artifacts/{NAME}/summary.csv)、[强制K=2诊断](artifacts/{NAME}/forced-K2-summary.csv)、[成对比较](artifacts/{NAME}/paired-comparison.json)及[结果](artifacts/{NAME}/report.json)保留。之前431个单元及图像不变，新增六幅固定编号图。
''')
    code('''display(pd.DataFrame(risk_report["decisions"]))
risk_checks = json.loads((risk_record / "science-validation.json").read_text())
print("核验通过：", risk_checks["passed"], "；评价帧：", risk_checks["evaluation_frames"], "；主评分：", risk_checks["primary_scores"], "；强制K=2：", risk_checks["forced_K2_scores"])
print("原生基线逐项一致：", risk_checks["baseline_previous_exact"], "；推断禁读真值：", risk_checks["labels_blocked_in_inference"])
print("网络重放最大差：", risk_checks["maximum_network_replay_error"], "；NLL重放最大相对差：", risk_checks["maximum_cost_replay_relative_error"])
''')
    for i,c in enumerate(cells):c['id']=f'risk-20261027-{i:02d}'
    InteractiveShell.instance();namespace={};executed=0;cwd=os.getcwd()
    try:
        os.chdir(PROJECT/'jupyter')
        for cell in cells:
            if cell.cell_type!='code':continue
            with capture_output(stdout=True,stderr=True,display=True) as capture:
                exec(compile(cell.source,'<risk-notebook>','exec'),namespace)
            assert not capture.stderr,capture.stderr
            cell.outputs=[]
            if capture.stdout:cell.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=capture.stdout))
            for out in capture.outputs:cell.outputs.append(nbformat.v4.new_output('display_data',data=out.data,metadata=out.metadata))
            executed+=1;cell.execution_count=executed
    finally:os.chdir(cwd)
    notebook.cells.extend(cells);assert notebook.cells[:431]==original
    anchors=re.findall(r'id="codex-figure-([^\"]+)"','\n'.join(c.source for c in notebook.cells if c.cell_type=='markdown'))
    assert len(anchors)==len(set(anchors))
    for c in cells:
        if c.cell_type=='markdown':assert not re.search(r'\]\([^)]*\.(?:pdf|ipynb|md)(?:[#?)][^)]*)?\)',c.source,re.I)
    images=[o['data']['image/png'] for c in cells if c.cell_type=='code' for o in c.outputs if 'image/png' in o.get('data',{})]
    expected=load(RECORD/'figure-delivery.json')['figures'];assert len(images)==len(expected)==6
    assert {hashlib.sha256(base64.b64decode(s)).hexdigest() for s in images}==set(expected.values())
    nbformat.validate(notebook);encoded=(nbformat.writes(notebook)+'\n').encode()
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=path.stem+'.',suffix='.partial',delete=False) as stream:
        stream.write(encoded);stream.flush();os.fsync(stream.fileno());temporary=Path(stream.name)
    assert temporary.read_bytes()==encoded
    nbformat.validate(nbformat.reads(encoded.decode(),as_version=4))
    assert path.read_bytes()==previous,'Notebook changed during delivery; preserve user edits'
    temporary.replace(path);assert path.read_bytes()==encoded
    save(RECORD/'notebook-delivery.json',dict(passed=True,previous_cells=431,new_cells=len(cells),executed_cells=executed,
        history_preserved=True,figures=expected,notebook_sha256=digest(path),source_sha256=digest(__file__)))
    print('RISK_NOTEBOOK_DELIVERED',len(cells),executed,flush=True)


if __name__=='__main__':deliver()
