"""Record the completed prototype without changing earlier notebook outputs."""
from pathlib import Path
import json
import re
import shutil
import sys

import nbformat as nbf
from IPython.core.interactiveshell import InteractiveShell
from IPython.utils.capture import capture_output

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from pileup.analysis.alpha_pair_recovery_energy_study import NAME, EXPERIMENT, OUTPUT, RECORD, load, save, digest


def pct(value): return f'{100*value:.2f}%'
def pp(value): return f'{100*value:+.2f}'


def main():
    r = load(OUTPUT/'audited-report.json')
    assert load(RECORD/'report-validation.json')['passed']
    a = r['aggregate']; names = ['原NLL＋静态能标', '只改耦合电荷换算', '孤立E/t波形重拟合', '恢复感知E/t波形重拟合']
    methods = ['static_nll', 'coupled_charge_only', 'isolated_energy_refit', 'recovery_energy_refit']
    table = '\n'.join(f"| {name} | {pct(a[m]['weak'])} | {pct(a[m]['strong'])} | {a[m]['boundary_frames'] if m != 'static_nll' else '未复核'} | {a[m]['not_converged'] if m != 'static_nll' else '未复核'} |" for m, name in zip(methods, names))
    contrasts = [next(x for x in r['paired_comparisons'] if x['method'] == 'recovery_energy_refit' and x['baseline'] == 'static_nll' and x['group'] == g) for g in ['weak', 'strong']]
    contrast_text = '；'.join(f"{name}{pp(x['difference'])}个百分点，条件95%区间[{pp(x['low95'])}, {pp(x['high95'])}]" for name, x in zip(['含弱组合', '不含弱组合'], contrasts))
    counts = [next(x for x in r['event_totals'] if x['method'] == m and x['kind'] == 'pair') for m in methods]
    count_table = '\n'.join(f"| {name} | {x['output_count']} | {x['missed']} | {x['extra']} | {x['temporal_split']} | {x['temporal_merge']} | {x['strict_10']} | {x['frame_success_10']} |" for name, x in zip(names, counts))
    failure_table = '\n'.join(f"| {names[methods.index(x['method'])]} | {'含弱' if x['group']=='weak' else '不含弱'} | {x['wrong_count']} | {x['strict_time']} | {x['energy']} | {x['success']} |" for x in r['failures'] if x['method'] in ['static_nll', 'recovery_energy_refit'])
    baseline, current = counts[0], counts[-1]
    p = r['heldout_response']['numerical_response_checks']
    md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
    setup = f'''from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Image
recovery_project = next(p for p in [Path.cwd(), *Path.cwd().parents]
                        if (p / "pileup/experiments/{NAME}/protocol.json").exists())
recovery_artifact = recovery_project / "jupyter/artifacts/{NAME}"
recovery_report = json.loads((recovery_artifact / "audited-report.json").read_text())
print("历史重放；沿用原盲算法选择的事件数；未替换默认方法。")
'''
    def image_code(filename, first=False):
        return code((setup if first else '')+f'display(Image(filename=str(recovery_artifact / "{filename}.png")))\n')
    new = [md(r'''<a id="recovery-energy-20261029"></a>
## 21. 从拟合电荷转向含恢复历史的 E/t 联合拟合（2026-10-29）

第15–20节发现：共享微胞实际Q低于同母事件孤立Q，静态能标可能与拟合偏差相互抵消。本轮实现新的能量条件平均响应，直接拟合沉积能量E与到达时间t；不是仅改已有电荷的单位或换一个优化器。仍固定95/351 ns双发光衰减，2 ns采样、完整2.048 μs电压、足够动态范围，纯α、中央入射、局部K≤2。

冻结孤立能标得到 $q_0(E)$，沿用旧孤立电荷守恒模板 $u$ 与电子学AR3算子 $A$。对分量j：

$$z_j=A[q_0(E_j)u(t-t_j;q_0(E_j))],\quad
h_j[n]=\frac{1}{14400}\sum_{m\le n}\frac{\max(z_{1-j}[m],0)}{c_A}\exp[-(n-m)2/29],$$
$$v_j=A^{-1}\{z_j\exp(-a h_j)\},\qquad Q_j=\frac{\sum_nz_j[n]\exp(-a h_j[n])}{c_A}.$$

这里 $c_A$ 为固定增益/采样间隔乘AR系数和；h是等效电荷活动代理，不是真实微胞占据率或实际光子计数。乘法衰减改变波形形状和面积，**不把共享Q重新归一回孤立Q**。只拟合一个无量纲强度a，结果2.515148；29 ns及14400胞沿用现有配置，没有改变模拟器。全响应计算6.144 μs，拟合只用观测的2.048 μs。

Rosado2019支持“平均响应依赖光子时间分布与恢复历史”的方向；Fixsen2004支持检验能量条件均值及加权残差。这里是项目经验代理，未复现Rosado式19、Fixsen式2或微观恢复状态推断，也不宣称统计最优。两篇均复用共享库AI待审笔记；器件与验证条件见第27节。

**图 0007-21-1：先弱后强40 ns时的因果历史代理和弱分量响应。** (a)作用于先到弱分量的邻居活动；(b)已知真实E/t时，孤立均值、恢复均值与实际共享分量。此图为评价诊断，真实E/t不输入盲重拟合。横轴均为窗口起点后的0–2048 ns；旧观测窗含独立触发抖动，窗口起点不等于事件起点。
'''), image_code(r['figures'][0], True), md('''## 22. 开发划分、参数冻结与前向响应检查

原开发库280个fold0母事件用于本轮强度拟合，另280个fold1母事件用于检查；历史重放的320个母事件与两者分开。强度仅用fold0的900对、8/40/200 ns、按条件/间隔/分量等权的分数Q残差选择，随后冻结，再评价fold1的900对与280单事件。原模板/能标来自此前工作；这些资料已在项目中使用过，不能称作新的独立确认。

含弱组合的两个分量一起统计时，fold1分量Q相对RMSE从7.05%降至3.06%；不含弱组合从8.84%降至1.78%。“含弱组”并非只统计弱分量。相应归一化创新残差RMSE从0.344降至0.335、0.208降至0.168；Q均值准确并不保证逐事件波形分解准确。

280个单事件的零邻居极限严格保持原孤立响应。实际模板检查发现：线性插值使模板在名义t前0.25 ns存在极小支撑，最大提前修正约0.000129 mV。因此原“名义t之前绝对零影响”断言没有通过；按模板实际支撑检查，之前的影响为零。保留失败记录，不把离散滤波因果性称为严格连续时间物理因果性，未改模型或重跑推断来消除该结果。

**图 0007-22-1：留出开发集的分量电荷相对RMSE。** 每格包含两个分量，(a)孤立响应，(b)恢复响应；含弱组合每格20对，不含弱组合每格40对，共900对。评价使用已知E/t，仅检验前向响应。
'''), image_code(r['figures'][1]), md('''## 23. 四种方法的共同波形比较

1600帧包括1200对、320单事件和80空窗。保留原联合NLL根据观测电压选出的K；三个新流程都不能补回原来没有输出的事件。原NLL＋静态能标作基线；“只改耦合电荷换算”保持原Q/t/波形不变；两种E/t重拟合使用同一固定带状权重，差别为恢复强度a=0或2.515148。

权重固定在原盲Q/t处，关闭原条件电荷总量的κ项（κ=0），因为本轮均值条件是E。它是GLS加权消融，不是重新校准的完整似然或置信区间。时间限于原起点±24 ns，L-BFGS-B最多60次迭代，能量限于冻结能标支持域；边界及不收敛结果全部保留。

严格整对成功要求：两个真实事件Q≥500 pe均匹配，时间误差≤min(10 ns, Δt/4)，各沉积能量误差≤10%，没有额外输出。含弱/不含弱分别有600对、15个条件/间隔，各格40对；下表为15格等权均值。边界及不收敛数覆盖全部1600帧。

| 方法 | 含弱组合 | 不含弱组合 | 边界帧 | 未收敛帧 |
| --- | ---: | ---: | ---: | ---: |
'''+table+'\n\n恢复重拟合相对原NLL：'+contrast_text+'''。区间为后验描述性配对bootstrap：40母事件块分两层各20、2000次，所有方法/条件共用抽样；不含能标、强度学习、选模及物理参数不确定度。

**图 0007-23-1：四方法的严格整对能量恢复数。** 每格分母40，颜色统一；(a)原NLL＋静态能标，(b)只改电荷换算，(c)孤立E/t重拟合，(d)恢复E/t重拟合。不要把前向Q预测RMSE改善当作本图恢复率改善。
'''), image_code(r['figures'][2]), md('''## 24. 孤立/空窗控制、漏检与额外输出

四方法的弱单事件严格能量恢复均为65/80，其余低/中/强均80/80；80空窗均0输出。中能单事件有1个额外输出，沿用原盲K，未删除。单事件匹配后误差分布另报中位数与sigma68=(p84−p16)/2，不把sigma68当FWHM或均值置信区间。

**图 0007-24-1：历史单事件控制。** (a)匹配后能量误差中位数±sigma68；(b)严格能量事件恢复率，分母始终为80个真事件。

双脉冲汇总（1200对、2400真分量）：

| 方法 | 输出 | 普通10ns漏检 | 额外输出 | 时间误拆标志 | 时间误并标志 | 严格E事件 | 严格整对 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
'''+count_table+'''

时间误拆/误并沿用原10 ns邻近图定义：未匹配真事件靠近已分配候选标为误并，额外候选靠近已匹配真事件标为误拆。初版评分中的split/merge只是计数多/少标志，已在独立审核表重命名，保留原文件；两定义不可混用。静态基线逐帧复现旧评分。输出总数不变，时间重拟合可改变匹配、漏检、额外输出和时间图标志。
'''), image_code(r['figures'][3]), md('''## 25. 查看完整堆积波形及恢复分量

**图 0007-25-1：完整0–2048 ns窗口中的观测电压、总拟合和分量。** (a)(c)high→high、8 ns；(b)(d)weak→high、40 ns。均取协议先导条件的首例，未按拟合好坏挑图。上排总波形；下排先到分量，与真实模拟分量作评价对照。下排采用普通10 ns匹配，不能据曲线接近就认定满足更严格的8 ns整对端点。强脉冲下小分量的分配仍可能有偏差。
'''), image_code(r['figures'][4]), md('''## 26. 还没有解决什么；下一步的准入条件

恢复响应明显改善不含弱事件的能量恢复，但含弱组没有净改善；因此保留为开发原型，不替换默认方法、不宣称提高事件检出数、不据此启动大规模高率生产。

按顺序把失败归为：K不符；K正确但严格时间不全；K与严格时间都正确但能量不全；全部成功。每个组合组均600对，这是一种互斥记账，不能作为失败的唯一物理成因。

| 方法 | 组合 | K不符 | 严格时间不足 | 能量不足 | 成功 |
| --- | --- | ---: | ---: | ---: | ---: |
'''+failure_table+'''

下一步先在开发资料中分开检验两个瓶颈：原K已经漏掉的弱事件，需要能量条件模型参与K=1/K=2选择并控制空/单误报；已经输出两事件却E偏差大的组合，需要检查低光量条件的波形均值/协方差及E/t相关性。先补实际模板起点约定和数值收敛控制，不能仅由电荷均值拟合好就认定响应充分。

进一步方法应另开冻结实验，用开发集确定候选与门槛，确认样本不再调参；待弱组和控制同时通过后，才用新母事件确认，并扩展0.1/1/10 MHz完整窗口未知K≥3与绝对计数谱。α/γ同等恢复、free-τ、位置/器件参数失配及实物校准仍未完成。
'''), md(f'''## 27. 依据、复现与完成验证

文献只读复用共享库：Rosado2019，DOI 10.1109/JSEN.2019.2938018，PDF2–4恢复/平均响应、PDF5–8验证；S13360-1325CS/1350CS，LYSO约300–2100 keV γ、400 ns积分及连续绿色LED，采样率/样本数未完整报告。适中非线性近似及连续光平台约20%差异保留，29 ns仍是同系列代理常数，并非当前6050PE实测。其上升/下降双指数差不等于项目95/351 ns两发光衰减。

Fixsen2004，DOI 10.1016/j.nima.2003.11.313，PDF2式1–2；模型TES约0.5–5 eV，每能点5000噪声抽样，非GAGG高率实测或堆积验证，采样率/完整窗口点数未列全。式1二次型不含一般似然的logdet，本轮只借鉴能量条件均值/权重的研究方向。两篇保留AI待审，没有人工采用确认；未回查原PDF或写入共享文献库。

实际研究路径`pileup/experiments/{NAME}/`；本地结果`jupyter/artifacts/{NAME}/`；运行记录`jupyter/run_records/{NAME}/`。模型/数据/源码和逐帧预测均有冻结哈希，模型冻结先于历史评价。推断子进程只接收电压、原盲起点和固定模型，truth/labels/ROOT三类读取探针均拒绝。7项数值检查、旧均值/权重兼容性、1180个留出响应检查、6400条评分和原基线复现已完成；原严格名义起点因果性检查的失败单独保留。

两个进程完成本批推断；6帧先导后处理其余1594帧。没有新G4/ROOT或神经网络训练，没有改变物理模型。35个恢复重拟合未收敛结果保留，不能当作全局最优。只看本笔记图表不触发训练或模拟。

从项目根目录可重建评分与图表；输入缓存、冻结模板和Python环境须可用：

```bash
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py evaluate
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run_report_v3.py
```

首次报告检查因旧模板的0.25 ns插值支撑触发断言；原报告源码及冻结记录保留，v2显式保留名义起点检查不通过，v3另将未保留的原基线优化器状态标为未复核，而非零失败。拟合输出完全未变。

本轮开始前已本地提交前一阶段：6779620，英文记录“Add shared-cell response diagnostics for SiPM energy reconstruction”。本轮新研究尚未提交或上传。
'''), code('''display(pd.read_csv(recovery_artifact / "failure-stages.csv"))
print(json.dumps(recovery_report["heldout_response"], ensure_ascii=False, indent=2))
''')]
    nbpath = ROOT/'jupyter/codex_0007_SiPM_Energy_Calibration.ipynb'
    before = nbf.read(RECORD/(nbpath.name+'.before'), as_version=4)
    nb = nbf.read(nbpath, as_version=4)
    assert nb.cells == before.cells and len(nb.cells) == 42
    old_intro = nb.cells[0].source
    nb.cells[0].source = re.sub(r'\*\*2026-10-29更新\*\*：.*?\n\n',
        '**2026-10-29更新**：恢复历史进入E/t联合拟合，最新结果见[第21–27节](#recovery-energy-20261029)。不含弱事件的严格整对能量恢复从53.33%提高到82.33%，含弱组合38.17%→38.00%，尚未改善。保留开发原型，默认算法未替换；第15–20节为响应诊断，第8–14节为能标确认。\n\n', old_intro, count=1, flags=re.S)
    shell = InteractiveShell.instance(); count = max(c.execution_count or 0 for c in nb.cells if c.cell_type == 'code')
    for cell in new:
        if cell.cell_type != 'code': continue
        count += 1
        with capture_output(stdout=True, stderr=True, display=True) as out: result = shell.run_cell(cell.source, store_history=False)
        if result.error_in_exec or result.error_before_exec: raise RuntimeError(str(result.error_in_exec or result.error_before_exec))
        cell.execution_count = count; cell.outputs = []
        for channel in ['stdout', 'stderr']:
            if getattr(out, channel): cell.outputs.append(nbf.v4.new_output('stream', name=channel, text=getattr(out, channel)))
        for item in out.outputs: cell.outputs.append(nbf.v4.new_output('display_data', data=item.data, metadata=item.metadata))
    nb.cells.extend(new); nbf.validate(nb)
    assert nb.cells[1:42] == before.cells[1:]
    for cell in new:
        if cell.cell_type == 'markdown':
            assert not any(re.search(r'\.(pdf|md|ipynb)(#.*)?$', s, re.I) for s in re.findall(r'\]\(([^)]+)\)', cell.source))
    nbf.write(nb, nbpath)
    planpath = ROOT/'jupyter/Research_Plan_20260921.ipynb'; plan = nbf.read(planpath, as_version=4)
    oldplan = nbf.read(RECORD/(planpath.name+'.before'), as_version=4)
    assert plan.cells == oldplan.cells and len(plan.cells) == 82
    plan.cells[0].source = plan.cells[0].source.replace('[第17节](#sipm-energy-history-20261029)', '[第18节](#recovery-energy-20261029)')
    plan.cells.append(md(f'''<a id="recovery-energy-20261029"></a>
## 18. 恢复感知 E/t 联合拟合原型（2026-10-29）

前一阶段共享微胞诊断已本地提交6779620，未上传。本轮把恢复历史引入能量条件波形，单强度a=2.515148由900开发对选择后冻结，另900对与280单事件检查。含弱组合两个分量的Q相对RMSE7.05%→3.06%，不含弱组合8.84%→1.78%。

1600帧历史重放比较原NLL＋静态能标、只改电荷换算、孤立E/t重拟合、恢复E/t重拟合，沿用原盲K。严格整对E±10%成功：不含弱组53.33%→82.33%，含弱组38.17%→38.00%；条件差异不能被总体改善掩盖。四方法单事件控制均为弱65/80、其余各80/80，80空窗无输出，中能单事件沿用1个额外输出。恢复重拟合35帧未收敛，保留全部结果。

真实模板审核发现0.25 ns起点前插值支撑，最大提前修正0.000129 mV；不作严格连续时间因果性声明。原名义t断言失败及报告修订都保留，没有事后修改模型。7项数值检查、原基线复现、预测冻结和真值读取隔离完成。无新ROOT，全部母事件已有历史使用，不能称新独立确认。

本地codex_0007追加第21–27节及五幅固定图0007-21-1至25-1；实验`pileup/experiments/{NAME}/`，采用记录L75。本轮新研究未提交/上传，默认算法不变。

下一步在开发集分别处理K漏检与已检出弱事件E/t偏差：统一模板实际起点、检查收敛；用能量条件响应参与K=1/K=2比较并控制空/单误报，同时检查弱光量波形均值/协方差。弱组与控制通过后再冻结新母事件确认，然后推进0.1/1/10 MHz全窗未知K≥3、绝对计数谱、α/γ同等恢复与模型失配。
'''))
    nbf.validate(plan); nbf.write(plan, planpath)
    tabledir = EXPERIMENT/'tables'; tabledir.mkdir(exist_ok=True)
    for path in OUTPUT.glob('*.csv'): shutil.copyfile(path, tabledir/path.name)
    save(EXPERIMENT/'result-summary.json', r)
    for name in ['report-validation.json', 'heldout-response-validation.json', 'initial-response-audit-failure.json',
                 'pilot-validation.json', 'inference-validation.json', 'compatibility.json', 'numerical-tests.json',
                 'environment.json', 'plot-inputs.json']:
        shutil.copyfile(RECORD/name, EXPERIMENT/name)
    save(RECORD/'notebook-validation.json', dict(passed=True, old_cells=42, new_cells=len(nb.cells),
        old_calculations_outputs_preserved=True, appended_executed_cells=6, figures=5,
        notebook_sha256=digest(nbpath), plan_old_cells=82, plan_new_cells=len(plan.cells)))
    print('RECOVERY_DELIVERED', len(nb.cells), len(plan.cells), flush=True)


if __name__ == '__main__': main()
