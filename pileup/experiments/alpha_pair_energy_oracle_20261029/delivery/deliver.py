"""Append the completed energy/count/timing diagnostic to Notebook0007."""
from pathlib import Path
import json
import re
import shutil
import sys

import nbformat as nbf
from IPython.core.interactiveshell import InteractiveShell
from IPython.utils.capture import capture_output

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'pileup/analysis').is_dir())
sys.path.insert(0, str(ROOT))
from pileup.analysis.alpha_pair_energy_oracle_study import NAME, EXPERIMENT, OUTPUT, RECORD, METHODS, load, save, digest

NAMES = ['未知K的能量盲选择', '已知K，波形搜索起点', '已知K与固定真实时间', '真实时间起点后释放E/t']


def pct(x): return f'{100*x:.2f}%'
def pp(x): return f'{100*x:+.2f}'


def main():
    report = load(OUTPUT/'audited-report.json')
    assert load(RECORD/'report-validation.json')['passed']
    aggregate = {(r['method'], r['group']): r for r in report['aggregate'] if r['fold'] == 1}
    table = '\n'.join(f"| {name} | {aggregate[m,'weak']['successful_pairs']}/300（{pct(aggregate[m,'weak']['strict_pair_success'])}） | "
        f"{aggregate[m,'strong']['successful_pairs']}/600（{pct(aggregate[m,'strong']['strict_pair_success'])}） | "
        f"{aggregate[m,'weak']['not_converged']} / {aggregate[m,'strong']['not_converged']} |" for name,m in zip(NAMES,METHODS))
    diffs = report['paired_intervals']
    def delta(after,before,group='weak'):
        x=next(x for x in diffs if (x['after'],x['before'],x['group'])==(after,before,group))
        return f"{pp(x['difference'])}个百分点（条件95%区间[{pp(x['low95'])}, {pp(x['high95'])}]；增加{x['gained']}对、减少{x['lost']}对）"
    weak_fixed=aggregate['fixed_true_time','weak'];weak_free=aggregate['free_from_true_time','weak']
    incomplete=300-weak_fixed['successful_pairs']
    searchweak=[s for s in report['search_summary'] if 'weak' in s['condition']]
    searchgap=sum(s['lower_cost_from_true_start'] for s in searchweak)
    control=[]
    for method,name in zip(METHODS,NAMES):
        rr=[r for r in report['summary'] if r['method']==method and r['fold']==1 and r['kind']=='single']
        empty=[r for r in report['summary'] if r['method']==method and r['fold']==1 and r['kind']=='empty']
        control.append(f"| {name} | {sum(r['strict_10'] for r in rr)}/280 | {sum(r['missed'] for r in rr)} | {sum(r['extra'] for r in rr)} | {sum(r['output_count'] for r in empty)}/80窗 |")
    countrows=[]
    for method,name in zip(METHODS,NAMES):
        rr=[r for r in report['summary'] if r['method']==method and r['fold']==1 and r['kind']=='pair']
        countrows.append(f"| {name} | {sum(r['output_count'] for r in rr)} | {sum(r['missed'] for r in rr)} | {sum(r['extra'] for r in rr)} | {sum(r['temporal_split'] for r in rr)} | {sum(r['temporal_merge'] for r in rr)} |")
    conclusion=f'''固定真实时间后，含弱组合仍有{incomplete}/300对未达到严格E±10%端点。该对照排除了事件数和时间估计错误，但仍包含噪声、平均响应失配、有限能量域和有限优化等因素；不能直接称物理不可辨识或性能上限。

真实时间起点后释放E/t得到{weak_free['successful_pairs']}/300。与固定真实时间相比{delta('free_from_true_time','fixed_true_time')}。释放时间始终保留原解作为可行起点，代价不增；准确度与代价仍须分别评价。

300个含弱验证对中，{searchgap}对的真实时间起点能找到比波形起点搜索低超过1的加权残差平方和。这里1只是预先声明的描述性差值，不是显著性阈值；低代价也不保证能量正确。完整逐条件差值及成败转换保留。

本批结果将后续优先级指向能量估计的噪声与响应模型检查：固定真实时间组全部收敛，却仍有208/300个含弱组合失败；仅提供K只增加2个成功对，真实时间起点后释放E/t相对已知K的波形搜索也仅净增3对。它们没有显示足以解决当前弱分量问题的大幅搜索收益。这是对本批开发结果的判断，不是全局最优性或物理不可分辨的证明。
'''
    reference_text = '\n'.join(f"| {'含弱' if r['group']=='weak' else '不含弱'} | {r['isolated_both_energy10']}/{r['pairs']} | {r['shared_fixed_time_pair_energy10']}/{r['pairs']} | {r['lost_after_pairing']} | {r['gained_after_pairing']} |" for r in report['isolated_reference']['totals'])
    md,code=nbf.v4.new_markdown_cell,nbf.v4.new_code_cell
    setup=f'''from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Image
oracle_project = next(p for p in [Path.cwd(), *Path.cwd().parents]
                      if (p / "pileup/experiments/{NAME}/protocol.json").exists())
oracle_artifact = oracle_project / "jupyter/artifacts/{NAME}"
oracle_report = json.loads((oracle_artifact / "audited-report.json").read_text())
print("1540帧：1260帧留出验证，280帧开发单事件阈值校准；无新ROOT。")
print("仅blind_energy_gate是盲输出；其余三组明确使用真实K或真实时间。")
'''
    new=[md('''<a id="energy-oracle-20261029"></a>
## 28. 用已知K/时间对照定位弱事件误差（2026-10-29）

上一轮的含弱组合未见改善。本轮按相同电压、相同恢复均值、相同逐帧权重，分别提供事件数或时间信息，检查能量分配、E/t耦合、有限搜索和事件数选择。K表示窗口内的物理事件数。真值辅助结果不计作盲效率，也不是严格信息上限。

| 流程 | 已知信息 | 时间如何获得 | 用途 |
| --- | --- | --- | --- |
| 未知K的能量盲选择 | 只有观测电压及冻结模型 | 从电压搜索K=0/1/2 | 本轮盲对照 |
| 已知K，波形搜索起点 | 真实K | 直接选择已冻结的对应K候选，不读真实时间 | 分离K选择损失 |
| 已知K与固定真实时间 | 真实K、真实时间 | 时间固定，仅拟合E | 检查时间已知后的能量误差 |
| 真实时间起点后释放E/t | 真实K、真实时间初始化 | 从上一组出发，时间可在整个观测窗内变化 | 检查E/t耦合及有限搜索 |

所有组均不以真实E或真实Q作起点/界限。固定时间组由电压线性解及5种电荷份额起点拟合E；自由时间组用白化残差least_squares。每帧权重由电压变量投影的单事件Q/t固定，所有K与真值对照共用同一权重。它与上一轮原NLL起点处的权重不同，因此本轮开发数据百分比不能直接与上一轮历史38%作配对增益。

继承原能量响应、强度2.515148、29 ns代理、14400微胞、固定95/351 ns、完整2.048 μs/2 ns观测和足够动态范围，不重新训练均值、能标或网络。旧模板在名义t前0.25 ns的插值支撑继续显式保留；本轮不改变前向模型，以免混淆诊断因素。

原计划缓存2520帧。6帧先导的保守两进程成本估计4201秒超过3600秒预算后，在完整推断/评分前缩为1540帧：fold1全部900对＋280单＋80空，以及fold0的280单事件用于门槛校准。新估计2567秒；未删原始缓存或被排除的4个先导结果。fold0/fold1各280光学母ROOT分离，全部已有历史使用，无新独立输运确认。

盲候选先完成并冻结，随后才在独立阶段传入真实时间。子进程禁止读取truth/labels/ROOT文件；真值诊断阶段只通过显式参数接收K/t，不接收E/Q。280个fold0单事件用于4个拟合Q区间的经验增峰秩门，阈值0.025，随后冻结再评价fold1；由于响应也用过fold0母事件，不声称其为严格校准的误报概率。
'''),code(setup),md('''## 29. 相同验证数据上的四组结果

留出验证共900对：含弱组合300对（15格、每格20对），不含弱组合600对（15格、每格40对）。严格整对成功要求所有真实Q≥500 pe分量的E_dep误差≤10%、时间误差≤min(10 ns, Δt/4)，且没有额外输出。下表的未收敛分含弱/不含弱，只统计双事件；所有未收敛输出均保留。

| 方法 | 含弱组合严格整对 | 不含弱组合严格整对 | 未收敛（含弱/不含弱） |
| --- | ---: | ---: | ---: |
'''+table+'''

**图 0007-29-1：已知事件数和时间对严格整对能量准确度的影响。** (a)本轮未知K盲选择；(b)真实K选择已有盲候选；(c)真实K/t、只拟合E；(d)真实时间初始化后释放E/t。(b)–(d)都不是盲成绩，每格明确显示实际分母。
'''),code(f'display(Image(filename=str(oracle_artifact / "{report["figures"][0]}.png")))\n'),md('''## 30. 把事件数、时间自由度和起点搜索分开

含弱组合：已知K相对盲选择'''+delta('true_count_blind_start','blind_energy_gate')+'；固定真实时间相对已知K的波形搜索'+delta('fixed_true_time','true_count_blind_start')+'。\n\n'+conclusion+f'''
**图 0007-30-1：配对信息对照。** (a)四组在相同验证对上的严格成功率；(b)配对差值与条件95%区间。按同一ROOT关联的连通母事件组重采样，得到{report['bootstrap']['blocks']}组，按相同实验组成分{len(report['bootstrap']['stratum_sizes'])}层，2000次、种子202610297301；重复数字化和跨历史数据集的同母事件归入同组。区间条件于当前已学习模型与门槛，不包含训练/选模和器件参数不确定度。
'''),code(f'display(Image(filename=str(oracle_artifact / "{report["figures"][1]}.png")))\n'),md('''## 31. 弱分量能量、控制样本及计数错误

**图 0007-31-1：弱分量误差及保留分母的能量恢复。** (a)普通10 ns时间匹配后的E误差中位数±sigma68；sigma68=(p84−p16)/2，非FWHM或均值置信区间。(b)普通时间匹配且弱E在±10%的数量，分母为全部真实弱分量；这是分量端点，不是严格整对端点。漏检不从分母中删掉。

留出单/空控制如下。真值辅助组已被告知物理事件数，其计数控制不是新盲误报性能证据。

| 方法 | 单事件严格E恢复 | 单事件漏检 | 单事件额外输出 | 空窗输出 |
| --- | ---: | ---: | ---: | ---: |
'''+ '\n'.join(control)+'''

900对/1800个真分量的普通10 ns计数核对：

| 方法 | 输出 | 漏检 | 额外输出 | 时间误拆标志 | 时间误并标志 |
| --- | ---: | ---: | ---: | ---: | ---: |
'''+ '\n'.join(countrows)+'''

误并表示未匹配真事件靠近已分配候选，误拆表示额外候选靠近已匹配真事件；它们不是简单的计数不足/过多。E±5/10/20%、偏差、sigma68、未收敛、边界及全部条件表均保留。
'''),code(f'display(Image(filename=str(oracle_artifact / "{report["figures"][2]}.png")))\n'),md('''## 32. 完整窗口中的总波形与真实分量

**图 0007-32-1：先强后弱与先弱后强40 ns的完整波形。** 左列/右列分别是fold1对应条件首例，不按结果挑选；(a)(b)总电压及拟合，(c)(d)弱分量，(e)(f)强分量。横轴为观测窗口起点后的0–2048 ns，保留前沿和全部观测尾部。真实分量仅用于评价画图，下排用普通10 ns匹配选择拟合分量。若没有匹配分量，不画出伪造的对应曲线。
'''),code(f'display(Image(filename=str(oracle_artifact / "{report["figures"][3]}.png")))\n'),md('''## 33. 解释边界与后续分支

真值辅助组分别改变可用信息，不能把它们的成功率写成新算法的盲效率。固定真实时间仍不准，也不能直接断言“脉冲一定不能分开”：当前均值是经验代理，权重未针对E重新校准，E有能标支持域，优化搜索有限，光子/微胞噪声仍然存在。

另找回每对的两个同ROOT光学母事件的既有孤立波形，两侧都使用已知K/t、只拟合E。逐对直接检查两个孤立估计是否同时在±10%，不假设两条成功率相乘：

| 组合 | 两个孤立波形均正确 | 共享波形固定真实t后正确 | 孤立正确/共享失败 | 孤立失败/共享正确 |
| --- | ---: | ---: | ---: | ---: |
'''+reference_text+'''

这只是同母事件参考：孤立波形来自原生数字化，共享波形来自离线数字化，窗口与随机实现也不同，不能把全部差值归因于恢复历史，也不是严格性能上限。280母ROOT的E_dep身份逐一核对，额外对照没有新拟合或模拟。

对于固定真实时间后仍失败的弱事件，后续先做同模型噪声闭合或固定母事件重复数字化对照，区分统计涨落、均值失配与权重失配；再决定如何改弱光量均值/协方差。对于真实时间起点明显优于波形起点的条件，再改候选搜索；对于已知K优于盲选择的条件，再改K=1/K=2判别并控制单事件误拆。不会把所有失败统一归为优化器未收敛。

本轮所有E均受冻结能标域约束，这是模型先验，不是ADC动态范围限制。特别是K=1拟合一对高能事件时可能达到能标上界，影响K选择；因此该盲对照不能外推到能标域外单事件或任意能谱。误报控制仅覆盖本批单事件能区。

保留默认方法。本轮没有生成新确认母事件；尚未完成γ、free-τ、K≥3高率全窗、位置/参数失配或实物验证。开发控制通过后才能冻结新方法作独立确认，再推进0.1/1/10 MHz完整窗口与α/γ同等恢复。
'''),md(f'''## 34. 完成验证、来源与重放

5项数值检查通过：白化与带状二次型一致、无噪声已知时间能量反演、释放时间代价不增、真实K仅选择已冻结盲候选、空窗与经验门的并列秩。另在实际冻结模板的4种能量组合×3间隔完成12项无噪声E反演，最大相对误差约1.19×10⁻⁸。它们是数值正确性检查，不是有噪声性能保证。

本轮继承此前只读核对的Rosado2019（DOI 10.1109/JSEN.2019.2938018，PDF2–4/5–8）及Fixsen2004（DOI 10.1016/j.nima.2003.11.313，PDF2式1–2）的条件性方法依据。前者为S13360-1325CS/1350CS、LYSO γ及连续LED的非线性响应，后者为非线性TES模型的能量估计；都不是当前GAGG逐事件性能证据。两篇AI待审状态保留，本轮未新增精读、回查原PDF或修改共享库；四组对照设计与所有数值均为本项目工作。

实际源码/协议在`pileup/experiments/{NAME}/`，输出在`jupyter/artifacts/{NAME}/`，记录在`jupyter/run_records/{NAME}/`。freeze保留原2520帧准备，execution-plan及v2运行器记录成本驱动的1540帧执行；原始计划/失败成本检查不覆盖。两进程、BLAS单线程，模型与输入SHA256保留；盲候选、门槛、真值诊断依次冻结。

从项目根目录、保留本地缓存与模型的环境中重建评分/图：

```bash
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run_v2.py evaluate
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run_all_reports.py
```

Notebook图表只读已生成结果，不自动启动拟合。未提交、未推送；前一阶段6779620及其后的现有改动均保留。
'''),code('''display(pd.DataFrame(oracle_report["aggregate"]))
print(json.dumps(oracle_report["bootstrap"], ensure_ascii=False, indent=2))
''')]
    path=ROOT/'jupyter/codex_0007_SiPM_Energy_Calibration.ipynb';nb=nbf.read(path,as_version=4)
    before=nbf.read(RECORD/(path.name+'.before'),as_version=4)
    assert nb.cells==before.cells and len(nb.cells)==55
    nb.cells[0].source=re.sub(r'\*\*2026-10-29更新\*\*：.*?\n\n',
        '**2026-10-29更新**：最新已知K/时间诊断见[第28–34节](#energy-oracle-20261029)。固定真实时间、真实K的波形搜索和真实时间初始化后释放E/t分别评价；真值辅助组不计作盲效率。第21–27节保留恢复响应原型，第15–20节保留共享微胞诊断。\n\n',nb.cells[0].source,count=1,flags=re.S)
    shell=InteractiveShell.instance();count=max(c.execution_count or 0 for c in nb.cells if c.cell_type=='code')
    for cell in new:
        if cell.cell_type!='code':continue
        count+=1
        with capture_output(stdout=True,stderr=True,display=True) as out:result=shell.run_cell(cell.source,store_history=False)
        if result.error_in_exec or result.error_before_exec:raise RuntimeError(str(result.error_in_exec or result.error_before_exec))
        cell.execution_count=count;cell.outputs=[]
        for channel in ['stdout','stderr']:
            if getattr(out,channel):cell.outputs.append(nbf.v4.new_output('stream',name=channel,text=getattr(out,channel)))
        for item in out.outputs:cell.outputs.append(nbf.v4.new_output('display_data',data=item.data,metadata=item.metadata))
    nb.cells.extend(new);nbf.validate(nb);assert nb.cells[1:55]==before.cells[1:]
    for cell in new:
        if cell.cell_type=='markdown':assert not any(re.search(r'\.(pdf|md|ipynb)(#.*)?$',s,re.I) for s in re.findall(r'\]\(([^)]+)\)',cell.source))
    nbf.write(nb,path)
    planpath=ROOT/'jupyter/Research_Plan_20260921.ipynb';plan=nbf.read(planpath,as_version=4)
    oldplan=nbf.read(RECORD/(planpath.name+'.before'),as_version=4);assert plan.cells==oldplan.cells and len(plan.cells)==83
    plan.cells[0].source=plan.cells[0].source.replace('[第18节](#recovery-energy-20261029)','[第19节](#energy-oracle-20261029)')
    plan.cells.append(md(f'''<a id="energy-oracle-20261029"></a>
## 19. 事件数、时间和能量的真值辅助诊断（2026-10-29）

在独立实验`{NAME}`继承冻结恢复均值与能标，对1540帧作四组对照：280个fold0单事件标定经验门，fold1的900对＋280单＋80空用于验证。原2520帧先导成本超预算，完整推断前缩减执行范围；全部输入/先导保留，未据准确度挑样本。fold0/fold1各280母事件分离，本轮仍是开发资料，无新ROOT。

| 方法 | 含弱组合 | 不含弱组合 | 未收敛（含弱/不含弱） |
| --- | ---: | ---: | ---: |
'''+table+'\n\n'+conclusion+f'''
各组共用观测导出的逐帧权重；与上一轮历史重放的母样本、权重及K规则不同，不能直接作38%基线的配对改善。只有未知K流程为盲输出，其余均标明真实K/t帮助。ROOT关联形成110个连通组，分9种设计组成分层bootstrap，避免跨数据集重复母事件被当成独立重复。

5项数值检查和实际模板12个无噪声反演通过；编号图0007-29-1至32-1及完整波形在codex_0007第28–34节。默认方法不变，新研究未提交/上传。

后续按诊断区分：固定真实时间仍失败时先作同模型噪声闭合/重复数字化，检查均值与权重；真实时间起点找到更好解时检查有限搜索；已知K改善时检查增峰判别与单事件误报。开发控制通过后才冻结新母事件确认，再扩高率K≥3、绝对能谱与α/γ。
'''))
    nbf.validate(plan);nbf.write(plan,planpath)
    tabledir=EXPERIMENT/'tables';tabledir.mkdir(exist_ok=True)
    for p in OUTPUT.glob('*.csv'):shutil.copyfile(p,tabledir/p.name)
    save(EXPERIMENT/'result-summary.json',report)
    for name in ['numerical-validation.json','pilot-validation.json','execution-plan-validation.json','blind-validation.json',
                 'oracle-validation.json','report-validation.json','environment.json','plot-provenance.json','mother-block-design.json']:
        shutil.copyfile(RECORD/name,EXPERIMENT/name)
    save(RECORD/'notebook-validation.json',dict(passed=True,old_cells=55,new_cells=len(nb.cells),
        appended_executed_cells=6,new_figures=4,original_calculations_outputs_preserved=True,
        notebook_sha256=digest(path),plan_old_cells=83,plan_new_cells=len(plan.cells)))
    print('ORACLE_NOTEBOOK_WRITTEN',len(nb.cells),len(plan.cells),flush=True)


if __name__=='__main__':main()
