"""Write the experiment README and append the project adoption record."""
from pathlib import Path
import json

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'pileup/analysis').is_dir())
NAME = 'alpha_pair_energy_oracle_20261029'
EXP, OUT, REC = [ROOT/p/NAME for p in ['pileup/experiments', 'jupyter/artifacts', 'jupyter/run_records']]


def load(p): return json.loads(p.read_text())
def pct(x): return f'{100*x:.2f}%'
def signed(x): return f'{100*x:+.2f}'


def main():
    report = load(OUT/'audited-report.json')
    methods = report['methods']; names = ['未知K能量盲选择', '已知K，波形搜索起点', '已知K与固定真实时间', '真实时间起点后释放E/t']
    aggregate = {(r['method'], r['group']): r for r in report['aggregate']}
    table = '\n'.join(f"| {name} | {aggregate[m,'weak']['successful_pairs']}/300（{pct(aggregate[m,'weak']['strict_pair_success'])}） | {aggregate[m,'strong']['successful_pairs']}/600（{pct(aggregate[m,'strong']['strict_pair_success'])}） | {aggregate[m,'weak']['not_converged']} / {aggregate[m,'strong']['not_converged']} |" for m,name in zip(methods,names))
    fixed = aggregate['fixed_true_time','weak']; free = aggregate['free_from_true_time','weak']
    effects = []
    for a,b,label in [('true_count_blind_start','blind_energy_gate','已知K减盲选择'),('fixed_true_time','true_count_blind_start','固定真实t减已知K波形搜索'),('free_from_true_time','fixed_true_time','释放t减固定真实t')]:
        r=next(x for x in report['paired_intervals'] if (x['after'],x['before'],x['group'])==(a,b,'weak'))
        effects.append(f"{label}：{signed(r['difference'])}个百分点，条件95%区间[{signed(r['low95'])},{signed(r['high95'])}]（增加{r['gained']}、减少{r['lost']}对）")
    reference='；'.join(f"{'含弱' if x['group']=='weak' else '不含弱'}组，两个同母事件孤立估计均合格{x['isolated_both_energy10']}/{x['pairs']}，共享波形固定真实时间后合格{x['shared_fixed_time_pair_energy10']}/{x['pairs']}" for x in report['isolated_reference']['totals'])
    blind=load(REC/'blind-validation.json');oracle=load(REC/'oracle-validation.json')
    controls=[]
    for method,name in zip(methods,names):
        ss=[r for r in report['summary'] if r['method']==method and r['fold']==1 and r['kind']=='single']
        ee=[r for r in report['summary'] if r['method']==method and r['fold']==1 and r['kind']=='empty']
        controls.append(f"| {name} | {sum(x['strict_10'] for x in ss)}/280 | {sum(x['missed'] for x in ss)} | {sum(x['extra'] for x in ss)} | {sum(x['output_count'] for x in ee)}/80窗 |")
    core=f'''已完成1540帧的共同对照：280个fold0单事件校准经验增峰门，fold1全部900对、280单事件、80空窗验证；fold0/fold1各280个光学母ROOT分离，均为既有开发资料，无新独立输运确认。

| 方法 | 含弱组合严格整对 | 不含弱组合严格整对 | 未收敛（含弱/不含弱双事件） |
| --- | ---: | ---: | ---: |
{table}

严格端点为两个真实Q≥500pe事件的E_dep±10%、时间误差≤min(10ns,间隔/4)、无额外输出。验证双脉冲的1800个真实分量全部超过500pe；15个含弱条件/间隔各20对，15个不含弱条件/间隔各40对，各组等权。只有第一行为盲输出，其余都明确使用真实K或真实时间，不是盲效率或严格上限。

{'；'.join(effects)}。按ROOT连通的110个母事件组、9种相同实验组成分层，2000次配对bootstrap，种子202610297301；同一母事件的跨数据集重复数字化归入同组。区间条件于已学习的模型与门槛，不含参数学习/选择和器件不确定度。

固定真实时间后，含弱组合仍有{300-fixed['successful_pairs']}/300对不满足能量端点；该部分不能归为事件数或时间估计错误，但仍包含噪声、均值/权重失配、能量界限和有限优化等因素。释放E/t得到{free['successful_pairs']}/300，代价不增并不保证能量更准确。

本批固定真实时间拟合全部收敛；已知K仅比盲选择增加2个含弱成功对，真实时间初始化后释放E/t比已知K的波形搜索净增3对。因此下一步优先区分噪声与响应/权重失配，候选搜索不是本批结果支持的主要改进方向；不据此声称搜索已达到全局最优。新盲选择在280个留出单事件中产生9个额外输出，亦不能视为已获得可靠误报控制。

同ROOT孤立参考：{reference}。两侧都使用已知K/t、只拟合E，直接逐对检查两份孤立结果，不假设成功概率相乘。孤立来自原生、共享来自离线数字化，窗口和随机实现不同，不能把全部差值归因于恢复或当作严格上限；母事件E_dep一致性核对通过。

| 方法 | 留出单事件严格E恢复 | 单事件漏检 | 单事件额外输出 | 空窗输出 |
| --- | ---: | ---: | ---: | ---: |
{chr(10).join(controls)}

真值辅助组的计数控制受已知K帮助，不作为新盲误报证据。全部条件的E偏差、sigma68=(p84−p16)/2、E±5/10/20%、漏检、额外输出、时间误拆/误并、边界与不收敛保留于tables。sigma68非FWHM；匹配子集误差不能替代全体真分量效率。
'''
    method='''继承上一轮冻结的恢复均值、孤立能标、a=2.515148、29ns代理和14400微胞，不改物理模型、不重新训练网络或能标。范围为中央纯α、固定95/351ns双衰减、2ns×1024完整电压、足够动态范围、局部K≤2。源能量不等于实际E_dep，评价按母事件实际沉积能量。

四组共用由观测电压变量投影单事件Q/t得到的固定带状权重。白化平方和等于rᵀC⁻¹r，所有K与真值对照用同一C；是固定GLS指标，非新完整似然或区间标定。权重、母样本和K选择规则与上一轮历史NLL重放不同，本轮百分比不能直接与此前38%作配对增益。

盲候选适配现有PairFitter起点逻辑，另加入8种间隔、两个方向和5种份额的恢复均值初值评分；保留原二事件起点与排名前4个，24次评估筛选，最优2个最多120次评估精化。自由时间界限为整个观测窗，E限冻结能标域；真实E/Q从不作为初值或界限。真实K组只在已冻结候选中选择对应K，不回填真实时间起点的更好解。固定真实时间组用电压线性初值＋5种份额，仅拟合E；释放时间组从此出发，是额外真值初始化诊断。

盲K选择保留峰高>1mV触发和预测Q≥500pe，使用共同权重下的正代价改进及单事件加一秩≤0.025，Q分箱1500/5000/12500pe。门只用280个fold0单事件（各箱40/80/80/80）校准并冻结，再连接fold1评价。响应亦用过fold0母事件，因此此经验门不具有本文已证明的误报率保证。

E支持域也是模型先验：K=1拟合一对高能事件可能达到能标上界，从而影响K选择；本轮没有验证能标域外单事件，不能推广到任意能谱。旧模板在名义t前0.25ns的插值支撑保留，不宣称严格连续时间物理因果性。没有ADC触顶整窗拒绝。
'''
    validation=f'''5项数值测试通过：白化二次型一致、无噪声固定t反演、释放t代价不增、已知K不替换盲参数、空窗与经验门并列秩。实际模板4能量组合×3间隔的12个无噪声反演最大相对E误差1.19×10⁻⁸；只是数值正确性，不是有噪声性能保证。

原2520帧计划的6帧先导保守成本4201秒超过3600秒预算，完整推断/评分前缩为1540帧，新预测2567秒；原输入、失败成本检查、源码及4个未纳入完整评价的先导结果保留。execution-plan、execution-freeze和v2运行器记录范围变化，模型/搜索配置未据准确度调整。盲候选阶段处理{blind['new_frames']}个新帧，用时{blind['seconds']:.2f}秒；真值诊断阶段处理{oracle['new_frames']}个新帧，用时{oracle['seconds']:.2f}秒；均最多2进程、BLAS单线程。

盲候选先冻结，再冻结经验门，最后运行并冻结真值诊断。子进程truth/labels/ROOT读取探针拒绝；诊断阶段经显式参数仅提供K/t，真实E/Q仅用于最后评分。原模型、输入、代码与逐帧输出哈希保存，无新ROOT或默认方法替换。
'''
    follow='''后续先对“固定真实时间仍失败”的弱事件做同模型噪声闭合/固定母事件重复数字化，分开噪声、均值与权重失配；对真实时间起点更优的条件检查有限搜索，对已知K改善的条件检查增峰规则及单事件误拆。开发控制通过后再冻结新母事件确认，之后扩0.1/1/10MHz完整窗口K≥3、绝对能谱、α/γ同等恢复、free-τ、位置/器件失配与实物校准。本轮不将真值辅助成绩用于默认算法升级。'''
    readme=f'''# 弱α能量误差：已知事件数与时间的对照

{core}
## 方法与信息边界

{method}
## 运行与验证

{validation}
源码和配置在source及各freeze清单，结果与表在result-summary.json和tables。原run.py保留最初完整计划；实际执行使用run_v2.py。保留本地波形、模型、Python环境后，从项目根目录运行：

```bash
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run_v2.py evaluate
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run_all_reports.py
```

第二条顺序生成主报告与同母孤立参考，均不重新拟合。完整输入缓存2520帧，实际评分只使用execution-plan中的1540帧；不能将残留先导文件数当正式样本数。数据、模型和Notebook本地维护，不随Git上传。

本地codex_0007追加第28–34节及图0007-29-1至32-1，保留原55单元计算/输出；Research_Plan追加第19节。文献采用条件承接docs/literature-use.md的L75，新增诊断记录L76。此前已有未提交改动保留，本轮未提交/上传。

## 下一步

{follow}
'''
    (EXP/'README.md').write_text(readme)
    doc=ROOT/'docs/literature-use.md';text=doc.read_text();assert '## L76 ' not in text
    evidence='''**文献依据层级**：承接L75，复用此前已核对的[Rosado2019笔记](../../papers/sources/2019_Rosado_IEEESensorsJ_Modeling-the-Nonlinear-Response-of-Silicon-Photomultipliers.md)（[DOI](https://doi.org/10.1109/JSEN.2019.2938018)，PDF2–4、5–8）及[Fixsen2004笔记](../../papers/sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md)（[DOI](https://doi.org/10.1016/j.nima.2003.11.313)，PDF2式1–2）。前者为S13360-1325CS/1350CS、LYSO约300–2100keV γ的400ns积分及连续LED，采样率/样本数未完整报告，近似适用域及1350连续平台约20%偏差保留；后者为非线性TES模型、约0.5–5eV、每能点5000次噪声抽样，非GAGG高率实测，完整采样率/点数未列全。

两篇均AI待审；未回查原PDF、未新增人工采用或修改共享库。本轮复核主题导航与既有笔记条件，没有新增论文性能主张。真实K/t对照、候选搜索和所有结果均为项目方案，不是Rosado式19、Fixsen式2或其最优性结论的复现。
'''
    text+='\n\n<a id="energy-oracle-20261029"></a>\n## L76 — 2026-10-29：已知K与时间的弱α能量误差诊断\n\n'+evidence+'\n**本项目计算**：\n\n'+core+'\n**方法和适用范围**：\n\n'+method+'\n**验证与复现**：\n\n'+validation+f'\n实验为`pileup/experiments/{NAME}/`，本地结果/记录为同名jupyter/artifacts、run_records目录。codex_0007新增第28–34节、图0007-29-1至32-1，Research_Plan新增第19节；既有计算保留，无新ROOT、未提交/上传。\n\n'+follow+'\n'
    doc.write_text(text)
    path=ROOT/'README.md';text=path.read_text();assert f'({NAME})' not in text and f'pileup/experiments/{NAME}/README.md' not in text
    values=[aggregate[m,'weak']['strict_pair_success'] for m in methods]
    text+=f'''\n随后完成[弱α的已知K/时间对照](pileup/experiments/{NAME}/README.md)：1540帧开发数据，包括900对留出双脉冲、280留出单事件、80空窗及280门槛校准单事件。含弱组合的严格整对E±10%成功依次为未知K盲选择{pct(values[0])}、已知K波形搜索{pct(values[1])}、固定真实时间{pct(values[2])}、真实时间起点后释放E/t{pct(values[3])}。后三项均有真值帮助；本轮样本/权重不同，不能与上一轮历史38%直接作配对增益。完整能量偏差、误拆漏检、同母孤立参考及完整波形见本地codex_0007第28–34节和[采用记录L76](docs/literature-use.md#energy-oracle-20261029)。原模型与默认方法保持，无新ROOT；下一步按诊断分开检查噪声/均值/权重与有限搜索。\n'''
    path.write_text(text)
    print('ORACLE_DOCS_WRITTEN', flush=True)


if __name__=='__main__':main()
