"""Write checked result summaries for the allocation experiment and L60."""
from __future__ import annotations

from .alpha_pair_allocation_study import PROJECT,EXPERIMENT,OUTPUT,RECORD
from .alpha_pair_information_study import load,save,digest
from .alpha_pair_charge_allocation import VARIANTS


def deliver():
    checks=load(RECORD/'science-validation.json');notebook=load(RECORD/'notebook-delivery.json')
    assert checks['passed'] and notebook['passed'] and load(RECORD/'figure-visual-review.json')['passed']
    report=load(OUTPUT/'report.json');summary=report['summary']
    lookup={(r['domain'],r['method'],r['condition'],r['gap_ns']):r for r in summary}
    def get(domain,method,condition,gap):return lookup[domain,method,condition,gap]
    lines=['| 输入／方法 | 先强后弱8 ns | 40 ns | 200 ns | 先弱后强8 ns |',
           '| --- | ---: | ---: | ---: | ---: |']
    for domain in ('native','matched'):
        for method in ('baseline',*VARIANTS):
            counts=[get(domain,method,'high_weak',g)['success'] for g in (8.,40.,200.)]
            counts.append(get(domain,method,'weak_high',8.)['success'])
            lines.append('| '+domain+' / '+method+' | '+' | '.join(f'{n}/40' for n in counts)+' |')
    results='\n'.join(lines)
    controls=[]
    for method in ('baseline',*VARIANTS):
        rr=[get('native',method,c,None) for c in ('all','weak','low','mid','high')]
        controls.append(f"- `{method}`：原生空／weak／low／mid／high额外输出帧依次为"+'/'.join(str(r['extra_frames']) for r in rr)+'（每组80帧）。')
    chosen=report['chosen']
    decision=f'按预定原生开发约束暂选`{chosen}`，后续仍需连续Q／失配和新母样本确认，不替换默认流程。' if chosen else '四种变体均未同时通过预定原生开发约束，不选择新默认流程。'
    baseline40=get('native','baseline','high_weak',40.)['success'];new40=get('native','allocation_profile','high_weak',40.)['success']
    baseline200=get('native','baseline','high_weak',200.)['success'];new200=get('native','allocation_profile','high_weak',200.)['success']
    diag=next(r for r in report['diagnostics'] if (r['domain'],r['method'],r['condition'],r['gap_ns'])==('native','allocation_profile','high_weak',40.))
    timing=[r for r in report['timing'] if r['domain']=='native']
    timing_text='；'.join(f"{r['method']} {r['incremental_seconds_mean']:.3f} s、{r['objective_evaluations_mean']:.0f}次函数调用" for r in timing)
    runtime=load(RECORD/'inference.json');pilot=load(RECORD/'pilot-inference.json')
    status_lines=['| 变体 | 开发要求全部通过 | 未通过项 |','| --- | --- | --- |']
    for row in report['decisions']:
        status_lines.append('| '+row['method']+' | '+str(row['passed'])+' | '+('、'.join(k for k,v in row['checks'].items() if not v) or '无')+' |')
    statuses='\n'.join(status_lines)

    readme=f'''# 总电荷与强弱分配的盲搜索对照

2026-10-27执行，承接局部信息诊断。本轮改变候选搜索和数值坐标，检验强弱Q反相关方向的有限优化；均值、协方差、logdet及物理模型保持。{decision}

## 数据与预定方法

280个独立原生单事件用于门槛，760帧已公开原生数据与760帧同Q/t高斯draw0用于开发比较，共1800帧；高斯draw1本轮未运行。每种输入包含三种组合×8/40/200 ns×40对，加空／四档单事件各80帧。训练280、门槛280和评价320个母事件互不重叠；跨条件复用母事件。没有新ROOT或网络训练。

保持纯α、固定95/351 ns双τ、2 ns×1024完整模拟电压、充分动态范围、未知局部K≤2。原生Q为实现雪崩等效pe，高斯Q为相同数值的条件参数，不作沉积能量；t为有效面首光子。

四种变体为`component_local`、`allocation_local`、`component_profile`、`allocation_profile`，均保留基线K0/K1和原最佳K2。local对原K2候选增加最多70次迭代；profile另外使用同一组最多8个起点，各筛选20次，最好的两个精修70次，并复用local结果。

新坐标为`t1/10,t2/10,100*ln(S),10*r`，S=Q1+Q2，Q1=L+(U−L)r，Q2=S−Q1；L=max(50,S−50000)，U=min(50000,S−50)，r∈[0,1]。完整保持每个Q在50–50000 pe的原域。时间范围仍为起点±24 ns，L-BFGS-B设置相同。数值缩放和有限差分的物理步长也改变；不声称增加信息，不给MLE目标加入Jacobian。

共同profile以观测单事件时刻为锚，前／后4/8/12/20/40/80/160/320 ns配第二候选；第一事件份额0.01/0.03/0.05/0.10/0.15/0.5/0.85/0.90/0.95/0.97/0.99。每点在可行域内数值优化lnS，最多24次、容差10⁻⁴，并保留原值。每方向选4个间隔的最佳份额作共同起点，不读取真实时刻、Q或K。有限剖面与局部优化均无全局最优保证。

各方法在相同280个原生单事件上重新标定K2增峰分布；原单Q分档1500/5000/12500 pe、加一／并列上尾秩p≤0.025和Q≥500保护保持。另存直接沿用原数值门槛的敏感性表。推断worker阻止打开truth、case、ROOT；所有失败和未收敛候选保留。

## 实际结果

严格整对成功要求两真实Q≥500、各Q误差≤20%、各时刻误差≤min(10 ns,间隔/4)、无额外输出。每点40对：

{results}

原生40 ns、allocation_profile的K2候选中，{diag['positive_cost_gain']}/40代价下降，但强制K2严格成功{diag['known_K_success']}/40，盲结果{new40}/40。不能以似然下降冒充逐事件恢复。

{chr(10).join(controls)}

全部条件的Q偏差、sigma68误差宽度、时间误差、漏检、额外输出、误拆／误并、K2输出和时间匹配，以及原门槛对照均保存在CSV。sigma68为匹配后Q相对误差的(p84−p16)/2，未标定沉积能量，不能当能谱FWHM。

预定原生约束包括空／单额外输出≤5%、单事件严格成功≥95%、单Q偏差中位绝对值≤5%；先强后弱40 ns至少提高10个百分点、200 ns不下降；先弱后强8 ns最多下降5个百分点，中／中各点≥95%，六个强弱组平均不下降。

{statuses}

{decision}

差值区间用40个母事件块、前20／后20分层、2000次共同权重bootstrap。中／中每组涉及40块，强弱各涉及前20块、每块2对，不能把全部帧独立重抽。它是已公开开发样本的条件性描述；退化[0,0]不表示总体差异确定，高斯draw0不是新光学母事件。

## 成本、验证和边界

原生每帧平均增量：{timing_text}。profile单独使用时包含共同剖面成本；本次共同工作只计算一次。数字排除已缓存基线，逐帧并发wall time不是实时吞吐保证，迭代上限相同不等于函数调用数相同。

7项测试通过：完整电荷域往返、总Q守恒、目标不变、可辨认相关电荷玩具问题、剖面保留原值、空候选保留、真值访问阻断。32帧先导后4进程运行，正式新增{runtime['new_frames']}帧用时{runtime['seconds']:.1f}s，先导{pilot['seconds']:.1f}s。实际环境与源码哈希保留；源HEAD为`{load(RECORD/'environment.json')['git_head']}`加冻结新增代码，上一轮信息分析改动完整保留。

1800帧／7200组变体候选，7600条主评分／13680条含原门槛评分，252行汇总，3600条强制K2诊断核验通过；基线逐条与原缓存一致，K0/K1不变，原K2代价不被劣化，重放最大相对差{checks['maximum_cost_replay_error']:.3g}。标量剖面{checks['scalar_profile_points']}点，其中{checks['scalar_profile_not_converged']}点优化器未报成功，仍保留最佳已算值；完整局部拟合未收敛统计另列。

不读取真实Q/t的四种方法都属于盲推断；强制K2图和实际分量仅用于诊断／显示。图85-1为预定第0对、40和200 ns、全长0–2048 ns，分别画总波形、强分量和弱分量；若最终只选一个，弱输出画零，绝不拿被拒绝K2候选冒充成功解堆积。

既有Fixsen2004支持条件均值／协方差原则，Liddick2012支持检测／估计分别评价；据AI精读稿待人工审核，本轮未重查PDF。坐标和剖面为项目实现，非论文已验证的GAGG改进；详见`docs/literature-use.md` L60。

局部K≤2、固定双τ、已公开三种组合和三档间隔限制结论。没有验证连续Q／位置失配、混合α/γ、free-τ、随机高率K≥3或实物性能。若同一似然搜索收益不足，后续应检验估计目标、响应／噪声模型和来自训练集的先验；禁止用四档真实能量约束测试答案。

## 文件与重放

本地Notebook为`jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb`第82–86节；旧410个单元及输出保留，新增{notebook['new_cells']}单元，{notebook['executed_cells']}代码单元已执行。固定图0004-82-1、83-1、83-2、84-1、85-1。

`jupyter/artifacts/alpha_pair_charge_allocation_20261027/`保存全部候选、剖面、评分、CSV和图；`jupyter/run_records/alpha_pair_charge_allocation_20261027/`保存输入来源、环境、验证与备份。`freeze.json`保存28个计算文件；报告源码与运行入口另由`report-freeze.json`保护。方法门槛是实际280校准帧导出的`gate-calibrations.json`。

从项目根目录运行，先保留全部本地输入／模型和运行记录：

```bash
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_charge_allocation_20261027/run.py infer --workers 4
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_charge_allocation_20261027/report.py analyze
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_charge_allocation_20261027/report.py figures
bash tools/run_gagg_python.sh pileup/experiments/alpha_pair_charge_allocation_20261027/report.py notebook
```

计算按index续跑，正式运行要求先导验证；不要在已冻结目录再次prepare。报告用冻结源码，无需切换Git。输入在`input-provenance.json`中逐一记录，图中真实分量另见`plot-truth-provenance.json`；仅克隆仓库不能替代本地模型NPZ、数据数组和运行记录。Jupyter／Notebook／模型数组按项目约定不上传。Notebook重放保留后续用户修改；文档生成另检查已交付文件哈希，避免覆盖新内容。
'''
    root_marker='<a id="alpha-pair-charge-allocation-20261027"></a>'
    root=f'''{root_marker}
2026-10-27完成[总电荷／分配坐标与共同剖面搜索对照](pileup/experiments/alpha_pair_charge_allocation_20261027/README.md)：四种变体在1800帧上比较，保持同一物理似然与电荷域，以相同280个原生单事件分别校准增峰门槛。原生先强后弱40 ns严格恢复由基线{baseline40}/40变为allocation_profile的{new40}/40，200 ns由{baseline200}/40变为{new200}/40。{decision}七项数值测试、13680条评分和3600条已知K诊断核验通过；五幅固定编号图、全长分量波形及完整控制记录在本地Notebook0004第82–86节与[L60](docs/literature-use.md#alpha-pair-charge-allocation-20261027)。
'''
    literature=f'''{root_marker}
### L60：2026-10-27——总电荷／分配坐标与共同剖面搜索的实际比较

**证据与审核边界**：本轮只读共享库入口、主题／审核索引、相关方法页及Fixsen2004完整既有笔记，并复用L59已核的Liddick2012记录，没有联网、重查原PDF或修改文献库。Fixsen等，*Pulse estimation in nonlinear detectors with nonstationary noise*，DOI [10.1016/j.nima.2003.11.313](https://doi.org/10.1016/j.nima.2003.11.313)，来源`papers/sources/2004_Fixsen_NIMA_Pulse-estimation-in-nonlinear-detectors-with-nonstationary-noise.md`，PDF2式1–2、PDF3；TES模拟0.5–5 eV、每能量5000次，采用条件均值／协方差原则，原文WLS不等于当前完整高斯NLL。Liddick、Darby、Grzywacz，2012，*Algorithms for pulse shape analysis using silicon detectors*，DOI [10.1016/j.nima.2011.12.034](https://doi.org/10.1016/j.nima.2011.12.034)，来源`papers/sources/2012_Liddick_NIMA_Algorithms-for-pulse-shape-analysis-using-silicon-detectors.md`，PDF3–8、图3–13；约65 µm DSSD、40 MSPS／12 bit、25 µs记录及实测单脉冲复制／缩放／平移合成双脉冲，只采用检测与估计分开评价。两篇科学内容仍为AI草稿待人工审核；B23-022没有新增人工科学确认，门选异文保留。当前新坐标、缩放及数值profile是项目方案，不声称是文献原算法或已获GAGG证实。

**预定设计**：280个原生独立单事件门槛帧＋760已公开原生评价＋760匹配高斯draw0评价，共1800帧。每种输入三组合×8/40/200 ns×40对，加空和四档单事件各80帧；draw1未运行。训练280／门槛280／评价320母事件互不交叉，条件间仍复用，不算新独立确认。纯α、固定95/351 ns双τ、2 ns×1024完整电压、充分动态范围、未知K≤2。原生Q为实现雪崩等效pe，高斯Q为相同数值条件参数；t为有效面首光子，不作沉积能量。没有新ROOT、训练或电荷先验，worker阻止打开truth／case／ROOT。

保持均值、协方差、logdet及K0/K1；所有变体保留原最佳K2。component_local与allocation_local都在原候选上最多70次局部迭代；后者用`t1/10,t2/10,100lnS,10r`，S=Q1+Q2，Q1=L+(U−L)r，L=max(50,S−50000)，U=min(50000,S−50)，保持各Q为50–50000 pe。MLE目标不加Jacobian，时刻范围仍起点±24 ns；数值缩放和有限差分的物理步长也改变，不宣称增加信息。

两个profile变体共享完全相同起点：观测单事件时刻前／后4/8/12/20/40/80/160/320 ns，第一Q份额0.01/0.03/0.05/0.10/0.15/0.50/0.85/0.90/0.95/0.97/0.99；每点完整可行域中数值优化lnS，最多24次、容差10⁻⁴，保留原值。各方向选4个间隔的最佳份额，共最多8个起点，各20次筛选后两个候选各70次精修；复用local结果。有限标量与局部搜索不是全局极小值证明。各新搜索用相同280单事件重新校准NLL增峰参考分布，Q分档与加一／并列上尾秩p≤0.025保持；原数值门槛另作敏感性对照。

**实际结果**：各Q误差≤20%、各t误差≤min(10 ns,间隔/4)、无额外输出，真Q均≥500，每点40对：

{results}

原生40 ns的allocation_profile有{diag['positive_cost_gain']}/40候选代价下降，强制K2成功{diag['known_K_success']}/40，实际盲成功{new40}/40，分别记录。全部条件和输入的漏检、额外输出、误拆／误并、Q偏差／sigma68误差宽度、t误差、K2与时间匹配、阈下候选及未收敛结果保留。sigma68为匹配条件下相对Q误差宽度，没有沉积能标，不作能谱FWHM。

{chr(10).join(controls)}

**选择要求与决定**：原生空／单额外输出≤5%，单事件严格成功≥95%、Q偏差中位绝对值≤5%；先强后弱40 ns至少增10个百分点、200 ns不减；先弱后强8 ns最多降5个百分点，中／中每点≥95%，六个强弱组平均不下降。{decision}

{statuses}

**不确定性与成本**：40母事件块按前20／后20分层，2000次共同权重重抽保留方法／条件／间隔／输入复用。中／中覆盖40块，强弱覆盖前20块、每块2对；不是逐帧独立区间。全成／全败的[0,0]为样本退化，不是总体保证。原生每帧增量：{timing_text}；profile包含独立使用时的共享剖面成本，本次实际只计算一次。时间排除缓存基线，并发逐帧wall time不作实时吞吐指标。

**完成的验证**：7项数值检查和32帧先导通过；1800帧／7200变体候选组、7600主评分／13680含门槛敏感性评分、252汇总、3600强制K2诊断核验。原基线逐条一致，K0/K1不变，K2最低代价不劣化，最大重放相对差{checks['maximum_cost_replay_error']:.3g}；{checks['scalar_profile_points']}个标量profile点中{checks['scalar_profile_not_converged']}个未报收敛，全部保留。源码、模型、配置和报告冻结，未改默认方法。实验为[alpha_pair_charge_allocation_20261027](../pileup/experiments/alpha_pair_charge_allocation_20261027/README.md)。本地`jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb`第82–86节，图0004-82-1、83-1、83-2、84-1、85-1；旧410单元保留。预定40／200 ns第0对显示完整0–2048 ns实际原生分量、条件均值与最终盲输出，选K=1时弱输出画零。

**待验证**：结果限于已公开三组合／三间隔、固定双τ、局部K≤2。连续Q／幅比外推、位置和噪声失配、α/γ、free-τ、高率K≥3及实物均未验证。若同一NLL搜索收益不足，优先检验估计目标、响应／噪声模型和训练得到的先验；不将四档真实能量用作测试先验。符合开发要求的方案仍需扩大条件后冻结、新母样本确认。有限样本估计与非规则模型选择的一手统计论证、Liddick门选异文等继续列为共享库单独核查待办，不从本项目修改文献库。
'''
    previous=load(RECORD/'document-delivery.json') if (RECORD/'document-delivery.json').exists() else None
    outputs={PROJECT/'README.md':(root_marker,root),PROJECT/'docs/literature-use.md':(root_marker,literature),EXPERIMENT/'README.md':(None,readme)}
    hashes={}
    for path,(marker,text) in outputs.items():
        original=path.read_text() if path.exists() else ''
        if marker:
            if marker in original:
                assert previous and digest(path)==previous['files'][str(path.relative_to(PROJECT))],'Retain subsequent user changes'
                original=original.split(marker)[0]
            content=original.rstrip()+'\n\n'+text.strip()+'\n'
        else:
            assert not original or previous and digest(path)==previous['files'][str(path.relative_to(PROJECT))],'Retain existing experiment README'
            content=text
        path.write_text(content)
        hashes[str(path.relative_to(PROJECT))]=digest(path)
    save(RECORD/'document-delivery.json',dict(passed=True,files=hashes,source_sha256=digest(__file__)))
    print('ALLOCATION_DOCUMENTS_DELIVERED',len(hashes),flush=True)


if __name__=='__main__':deliver()
