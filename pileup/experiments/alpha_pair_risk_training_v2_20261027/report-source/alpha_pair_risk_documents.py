"""Checked experiment/README/literature summaries for continuous-Q learning."""
from __future__ import annotations

from .alpha_pair_risk_study import PROJECT,EXPERIMENT,OUTPUT,RECORD,NAME
from .alpha_pair_information_study import load,save,digest


def deliver():
    for name in ('science-validation','notebook-delivery','figure-visual-review'):
        assert load(RECORD/(name+'.json'))['passed'],name
    report=load(OUTPUT/'report.json');check=load(RECORD/'science-validation.json');note=load(RECORD/'notebook-delivery.json')
    methods=('baseline','nn_absolute','nn_logarithmic','nn_relative','mlp_nll','mlp_regularized')
    aggregates={(r['domain'],r['method'],r['part']):r for r in report['aggregates']}
    summary={(r['domain'],r['method'],r['condition'],r['band']):r for r in report['summary']}
    lines=['| 方法 | 新高斯600对 | 其中分布改变240对 | 已公开原生360对 |','| --- | ---: | ---: | ---: |']
    for method in methods:
        rr=[aggregates[d,method,part] for d,part in (('continuous','all'),('continuous','prior_shift'),('native','all'))]
        lines.append('| '+method+' | '+' | '.join(f"{r['success']}/{r['frames']} ({100*r['efficiency']:.1f}%)" for r in rr)+' |')
    results='\n'.join(lines)
    lines=['| 方法 | 原生先强后弱8 ns | 40 ns | 200 ns |','| --- | ---: | ---: | ---: |']
    for method in methods:
        lines.append('| '+method+' | '+' | '.join(str(summary['native',method,'high_weak',band]['success'])+'/40' for band in ('8','40','200'))+' |')
    highweak='\n'.join(lines)
    lines=['| 方法 | 开发要求全部通过 | 未通过项 |','| --- | --- | --- |']
    for row in report['decisions']:
        lines.append('| '+row['method']+' | '+str(row['passed'])+' | '+('、'.join(k for k,v in row['checks'].items() if not v) or '无')+' |')
    statuses='\n'.join(lines)
    promising=[r['method'] for r in report['decisions'] if r['passed']]
    decision='通过预定开发要求的方案：'+('、'.join(promising) if promising else '无')+'；不自动替换默认方法。'
    selected=report['selected_training_risk']
    train=[]
    for kind in ('absolute','logarithmic','relative'):
        r=load(EXPERIMENT/'models'/kind/'selection.json');runtime=load(RECORD/('training-'+kind+'.json'))
        train.append(f"{kind}第{r['step']}步，验证强制K2成功{100*r['forced_pair_success']:.2f}%，总训练{runtime['seconds']:.1f}s")
    text=f'''# 连续电荷监督学习与物理拟合

2026-10-27执行，承接同一似然下电荷分配搜索未明显改善的结果。比较三种电荷损失的小型MLP，并实际运行MLP＋原NLL和MLP＋学习型Q惩罚NLL。{decision}

## 数据、物理量与冻结顺序

纯α，固定95/351 ns双τ，2 ns×1024完整模拟电压，充分动态范围，局部未知K≤2。Q使用pe，原生数据为实现雪崩等效电荷；新高斯Q是均值／协方差的条件参数，不作新实现雪崩电荷或沉积能量。t为有效面首光子。本轮无新ROOT或新光学母事件，未改变探测器参数。

9600训练＋960验证，事件数概率0.1/0.3/0.6。单Q为500–25000 pe对数均匀，双事件混合独立对数均匀、强弱随机先后和近等强，个体Q同域、总Q≤48000。间隔4–320 ns的均匀／对数均匀各一半。它是声明的模拟设计分布，非实测源能谱或旧四档测试答案。使用冻结的ConservingResponse与JointNoise条件高斯抽样；FIR前三点的残差前缀设零，原物理NLL评价完整1021个创新。

训练／验证／评价参数种子基数依次270310001/270310002/270310003，每帧seed=base×100000+index；噪声使用SeedSequence([seed,993])独立流。网络初始化270310010，小批次270310011；每种损失只有一个训练初始化。每500步按公共验证指标选择：强制K2严格整对正确最高，其次Q绝对ln误差中位最低，再次count BCE最低。各模型及混合所选损失/宽度先冻结，随后生成新评价。

新连续评价1000帧：360对为近等强／先强后弱／先弱后强×6–12/24–64/128–256 ns×40对，弱600–1500、强12000–24000 pe；另240对两Q独立均匀500–23000、间隔均匀4–320 ns，检验训练分布改变；80空＋四档连续单Q各80。它是模型内新参数／噪声样本，不是未见物理确认。

原生760帧沿用已公开确认批（三组合×8/40/200 ns×40对＋400控制），另280个原生单事件校准门槛。母事件分组与旧模型标定／门槛／评价保持不交叉，条件间复用。旧评价本轮只作公开迁移诊断。所有方法同一评价帧，基线逐项重放核对。

## 实现与先导修订

观测单事件初拟合生成1408维特征：完整1021个白化FIR残差＋385个观测初拟合对齐残差＋lnQ与时刻；asinh(z/3)、训练集均值／标准差归一化后截到±8。这是网络数值预处理，不是ADC裁剪。MLP隐藏层128/96、SiLU，输出两个时间槽（−384至384 ns/2 ns网格、峰附近五点平均）、总Q和分配份额，以及K2分数。最终Q候选clip到50–50000 pe，盲输出仍要求每Q≥500。训练目标为模拟事件数/Q/t参数，没有将测试真值或实际分量波形送入推断。

三损失为固定尺度绝对平方((Qhat−Q)/5000)²、对数平方ln(Qhat/Q)²、相对平方(Qhat/Q−1)²。总损失count BCE＋3×双事件平均Q损失＋0.4×双事件时间CE，时间软标签σ=2 ns；3000步、batch64、Adam0.0007、相同初始化及小批次流。NumPy实现，未新增ML依赖。

初版`alpha_pair_risk_training_20261027`完整保留。24帧先导及训练输入发现单事件初拟合总Q偏差：5757个训练双事件有740个所需|ln(Qsum/Qpilot)|>0.1，范围−0.758758至0.425039。因此正式训练及评价前另建v2，把总量头从Qpilot×exp(0.1tanh(z))改成Qpilot×exp(tanh(z))。生成函数AST相同、数组逐字节哈希相同，v2复用初版训练/验证数据，并重做梯度、过拟合先导和物理混合检查；没有重抽样选容易样本。`preflight-head-amendment.json`记录依据与顺序。

验证选择：{'；'.join(train)}。两个混合流程都用`{selected}`，lnQ惩罚宽度为{report['regularization_sigma']}。网络和原最佳K2两起点，各最多70次L-BFGS-B、时刻起点±24 ns，保留最佳初值；regularized加入Σ[(lnQ−lnQnetwork)/σ]²，σ=max(0.2,验证全部双事件的68%绝对lnQ误差)，包含错时刻。网络与拟合看同一波形，所以是数据依赖正则化，不是独立先验或新增观测证据。

网络分数不当似然、不再乘事件率先验。各输出使用同一280单事件、拟合单Q边界1500/5000/12500 pe、加一并保留并列的上尾秩≤0.025，各自校准增峰；保留原K0/K1及每Q≥500。整体比较同时改变候选和检测分数，额外报告强制K2诊断，不能把所有差异归因于拟合。推断禁止打开truth/case/ROOT；失败与未收敛输出均计入。

## 实际结果与控制

严格整对正确：两个真实Q≥500，每个Q误差≤20%，每个时刻误差≤min(10 ns,gap/4)，无额外输出。合计数的条件配比是人为设计，不能解释为某注入率的效率。

{results}

{highweak}

各组空／单控制、Q偏差/σ68、时间误差、漏检、额外输出、误拆／误并、K2输出及强制K2候选均记录在CSV。误差宽度是时间匹配条件下的(p84−p16)/2，没有沉积能标，不能当能谱FWHM。额外输出可能为有信号区域的错时刻或误拆。

预定开发要求：新连续双事件总成功率至少+5个百分点；连续及原生各空／单组额外输出≤5%，单组成功≥95%，匹配Q偏差中位绝对值≤5%；原生双事件总数不下降。

{statuses}

{decision} 成对区间2000次bootstrap（种子270310099）：连续按独立帧，原生每条件按母事件块（强弱20块×2对，中／中40块）；没有把原生合并数按独立帧计算区间。区间限当前模型/公开资料，退化区间不等于物理确定性。无训练多种子、连续能量新母样本、位置失配或实物确认，不能推断普适间隔极限、最高率、α/γ、free-τ或高率K≥3。

## 验证、输出与重放

6项数值测试通过：三风险标量最优值、三损失输出梯度、网络梯度、时间解码/总Q、门槛并列与阈下保护、保存/Adam下降。24帧先导后完成三次正式训练；四核上BLAS单线程。1760评价帧×六方法＝10560主评分，960双事件×六方法＝5760强制K2诊断；原生基线逐项一致，K0/K1不变，NLL混合保留不差于原K2的值。网络重放最大差{check['maximum_network_replay_error']:.3g}、NLL最大相对差{check['maximum_cost_replay_relative_error']:.3g}。网络耗时表只包括batch forward/解码，排除了初拟合特征和基线；混合耗时也是增量，不声称实时吞吐。

本地Notebook0004第87–92节新增{note['new_cells']}单元、执行{note['executed_cells']}代码单元；原431单元及输出保持。图0004-87-1、88-1、89-1、90-1、90-2、91-1固定编号。预定native先强后弱40/200 ns第0对展示全0–2048 ns、原生真实分量与条件均值；最终只选K1时第二输出画零。

源码在`inference-source/`及`evaluation-source/`，`freeze.json`、`evaluation-freeze.json`、`training-freeze.json`与`report-freeze.json`分别保护计算、推断、权重与报告。首次生成基线发现冻结包缺少既有`alpha_pair_residual`依赖，尚未写出任何基线预测；另加`dependency-source/`与`dependency-freeze.json`、`evaluate_complete.py`补齐，来自旧确认实验的相同源码，未改算法、权重或测试配置，原快照保留。实际HEAD为`{load(RECORD/'environment.json')['git_head']}`加冻结代码；之前信息/分配研究未提交改动保留。

从项目根目录运行，保留本地输入、模型及运行记录，无需切换Git：

```bash
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/evaluate_complete.py baseline --workers 4
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/evaluate_complete.py predict --domain continuous --workers 4
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/evaluate_complete.py predict --domain native --workers 4
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/report.py analyze
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/report.py figures
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/report.py notebook
```

推断按index续跑。数据生成/训练不覆盖现有目录；`run.py generate --role train|validation|evaluation`、`train --risk absolute|logarithmic|relative`与`freeze-training`供新独立输出路径使用，先做对应prepare/先导并保留冻结顺序。v2数据来源、生成记录、先导修订及全部文件哈希在`jupyter/run_records/{NAME}/`；数组、候选、评分/CSV/图在`jupyter/artifacts/{NAME}/`。Notebook、模型和数据按约定仅本地维护，单独克隆代码不包含输入或训练权重。

文献来源和采用边界见`docs/literature-use.md` L61：Morad2024的分阶段评价和Kim2023的监督/弱信号边界仅作为设计依据；本轮MLP/三风险/惩罚不是其原算法复现。
'''
    marker='<a id="alpha-pair-risk-training-20261027"></a>'
    best=aggregates['continuous','mlp_regularized','all'];base=aggregates['continuous','baseline','all']
    native=aggregates['native','mlp_regularized','all'];nativebase=aggregates['native','baseline','all']
    root=f'''{marker}
2026-10-27完成[连续Q监督学习与物理拟合](pileup/experiments/{NAME}/README.md)：9600训练、960验证，比较绝对／对数／相对Q三损失；按验证选择`{selected}`用于两种混合拟合，冻结后另生成1000个连续Q高斯评价样本。带学习型Q惩罚的NLL相对原物理基线，严格双事件成功为新高斯{base['success']}/600→{best['success']}/600、已公开原生迁移{nativebase['success']}/360→{native['success']}/360；完整六流程及分条件结果保留。{decision}没有新光学母事件，不能当物理独立确认。6项数值检查、10560主评分与5760已知K诊断核验通过；Notebook0004第87–92节、六幅固定编号图及[L61](docs/literature-use.md#alpha-pair-risk-training-20261027)记录训练、控制和全窗分量。
'''
    literature=f'''{marker}
### L61：2026-10-27——连续Q监督目标与物理拟合的实际比较

**既有证据与审核状态**：只读共享文献库的主题/审核入口、数字脉冲/PSD/光产额页及两篇既有sources记录，本轮未联网、未重查PDF、未修改库。Morad等，2024，*Model-Based Deep Learning Algorithm for Detection and Classification at High Event Rates*，DOI [10.1109/TNS.2024.3371573](https://doi.org/10.1109/TNS.2024.3371573)，来源`papers/sources/2024_Morad_TNS_Model-Based-Deep-Learning-Algorithm-for-Detection-and-Classification-at-High-Event-Rates.md`，PDF3–5式4–9、算法1–2，PDF5–7及PDF10表II。25 mm CLYC、PMT800 V、1 GS/s、17 mV阈值，Poisson线性叠加+AWGN、K≤3；采用检测/到达时刻/分类分开评价，不将6 Mcps检出81.8%、条件分类89.1%或联合72.9%误作本项目电荷恢复率。原文误报14.1%的分母仍有解释缺口。人工审核第553项仅限定确认联合指标区别及适用条件，未全面确认方法；原笔记关于NoA posterior/likelihood及再乘Poisson项的疑问保持。当前K2 logit不当似然，不重复乘源率先验。

Kim等，2023，*Restoring original signals from pile-up using deep learning*，NIMA1055/168492，DOI [10.1016/j.nima.2023.168492](https://doi.org/10.1016/j.nima.2023.168492)，来源`papers/sources/2023_Kim_NIMA_Restoring-original-signals-from-pile-up-using-deep-learning.md`，PDF2–6、表1–2。TexAT Micromegas，原512×40 ns、100点峰对齐输入，约40000组合样本60/20/20，CNN波形MSE；弱Type2幅度误差15.2%与强Type1约5%说明弱目标单列的必要性，不移用为GAGG性能。时间表1.245×40=49.8 ns与正文48 ns差异保留；其第二信号常作为噪声，不等于本项目两个α都要恢复。科学内容仍为AI精读稿待人工审核；关联书目更新不提升审核级别。

**本项目方案及实际执行**：NumPy MLP1408→128→96→773，观测单事件初拟合创新作为输入；只在训练使用模拟K/Q/t。三组相同初始化/小批次，绝对((Qhat−Q)/5000)²、ln(Qhat/Q)²、相对(Qhat/Q−1)²。总损失BCE+3Q+0.4时间CE，3000步、batch64、lr0.0007，初始化270310010、抽样270310011；不是上述文献的原网络或性能复现。9600训练/960验证来自连续Q设计分布，非四档测试真值先验。初版总Q头修正过窄，先导/训练输入核查后、正式训练和评价前另建v2，将exp(±0.1)改成exp(±1)；数据及生成函数一致并有哈希、初版保留。

每500步按验证强制K2严格成功、Q绝对ln误差中位、count BCE选择。选择结果：{'；'.join(train)}。混合方案选择`{selected}`，两种物理精修共用原K2/网络两起点、各70次迭代及±24 ns范围；带惩罚者增加Σ[(lnQ−lnQnet)/σ]²，σ取验证全部双事件68%绝对lnQ误差、下限0.2，当前为{report['regularization_sigma']}。这是同一波形导出的数据依赖正则，不是独立贝叶斯先验。网络分数及各输出用同一280原生单事件分档校准加一上尾秩≤.025，K0/K1沿用基线，每候选Q≥500；所有失败保留，worker禁读truth/case/ROOT。

纯α、95/351 ns固定双τ、2 ns×1024完整模拟电压、充分动态范围、未知局部K≤2。原生Q为实现雪崩等效pe，高斯Q为条件生成参数，t为有效面首光子，不是沉积能量。训练/验证/评价种子基数270310001/2/3，每帧base×100000+index，噪声SeedSequence([seed,993])。冻结权重后新高斯1000帧含600对和400控制；已有原生760帧是公开迁移诊断，另280门槛。高斯沿用同一已标定模型，不是新G4/物理确认。

严格整对正确为每Q±20%、每t±min(10 ns,gap/4)、无额外输出。实际结果：

{results}

{highweak}

{statuses}

{decision} 合计条件配比由本轮设计，不代表源率。完整误差偏差/σ68、漏检/额外输出/误拆/误并、强制K2、母事件块配对区间均保留；σ68是时间匹配条件下相对Q误差宽度，非能谱FWHM。6项数值检查、24帧先导、模型/代价重放及10560主评分/5760诊断通过；不依残差下降或训练loss宣称恢复改善。

**记录与待验证**：`jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb`第87–92节，六幅图0004-87-1、88-1、89-1、90-1、90-2、91-1；冻结实验`pileup/experiments/{NAME}/`，本地结果/记录在同名artifacts/run_records。旧431单元及此前研究保留。仅一个训练种子，尚需训练随机性、连续能量新光学母事件、位置/噪声失配、完整随机高率和α/γ验证。后续先按时间候选、Q偏差和计数门槛分解失败；如采用学习型先验，须单独校准覆盖/偏差，不使用旧测试四档Q作答案。需要新统计/后验校准文献整理时交由文献库处理，本项目不修改共享库。
'''
    updates={EXPERIMENT/'README.md':text,PROJECT/'README.md':root,PROJECT/'docs/literature-use.md':literature}
    record=load(RECORD/'document-delivery.json') if (RECORD/'document-delivery.json').exists() else None
    for path,content in updates.items():
        relative=str(path.relative_to(PROJECT))
        if path==EXPERIMENT/'README.md':
            if path.exists():assert record and digest(path)==record['files'][relative]
            path.write_text(content)
        else:
            old=path.read_text()
            if marker in old:
                assert record and digest(path)==record['files'][relative]
                old=old[:old.index(marker)].rstrip()+'\n'
            backup=RECORD/('README.before.md' if path==PROJECT/'README.md' else 'literature.before.md')
            if not backup.exists():backup.write_text(old)
            path.write_text(old.rstrip()+'\n\n'+content)
    save(RECORD/'document-delivery.json',dict(passed=True,files={str(p.relative_to(PROJECT)):digest(p) for p in updates},source_sha256=digest(__file__)))
    print('RISK_DOCUMENTS_DELIVERED',flush=True)


if __name__=='__main__':deliver()
