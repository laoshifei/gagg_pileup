"""Figures and append-only notes for conditional grouped mean-shape stability."""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR','/tmp/gagg_shape_stability_mpl')
os.environ.setdefault('IPYTHONDIR','/tmp/gagg_shape_stability_ipython')
os.environ.setdefault('JUPYTER_RUNTIME_DIR','/tmp/gagg_shape_stability_jupyter')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nbformat
from nbclient import NotebookClient
import numpy as np

from .alpha_mean_shape_stability_study import ROOT,EXP,ART,REC,NAME,read,save,sha,verify


def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(map(str,r))+' |' for r in rows])


def write_csv(name,rows):
    (EXP/'tables').mkdir(exist_ok=True)
    with (EXP/'tables'/f'{name}.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader()
        writer.writerows({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v
                         for k,v in r.items()} for r in rows)


def pooled_mean_control():
    folds=np.repeat(np.arange(4),8);values=np.repeat([-1.,-.5,.5,1.],8)
    residual=np.array([values[i]-values[folds!=folds[i]].mean() for i in range(32)])
    ratio=float(np.mean(residual**2)/np.mean(values**2))
    assert abs(residual.mean())<1e-14 and abs(ratio-16/9)<1e-14
    result=dict(description='Synthetic equal-fold mean subtraction;not an extra detector candidate',
        old_mean=float(values.mean()),corrected_pooled_mean=float(residual.mean()),
        corrected_to_old_mse=ratio,algebra='Each held-out residual subtracts the other three fold means;pooled mean is identically zero,while MSE can increase')
    save(EXP/'pooled-mean-control.json',result)
    return result


def figures(summary,partitions,bootstraps):
    plt.rcParams.update({'font.size':10,'figure.dpi':140,'axes.grid':True,'grid.alpha':.2})
    paths={};gaps=['8 ns','40 ns','200 ns'];x=np.arange(3)
    fig,axes=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    ratio=[r['mse_ratio'] for r in summary['projection']]
    ci=np.array([r['refit_bootstrap_ci95']['mse_ratio'] for r in summary['projection']])
    axes[0].vlines(x,ci[:,0],ci[:,1],color='tab:blue',lw=2)
    axes[0].scatter(x,ratio,color='tab:blue',s=40,zorder=3)
    axes[0].axhline(1,color='black',ls='--',lw=.8)
    axes[0].set(xticks=x,xticklabels=gaps,ylabel='Corrected / old projected MSE',
                title='(a) Primary 4-fold: refit bootstrap 95% range')
    for i,value in enumerate(ratio):axes[0].annotate(f'{value:.3f}',(i,value),xytext=(8,0),textcoords='offset points')
    wave=summary['primary']['waveform'];x=np.arange(4)
    axes[1].bar(x-.16,[r['mean_mse_ratio'] for r in wave],width=.32,label='Pooled mean waveform')
    axes[1].bar(x+.16,[r['mother_mse_ratio'] for r in wave],width=.32,label='Mean of mother MSEs')
    axes[1].axhline(1,color='black',ls='--',lw=.8)
    axes[1].set(xticks=x,xticklabels=['weak','low','mid','high'],ylabel='Corrected / old waveform MSE',
                title='(b) Distinct losses: mean bias and individual mothers')
    axes[1].legend(fontsize=8)
    fig.suptitle('Figure 0007-115-1: Development cross-fitting of the incremental correction')
    p=ART/'figures/fig_0007_115_1_crossfit.png';fig.savefig(p);plt.close(fig);paths['115']=p
    fig,axes=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    for j,label in enumerate(gaps):
        values=[r['projection'][j]['mse_ratio'] for r in partitions]
        axes[0].plot(np.arange(1,11),values,'o-',label=label)
    axes[0].axhline(1,color='black',ls='--',lw=.8);axes[0].axvline(1,color='grey',ls=':',lw=.8)
    axes[0].set(xlabel='Fixed partition (1 is primary)',ylabel='Corrected / old projected MSE',
                title='(a) Same 32 strong mothers in ten partitions',xticks=np.arange(1,11))
    axes[0].legend(fontsize=8)
    shifts=np.array([[m['nominal_shift_100logE'] for m in r['models']] for r in bootstraps]).reshape(-1,3)
    axes[1].boxplot([shifts[:,j] for j in range(3)],labels=gaps,whis=(2.5,97.5),showfliers=False)
    axes[1].axhline(0,color='black',lw=.8)
    axes[1].set(ylabel='Projected correction / 100 delta ln E',
                title='(b) Refit fold correction at nominal 7.5 MeV')
    fig.suptitle('Figure 0007-116-1: Split sensitivity and training uncertainty (dependent resamples)')
    p=ART/'figures/fig_0007_116_1_training_uncertainty.png';fig.savefig(p);plt.close(fig);paths['116']=p
    fig,axes=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    rows=summary['variance'];x=np.arange(3)
    between=np.array([r['between_variance_unclipped'] for r in rows]);within=np.array([r['within_variance']/8 for r in rows])
    axes[0].bar(x,between,label='Between mothers, B')
    axes[0].bar(x,within,bottom=between,label='Within / 8')
    axes[0].set(xticks=x,xticklabels=gaps,ylabel='Variance / (100 delta ln E)$^2$',
                title='(a) Old residual: variance of 8-repeat means')
    axes[0].legend(fontsize=8)
    repeats=np.array([1,2,4,8,16,32,64])
    for r,label in zip(rows,gaps):
        value=np.ceil(1.96**2*(r['between_variance_unclipped']+r['within_variance']/repeats)/100)
        line=axes[1].plot(repeats,value,'o-',label=label)[0]
        axes[1].axhline(r['planning_mothers_halfwidth10_infinite_repeats'],color=line.get_color(),ls=':',lw=.8)
    axes[1].set_xscale('log',base=2);axes[1].set_xticks(repeats,labels=repeats)
    axes[1].set(xlabel='Digitizations per mother',ylabel='Estimated total independent strong mothers',
                title='(b) Normal-approximation planning, half-width 10')
    axes[1].legend(fontsize=8)
    fig.suptitle('Figure 0007-117-1: More repeats cannot remove between-mother variation')
    p=ART/'figures/fig_0007_117_1_mother_budget.png';fig.savefig(p);plt.close(fig);paths['117']=p
    return paths


def main():
    verify();assert read(EXP/'replay.json')['passed']
    result_dir=EXP/'results/primary'
    for p,d in read(result_dir/'output-freeze.json').items():assert sha(ROOT/p)==d,p
    summary=read(result_dir/'summary.json');parts=read(result_dir/'partitions.json');boots=read(result_dir/'bootstrap.json')
    plots=figures(summary,parts,boots);algebra=pooled_mean_control()
    write_csv('primary_projection',summary['projection']);write_csv('primary_waveform',summary['primary']['waveform'])
    write_csv('variance_planning',summary['variance'])
    write_csv('partition_projection',[dict(partition=i,**r) for i,p in enumerate(parts) for r in p['projection']])
    write_csv('bootstrap_projection',[dict(draw=i,**r) for i,p in enumerate(boots) for r in p['projection']])
    projection_table=table(['间隔','旧→修正平均投影','逐母投影MSE比及重训练95%分位范围','改善划分/10'],[
        [f"{r['gap_ns']:g} ns",f"{r['old_mean']:+.2f}→{r['corrected_mean']:+.2f}",
         f"{r['mse_ratio']:.4f} [{r['refit_bootstrap_ci95']['mse_ratio'][0]:.4f}, {r['refit_bootstrap_ci95']['mse_ratio'][1]:.4f}]",
         r['improving_partitions']] for r in summary['projection']])
    wave_table=table(['档位','组平均波形MSE比','逐母波形MSE均值比'],[
        [('弱','低','中','强')[r['group']],f"{r['mean_mse_ratio']:.4f}",f"{r['mother_mse_ratio']:.4f}"]
        for r in summary['primary']['waveform']])
    variance_table=table(['间隔','W','B（未截零）','8次平均母间份额','当前均值标准误','8次/母的规划总母数','无限重复的规划下限'],[
        [f"{r['gap_ns']:g} ns",f"{r['within_variance']:.2f}",f"{r['between_variance_unclipped']:.2f}",
         f"{100*r['between_fraction_mean8']:.2f}%",f"{r['mean_standard_error']:.2f}",
         r['planning_mothers_halfwidth10_repeat8'],r['planning_mothers_halfwidth10_infinite_repeats']]
        for r in summary['variance']])
    shifts=np.array([[m['nominal_shift_100logE'] for m in r['models']] for r in boots]).reshape(-1,3)
    shift_quantiles=np.percentile(shifts,[2.5,50,97.5],axis=0)
    shift_text='；'.join(f'{gap} ns为[{shift_quantiles[0,j]:+.2f}, {shift_quantiles[2,j]:+.2f}]'
                         for j,gap in enumerate((8,40,200)))
    key=(f"128个历史开发母、十组按母分层四折和1000次重新训练重采样后，8/40 ns逐母投影MSE比为"
         f"{summary['projection'][0]['mse_ratio']:.4f}/{summary['projection'][1]['mse_ratio']:.4f}，"
         f"重训练95%分位范围均跨1；十种划分分别仅{summary['projection'][0]['improving_partitions']}/10和"
         f"{summary['projection'][1]['improving_partitions']}/10改善，均未通过预设开发稳定性门。")
    caution=('四折的留出残差不是独立估计：同一母事件会进入其他折的训练。对于等规模四折，'
             '若修正只是“其他三折的均值”，合并留出残差的均值就恒为零，即使逐母MSE增大。'
             '因此合并均值下降及其窄重采样范围不能单独证明泛化；本轮主要终点为逐母投影MSE。')
    limitations=('只对新增修正交叉拟合；旧能量响应/空间先验曾用过这些128母，且本轮方法在观察L93失败后选定，'
                 '不能当作端到端独立验证。原64个confirmation母和L92的32个强noise-fold母没有进入本轮拟合/评分。'
                 '固定四折层的重采样每层仅8母，包含本批母均值变化和修正重训练，未重新生成数字化噪声，'
                 '也未传播旧模型、投影算子、不同划分和器件系统误差；95%分位范围是条件敏感性描述，不保证新母池总体覆盖率。'
                 '十种划分、四折以及1000次重采样都不是新增独立实验。')
    next_text=('下一步应先冻结较大且未参与旧均值/能标训练的强母事件确认设计，保持本轮修正规则固定，'
               '以逐母投影MSE、均值不确定度及孤立能量控制为主要检查。样本预算优先增加独立光学母事件；'
               '8 ns当前方差给出的粗略量级为334个强母、每母8次，尚非所需新增数量或保证性样本量。'
               '不要在既有confirmation结果上选择修正幅度、平滑宽度或最佳划分。孤立和共享控制通过前，不推进双脉冲能量或高率扫描。')
    sections=[('114',f'''## 114. 条件均值在开发母事件间是否稳定（2026-11-02）

L93的固定均值修正通过面积/孤立能标，却未通过8 ns弱方向迁移门。本轮研究问题是该修正对开发母事件选择有多敏感，以及增加同母重复能解决多少不确定性。只复用L93训练128母的1024次既有数字化（四档各32母、每母8次）；不新增ROOT、输运、数字化或能量拟合。实际源能量1.625/2.5/5.486/8 MeV，E_dep、原ROOT/宏/RunInfo和种子随实验保留。

每档按固定路径哈希分成四折，每折8母/档。每个修正在96母（每档24）上训练，再对其余32母（每档8）的8次均值评分；一个母的全部重复始终同行。十个固定哈希种子737000–737009对应十种划分，第一组为主要结果，其他为划分敏感性。平滑4 ns、因果保护2 ns、起始8 ns、末端160 ns渐消、面积投影、log q0插值及非负尺度规则与L93完全相同，没有超参数搜索。

再在主要划分的“能区×折”16层内各抽取8母（有放回），每次重训四个修正，共1000次。重复抽到的同一母仍保留原折，训练/留出集合不共享母身份。种子737100000+重采样编号。评价是真实E处的祖先残差诊断；虽然新增修正不读取留出信号，真实E仍用于标定和评分，不是盲算法。

开发稳定性门预先冻结：8/40 ns的主要逐母投影MSE比<1、重训练95%分位上界<1、至少8/10划分改善；四档主要组平均波形MSE比≤1.1；所有运行的非负尺度≥0.5。门仅决定是否值得开展新孤立确认，不授权直接进入双脉冲。沿用L91已有精读的有限方法背景；本轮分组/重采样设计为项目分析，未新增论文性能主张或原PDF核查。
'''),('115',f'''## 115. 逐母残差与合并均值必须区分

{key}

{projection_table}

投影方向仍为固定旧模型、名义强7.5/弱0.25 MeV、弱间隔8/40/200 ns的联合logE影响算子。均值单位100δlnE，MSE单位(100δlnE)²，比值无量纲；不是实际能量百分比、非线性双脉冲RMSE或检出效率。主要MSE先对32个强母各自8次均值求投影平方，再平均；没有把四折平均投影的平方冒充该指标。

{caution} 数值反例取四折常数−1、−0.5、0.5、1，各8样本；其合并均值仍为0，而留出MSE/原MSE={algebra['corrected_to_old_mse']:.6f}=16/9。此为代数说明，没有作为额外探测器候选。

{wave_table}

组平均残差波形均明显减小，但弱、低、强档逐母波形MSE并未因此降低；随机母差异和学习误差不能由组均值曲线代表。组平均波形MSE与逐母波形MSE用同一2048 ns/1024点均匀时间权重，二者都不同于GLS能量损失。

图 0007-115-1：(a)主要四折逐母投影MSE比和重训练95%分位范围；(b)四能区组平均波形与逐母波形两种损失的比值。第一幅范围条件于已暴露开发池和固定划分。
'''),('116',f'''## 116. 划分敏感性与训练不确定度

十个划分的8/40/200 ns逐母投影MSE比范围依次为{'; '.join(f"[{r['partition_mse_ratio_range'][0]:.4f}, {r['partition_mse_ratio_range'][1]:.4f}]" for r in summary['projection'])}。同一32个强母反复评价，不能把改善7/10或3/10当作独立成功概率、p值或十次独立复现；也不挑选最好划分的结果替代主要划分。

将1000次重采样中的四个训练折修正统一作用于名义7.5 MeV，投影修正量的合并2.5–97.5百分位为：{shift_text}，单位100δlnE。该量是“将扣掉多少旧残差”，不是修正后能量误差。这4000个模型互相依赖，合并分位只展示不同训练折/重采样的边际散布；8/40 ns包含零，说明短间隔修正方向对训练样本并不稳固。

{limitations}

全部4040个正式训练模型的全局非负尺度最小为{summary['minimum_scale']:.1f}，没有因非负性而缩小或删除模型。未通过的两个短间隔稳定性门和所有划分/重采样均完整保留。

图 0007-116-1：(a)十组固定划分的投影MSE比，第一组为主要划分；(b)1000次重采样×四训练折的名义强脉冲修正投影箱线图，箱为25–75%、须为2.5–97.5%，不把箱内4000模型视作独立样本。
'''),('117',f'''## 117. 母间差异主导：样本预算粗估

对32个强开发母的旧残差投影，W为母内8次无偏方差的母均值，B=母均值样本方差−W/8，保留负B估计。其估计仍可含未去净的能量趋势、几何和固定模型失配，不把B全部归为微观光子统计。

{variance_table}

W、B单位(100δlnE)²；标准误单位100δlnE。母间份额为B/(B+W/8)，不包含平均偏差平方。8/40 ns分别约82.6%/81.4%，所以仅增加同母数字化难以显著缩小总误差。这里分解的是旧开发残差，不能直接当作修正后新母事件的方差。

规划取正态近似的95%均值半宽10（单位100δlnE，不是实际E±10%），总母数N≈ceil(1.96²×(B+W/R)/10²)。R=8得到8/40 ns约334/68个独立强母；R→∞的估计下限仍约276/55。数字是总量，不是已经批准或必需新增的母数；方差只由32母估计，未给规划数字的区间、未含标定偏差及总体改变，不能宣称此预算保证通过验证。200 ns规划数字12/10也不构成减少现有数据的理由。

图 0007-117-1：(a)八次平均的母间B与母内W/8；(b)增加每母重复次数时的规划总母数，虚线为无限重复下限。分析/绘图均使用未删选的32母；仅规划无限重复项为防负数才对B取max(B,0)，原估计保留。
'''),('118',f'''## 118. 核验、限制和下一步

四项针对性测试检查母身份唯一性/分层平衡、重采样重复身份不跨折、矢量化电压与标量及非负尺度一致、留出信号改为NaN不影响训练、有限重复修正保留负B（后两者分属测试用例）。完整128母修正重新计算后，形状和电荷结点与L93逐项相同；缓存非负约束与原逐点计算一致，最大矢量化电压差{read(EXP/'numerical.json')['max_vectorized_difference_mV']:.3g} mV。

正式40个划分模型＋4000个重采样模型，共4040次修正标定；成本先导32、完整规则复算1、局部重放20另列，共4093次实验标定调用（不含单元测试中的合成例）。无新数字化、ROOT或E拟合。先导前8个重采样与正式结果相同；重放划分0/9和重采样8/511/999，5份摘要及12个数组逐项一致。工作进程禁止读取原confirmation数据与ROOT。所有输入、源码、规则与种子在当前结果产生前冻结，原母宏/RunInfo及已知历史暴露保留。

准备阶段原C++编译命令曾发生一次编译器内部崩溃；同命令重试与随后构建成功，未改C++或编译选项。失败源/协议副本与处理记录在本地run_records中，发生于设计冻结和研究计算之前，不删除或计作物理样本。冻结后数值控制与计算完成。

{next_text}

本阶段没有验证孤立/双脉冲E偏差、峰分辨率、漏检、误拆分/误合并、重建效率或能谱选择偏差，不将未测指标填写为零。更换默认算法与高率扩展仍无证据。研究对象继续固定纯α、95/351 ns（快比例0.3）、29 ns恢复代理、14400微胞和12×12开发先验；相位0只用于本轮统计，四相位只作冻结模型非负约束。

源码、规则、实际宏/RunInfo、输入与输出哈希及表格在pileup/experiments/{NAME}；图和数组在jupyter/artifacts/{NAME}，日志/旧笔记副本在jupyter/run_records/{NAME}。旧笔记计算和图像保留，L94与研究计划第35节同步。
''')]
    append_notebooks(sections,plots,key,caution,next_text)
    write_documents(key,projection_table,wave_table,variance_table,caution,limitations,next_text)
    outputs=[*plots.values(),EXP/'pooled-mean-control.json',EXP/'README.md',ROOT/'README.md',ROOT/'docs/literature-use.md',
        *[ROOT/'jupyter'/x for x in ('codex_0007_SiPM_Energy_Calibration.ipynb','Research_Plan_20260921.ipynb','Run_Catalog.ipynb')],
        *sorted((EXP/'tables').glob('*.csv'))]
    save(EXP/'report-output-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in outputs})
    print('Wrote3figures,sections114–118,plan35,Run_Catalog,L94;new displays executed',flush=True)


def append_notebooks(sections,plots,key,caution,next_text):
    names=('codex_0007_SiPM_Energy_Calibration.ipynb','Research_Plan_20260921.ipynb','Run_Catalog.ipynb')
    notebooks={}
    for name in names:
        nb=read(ROOT/'jupyter'/name);before=read(REC/(name+'.before'))
        assert nb==before,name;notebooks[name]=nb
    assert len(notebooks[names[0]]['cells'])==185
    cells=[]
    for number,body in sections:
        cells.append(nbformat.v4.new_markdown_cell(body,id=f'shape-stability-{number}'))
        if number in plots:
            path=plots[number].relative_to(ROOT)
            source=("from pathlib import Path\nfrom IPython.display import display, Image\n"
                    "root = Path.cwd()\nif root.name == 'jupyter':\n    root = root.parent\n"
                    f"display(Image(filename=str(root / '{path.as_posix()}')))\n")
            cells.append(nbformat.v4.new_code_cell(source,id=f'shape-stability-figure-{number}'))
    display_nb=nbformat.v4.new_notebook(cells=cells)
    NotebookClient(display_nb,timeout=120,kernel_name='python3',resources={'metadata':{'path':str(ROOT)}}).execute()
    nb=notebooks[names[0]];nb['cells'].extend(display_nb['cells'])
    nbformat.validate(nbformat.from_dict(nb));(ROOT/'jupyter'/names[0]).write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n')
    for name,heading in [(names[1],'## 35. 均值修正对母事件选择不稳定：先增加独立母事件证据（2026-11-02）'),
                         (names[2],'## 条件均值的分组交叉验证与重训练重采样（2026-11-02）')]:
        body=heading+'\n\n'+key+'\n\n'+caution+'\n\n'+next_text
        body+='\n\n本轮没有新输运、数字化、ROOT或能量拟合；详见codex_0007第114–118节、图0007-115-1至117-1及L94。'
        if name==names[2]:
            body+=(f'\n\n实验pileup/experiments/{NAME}。只复用alpha_mean_shape_20261102训练128母的1024条既有数字化，'
                   '每母8次、相位0；源能量四档1.625/2.5/5.486/8 MeV，E_dep逐母保留。'
                   '10划分×4折＋1000重采样×4折=4040次正式修正标定；先导32、完整规则复算1、重放20，共4093次实验标定调用。'
                   '复用原种子736000000+母ID×100+重复ID，读出另加10000000；本轮划分种子737000–737009，重采样737100000+编号。'
                   '原ROOT/宏/RunInfo在mothers.json、original-metadata/、original-runinfo/，原数字化任务在original-digitization-plan.json。'
                   '4096 ns原数字化、前2048 ns观测；最多8进程、BLAS单线程。不存在新的源注入率或生产ROOT映射。')
        nb=notebooks[name];nb['cells'].append(nbformat.v4.new_markdown_cell(body,id='shape-stability-summary-20261102'))
        nbformat.validate(nbformat.from_dict(nb));(ROOT/'jupyter'/name).write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n')


def write_documents(key,projection_table,wave_table,variance_table,caution,limitations,next_text):
    docs=f'''\n\n<a id="mean-shape-stability-20261102"></a>
## L94 — 2026-11-02：分组重训练未支持短间隔均值修正稳定推广

**项目结果。** {key}

{projection_table}

均值单位100δlnE，MSE单位(100δlnE)²；不是E百分比或双脉冲能量误差。200 ns为次要结果。主要逐母MSE对32个强母的8次残差均值分别平方再平均；四折修正都没有使用自身留出母的信号。

**合并均值的限制。** {caution} 数值反例的合并均值仍为0，MSE却扩大16/9。四档组平均/逐母波形损失也显示不同结论：

{wave_table}

**冻结方法与输入。** 仅用L93的128开发母（四档各32）、1024条旧数字化；源能量1.625/2.5/5.486/8 MeV，实际E_dep分开记录。十组固定哈希四折，每模型96训练/32留出；1000次在主要划分的16个能区×折层内整母重采样，重复身份保留原折，重训四模型。保留L93的4 ns平滑、2 ns因果保护、面积投影与log q0插值，未搜索幅度、平滑或最佳划分。所有4040个正式模型的非负尺度均为1。固定原响应/空间先验与弱方向算子，不在confirmation结果上选择模型。

**方差与规划。**

{variance_table}

W为母内方差均值，B=母均值方差−W/8且不截零，仍可能包含能量趋势/几何/模型失配；份额为B/(B+W/8)，不含均值偏差平方。规划以正态近似95%均值半宽10个100δlnE单位，N≈ceil(1.96²(B+W/R)/100)，不是实际能量±10%。334/68是8/40 ns条件下的总强母数量粗估，不是必需新增或保证通过的样本量；参数只由32母估计，未传播规划不确定度。无限重复时的276/55下限说明重复不能消除母间变化。

**证据层级与核验。** {limitations} 沿用L91既有精读及审核边界，无新原文核查或论文性能主张。分组、重采样、门限和样本规划为项目分析。4项针对性测试通过；完整128母修正规则精确复现L93，缓存非负尺度/插值与原标量一致；先导8次重采样摘要与正式相同；重放2组划分及3次重采样，5份摘要与12个数组逐项一致。正式4040、先导32、规则复算1、重放20，共4093次实验标定调用，无新增数字化、ROOT或E拟合。未测E偏差/峰分辨率、盲漏检、误拆并、效率或能谱选择偏差，不填零。

**下一步。** {next_text}

位置：`pileup/experiments/{NAME}/`；可复用模块`pileup/analysis/alpha_mean_shape_stability.py`、`alpha_mean_shape_stability_study.py`、`alpha_mean_shape_stability_report.py`。本地codex_0007第114–118节、图0007-115-1至117-1，研究计划第35节与Run_Catalog同步。全部旧笔记/图保留，本轮没有默认升级或远程上传。
'''
    p=ROOT/'docs/literature-use.md';old=(REC/'literature-use.md.before').read_text()
    assert p.read_text()==old and '## L94 ' not in old;p.write_text(old+docs)
    note=f'''# 条件均值的母事件分组稳定性

{key}

{projection_table}

{caution}

仅复用上轮训练128母的1024条数字化，每母8次，相位0。4能区各32母；每次三折96母标定、32母留出，10个固定划分；1000次按能区×折整母重采样并重训。旧响应/空间先验固定且已有历史暴露，仅新增修正交叉拟合。原64个confirmation母不进入本轮拟合/评分。没有新ROOT、数字化、能量拟合或默认算法替换。

{variance_table}

方差单位(100δlnE)²，均值标准误单位100δlnE。规划总母数使用正态95%半宽10的粗估，不是实际能量±10%、新增母数量或保证性预算。{next_text}

{limitations}

项目根目录下可执行：

```bash
python3 -B pileup/experiments/{NAME}/run.py verify
python3 -B pileup/experiments/{NAME}/validate.py
```

run.py直接导入本目录源码副本，不依赖切换Git版本。模型/原宏/RunInfo、折分配、种子、输入/源码哈希见protocol.json、design-freeze.json、folds.json、mothers.json、original-metadata/、original-runinfo/。results/primary保存所有划分、重采样与摘要，失败和不改善的划分均保留；不按结果选最好划分。统计输入仅是jupyter/artifacts/{NAME}/inputs.npz与positivity-cache.npz；二者已在本地自包含复制，不要求再读取历史confirmation数组。数组/模型/ROOT和Notebook不随Git发布。

如需完整复算，选一个尚不存在的label；输出写新目录，比较摘要、全部划分和重采样逐值一致。此轮实际只运行了成本先导和预定局部重放，并未额外完整复算：

```bash
python3 -B pileup/experiments/{NAME}/reproduce.py --label another_replay
```

正式4040次修正标定，成本先导32、完整旧规则复算1、局部重放20，共4093次实验标定调用；没有新数字化或E拟合。4项针对性测试、完整规则一致、分层身份隔离和12数组精确重放通过。最多8进程、BLAS单线程。准备阶段原C++编译器一次内部崩溃后同命令重试成功，发生在冻结/研究计算前；失败快照与处理记录留存在run_records，冻结后数值控制通过。

图与数组在jupyter/artifacts/{NAME}，运行记录在jupyter/run_records/{NAME}。本地codex_0007第114–118节含三幅固定编号图，研究计划第35节、Run_Catalog和L94同步。报告源码与输出单独冻结，run_report.py只执行新增显示单元并检查原笔记一致，避免重复追加。
'''
    (EXP/'README.md').write_text(note)
    p=ROOT/'README.md';old=(REC/'README.md.before').read_text();assert p.read_text()==old
    p.write_text(old+f'''\n\n2026-11-02完成[均值修正的分组稳定性诊断](pileup/experiments/{NAME}/README.md)：128个开发母、十组四折和1000次整母重采样重训练，8/40 ns逐母残差投影MSE比0.9497/0.9972，95%重训练分位范围均跨1；十种划分仅7/10与3/10改善，未通过稳定性门。合并留出均值接近零可能部分来自均值扣除的代数结构，不能单独证明泛化。8/40 ns八次平均的投影方差约82.6%/81.4%来自母间差异，下一步应优先增加独立母事件证据并冻结新确认设计。没有新数字化、ROOT、E拟合或默认升级；四项测试及12数组重放通过，结果见[采用记录L94](docs/literature-use.md#mean-shape-stability-20261102)与本地codex_0007第114–118节、三幅固定编号图。\n''')


if __name__=='__main__':main()
