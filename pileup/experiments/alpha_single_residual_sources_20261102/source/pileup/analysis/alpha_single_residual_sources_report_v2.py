"""Figures and append-only local notes for the grouped residual diagnostic."""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import re
import shutil

os.environ.setdefault('MPLCONFIGDIR', '/tmp/gagg_residual_sources_mpl')
os.environ.setdefault('IPYTHONDIR', '/tmp/gagg_residual_sources_ipython')
os.environ.setdefault('JUPYTER_RUNTIME_DIR', '/tmp/gagg_residual_sources_jupyter')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import nbformat
from nbclient import NotebookClient

from .alpha_single_residual_sources_study import ROOT, EXP, ART, REC, NAME, read, save, sha, verify


def table(headers, rows):
    return '\n'.join(['| '+' | '.join(headers)+' |', '| '+' | '.join(['---']*len(headers))+' |']+
        ['| '+' | '.join(map(str,row))+' |' for row in rows])


def figures(summary, mothers, a):
    plt.rcParams.update({'font.size':10,'figure.dpi':130,'axes.grid':True,'grid.alpha':.2})
    outputs = []
    t = a['time_ns']-128.
    visible = (t>=-32)&(t<=600)
    fig, axes = plt.subplots(2,2,figsize=(12,7),constrained_layout=True)
    for j,group in enumerate(('high','weak')):
        ax = axes[0,j]
        for role,color in [('train','tab:blue'),('validation','tab:orange')]:
            idx = [r['id'] for r in mothers if r['group']==group and r['role']==role]
            values = a['signal_residual_mV'][idx].mean(1)
            mean = values.mean(0); se=values.std(0,ddof=1)/np.sqrt(len(values))
            ax.plot(t[visible],mean[visible],color=color,label=role+' (16 mothers)')
            ax.fill_between(t[visible],(mean-se)[visible],(mean+se)[visible],color=color,alpha=.15)
        ax.axhline(0,color='black',lw=.6)
        ax.set(title=f'({chr(97+j)}) {group}: mean residual +/- 1 SE',ylabel='Residual / mV',xlabel='Time after onset / ns')
        ax.legend(fontsize=8)
        ax = axes[1,j]
        ax.plot(t[visible],a[group+'_event_signal_between'][visible],label='Between mothers (unclipped)')
        ax.plot(t[visible],a[group+'_event_signal_within'][visible]/8,label='Within / 8')
        ax.axhline(0,color='black',lw=.6)
        ax.set(title=f'({chr(99+j)}) {group}: estimated variance of 8-repeat mean',ylabel='Variance / mV$^2$',xlabel='Time after onset / ns')
        ax.legend(fontsize=8)
    fig.suptitle('Figure 0007-104-1: Isolated residual means and variance sources')
    p=ART/'figures/fig_0007_104_1_residual_sources.png';fig.savefig(p);plt.close(fig);outputs.append(p)
    rows=[r for r in summary['variance'] if r['group']=='high' and r['kind']=='event_signal' and r['projection'].startswith('weak_')]
    x=np.arange(len(rows));labels=['8 ns','40 ns','200 ns']
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    between=np.array([r['between_variance_unclipped'] for r in rows]);within=np.array([r['within_variance']/8 for r in rows])
    axes[0].bar(x,between,label='Between mothers');axes[0].bar(x,within,bottom=between,label='Within / 8')
    axes[0].set(xticks=x,xticklabels=labels,ylabel='Variance / (100 delta ln E)$^2$',title='(a) Estimated variance after 8 repeats')
    axes[0].legend(fontsize=9)
    mean=np.array([r['mean'] for r in rows]);ci=np.array([r['mean_ci95'] for r in rows])
    axes[1].errorbar(x,mean,yerr=np.stack([mean-ci[:,0],ci[:,1]-mean]),fmt='o',capsize=5)
    axes[1].axhline(0,color='black',lw=.7)
    axes[1].set(xticks=x,xticklabels=labels,ylabel='100 delta ln E',title='(b) Mean projection, mother-bootstrap 95% CI')
    fig.suptitle('Figure 0007-105-1: Strong isolated residual along fixed weak-energy directions')
    p=ART/'figures/fig_0007_105_1_projected_sources.png';fig.savefig(p);plt.close(fig);outputs.append(p)
    fig,axes=plt.subplots(1,3,figsize=(14,4),constrained_layout=True)
    primary=[r for r in summary['prediction'] if r['accessible_before_weak']]
    full=[r for r in summary['prediction'] if not r['accessible_before_weak']]
    x=np.arange(3)
    for offset,values,key,label,color in [(-.25,primary,'energy_vs_mean','Prefix: energy only','tab:blue'),
        (0.,primary,'shape_vs_mean','Prefix: energy + shape','tab:orange'),
        (.25,full,'shape_vs_mean','Full isolated: energy + shape','tab:green')]:
        ratio=np.array([r[key]['ratio'] for r in values]);ci=np.array([r[key]['ci95'] for r in values])
        axes[0].bar(x+offset,ratio,width=.24,color=color,label=label)
        axes[0].errorbar(x+offset,ratio,yerr=np.stack([ratio-ci[:,0],ci[:,1]-ratio]),fmt='none',ecolor='black',capsize=2)
    axes[0].axhline(1,color='black',ls='--',lw=.8)
    axes[0].set(xticks=x,xticklabels=labels,ylabel='MSE / training-mean predictor MSE',title='(a) Validation: 16 strong mothers')
    axes[0].legend(fontsize=7)
    for ax,gap,letter in zip(axes[1:],(8,40),('b','c')):
        key=f'gap{gap}_prefix{gap}'
        target=a[key+'_target'].mean(1);prediction=a[key+'_voltage_energy_plus_shape'].mean(1)
        ax.scatter(target,prediction,s=28)
        limits=[min(target.min(),prediction.min()),max(target.max(),prediction.max())]
        ax.plot(limits,limits,color='black',ls='--',lw=.7)
        ax.set(xlabel='Target / 100 delta ln E',ylabel='Prediction / 100 delta ln E',
               title=f'({letter}) {gap} ns: mother averages for display')
    fig.suptitle('Figure 0007-106-1: Opposite-half prediction; plotted averages are not individual fits')
    p=ART/'figures/fig_0007_106_1_predictability.png';fig.savefig(p);plt.close(fig);outputs.append(p)
    return outputs


def main():
    verify()
    summary=read(EXP/'result-summary.json');mothers=read(EXP/'mothers.json')
    assert read(EXP/'reproductions/complete_replay/reproduction-equality.json')['all_summary_values_equal']
    for p,d in read(EXP/'analysis-freeze.json').items():assert sha(ROOT/p)==d,p
    with np.load(ART/'analysis.npz') as f:a={k:f[k].copy() for k in f.files}
    plots=figures(summary,mothers,a)
    (EXP/'tables').mkdir(exist_ok=True)
    for name,key in [('variance_sources','variance'),('prediction','prediction'),('isolated_energy','isolated_energy')]:
        rows=summary[key]
        with (EXP/'tables'/f'{name}.csv').open('w',newline='') as out:
            writer=csv.DictWriter(out,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader()
            writer.writerows({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in rows)
    projected=[r for r in summary['variance'] if r['group']=='high' and r['kind']=='event_signal' and r['projection'].startswith('weak_')]
    variance_table=table(['弱方向间隔','平均投影及95%区间','母内方差 W','校正后母间方差 B','8次平均的母间方差占比'],[
        [r['projection'].replace('weak_logE_gap','')+' ns',f"{r['mean']:+.2f} [{r['mean_ci95'][0]:.2f}, {r['mean_ci95'][1]:.2f}]",
         f"{r['within_variance']:.2f}",f"{r['between_variance_unclipped']:.2f}",f"{100*r['between_fraction_mean8']:.2f}%"] for r in projected])
    prediction_table=table(['间隔','仅能量/均值 MSE比','能量＋形状/均值 MSE比及95%区间','能量＋形状/仅能量 MSE比','通过'],[
        [f"{r['gap_ns']:g} ns",f"{r['energy_vs_mean']['ratio']:.4f}",
         f"{r['shape_vs_mean']['ratio']:.4f} [{r['shape_vs_mean']['ci95'][0]:.4f}, {r['shape_vs_mean']['ci95'][1]:.4f}]",
         f"{r['shape_vs_energy']['ratio']:.4f}",'是' if r['gate_passed'] else '否']
        for r in summary['prediction'] if r['accessible_before_weak']])
    singles=[r for r in summary['isolated_energy'] if r['role']=='validation' and r['horizon_ns']==1920.]
    single_table=table(['档位','E偏差/%','E RMSE/%','E误差SD/%','E±10%','边界','未收敛'],[
        [r['group'],f"{r['bias_percent']:+.3f}",f"{r['rmse_percent']:.3f}",f"{r['error_sd_percent']:.3f}",
         f"{r['within10_count']}/{r['repeat_fits']}",r['boundary'],r['nonconverged']] for r in singles])
    ranges={g:(min(r['edep_MeV'] for r in mothers if r['group']==g),max(r['edep_MeV'] for r in mothers if r['group']==g)) for g in ('weak','high')}
    counts=read(REC/'features.json')
    total_boundary=sum(r['boundary'] for r in counts['rows']);total_unconverged=sum(r['fits']-r['converged'] for r in counts['rows'])
    result_text=(f"强孤立残差沿8/40 ns弱能量方向的平均投影为{projected[0]['mean']:+.2f}/{projected[1]['mean']:+.2f}（单位100δlnE），"
        f"母事件重采样95%区间均在零以上。按B/(B+W/8)估计，8次数字化平均时母间差异仍占投影方差的"
        f"{100*projected[0]['between_fraction_mean8']:.2f}%/{100*projected[1]['between_fraction_mean8']:.2f}%。"
        "这同时支持存在平均响应偏差和不能靠同母重复平均消除的差异；B还可能包含未去净的能量趋势/几何效应，不能全部归为纯光子统计。")
    prediction_text=("8/40 ns可用前段中，加入四个形状特征相对仅电压能量代理的MSE比为"
        f"{summary['prediction'][0]['shape_vs_energy']['ratio']:.4f}/{summary['prediction'][2]['shape_vs_energy']['ratio']:.4f}，"
        "相对训练均值基线的配对95%区间均跨1；两项主要条件均未通过事先冻结的推进门。"
        "没有证据据此增加逐事件形状自由度或替换默认重建，也没有证明其他特征/非线性模型必然无效。")
    sections=[
        ('103',f'''## 103. 强孤立残差来源与可预测性：输入和冻结协议（2026-11-02）

研究问题是平均模板偏差、光学母事件差异和数字化涨落各贡献多少，以及观测电压能否预测稳定的强孤立残差。沿用第99–102节/L91的文献限定范围，本轮没有新增论文结论或重读共享库。

从旧noise fold1按固定路径哈希选取强、弱各32个母ROOT；每档16训练、16保留评价。64母与旧280均值训练母、空间先验128母和最近四对诊断母不重叠，但64母已用于历史噪声标定，不能称新未见输运确认。源能量为8/1.625 MeV；实际E_dep范围分别为{ranges['high'][0]:.6f}–{ranges['high'][1]:.6f}和{ranges['weak'][0]:.6f}–{ranges['weak'][1]:.6f} MeV。逐母用真实E_dep计算诊断残差，不把源能量当沉积能量。

每母16次独立数字化，共1024次正式运行；先导8次、重放4次另列。完整数字化4096 ns，观测前2048 ns（1024点，2 ns周期、0.5 ns孔径）；第一枚有效面光子对齐128 ns，已知K=1和时刻，固定相位。传感器/读出、95/351 ns均值形状、29 ns恢复代理、14400微胞和12×12空间先验沿用旧版本。此处无Poisson源率扫描，没有新G4、ROOT或实物测量；实际原ROOT的RunInfo、宏和种子另存，不用当前宏推断历史。

事先冻结训练均值、仅电压能量代理、能量代理＋四个形状特征三种预测器。ridge只在16个强训练母内逐母留一选择；评价16母的重复始终一起保留。输入来自8/40/200 ns或完整孤立窗，仅接收完整孔径内观测电压；真实E、祖先分量、光子落点和源能量不作特征。目标是另一半8次数字化的强祖先残差投影，以排除同次噪声相关；目标本身仍有有限重复噪声，只测试稳定母响应，不覆盖全部逐次涨落。

8和40 ns分别要求形状预测/训练均值MSE比<0.90、该比95%区间上界<1、形状预测/仅能量MSE比<0.95，两间隔都通过才考虑另建双脉冲验证。这些是本项目工程门，不是文献保证。完整孤立窗只作单独信息对照，也被压缩为四个残差块，因此不是任意全波形算法的信息上界。
'''),
        ('104',f'''## 104. 母内与母间残差分解

以事件祖先电压减去其E_dep对应的冻结因果孤立模板，另保留完整含暗计数/电子学噪声的电压残差。统计模型为r_mk=b+u_m+ε_mk；W取母内无偏样本方差的母均值，B取母均值间样本方差减W/16。有限重复修正后的负B估计保留，不静默截零。均值区间和投影方差区间按母事件整组重采样4000次；重复不增加独立母数，区间未传播原模型/噪声标定不确定度。

逐母已减去实际E_dep模型，但弱档能量范围仍有差别，剩余的条件均值趋势、固定几何或模型失配都可能进入B；本轮没有唯一分离微观光子统计与全部系统误差。以下图中±1SE只是点态母均值标准误，不是同时置信带。

图 0007-104-1：(a)强、(b)弱档训练/保留母的平均祖先残差和±1SE；(c)(d)未截零的母间方差B与母内W/8，展示对8次平均波形的估计贡献，单位mV²。
'''),
        ('105',f'''## 105. 沿弱能量敏感方向的影响

固定名义E_dep为强7.5 MeV、弱0.25 MeV，使用现有空间模型、FIR及固定GLS度量构造联合logE导数的伪逆。所有母用同一投影算子；强孤立残差位于128 ns，弱事件名义时刻为其后8/40/200 ns。此处真实E只用于诊断残差标签，尚未拟合双脉冲。

{variance_table}

均值单位为100δlnE；方差单位为(100δlnE)²。它们不是相对E百分比、实际非线性双脉冲误差或器件分辨率。最后一列是从16次重复估计B、W后推算8次平均的B/(B+W/8)，不包括平均偏差b²。

{result_text}

图 0007-105-1：(a)8次平均的投影方差来源；(b)平均投影及按32个强母整组重采样的95%区间。完整电压/祖先信号、孤立E方向及未截零方差的逐项结果见实验tables/variance_sources.csv。
'''),
        ('106',f'''## 106. 只用观测电压的可预测性检验

前段能量代理采用已知孤立时刻的单位电压权重拟合，不是默认完整GLS；随后把该前段残差分成四段并按拟合波形RMS归一化。标准化和ridge惩罚选择只使用训练母。所有可用孔径的结束时间不晚于弱脉冲名义到达时刻；单元测试以未来样点NaN干预验证未读取未来。独立特征进程的文件读取守卫禁止读取真实E、祖先信号和光子数据。

{prediction_table}

目标是另一半数字化的强残差弱E投影，并非真正的弱E拟合误差。MSE先在每个评价母的16次预测内取均值，再按16母等权汇总；4000次配对母重采样的区间条件于当前训练模型，不包含重训或标定不确定度。窄能区的能量代理可能利用能区先验，因此保留“仅能量”对照，不能直接推广至未知宽能谱。

{prediction_text}

图 0007-106-1：(a)前段能量/形状和完整孤立窗形状预测相对训练均值的MSE比；(b)(c)8/40 ns按母平均后的目标/预测散点，只为显示母差异，主要MSE仍逐次预测计算。完整窗四段的时间分辨粗于短前段，二者不是嵌套特征集；其结果不能作为可辨识性的普适上界。
'''),
        ('107',f'''## 107. 核验、孤立能量控制与下一步

四项针对性检查通过：已知随机效应的方差恢复及负B保留、对侧半组目标隔离、未来样点不可见、按母分组的回归。64母原truth文件中的E_dep/源能量、ROOT RunInfo中的种子与单事件配置逐项一致。先导8次与正式相应重复、独立目录重放4次与正式数组均逐项相同；残差代数闭合最大误差{summary['exact_voltage_closure_mV']:.3g} mV，固定投影算子作用于名义导数得到单位阵。全部4096次前段拟合保留，未收敛{total_unconverged}，边界{total_boundary}；每前段/档位完整统计另存。另用reproduce.py在新目录完整复算一遍：再做1036次数字化和4096次前段拟合，结果摘要逐值相同。合计实际2072次数字化，统计样本仍只用第一批1024次和64母，不把重放加入分母。

保留母上完整孤立观测窗的单位电压权重能量拟合如下；每档16母×16重复，256次不是256个独立母事件。

{single_table}

这是K=1、已知时刻的孤立控制；未测盲漏检、额外脉冲、误拆分、误合并、重建效率或能谱选择效率，这些指标不能填零。全套数据保留，没有按拟合质量删除母事件或重复。

下一步优先独立检验跨母稳定的条件平均响应修正，保持孤立能量标度和共享微胞响应一致；不要直接用本批逐事件真值残差修正重建。需要在新的开发划分上检查b(E,t)及时间对齐，覆盖多个沉积能量，冻结后才回到同波形、同预算的强弱双脉冲比较。随机母差异继续进入误差/协方差预算，避免与新增低秩项重复计数。此前L54的负结果和本轮形状预测未过门均保留。

本轮只完成来源与可预测性诊断，没有修改物理参数、平均模板、默认算法或开展新的双脉冲能量拟合。实际源码、固定协议、原RunInfo/宏、种子、输入/输出哈希及结果在pileup/experiments/{NAME}；数组/图在jupyter/artifacts/同名目录；日志和修改前Notebook副本在jupyter/run_records/同名目录。沿用L91据既有精读笔记的有限方法依据，不新增文献数值结论。
''')]
    plot_for={'104':plots[0],'105':plots[1],'106':plots[2]}
    newcells=[]
    for number,body in sections:
        newcells.append(nbformat.v4.new_markdown_cell(body,id=f'residual-sources-{number}'))
        if number in plot_for:
            p=plot_for[number].relative_to(ROOT)
            code=("from pathlib import Path\nfrom IPython.display import display, Image\n"
                  "root = Path.cwd()\nif root.name == 'jupyter':\n    root = root.parent\n"
                  f"display(Image(filename=str(root / '{p.as_posix()}')))\n")
            newcells.append(nbformat.v4.new_code_cell(code,id=f'residual-sources-figure-{number}'))
    nbpath=ROOT/'jupyter/codex_0007_SiPM_Energy_Calibration.ipynb'
    current=json.loads(nbpath.read_text())
    before=json.loads((REC/(nbpath.name+'.before')).read_text())
    assert current['cells']==before['cells'] and len(current['cells'])==168
    display_nb=nbformat.v4.new_notebook(cells=newcells)
    NotebookClient(display_nb,timeout=120,kernel_name='python3',resources={'metadata':{'path':str(ROOT)}}).execute()
    current['cells'].extend(display_nb['cells'])
    nbformat.validate(nbformat.from_dict(current))
    nbpath.write_text(json.dumps(current,ensure_ascii=False,indent=1)+'\n')
    short=f'''完成64个历史开发母（强/弱各32）的1024次正式独立数字化；新形状训练/评价各16个强母，数据与均值训练和最近诊断母不重叠，但曾用于历史噪声标定。{result_text}\n\n{prediction_text} 下一步优先做条件平均响应的独立修正控制，覆盖更多沉积能量并保持孤立能标/共享恢复一致，再决定是否进入双脉冲能量比较。当前没有新ROOT、默认替换或盲效率结论。完整记录在codex_0007第103–107节、图0007-104-1至106-1和L92。'''
    for filename,heading in [('Research_Plan_20260921.ipynb','## 33. 残差来源已分解：优先验证稳定均值偏差（2026-11-02）'),
                             ('Run_Catalog.ipynb','## 强孤立残差来源与可预测性：离线重复数字化（2026-11-02）')]:
        p=ROOT/'jupyter'/filename;nb=json.loads(p.read_text())
        prior=json.loads((REC/(filename+'.before')).read_text())
        assert nb['cells']==prior['cells']
        body=heading+'\n\n'+short
        if filename=='Run_Catalog.ipynb':
            body+=f'\n\n实验目录：pileup/experiments/{NAME}。先导8次、正式1024次、独立重放4次；另在reproductions/complete_replay完整重复该流程，合计2072次实际数字化和8192次前段特征拟合，统计只用首批1024次正式波形。两次摘要逐值相同。最多8进程、BLAS单线程。原母ROOT映射、种子、宏与RunInfo分别见mothers.json、original-metadata/、original-runinfo/；新数字化种子735000000+母ID×100+重复ID，读出种子另加10000000。无新ROOT，源注入率不适用于这批受控单脉冲；原运行条件按保存的RunInfo解释。'
        nb['cells'].append(nbformat.v4.new_markdown_cell(body,id='residual-sources-summary-20261102'))
        nbformat.validate(nbformat.from_dict(nb))
        p.write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n')
    docs=f'''\n\n<a id="single-residual-sources-20261102"></a>
## L92 — 2026-11-02：强孤立残差的母间差异与有限形状可预测性

**项目结果。** {result_text}

{variance_table}

平均投影单位100δlnE，方差单位(100δlnE)²；不是相对E百分比或实际双脉冲误差。B由母均值样本方差减W/16估计；母间份额是B/(B+W/8)，不包含平均偏差平方。负方差估计保留在完整表中。

**可预测性。** {prediction_text}

{prediction_table}

16强母训练、16强母保留评价，逐母留一选择ridge；同母重复不跨训练/评价。已知K=1/时刻，只从弱到达前完整孔径的电压拟合能量并提取四个残差块；另一半8重复均值是诊断目标，不共享同次数字化噪声。目标仍有有限重复误差，单位100δlnE，MSE比不是能量恢复效率。完整孤立窗为另一种压缩特征对照，不是可辨识性上界。上述失败只约束本轮特征和预测器，不证明所有形状方法无效。

**数据与核验。** 64母来自旧noise fold1（强/弱各32），与均值训练、空间先验及最近诊断母不重叠，但已用于历史噪声标定，不能称新输运确认。源能量8/1.625 MeV；E_dep分别{ranges['high'][0]:.6f}–{ranges['high'][1]:.6f}/{ranges['weak'][0]:.6f}–{ranges['weak'][1]:.6f} MeV。固定相位、每母16重复，1024正式数字化，先导8与重放4另列；完整4096 ns数字化，观测2048 ns/2 ns。四项针对性测试、原RunInfo/能量标签、重放逐数组相等、残差闭合及投影恒等式通过。4096次前段能量拟合未收敛{total_unconverged}、边界{total_boundary}；全部保留。另在新目录完整复算，摘要逐值相同；两遍合计2072次实际数字化，统计不重复计数。

**采用与下一步。** 沿用L91三篇既有精读及审核边界，没有新增外部论文主张。当前证据优先支持在新的开发划分检验条件均值b(E,t)与时间对齐，覆盖更多沉积能量、保持孤立能标和共享响应，再以同输入/预算验证双脉冲能量。母差异及数字化涨落保留在误差模型中，不与低秩项重复计数。本轮未修改物理/平均响应参数、未新增双脉冲能量拟合或默认算法；盲漏检、误拆并、效率与能谱选择偏差未测。

位置：`pileup/experiments/{NAME}/`及可复用模块`pileup/analysis/alpha_single_residual_sources.py`、`alpha_single_residual_sources_study.py`、`alpha_single_residual_sources_report.py`；本地codex_0007第103–107节、图0007-104-1至106-1，研究计划第33节与Run_Catalog同步。实际原宏/RunInfo、种子、来源哈希及冻结源码随实验记录。
'''
    docpath=ROOT/'docs/literature-use.md'
    old=(REC/'literature-use.md.before').read_text()
    assert docpath.read_text()==old and '## L92 ' not in old
    docpath.write_text(old+docs)
    readme=f'''# 强孤立残差来源与可预测性诊断

本实验完成64个历史光学母事件的1024次重复数字化，区分平均模板误差、母间差异和母内数字化涨落，并在母事件不重叠的保留数据上测试观测电压前段对稳定强孤立残差的预测。已知K=1和时刻，目标用真实E/祖先信号构造；预测特征只读电压。没有新G4/ROOT或盲恢复结果。

{result_text}

{prediction_text}

{prediction_table}

下一步优先验证跨母稳定的条件均值修正。原L54低秩协方差负结果、本轮失败条件和所有边界/未收敛统计均保留。结果仅适用于当前两个窄源能区、固定相位和器件近似；投影单位100δlnE，不能写成实际弱能量误差。

运行须在项目根目录。`run.py`优先导入本目录实际源码副本，`protocol.json`和`design-freeze.json`在新波形产生前冻结；`input-freeze.json`记录旧ROOT、原truth/metadata和派生光子缓存。`original-runinfo/`、`original-metadata/`保留实际历史宏/种子，`mothers.json`列出源能量与E_dep。`pilot-gate.json`为成本门，`result-summary.json`为预定预测门结果。

```bash
python3 -B pileup/experiments/{NAME}/run.py verify
```

可在新目录完整复算数字结果，无需切换Git版本；reproduce.py使用本目录冻结源码、复制配置与输入，在三个同名reproductions/子目录保存输出。已用complete_replay验证整个流程，摘要逐值相同。以下运行须选择未存在的label，适用于已验证的Linux/WSL环境：

```bash
python3 -B pileup/experiments/{NAME}/reproduce.py --label another_replay
```

若本机尚无age.so，先执行build.json列出的编译命令。每次数字化种子735000000+母ID×100+重复ID，读出种子另加10000000。最多8进程、BLAS单线程。正式1024、先导8、局部重放4，另一次完整复算重复上述1036次，合计2072次实际数字化；统计仍仅使用首批1024正式波形。绘图/笔记报告是独立入口run_report.py，追加笔记时检查原内容一致，避免重复追加。

原ROOT/数组及Notebook按仓库约定仅本地保存，克隆仓库后须另行提供数据。图与分析数组位于`jupyter/artifacts/{NAME}/`，日志与修改前笔记在`jupyter/run_records/{NAME}/`。本地codex_0007第103–107节保存直接可见的三幅固定编号图，研究计划与Run_Catalog同步。
'''
    (EXP/'README.md').write_text(readme)
    outputs=[*plots,nbpath,ROOT/'jupyter/Research_Plan_20260921.ipynb',ROOT/'jupyter/Run_Catalog.ipynb',docpath,EXP/'README.md',*sorted((EXP/'tables').glob('*.csv'))]
    save(EXP/'report-output-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in outputs})
    print('Wrote3figures,Notebook sections103–107,plan33,catalog,L92;new display cells executed',flush=True)


if __name__=='__main__':
    main()
