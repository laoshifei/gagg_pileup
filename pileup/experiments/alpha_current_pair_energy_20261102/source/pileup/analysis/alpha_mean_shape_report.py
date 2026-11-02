"""Reproducible figures and append-only notes for the mean-shape experiment."""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/gagg_mean_shape_mpl')
os.environ.setdefault('IPYTHONDIR', '/tmp/gagg_mean_shape_ipython')
os.environ.setdefault('JUPYTER_RUNTIME_DIR', '/tmp/gagg_mean_shape_jupyter')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nbformat
from nbclient import NotebookClient
import numpy as np

from .alpha_mean_shape_study import ROOT, EXP, ART, REC, NAME, read, save, sha, verify


def table(headers, rows):
    return '\n'.join(['| '+' | '.join(headers)+' |', '| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(map(str, r))+' |' for r in rows])


def csv_table(name, rows):
    (EXP/'tables').mkdir(exist_ok=True)
    with (EXP/'tables'/f'{name}.csv').open('w', newline='') as out:
        writer=csv.DictWriter(out, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows({k:json.dumps(v, ensure_ascii=False) if isinstance(v,(dict,list)) else v
                         for k,v in r.items()} for r in rows)


def energy_row(result, group, method, phase=0., kind='mean'):
    return next(r for r in result['energy'] if r['group']==group and r['method']==method
                and r['phase_ns']==phase and r['kind']==kind)


def figures(result, pair):
    plt.rcParams.update({'font.size':10, 'figure.dpi':140, 'axes.grid':True, 'grid.alpha':.2})
    output={}
    with np.load(EXP/'shape-correction.npz') as f:
        grid=f['grid']; profiles=f['profiles']*float(f['scale'])
        knots=np.exp(f['log_charge_knots'])
    fig,axes=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    for i,label in enumerate(('weak','low','mid','high')):
        axes[0].plot(grid,profiles[i]*1e3,label=label)
        axes[1].plot(grid,2*np.cumsum(profiles[i]*knots[i]),label=label)
    axes[0].set(xlim=(-8,600),xlabel='Time after onset / ns',ylabel='Correction / (mV per 1000 q0 pe)',
                title='(a) Training-only, area-projected profiles')
    axes[1].set(xlim=(-16,1620),xlabel='Time after onset / ns',ylabel='Cumulative correction / (mV ns)',
                title='(b) At each charge knot: complete area is zero')
    for ax in axes:ax.axhline(0,color='black',lw=.7);ax.legend(fontsize=8)
    fig.suptitle('Figure 0007-109-1: Frozen isolated mean-shape correction')
    p=ART/'figures/fig_0007_109_1_mean_shape.png';fig.savefig(p);plt.close(fig);output['109']=p
    fig,axes=plt.subplots(2,2,figsize=(12,7),constrained_layout=True)
    phase=np.array([0.,.5,1.,1.5])
    for group,ax in enumerate(axes.flat):
        for method,style in [('old','o-'),('corrected','s-')]:
            rows=[energy_row(result,group,method,p) for p in phase]
            ax.plot(phase,[r['bias_percent'] for r in rows],style,label=method)
        ax.axhline(0,color='black',lw=.7)
        ax.axhspan(-5,5,color='green',alpha=.05)
        ax.axhline(5,color='grey',ls=':',lw=.7);ax.axhline(-5,color='grey',ls=':',lw=.7)
        ax.set(xlabel='ADC phase shift / ns',ylabel='Mean E_dep error / %',
               title=f'({chr(97+group)}) '+('weak','low','mid','high')[group]+': 16 mothers, 8-repeat means')
        ax.legend(fontsize=8)
    fig.suptitle('Figure 0007-110-1: Known-onset isolated energy bias across ADC phases')
    p=ART/'figures/fig_0007_110_1_energy_phases.png';fig.savefig(p);plt.close(fig);output['110']=p
    fig,axes=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    x=np.arange(3)
    for method,offset in [('old',-.08),('corrected',.08)]:
        mean=np.array([r[method+'_mean'] for r in result['projection']])
        ci=np.array([r[method+'_ci95'] for r in result['projection']])
        axes[0].errorbar(x+offset,mean,yerr=np.stack([mean-ci[:,0],ci[:,1]-mean]),
                         fmt='o',capsize=4,label=method)
    axes[0].axhline(0,color='black',lw=.7)
    axes[0].set(xticks=x,xticklabels=['8 ns','40 ns','200 ns'],ylabel='100 delta ln E',
                title='(a) Fixed old weak-energy directions, 95% CI')
    axes[0].legend(fontsize=8)
    ratios=[r['ratio'] for r in result['waveform']]
    bars=axes[1].bar(np.arange(4),ratios,color=['tab:blue' if r<=1.1 else 'tab:red' for r in ratios])
    axes[1].bar_label(bars,fmt='%.3f',padding=3)
    axes[1].axhline(1,color='black',ls='--',lw=.7);axes[1].axhline(1.1,color='tab:red',ls=':',lw=.7)
    axes[1].set(xticks=np.arange(4),xticklabels=['weak','low','mid','high'],ylabel='Corrected / old MSE',
                title='(b) Group-mean waveform residual, phase 0')
    axes[1].set_ylim(0,max(1.3,max(ratios)*1.18))
    fig.suptitle('Figure 0007-111-1: Residual changes on 64 historical validation mothers')
    p=ART/'figures/fig_0007_111_1_residuals.png';fig.savefig(p);plt.close(fig);output['111']=p
    if pair is not None:
        fig,axes=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
        x=np.arange(5);labels=['H-W 8','H-W 40','H-W 200','H-W 1200','W-H 40']
        for method,offset in [('old',-.15),('corrected',.15)]:
            axes[0].bar(x+offset,[r[method]['error_pp'] for r in pair['groups']],width=.3,label=method)
        axes[0].axhspan(-5,5,color='green',alpha=.05)
        axes[0].set(xticks=x,xticklabels=labels,ylabel='Weak area-change error / pp',title='(a) Shared / isolated area')
        axes[0].tick_params(axis='x',labelrotation=20);axes[0].legend(fontsize=8)
        for component,offset in [('strong',-.15),('weak',.15)]:
            axes[1].bar(x+offset,[r['mse_ratios'][component+'_mse_mV2'] for r in pair['groups']],
                         width=.3,label=component)
        axes[1].axhline(1.1,color='tab:red',ls=':',lw=.8)
        axes[1].set(xticks=x,xticklabels=labels,ylabel='Corrected / old component MSE',title='(b) 8-repeat component means')
        axes[1].tick_params(axis='x',labelrotation=20);axes[1].legend(fontsize=8)
        fig.suptitle('Figure 0007-112-1: Truth-assisted shared response, 16 mother pairs')
        p=ART/'figures/fig_0007_112_1_shared.png';fig.savefig(p);plt.close(fig);output['112']=p
    return output


def main():
    verify(True)
    assert read(EXP/'replay-validation.json')['exact_arrays']
    result=read(EXP/'single-result.json'); gate=read(EXP/'single-gate.json')
    pair=read(EXP/'pair-result.json') if (EXP/'pair-result.json').exists() else None
    if not gate['passed']:assert pair is None and not (EXP/'pair-input-freeze.json').exists()
    for name in ('train-input-freeze.json','single-input-freeze.json','fit-output-freeze.json'):
        for p,d in read(EXP/name).items():assert sha(ROOT/p)==d,p
    plots=figures(result,pair)
    for key in ('energy','waveform','projection'):csv_table('single_'+key,result[key])
    csv_table('single_fit_details',read(EXP/'single-fit-details.json'))
    if pair is not None:csv_table('pair_groups',pair['groups']);csv_table('pair_details',pair['details'])
    mothers=read(EXP/'mothers.json'); numerical=read(EXP/'numerical.json')
    calibration=read(EXP/'calibration.json'); counts=read(REC/'fits.json')
    digitizations={s:read(REC/(s+'.json'))['digitizations'] for s in
                   ('pilot','train','single','pair_pilot','pair','replay') if (REC/(s+'.json')).exists()}
    actual_count=sum(digitizations.values())
    total_boundary=sum(r['boundary'] for r in result['energy'])
    total_nonconverged=sum(r['nonconverged'] for r in result['energy'])
    max_phase_span=max(np.ptp([r['bias_percent'] for r in result['energy'] if r['group']==g
                              and r['method']=='corrected' and r['kind']=='mean']) for g in range(4))
    max_rmse_ratio=max(energy_row(result,g,'corrected',p,'individual')['rmse_percent']/
                       energy_row(result,g,'old',p,'individual')['rmse_percent']
                       for g in range(4) for p in (0.,.5,1.,1.5))
    range_table=table(['档位','源能量/MeV','训练 E_dep/MeV（32母）','评价 E_dep/MeV（16母）'],[
        [('弱','低','中','强')[g],(1.625,2.5,5.486,8.)[g],
         *[f"{min(r['edep_MeV'] for r in mothers if r['group']==g and r['role']==role):.6f}–"
           f"{max(r['edep_MeV'] for r in mothers if r['group']==g and r['role']==role):.6f}"
           for role in ('train','validation')]] for g in range(4)])
    energy_table=table(['档位','均值波形E偏差/%：旧→修正','单次E RMSE/%：旧→修正',
                        '单次E误差SD/%：旧→修正','单次±10%：旧→修正'],[
        [('弱','低','中','强')[g],
         '→'.join(f"{energy_row(result,g,m)['bias_percent']:+.3f}" for m in ('old','corrected')),
         '→'.join(f"{energy_row(result,g,m,kind='individual')['rmse_percent']:.3f}" for m in ('old','corrected')),
         '→'.join(f"{energy_row(result,g,m,kind='individual')['error_sd_percent']:.3f}" for m in ('old','corrected')),
         '→'.join(f"{energy_row(result,g,m,kind='individual')['strict10_count']}/128" for m in ('old','corrected'))]
        for g in range(4)])
    projection_table=table(['名义弱间隔','旧平均及95%区间','修正平均及95%区间','绝对均值限值','通过'],[
        [f"{r['gap_ns']:g} ns",f"{r['old_mean']:+.2f} [{r['old_ci95'][0]:.2f}, {r['old_ci95'][1]:.2f}]",
         f"{r['corrected_mean']:+.2f} [{r['corrected_ci95'][0]:.2f}, {r['corrected_ci95'][1]:.2f}]",
         f"{r['absolute_mean_limit']:.2f}",'是' if r['passed'] else '否'] for r in result['projection']])
    wave_text='、'.join(f"{r['ratio']:.3f}" for r in result['waveform'])
    projection_text='；'.join(f"{r['gap_ns']:g} ns由{r['old_mean']:+.2f}变为{r['corrected_mean']:+.2f}"
                             for r in result['projection'][:2])
    verdict=('通过孤立门' if gate['passed'] else '未通过孤立阶段推进门，失败项：'+', '.join(gate['failures']))
    result_text=(f"面积保持的条件均值修正{verdict}。四档孤立能标、相位和RMSE门通过，"
                 f"组平均波形残差MSE比（修正/旧）依弱/低/中/强为{wave_text}，"
                 f"但沿固定旧弱能量方向的平均投影{projection_text}，单位100δlnE。")
    limits=('评价集仅64个已暴露历史母事件，强档16母；重复和四相位不增加独立输运数。'
            '训练128母也曾参与旧均值/空间先验。本轮按母隔离只约束新增修正，不能称全新确认。'
            '均值波形仍含有限数字化涨落。固定95/351 ns、快比例0.3、29 ns恢复代理、14400微胞、'
            '12×12空间先验和2 ns采样；未传播器件、能标、原模型及新均值训练误差。'
            '组均值通过不能保证每母无偏；误差SD不是单能峰实测FWHM。已知K=1和时刻，'
            '未测盲漏检、误拆分、误合并、能谱选择效率或高率性能。')
    next_text=('下一步先检查条件均值在不同光学母事件划分中的稳定性与训练不确定度，'
               '用开发母事件的分组交叉验证同时考察全波形和固定弱能量敏感方向；'
               '冻结后再作未参与调参的母事件确认。当前评价结果不可用于反复选择平滑宽度、修正幅度或方向。'
               '若主要误差来自母间涨落，应增加独立母事件而非只堆同母重复；'
               '通过孤立与共享门后再做同输入、同权重和搜索预算的双脉冲能量比较。')
    if pair is None:
        shared_text=('按冻结协议停在孤立门，没有运行共享/孤立双脉冲数字化、双脉冲能量拟合或更换默认重建。'
                     '预先保存的80个双脉冲任务只是计划，不能计为已运行。'
                     '数值上共享减孤立项保持一致，并不能代替实测共享响应的迁移验证。')
    else:
        shared_text=('孤立门通过后运行80个双脉冲条件、每条件8次共享及两次孤立参考，1920次正式数字化。'
                     f"共享门通过状态：{read(EXP/'pair-gate.json')['passed']}。仅为真实E/祖先分量辅助的前向响应检查，未拟合双脉冲能量。")
    sections=[('108',f'''## 108. 面积保持的条件平均响应：冻结设计（2026-11-02）

研究问题：在保留孤立能量标度和原共享恢复项的前提下，开发母事件得到的平均形状修正能否降低弱脉冲敏感方向的残差？沿用第99–107节/L91–L92的证据；据既有精读笔记，Fixsen等支持考察能量条件均值，但没有保证本项目修正或堆积性能。此处平滑、面积投影、插值及推进门都是项目假设，没有新原文核查或人工审核变更。

训练采用旧均值/空间先验的128母，每档32母、每母8次数字化，共1024次。评价从旧confirmation池按SHA256(736000:ROOT路径)排序，每档取16母；排除最近空间/因果控制32母，与本轮训练无交集。评价每母4相位×8次，共2048次；历史角色和暴露保留，不称新未见输运。

{range_table}

第一枚有效面光子移至128+phase ns，phase为0、0.5、1、1.5 ns。保留真实光子相对时间、位置与能量；数字化4096 ns，观测前2048 ns/1024点，2 ns周期、0.5 ns孔径，中心0.25+2n ns。源能量与沉积能量E_dep分开；真实E只供标定和评分，推断给定K=1/时刻且只读取观测电压与冻结模型。

单事件门：各组/相位均值波形能量偏差绝对值≤5%；单次E RMSE修正/旧≤1.1；各组相位偏差跨度≤1个百分点；相位0组平均残差MSE比≤1.1；8/40 ns弱方向平均投影绝对值≤max(10,0.7×旧绝对值)，单位100δlnE；未收敛为0。200 ns投影为次要结果。上述均为事先冻结的点估计工程门，不是统计显著性或文献阈值。通过才展开共享控制，本轮预先限制最多8进程，BLAS单线程。
'''),('109',f'''## 109. 训练修正与面积/因果控制

每母以8次祖先信号均值减旧孤立模型，再除以冻结q0(真实E_dep)，按能量组等母权平均。使用固定4 ns Gaussian平滑、事件后2 ns保护区与8 ns平滑起始、末端160 ns渐消至约1600 ns，并按旧模板支撑进行reference/(reference+0.01peak)缩减。随后沿正模板×窗口方向扣除面积模，使每个2 ns网格修正的完整离散积分为零。四个电荷结点在log q0上线性插值，结点外夹持，最终乘候选E的q0；没有逐事件形状自由参数，也没有用实现电荷归一化。

时间插值后不截零，以保持面积。训练模型上81能量点×4相位的非负性约束得到全局尺度{calibration['scale']:.3f}，≥0.5前置门通过。事件前修正严格为零；完整积分最大变化{numerical['max_complete_area_change_mV_ns']:.3g} mV·ns。原返回Q(E)、q0(E)、光量映射和shared−isolated项保持，数值差最大{max(r['shared_minus_isolated_difference_mV'] for r in numerical['shared_controls']):.3g} mV；零修正恢复旧模型。两分量logE梯度与中心差分最大相对差{max(numerical['gradient_relative_errors']):.3g}。

面积保持只约束完整积分，不能保证任意截短窗口的积分相同，也不能保证GLS能标不变，因此仍须实际拟合检查。这是经验孤立锚点修正，未修改微观闪烁发光时间分布、恢复参数或新数字化器。

图 0007-109-1：(a)四档训练得到的单位q0形状修正；(b)各电荷结点上的累积修正面积，完整积分回到零。图中只显示开发数据构造的模型，不是评价性能。
'''),('110',f'''## 110. 孤立能量与采样相位评价

每个母/相位分别拟合8条完整含暗信号/电子学噪声的电压和它们的均值；旧/修正使用完全相同的电压导出初值、固定旧模型GLS度量及标量有界求解器。256个条件×9波形×2方法，共{counts['fits']}次正式拟合；成本先导另{counts['pilot_fits']}次。拟合工作进程禁止读取ROOT、母标签、光子缓存和祖先真值。真实E只在后续评分时读取，已知K/t仍属真值帮助。

以下为相位0；每组16母，单次128波形，完整四相位的偏差、RMSE、误差SD、±10%计数、边界和未收敛保存在tables/single_energy.csv。均值列先逐母拟合8次均值，再平均相对误差，不是把16母拼成一个脉冲。

{energy_table}

修正后最大相位偏差跨度{max_phase_span:.5f}个百分点，所有组/相位单次RMSE比最大{max_rmse_ratio:.4f}；均通过冻结门，但单次RMSE没有因此一致改善。正式拟合边界共{total_boundary}、未收敛{total_nonconverged}，全部行保留。±10%是已知K/t孤立E条件计数，不是盲检出效率。

图 0007-110-1：(a)弱、(b)低、(c)中、(d)强档均值波形能量偏差随相位的变化，背景区域表示±5%工程门；同母不同相位并非独立样本。
'''),('111',f'''## 111. 全波形改善不保证弱方向改善

{result_text}

{projection_table}

投影使用L92冻结的旧模型名义强7.5 MeV/弱0.25 MeV联合logE影响算子；每个评价强母在真实E处计算孤立祖先残差，再投影到同一个算子，单位100δlnE。它不是相对E百分比、双脉冲能量RMSE或恢复效率。旧/修正没有各自换方向以追逐改善。95%区间来自16个强母整组4000次配对重采样，未重训修正或传播原标定误差；修正变化区间很窄只是说明同一固定修正对相近能量施加了接近的偏移。

L92另32个noise-fold母的旧8 ns均值为+41.65，而本批16个confirmation母为−4.83且区间跨零；两批不构成配对比较，不能将此差异直接归为某个物理机制或宣称总体均值变号。有限母样本与数字化误差仍很大；本轮不支持把上一批的平均方向当作跨母池普遍可迁移的修正。全波形MSE是均匀时间权重，弱能量方向是另一函数，二者并不等价。

图 0007-111-1：(a)旧/修正8、40、200 ns方向均值与母重采样95%区间；(b)各能区相位0组平均残差的全波形MSE比，1.1为预设上限。组平均MSE不是每母残差MSE平均，不能代表全部随机涨落。
'''),('112',f'''## 112. 推进门与停止位置

{shared_text}

仅失败8 ns方向门就停止，保留40/200 ns和孤立能标通过的结果，不从条件中挑选好结果宣布候选有效。当前门基于点估计；8 ns修正平均的区间跨零，因此失败不能解释为已经统计证明所有条件均值修正有害。

{limits}
'''),('113',f'''## 113. 复现记录与下一步

原ROOT的SHA256、实际truth E_dep、传感器/读出配置、RunInfo种子逐母核对，实际原宏随original-metadata/与original-runinfo/保留。四项针对性测试覆盖面积/分数相位、零修正、等母归一化与非法面积拒绝；数值面积、事件前零值、共享项不变和梯度检查完成。先导与正式相同输入的波形、拟合结果（除计时）一致；重放训练母1与评价条件197的前两次数字化，所有电压/祖先/电荷数组逐项相同。

正式训练1024、孤立评价2048；其他数字化逐项为{json.dumps({k:v for k,v in digitizations.items() if k not in ('train','single')},ensure_ascii=False)}，合计实际{actual_count}次。正式拟合{counts['fits']}、先导{counts['pilot_fits']}次；重放不增加统计样本。代码与协议在新波形前冻结，形状参数在读取评价结果前冻结。没有新G4或ROOT。

{next_text}

实际源码、宏/RunInfo、输入/输出哈希、种子与门结果在pileup/experiments/{NAME}；图和数组在jupyter/artifacts/{NAME}，运行日志与修改前笔记副本在jupyter/run_records/{NAME}。沿用L91既有精读和审核范围；新修正本身没有外部验证依据。旧笔记计算与图像全部保留，本阶段不自动上传。
''')]
    if pair is not None:
        sections[4]=(sections[4][0],sections[4][1]+'\n图 0007-112-1：(a)共享弱面积变化预测误差；(b)强/弱分量均值残差MSE比。\n')
    newcells=[]
    for number,body in sections:
        newcells.append(nbformat.v4.new_markdown_cell(body,id=f'mean-shape-{number}'))
        if number in plots:
            path=plots[number].relative_to(ROOT)
            code=("from pathlib import Path\nfrom IPython.display import display, Image\n"
                  "root = Path.cwd()\nif root.name == 'jupyter':\n    root = root.parent\n"
                  f"display(Image(filename=str(root / '{path.as_posix()}')))\n")
            newcells.append(nbformat.v4.new_code_cell(code,id=f'mean-shape-figure-{number}'))
    nbpath=ROOT/'jupyter/codex_0007_SiPM_Energy_Calibration.ipynb'
    old_notebooks={}
    for filename in (nbpath.name,'Research_Plan_20260921.ipynb','Run_Catalog.ipynb'):
        current=read(ROOT/'jupyter'/filename);before=read(REC/(filename+'.before'))
        assert current==before,filename
        old_notebooks[filename]=current
    assert len(old_notebooks[nbpath.name]['cells'])==176
    display_nb=nbformat.v4.new_notebook(cells=newcells)
    NotebookClient(display_nb,timeout=120,kernel_name='python3',resources={'metadata':{'path':str(ROOT)}}).execute()
    current=old_notebooks[nbpath.name];current['cells'].extend(display_nb['cells'])
    nbformat.validate(nbformat.from_dict(current));nbpath.write_text(json.dumps(current,ensure_ascii=False,indent=1)+'\n')
    short=result_text+'\n\n'+shared_text+'\n\n'+next_text
    for filename,heading in [('Research_Plan_20260921.ipynb','## 34. 条件均值修正未通过弱方向迁移检查（2026-11-02）'),
                             ('Run_Catalog.ipynb','## 面积保持的平均形状修正：离线孤立检查（2026-11-02）')]:
        body=heading+'\n\n'+short+'\n\n详见codex_0007第108–113节、图0007-109-1至111-1及L93。'
        if filename=='Run_Catalog.ipynb':
            body+=(f'\n\n实验pileup/experiments/{NAME}，本轮无新ROOT。128训练母＋64评价母，训练1024、评价2048、先导4、重放4，共{actual_count}次数字化；'
                   f'{counts["fits"]}次正式能量拟合＋{counts["pilot_fits"]}次先导。训练种子736000000+母ID×100+重复ID，'
                   '评价736100000+母ID×100+重复ID，读出另加10000000；评价相位复用种子但改变随机数消耗，并非相同雪崩。'
                   '计划双脉冲任务未运行，不能计入完成数。mothers.json保存旧ROOT映射，original-metadata/和original-runinfo/保存实际宏、原输运种子和条件；'
                   'train/single-plan.json保存实际新数字化任务。4096 ns数字化、前2048 ns观测，最多8进程、BLAS单线程；没有Poisson源率扫描。')
        nb=old_notebooks[filename];nb['cells'].append(nbformat.v4.new_markdown_cell(body,id='mean-shape-summary-20261102'))
        nbformat.validate(nbformat.from_dict(nb));(ROOT/'jupyter'/filename).write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n')
    write_documents(result_text,energy_table,projection_table,shared_text,limits,next_text,actual_count,
                    counts,total_boundary,total_nonconverged,numerical)
    outputs=[*plots.values(),nbpath,ROOT/'jupyter/Research_Plan_20260921.ipynb',ROOT/'jupyter/Run_Catalog.ipynb',
             ROOT/'docs/literature-use.md',ROOT/'README.md',EXP/'README.md',*sorted((EXP/'tables').glob('*.csv'))]
    save(EXP/'report-output-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in outputs})
    print('Wrote figures, sections108–113, plan34, Run_Catalog and L93;new display cells executed',flush=True)


def write_documents(result_text,energy_table,projection_table,shared_text,limits,next_text,
                    actual_count,counts,boundaries,nonconverged,numerical):
    docs=f'''\n\n<a id="mean-shape-20261102"></a>
## L93 — 2026-11-02：面积保持的均值修正未通过8 ns弱方向迁移门

**项目结果。** {result_text}

{projection_table}

投影单位100δlnE，使用L92固定旧模型名义强7.5/弱0.25 MeV算子，真实E与祖先信号仅用于残差诊断；不是双脉冲能量误差。区间为16个强母的4000次重采样，条件于冻结训练模型，未传播训练误差。8/40 ns门为绝对均值≤max(10,0.7×旧绝对值)，200 ns为次要结果。8 ns区间仍跨零，工程门失败不等于证明全部均值修正有害。

**采用及条件。** 128个旧均值/空间先验开发母（四档各32）各8次数字化训练修正；64个与训练不重叠的旧confirmation母（各16）按固定路径哈希选择，排除最近32个空间/因果控制母，4相位×8次评价。历史暴露完整保留。源能量1.625/2.5/5.486/8 MeV，评价E_dep范围分别0.208451–0.285529、1.389439–1.512326、4.917246–4.977028、7.565534–7.604877 MeV。

修正按真实E_dep的q0归一化、等母平均；固定4 ns平滑、2 ns因果保护、8 ns起始坡、约1600 ns支撑与末端160 ns渐消，扣除完整面积模后在log q0结点插值。全局非负性尺度1，未用评价调节。q0(E)、原返回Q(E)、光量逆映射与共享减孤立项保持，微观发光/恢复参数未改。4 ns等数值及门限是项目假设，没有外部性能保证。

**孤立E控制。** 相位0结果如下，旧/修正使用相同电压导出的初值与固定GLS权重；每组128条单次波形只有16个母事件，另有16条8次均值。

{energy_table}

完整四相位偏差、RMSE、误差SD、±10%计数、边界/未收敛均保留。正式{counts['fits']}次拟合，边界{boundaries}、未收敛{nonconverged}；工作进程禁止读真实E/祖先/光子。已知K/t，条件±10%不能作盲效率。{shared_text}

**验证与限制。** 原ROOT/RunInfo/实际宏、标签、种子与输入哈希核对；4项面积/相位/归一化测试、数值因果/共享不变及梯度检查通过。完整积分最大变化{numerical['max_complete_area_change_mV_ns']:.3g} mV·ns，事件前修正为零。正式3072次数字化＋先导4＋局部重放4，共{actual_count}次；先导/重放数组精确一致。另{counts['pilot_fits']}次先导拟合与正式除计时外相同。{limits}

**证据与下一步。** 沿用L91据既有精读笔记的Fixsen条件均值思想及Busch/Zhang适用范围，没有新原文核查、论文数值主张或审核变更。L92旧32强母的8 ns平均投影+41.65与本批−4.83并非配对估计；本批区间跨零，尚不能分离有限母样本、数字化误差与池间差异。全波形MSE与弱方向敏感度不同，本轮不支持直接迁移这一固定修正。{next_text}

位置：`pileup/experiments/{NAME}/`；可复用模块`pileup/analysis/alpha_mean_shape_correction.py`、`alpha_mean_shape_study.py`、`alpha_mean_shape_report.py`。本地codex_0007第108–113节、图0007-109-1至111-1，研究计划第34节与Run_Catalog同步。原计算/图像保留，新阶段不自动上传。
'''
    p=ROOT/'docs/literature-use.md';old=(REC/'literature-use.md.before').read_text()
    assert p.read_text()==old and '## L93 ' not in old;p.write_text(old+docs)
    note=f'''# 面积保持的条件均值修正

{result_text}

{projection_table}

{shared_text}

训练128个历史母事件，评价64个不同历史母事件，四个源能区、四个采样相位。训练1024、评价2048、先导4、重放4，共{actual_count}次离线数字化；正式{counts['fits']}、先导{counts['pilot_fits']}次孤立能量拟合。没有新ROOT。源能量与真实E_dep在mothers.json中分开记录，原ROOT的实际宏/种子保存在original-metadata/、original-runinfo/。

q0(E)、返回Q(E)与共享减孤立项保持；新形状为训练残差的固定面积投影和log q0插值，不是逐事件自由形状参数。真实E只用于训练/评分，拟合只读电压和冻结模型，K=1/时刻已知。全波形残差MSE降低仍不保证弱方向改善。投影单位100δlnE，不是实际双脉冲能量误差。

{limits}

{next_text}

运行目录为项目根目录。run.py使用本目录实际源码副本，不要求切换Git版本；protocol.json、design-freeze.json与input-freeze.json记录波形产生前的设计/输入，model-freeze.json记录评价前固定的修正。single-gate.json保留唯一失败项，pair-plan.json为未执行计划。数值数组、模型npz、ROOT和Notebook按项目约定只在本地保存，克隆仓库后须另行提供这些输入。

```bash
python3 -B pileup/experiments/{NAME}/run.py verify
python3 -B pileup/experiments/{NAME}/validate.py
```

如需完整重跑，reproduce.py在新目录镜像冻结输入并写独立结果，最多8进程；原结果不会被覆盖。新目录必须不存在。此轮仅实际执行了预先计划的4次局部重放，完整复算入口不表示本轮已额外完整运行：

```bash
python3 -B pileup/experiments/{NAME}/reproduce.py --label another_replay
```

新数字化种子：训练736000000+母ID×100+重复ID；评价736100000+母ID×100+重复ID，读出另加10000000。各相位复用种子并不保证相同雪崩。运行条件与环境见protocol.json、environment.json、build.json；需重建age.so时使用build.json的编译命令。report-source-freeze.json单独固定报告源码；run_report.py只执行新显示单元并追加笔记，检查旧笔记完全一致，禁止重复追加。

图与数组在jupyter/artifacts/{NAME}，运行记录/修改前副本在jupyter/run_records/{NAME}。本地codex_0007第108–113节含三幅固定编号图，研究计划第34节、Run_Catalog及docs/literature-use.md的L93同步。
'''
    (EXP/'README.md').write_text(note)
    p=ROOT/'README.md';old=(REC/'README.md.before').read_text();assert p.read_text()==old
    p.write_text(old+f'''\n\n2026-11-02完成[面积保持的条件均值修正验证](pileup/experiments/{NAME}/README.md)：128个历史母训练、64个不同历史母评价，共3072次正式数字化和4608次已知K/t的孤立能量拟合。四档能标/相位门通过，四档组平均波形MSE均降低，但8 ns弱方向平均投影由−4.83变为−29.05（100δlnE），未通过预定推进门，因而停止双脉冲扩展。该投影不是实际弱能量误差；评价只有16个强母且区间跨零，历史池差异与训练不确定度仍待检查。面积、因果、共享项不变及精确重放检查通过，没有新ROOT或默认替换。下一步先验证条件均值跨母划分的稳定性，避免用已暴露评价集调修正。详情见[采用记录L93](docs/literature-use.md#mean-shape-20261102)，本地codex_0007第108–113节与三幅固定编号图。\n''')


if __name__=='__main__':main()
