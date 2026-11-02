"""Append the independent four-energy/phase and conditional shared controls."""
from __future__ import annotations
import csv
import json
import os

os.environ.setdefault('MPLCONFIGDIR','/tmp/gagg_current_shared_mpl')
os.environ.setdefault('IPYTHONDIR','/tmp/gagg_current_shared_ipython')
os.environ.setdefault('JUPYTER_RUNTIME_DIR','/tmp/gagg_current_shared_jupyter')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nbformat
from nbclient import NotebookClient
import numpy as np

from . import alpha_current_shared_study as study
from .alpha_current_shared import ENERGIES,PHASES,CONDITIONS
ROOT,EXP,ART,REC=study.ROOT,study.EXP,study.ART,study.REC
read,save,sha=study.read,study.save,study.sha


def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(map(str,row))+' |' for row in rows])


def csv_table(name,rows):
    (EXP/'tables').mkdir(exist_ok=True)
    with (EXP/'tables'/f'{name}.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader()
        writer.writerows({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in row.items()} for row in rows)


def figures(single,pair):
    plt.rcParams.update({'font.size':10,'figure.dpi':140,'axes.grid':True,'grid.alpha':.2})
    paths={}
    fig,axes=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    for group,energy in enumerate(ENERGIES):
        color=f'C{group}'
        for method,style in [('old','--'),('corrected','-')]:
            for axis,kind,key in [(axes[0],'mean','bias_percent'),(axes[1],'individual','rmse_percent')]:
                rows=[r for r in single['energy'] if r['group']==group and r['method']==method and r['kind']==kind]
                axis.plot([r['phase_ns'] for r in rows],[r[key] for r in rows],style+'o',color=color,
                          label=f'{energy:g} MeV / {method}')
    axes[0].axhline(0,color='grey',lw=.7)
    axes[0].set(xlabel='ADC phase / ns',ylabel='Mean-waveform E bias / %',title='(a) 32 fresh mothers per source energy')
    axes[1].set(xlabel='ADC phase / ns',ylabel='Individual E RMSE / %',title='(b) 8 repeats per mother; known K=1 / time')
    axes[1].legend(fontsize=7,ncol=2,loc='upper center',bbox_to_anchor=(.5,.92))
    fig.suptitle('Figure 0007-131-1: Frozen isolated energy fits on independent mothers')
    p=ART/'figures/fig_0007_131_1_energy.png';fig.savefig(p);plt.close(fig);paths['131']=p
    with np.load(ART/'single-summary.npz') as data:
        t=data['time_ns'];fig,axes=plt.subplots(2,2,figsize=(12,7),constrained_layout=True)
        for col,group in enumerate((0,3)):
            actual=data[f'g{group}_p0.0_actual_mV'].mean(0)
            axes[0,col].plot(t,actual,color='black',label='Ancestor mean')
            for method,color in [('old','tab:blue'),('corrected','tab:orange')]:
                residual=data[f'g{group}_p0.0_{method}_residual_mV'].mean(0)
                axes[0,col].plot(t,actual-residual,color=color,lw=1,label=method+' at true E')
                axes[1,col].plot(t,residual,color=color,lw=1,label=method)
            axes[0,col].set(title=f'({chr(97+col)}) Source {ENERGIES[group]:g} MeV; phase 0',ylabel='Signal / mV')
            axes[1,col].set(title=f'({chr(99+col)}) Observed - predicted group mean',ylabel='Residual / mV')
            for axis in axes[:,col]:axis.set(xlim=(0,2048),xlabel='Aligned time / ns');axis.legend(fontsize=8)
        fig.suptitle('Figure 0007-132-1: Complete observed window; independent isolated response')
        p=ART/'figures/fig_0007_132_1_waveforms.png';fig.savefig(p);plt.close(fig);paths['132']=p
    if pair is None:return paths
    fig,axes=plt.subplots(1,2,figsize=(13,4.5),constrained_layout=True)
    for j,(order,gap) in enumerate(CONDITIONS):
        rows=[r for r in pair['groups'] if r['condition']==order and r['gap_ns']==gap]
        label=f'{order} {gap:g} ns'
        for method,style in [('old','--'),('corrected','-')]:
            axes[0].plot(PHASES,[r[method]['weak_area_error_pp'] for r in rows],style+'o',color=f'C{j}',
                         label=label if method=='corrected' else None)
        axes[1].plot(PHASES,[r['mse_ratios']['weak_mse_mV2'] for r in rows],'-o',color=f'C{j}')
        axes[1].plot(PHASES,[r['mse_ratios']['strong_mse_mV2'] for r in rows],'--s',color=f'C{j}')
    for y in (-5,5):axes[0].axhline(y,color='grey',ls=':',lw=1)
    axes[0].axhline(0,color='grey',lw=.7);axes[0].legend(fontsize=7)
    axes[0].set(xlabel='ADC phase / ns',ylabel='Predicted - actual weak window-area change / pp',
                title='(a) Solid: correction; dashed: old model')
    axes[1].axhline(1.1,color='grey',ls=':',lw=1);axes[1].axhline(1,color='grey',lw=.7)
    axes[1].set(xlabel='ADC phase / ns',ylabel='Corrected / old component residual MSE',
                title='(b) Solid: weak; dashed: strong component')
    fig.suptitle('Figure 0007-133-1: Shared forward response; 32 independent mother pairs')
    p=ART/'figures/fig_0007_133_1_shared_gates.png';fig.savefig(p);plt.close(fig);paths['133']=p
    with np.load(ART/'pair-examples.npz') as data:
        fig,axes=plt.subplots(2,2,figsize=(12,7),constrained_layout=True)
        for col,i in enumerate((0,4)):
            for row,component in enumerate((0,1)):
                axis=axes[row,col]
                for key,color,style in [('actual','black','-'),('isolated','grey',':'),('old','tab:blue','--'),('corrected','tab:orange','-')]:
                    axis.plot(t,data[f'case{i}_{key}_mV'][component],color=color,ls=style,lw=1,label=key)
                axis.set(xlim=(0,2048),xlabel='Aligned time / ns',ylabel='Component / mV',
                    title=f'({chr(97+row*2+col)}) {[8,40][col]} ns; '+['strong','weak'][row])
                axis.legend(fontsize=8)
        fig.suptitle('Figure 0007-134-1: First frozen mother pair; phase 0; true E/time forward inputs')
        p=ART/'figures/fig_0007_134_1_components.png';fig.savefig(p);plt.close(fig);paths['134']=p
    return paths


def main():
    study.verify();replay=read(EXP/'replay.json');assert replay['passed']
    single=read(EXP/'single-result.json');sg=read(EXP/'single-gate.json')
    pair=read(EXP/'pair-result.json') if (EXP/'pair-result.json').exists() else None
    pg=read(EXP/'pair-gate.json') if pair else None
    mothers=read(EXP/'mothers.json');paths=figures(single,pair)
    for key in ('energy','waveform'):csv_table('single_'+key,single[key])
    inventory=[{k:r[k] for k in ('id','role','group','root','macro','seed','sipm_seed','readout_seed',
                               'source_energy_MeV','edep_MeV','root_sha256','log')} for r in mothers]
    csv_table('root_catalog',inventory)
    def e(group,phase,method,kind):
        return next(r for r in single['energy'] if r['group']==group and r['phase_ns']==phase and
                    r['method']==method and r['kind']==kind)
    energy_table=table(['源能量/MeV','相位/ns','均值E偏差/%：旧→修正','单次E RMSE/%：旧→修正'],[
        [energy,phase,f"{e(g,phase,'old','mean')['bias_percent']:+.4f}→{e(g,phase,'corrected','mean')['bias_percent']:+.4f}",
         f"{e(g,phase,'old','individual')['rmse_percent']:.4f}→{e(g,phase,'corrected','individual')['rmse_percent']:.4f}"]
        for g,energy in enumerate(ENERGIES) for phase in PHASES])
    wave_table=table(['源能量/MeV','相位/ns','组均值残差MSE/mV²：旧→修正','修正/旧','逐母残差MSE比'],[
        [ENERGIES[r['group']],r['phase_ns'],f"{r['old_group_mean_mse_mV2']:.5f}→{r['corrected_group_mean_mse_mV2']:.5f}",
         f"{r['group_mean_mse_ratio']:.4f}",f"{r['corrected_per_mother_mse_mV2']/r['old_per_mother_mse_mV2']:.4f}"]
         for r in single['waveform']])
    ranges=table(['源能量/MeV','正式母数','实际沉积能量范围/MeV'],[
        [energy,32,f"{min(r['edep_MeV'] for r in mothers if r['role']=='formal' and r['group']==g):.6f}–"
                   f"{max(r['edep_MeV'] for r in mothers if r['role']=='formal' and r['group']==g):.6f}"]
        for g,energy in enumerate(ENERGIES)])
    spans=[float(np.ptp([e(g,p,'corrected','mean')['bias_percent'] for p in PHASES])) for g in range(4)]
    failures=sg['failures']+(pg['failures'] if pg else [])
    passed=sg['passed'] and bool(pg and pg['passed'])
    key=('新增136个当前环境单α ROOT（8先导、128正式，四档各32正式母），固定L93模型。'
         f"四档/四相位孤立门{'通过' if sg['passed'] else '未通过'}；"+
         (f"共享微胞前向门{'通过' if pg['passed'] else '未通过'}。" if pg else '按预设规则未启动共享双脉冲验证。'))
    next_text=('下一步冻结已知K=2/时刻的成对能量拟合比较，复用本轮未改动的共享电压、相同权重和搜索预算；'
               '先区分弱脉冲能量偏差、分辨率和条件能量成功率，仍不替换默认模型或启动盲高率重建。' if passed else
               '停止进入双脉冲能量拟合与高率推广；保留全部失败条件。下一步在另建开发设计中定位模型失配，'
               '区分组均值偏差与母间涨落。若调整均值修正、权重或门限，须另留新的确认样本；本批不得再作为独立确认。')
    total_digits=4096+256+replay['digitizations']+(19200+36 if pair else 0)
    counts=(f"原生136次SiPM数字化单列；离线孤立正式4096、先导256；"
            +(f"共享恢复开启15360、关闭3840、共享先导36；" if pair else '共享正式、关闭恢复和共享先导均为0；')
            +f"重放{replay['digitizations']}，合计{total_digits}次离线随机数字化。"
            '孤立能量拟合正式9216、先导72、重放36，合计9324；双脉冲能量拟合0。')
    pilot=read(EXP/'pilot-gate.json');audits=pilot['readout_audits']
    validation=(f"8个原生先导完成电荷/雪崩/读出核对，同雪崩Python信号最大差"
                f"{max(r['python_signal_max_error_mV'] for r in audits):.3g} mV，内部步长减半的峰值相对差最大"
                f"{100*max(r['half_dt_max_error_over_peak'] for r in audits):.6f}%。"
                f"局部精确重放{replay['exact_arrays']}个数组，36次孤立拟合除计时外一致。"
                '原39份物理源码与已验证二进制哈希未变；此次新增分析逻辑，不重新编译G4。')
    single_method=('128个正式母每母四相位0/0.5/1/1.5 ns、每相位8次数字化。按母聚类，512个条件不是512个独立输运样本。'
        '所有面光子保留，首有效光子对齐128+相位 ns；离线4096 ns，观测前2048 ns/1024点，2 ns采样、0.5 ns孔径。'
        '离线种子740300000+母ID×100+重复ID，读出加10000000；跨相位共用种子不保证相同随机雪崩。'
        '孤立拟合已知K=1/时刻，只读完整含噪电压构造相同初值/GLS权重/预算；进程禁止读取真实E、ROOT及祖先信号。'
        '真实E_dep仅用于评分；前向波形诊断则明确使用真实E/t。')
    pair_table='共享阶段未执行；预先冻结的640个开启恢复条件和160个关闭恢复条件不能计入实际结果。'
    pair_text=pair_table;off_text=pair_table
    if pair:
        csv_table('shared_groups',pair['groups']);csv_table('recovery_off',pair['recovery_off'])
        csv_table('shared_details',read(EXP/'pair-details.json')['on'])
        pair_table=table(['顺序/间隔ns（相位0）','实际弱面积变化/%','预测误差/百分点：旧→修正','修正误差95%区间','弱/强MSE比','弱MSE比95%区间','冻结门'],[
            [f"{r['condition']}/{r['gap_ns']:g}",f"{r['corrected']['actual_weak_area_change_percent']:+.3f}",
             f"{r['old']['weak_area_error_pp']:+.3f}→{r['corrected']['weak_area_error_pp']:+.3f}",
             f"[{r['corrected_area_error_ci95'][0]:.3f}, {r['corrected_area_error_ci95'][1]:.3f}]",
             f"{r['mse_ratios']['weak_mse_mV2']:.4f}/{r['mse_ratios']['strong_mse_mV2']:.4f}",
             f"[{r['mse_ratio_ci95']['weak_mse_mV2'][0]:.4f}, {r['mse_ratio_ci95']['weak_mse_mV2'][1]:.4f}]",
             '通过' if r['passed'] else '未通过'] for r in pair['groups'] if r['phase_ns']==0])
        weak_ratios=[r['mse_ratios']['weak_mse_mV2'] for r in pair['groups']]
        strong_ratios=[r['mse_ratios']['strong_mse_mV2'] for r in pair['groups']]
        pair_table+=(f"\n\n全部20条件中，修正弱面积平均误差最大绝对值"
            f"{max(abs(r['corrected']['weak_area_error_pp']) for r in pair['groups']):.4f}个百分点；"
            f"弱分量MSE比范围{min(weak_ratios):.5f}–{max(weak_ratios):.5f}，其中{sum(v>1 for v in weak_ratios)}/20略高于1；"
            f"强分量范围{min(strong_ratios):.5f}–{max(strong_ratios):.5f}。"
            '通过支持冻结容差内的响应一致性；弱分量误差大多略增，尚无双脉冲能量改善结论。')
        pair_text=('正式弱源/强源各32母一一配对，32对母之间无共享父母；每对覆盖强后弱8/40/200/1200 ns及弱后强40 ns、四相位。'
            '640个条件仍只有32对独立母，每条件8重复的共享＋两个孤立对照。源能量1.625/8 MeV；模型输入使用逐母实际沉积能量。'
            '种子740500000+母对ID×1000+重复ID，读出加10000000；条件间共用种子不等于同一雪崩。'
            '全部正式预测在共享数字化前保存、冻结。相同孤立Q(E)使修正模型共享项可复用旧项，6个先导与直接计算逐点核对小于1e-10 mV。\n\n'
            +pair_table+'\n\n弱面积是前2048 ns观测窗内祖先分量电压积分相对其孤立对照的变化，不是完整电荷。'
            '1200 ns条件的弱尾部部分落在窗外，不能与完整积分混用。波形MSE先对每对母8次祖先均值算残差平方，再平均32对；'
            '弱/强分别判定≤1.1，修正弱面积平均误差绝对值≤5个百分点。完整20条件和10000次整母对配对重采样区间在shared_groups表中，'
            '区间为逐终点描述性95%，条件于冻结模型，不传播模型训练误差、不宣称同时覆盖；门用冻结点估计。')
        off_text=table(['顺序/间隔ns','恢复关闭弱面积变化/%','整母对95%区间','冻结门'],[
            [f"{r['condition']}/{r['gap_ns']:g}",f"{r['mean_weak_area_change_percent']:+.4f}",
             f"[{r['ci95'][0]:.4f}, {r['ci95'][1]:.4f}]",'通过' if r['passed'] else '未通过'] for r in pair['recovery_off']])
        off_text+=('\n\n恢复关闭仅将离线SensorConfig.recovery_enabled从True改为False，其他设置保持；全部32对、五条件、相位0，'
                   '以弱面积均值绝对变化≤3个百分点检查控制。此变体不修改物理模型或用于拟合修正。'
                   f"全部预测分量最小值{pair['minimum_model_mV']:.6g} mV，作为形状诊断保留。")
    sections=[('130',f'''## 130. 当前环境四档与四相位独立检查（2026-11-02）

{key}

L96强源独立迁移检验通过后，冻结L93全部模型/权重/门限，新增8个成本先导及128个正式母。没有重训、平滑选择或按效果挑选样本。模型开发环境与当前Geant4 11.3.2/ROOT6.34.10不同，沿用跨环境检验表述。

{ranges}

一个20 μs ROOT仅含一个5 μs中央正入射α，无背景；不对应泊松注入率。G4内部串行、最多8独立任务、启动≥1.05 s。G4种子740000000+ID，原生雪崩/读出加100000/200000；ID0–7先导，ID8–135正式，逐能区32母。模型及39份物理源码与L96相同，95/351 ns光学常数、快比例0.3、29 ns恢复代理和14400微胞均未调整。实际代码基于aa0ec00及冻结时新驱动，完整源码、配置、宏、RunInfo及哈希保存在新快照。

{single_method}
'''),('131',f'''## 131. 独立孤立能量检查

{energy_table}

各能区相位偏差跨度依次为{', '.join(f'{v:.6f}' for v in spans)}个百分点。门要求旧/修正均值E偏差绝对值≤5%，修正/旧单次E RMSE≤1.1，修正相位偏差跨度≤1个百分点，未收敛0；并须同时通过下一节的波形门。完整能量表还保留误差SD、±10%计数、边界和未收敛数。误差SD不当作实测峰FWHM，±10%计数不是盲效率；重复观测仍按母聚类。

图 0007-131-1：(a)四档各32母的均值波形拟合能量偏差；(b)单次拟合RMSE。虚线旧模型，实线冻结修正；横轴为四个采样相位。
'''),('132',f'''## 132. 独立孤立前向波形检查

{wave_table}

组均值残差MSE先对32母的8重复均值残差平均，再在1024采样点平方平均；冻结比值门≤1.1。逐母残差MSE另列，两者不能互换。孤立门整体{'通过' if sg['passed'] else '未通过'}；失败项：{'；'.join(sg['failures']) or '无'}。全部条件原样保留，不用观察到的结果放宽门限。

图 0007-132-1：(a)(b)弱/强源相位0的完整观测窗平均信号和真实E处预测；(c)(d)对应组均值残差。32母、每母8重复，不是逐事件拟合后的残差。
'''),('133',f'''## 133. 共享微胞前向检查

{pair_text}
'''+('\n\n图 0007-133-1：(a)五种顺序/间隔的弱窗面积误差；(b)弱/强分量逐母残差MSE比。每个相位均使用同32对母，不把80个点视为独立试验。' if pair else '')),
    ('134',f'''## 134. 恢复关闭对照与分量解释

{off_text}
'''+('\n\n图 0007-134-1：(a)(b)第一对冻结母在8/40 ns条件的强分量；(c)(d)弱分量。图示为相位0，前向输入包含真实E/t；灰色孤立对照与共享实际祖先均值分开。完整窗展示，不依据示例挑选统计结论。' if pair else '')),
    ('135',f'''## 135. 运行、验证及下一步

{counts}

{validation}

{next_text}

未评价盲漏检、误拆分、误合并、背景、能谱选择效率或真实器件性能；均标为未测，不能填0。冻结门是项目设计，无文献性能保证。沿用L91/L93/L96已保存的文献条件和人工审核状态，没有新增论文检索、PDF核查或共享库修改。模型假设、前向真值诊断和已知K/t电压拟合分别解释。

位置：pileup/experiments/alpha_current_shared_20261102；模型/源宏/源码/运行器和freeze清单保留。ROOT统一g4/data，图和数组位于jupyter/artifacts同名目录，日志及旧笔记备份位于jupyter/run_records同名目录。完整ROOT与实际宏、三组种子及状态写入Run_Catalog；研究计划第38节及L97同步。原Notebook计算和图像保留。
''')]
    append_notebooks(sections,paths,key,next_text,inventory,counts)
    write_documents(key,next_text,energy_table,wave_table,pair_table,counts,validation,failures,sg,pg)
    outputs=[*paths.values(),EXP/'README.md',ROOT/'README.md',ROOT/'docs/literature-use.md',
             *sorted((EXP/'tables').glob('*.csv')),*[ROOT/'jupyter'/name for name in
             ('codex_0007_SiPM_Energy_Calibration.ipynb','Research_Plan_20260921.ipynb','Run_Catalog.ipynb')]]
    save(EXP/'report-output-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in outputs})
    print(f'Wrote {len(paths)} figures, sections130–135, plan38, L97 and136-ROOT catalog',flush=True)


def append_notebooks(sections,paths,key,next_text,inventory,counts):
    names=('codex_0007_SiPM_Energy_Calibration.ipynb','Research_Plan_20260921.ipynb','Run_Catalog.ipynb')
    notebooks={name:read(ROOT/'jupyter'/name) for name in names}
    for name in names:assert notebooks[name]==read(REC/(name+'.before')),name
    assert [len(notebooks[name]['cells']) for name in names]==[210,106,67]
    cells=[]
    for number,body in sections:
        cells.append(nbformat.v4.new_markdown_cell(body,id=f'current-shared-{number}'))
        if number in paths:
            source=("from pathlib import Path\nfrom IPython.display import display, Image\n"
                    "root = Path.cwd()\nif root.name == 'jupyter':\n    root = root.parent\n"
                    f"display(Image(filename=str(root / '{paths[number].relative_to(ROOT).as_posix()}')))\n")
            cells.append(nbformat.v4.new_code_cell(source,id=f'current-shared-figure-{number}'))
    new=nbformat.v4.new_notebook(cells=cells)
    NotebookClient(new,timeout=120,kernel_name='python3',resources={'metadata':{'path':str(ROOT)}}).execute()
    notebooks[names[0]]['cells'].extend(new['cells'])
    body=('## 38. 当前环境四档/四相位独立母与共享响应门（2026-11-02）\n\n'+key+'\n\n'+next_text+
          '\n\n'+counts+'\n\n已知K/t孤立拟合与真值前向诊断均不是盲重建。失败项完整保存在single-gate/pair-gate中；'
          '条件未执行时不把计划数计成实际输出。详见codex_0007第130–135节及L97。')
    notebooks[names[1]]['cells'].append(nbformat.v4.new_markdown_cell(body,id='current-shared-plan-20261102'))
    mapping=table(['ID/用途','ROOT','实际宏','G4/雪崩/读出种子','E_source/E_dep MeV'],[
        [f"{r['id']}/{r['role']}",r['root'],r['macro'],f"{r['seed']}/{r['sipm_seed']}/{r['readout_seed']}",
         f"{r['source_energy_MeV']:.6f}/{r['edep_MeV']:.9f}"] for r in inventory])
    body=('## 当前环境独立四档/四相位及条件共享响应（2026-11-02）\n\n'+key+'\n\n'
          '实验alpha_current_shared_20261102，全部136 ROOT完成，无失败、替换或效果筛选。先导ID0–7，正式ID8–135，四档各32母。'
          '每ROOT一个20 μs窗，5 μs单α中央正入射，无背景，不对应泊松注入率；G4串行、最多8任务，启动≥1.05 s，秒级文件名。'
          'Geant4 11.3.2/ROOT6.34.10/g++11.4.0，NumPy1.26.4/SciPy1.15.3。基于aa0ec00及冻结时未提交驱动，'
          '39份物理源码与已验证二进制哈希不变，actual宏、RunInfo、源代码、配置与design-freeze保留。\n\n'+mapping+'\n\n'
          '全部正式母四相位0/0.5/1/1.5 ns、每相位8重复，离线4096 ns、观测前2048 ns。'
          '单脉冲种子740300000+母ID×100+重复ID；共享计划种子740500000+母对ID×1000+重复ID，读出均加10000000。'
          '共享先导采用独立的先导母0/6、种子740700000+重复ID；各阶段实际完成情况由generation/gate文件记录。'
          '同32对母的多相位/多间隔不是新增独立母。恢复关闭仅离线diagnostic，参数值在protocol保存。\n\n'+counts+'\n\n'
          '所有ROOT/实际宏/日志/SHA256见tables/root_catalog.csv；8先导另做原生确定性读出与步长减半核对，不增加随机数字化数。'
          '数组/图在jupyter/artifacts/alpha_current_shared_20261102，日志与旧Notebook备份在jupyter/run_records同名目录。\n\n'+next_text)
    notebooks[names[2]]['cells'].append(nbformat.v4.new_markdown_cell(body,id='current-shared-runs-20261102'))
    for name,nb in notebooks.items():
        nbformat.validate(nbformat.from_dict(nb));(ROOT/'jupyter'/name).write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n')


def write_documents(key,next_text,energy_table,wave_table,pair_table,counts,validation,failures,sg,pg):
    failure_text='；'.join(failures) or '无'
    (EXP/'README.md').write_text(f'''# 当前环境独立四档、四相位及条件共享响应检查

{key}

{energy_table}

{wave_table}

孤立能量门：旧/修正均值E偏差绝对值≤5%，单次E RMSE修正/旧≤1.1，修正相位偏差跨度≤1个百分点，未收敛0；同时要求组均值波形残差MSE比≤1.1。失败项：{failure_text}。统计和门限在新数据前冻结，没有重新训练或按结果改门。

{pair_table}

{'共享阶段已执行，' if pg else '尚未执行的共享预案：'}仅比较真实E/t前向预测与祖先分量。弱面积限于前2048 ns观测窗，1200 ns晚脉冲尾部未完整入窗；MSE先算每对母8次均值残差，再平均32对。共享门为修正弱窗面积误差≤5个百分点、弱/强MSE比均≤1.1；恢复关闭在相位0检查绝对面积变化≤3个百分点。10000次整母对配对重采样区间是逐终点描述性区间，不宣称联合覆盖。已知K=1/t孤立能量拟合只读电压，真实E用于评分；没有双脉冲能量拟合或默认升级。

全部136 ROOT保留：8先导、128正式；四源能量1.625/2.5/5.486/8 MeV各32正式母，实际沉积能量逐母记录。单α在20 μs窗的5 μs中央正入射，无背景、无泊松注入率。G4种子740000000+ID，原生雪崩/读出加100000/200000。串行G4、最多8独立任务、启动≥1.05 s，ROOT统一位于g4/data，源宏/种子/RunInfo/日志映射在tables/root_catalog.csv。正式样本未筛选或替换。

离线4096 ns、观测前2048 ns，首有效面光子128+相位 ns，2 ns采样/0.5 ns孔径。每母或母对每条件8重复；共享32对固定母复用于五顺序/间隔、四相位，不能计为640对独立母。完整设计及所有随机种子公式在protocol.json和各plan.json。恢复关闭只改离线SensorConfig.recovery_enabled=True→False，不修改原G4或学习新模型。

{counts}

{validation}

项目根目录核对已有输出：

```bash
python3 -B pileup/experiments/{study.NAME}/run.py verify
python3 -B pileup/experiments/{study.NAME}/validate.py
```

同环境数值重跑可先执行`python3 -B pileup/experiments/{study.NAME}/reproduce.py --target-root /path/to/new-project`。该目录必须不存在；脚本仅复制冻结输入与二进制，建立空输出和运行清单，不启动生产、不覆盖旧结果。进入新项目后按protocol顺序执行复制的run.py各阶段，并遵守各门；可运行至replay而不依赖Git切换。自动报告仍需要本地既有Notebook上下文。空项目的输入哈希验证及拒绝覆盖已实测通过，此检查没有生成新ROOT。

run.py使用本目录实际source副本；完整源码/模型/宏/运行器和freeze清单保留，重跑无需切换Git提交。既有生产阶段拒绝覆盖；需要新运行时另建实验和输出。当前环境Geant4 11.3.2/ROOT6.34.10；39份物理源码及L93模型保持，旧模型开发环境不同，因此仍为跨环境验证。build.json记录已验证二进制来源/编译方式，transport-environment.json记录数据集；新机器可用source/g4及年龄模型C++源码重新构建。ROOT、模型NPZ、二进制、输出数组、Notebook仅本地保留，不随Git分发，跨机器复现需要另提供这些文件。

图及数组：jupyter/artifacts/{study.NAME}；日志、旧笔记备份：jupyter/run_records/{study.NAME}。Notebook第130–135节、研究计划第38节、Run_Catalog和L97同步。验证脚本只读核对，不重复生产。

{next_text}
''')
    literature=ROOT/'docs/literature-use.md';before=(REC/'literature-use.md.before').read_text()
    assert literature.read_text()==before and '## L97 ' not in before
    literature.write_text(before+f'''

<a id="current-shared-20261102"></a>
## L97 — 2026-11-02：当前环境四档/四相位新母与条件共享响应

**项目结果。** {key}

{energy_table}

{wave_table}

**设计与证据层级。** L96强源检验后冻结模型，独立新增136个ROOT：8先导与128正式分开，1.625/2.5/5.486/8 MeV各32正式母；实际E_dep另列。中央正入射、单α在20 μs窗内5 μs，无背景和泊松注入率；G4种子740000000+ID，原生雪崩/读出加100000/200000。39份物理源码与L93模型不变，包括95/351 ns、快比例0.3、29 ns恢复代理、14400微胞和空间12×12先验。当前Geant4 11.3.2/ROOT6.34.10与旧模型开发环境不同，保留迁移验证边界。原生/离线数字化分开记账。

**孤立方法与门。** 每母四相位0/0.5/1/1.5 ns，每相位8重复，4096 ns离线、前2048 ns观测；首有效光子对齐128+相位 ns。已知K=1/t，电压共同初值/GLS权重/预算，拟合进程禁止读取真值；E_dep只在评分阶段使用。能量要求两模型均值E偏差≤5%、修正/旧单次RMSE≤1.1、修正相位跨度≤1个百分点、未收敛0；同时要求组均值前向残差MSE比≤1.1。波形门用先平均母再平方，逐母MSE单独保留，不能互换。失败项：{failure_text}。

**共享条件与实际范围。**

{pair_table}

{'共享阶段使用' if pg else '共享预案为'}32对独立弱/强母、强后弱8/40/200/1200 ns及弱后强40 ns、四相位；每条件8次共享及两个孤立对照，恢复关闭仅五条件相位0。{'本轮在孤立门通过后执行全部条件，先导母未进入正式统计。' if pg else '仅当孤立门通过才执行，未执行的计划不计实际结果。'}前向诊断明确使用真实E/t；全部预测在正式共享输出前冻结。弱窗面积只积分前2048 ns，1200 ns晚脉冲尾部部分窗外，不称完整电荷。共享门要求修正弱面积误差≤5个百分点、逐母弱/强分量MSE比均≤1.1；恢复关闭弱面积均值变化≤3个百分点。10000次整母对配对区间逐终点描述，不传播旧模型训练误差、不宣称同时覆盖；门均为项目冻结设计，不是文献保证。

**计数与验证。** {counts} {validation} 两个计划/共享项复用测试通过。全部ROOT实际宏、RunInfo、三组种子、单窗、哈希及光子缓存核对，启动间隔≥1 s；源码、模型和输出分别冻结。报告只追加已有主题Notebook，保留旧计算/图像。

**采用与下一步。** {next_text} 盲漏检、误拆分/误合并、能谱选择效率、背景和真实器件性能未测；孤立误差SD不当作实测FWHM，±10%条件能量计数不当作盲效率。

**文献来源范围。** 本轮复用L91/L93/L96已记录的条件均值思想与方法适用条件，未新检索论文、未回查原PDF、未更改原人工审核状态或共享库。新增结果属于项目计算，不作为外部论文结论。

位置：`pileup/experiments/{study.NAME}/`；`pileup/analysis/alpha_current_shared.py`、`alpha_current_shared_study.py`、`alpha_current_shared_report.py`。本地codex_0007第130–135节，固定图号0007-131-1/132-1{'/133-1/134-1' if pg else ''}；研究计划第38节及Run_Catalog逐ROOT清单。
''')
    readme=ROOT/'README.md';before=(REC/'README.md.before').read_text();assert readme.read_text()==before
    readme.write_text(before+f'\n\n2026-11-02完成[当前环境独立四档/四相位响应检查](pileup/experiments/{study.NAME}/README.md)：'
        +key+' '+next_text+' 136个ROOT与全部条件保留，未升级默认重建。结果见本地codex_0007第130–135节及'
        '[采用记录L97](docs/literature-use.md#current-shared-20261102)。\n')


if __name__=='__main__':main()
