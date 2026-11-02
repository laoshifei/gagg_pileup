"""Append the fixed known-K/time energy pilot without rewriting old cells."""
from __future__ import annotations
import csv
import json
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/gagg_current_pair_energy_mpl')
os.environ.setdefault('IPYTHONDIR','/tmp/gagg_current_pair_energy_ipython')
os.environ.setdefault('JUPYTER_RUNTIME_DIR','/tmp/gagg_current_pair_energy_jupyter')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nbformat
from nbclient import NotebookClient
import numpy as np
from . import alpha_current_pair_energy_study as study
from .alpha_current_pair_energy import CONDITIONS
ROOT,EXP,ART,REC=study.ROOT,study.EXP,study.ART,study.REC
read,save,sha=study.read,study.save,study.sha


def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(map(str,row))+' |' for row in rows])


def csv_table(name,rows):
    (EXP/'tables').mkdir(exist_ok=True)
    with (EXP/'tables'/f'{name}.csv').open('w',newline='') as stream:
        fields=list(dict.fromkeys(key for row in rows for key in row))
        writer=csv.DictWriter(stream,fieldnames=fields,lineterminator='\n');writer.writeheader()
        writer.writerows({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in row.items()} for row in rows)


def figures(summary,details):
    plt.rcParams.update({'font.size':10,'figure.dpi':140,'axes.grid':True,'grid.alpha':.2})
    paths={};groups=[r for r in summary['groups'] if r['kind']=='individual']
    x=np.arange(5);labels=['H-W 8','H-W 40','H-W 200','H-W 1200','W-H 40']
    fig,axes=plt.subplots(2,2,figsize=(12,7),constrained_layout=True)
    for method,color,offset in [('old','tab:blue',-.12),('corrected','tab:orange',.12)]:
        for row,component in enumerate(('weak','strong')):
            for col,key in enumerate(('bias_percent','rmse_percent')):
                axis=axes[row,col]
                axis.bar(x+offset,[r[method][component+'_'+key] for r in groups],width=.24,color=color,label=method)
                axis.set(xticks=x,xticklabels=labels,xlabel='Order / gap (ns)',ylabel=component+' E '+('bias / %' if col==0 else 'RMSE / %'),
                         title=f'({chr(97+row*2+col)}) '+component+'; 8 mother pairs per condition')
                axis.axhline(0,color='grey',lw=.7)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Figure 0007-137-1: Individual noisy pair energies; known K=2 / time')
    p=ART/'figures/fig_0007_137_1_energy.png';fig.savefig(p);plt.close(fig);paths['137']=p
    fig,axes=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    for method,color,offset in [('old','tab:blue',-.12),('corrected','tab:orange',.12)]:
        axes[0].bar(x+offset,[r[method]['strict_pair_percent'] for r in groups],width=.24,label=method,color=color)
    axes[0].set(xticks=x,xticklabels=labels,xlabel='Order / gap (ns)',ylabel='Both E within 10% and converged / %',
                ylim=(0,105),title='(a) 8 individual noisy observations per condition')
    axes[0].legend(fontsize=8)
    for j,r in enumerate(summary['short_pooled']):
        ci=r['comparison']['weak_rmse_ratio_ci95'];mean=r['comparison']['weak_rmse_ratio']
        axes[1].vlines(j,*ci,color=f'C{j}',lw=2);axes[1].scatter(j,mean,color=f'C{j}',zorder=3)
    axes[1].axhline(1,color='grey',lw=.8);axes[1].axhline(.95,color='grey',ls=':',lw=.8)
    axes[1].set(xticks=[0,1],xticklabels=['Individual','8-repeat mean'],ylabel='Corrected / old weak E RMSE',
                title='(b) Pooled H-W 8/40 ns; paired mother 95% intervals')
    fig.suptitle('Figure 0007-138-1: Conditional energy success and paired pilot uncertainty')
    p=ART/'figures/fig_0007_138_1_success.png';fig.savefig(p);plt.close(fig);paths['138']=p
    with np.load(ART/'examples.npz') as data:
        t=data['time_ns'];background=data['background_mV']
        fig,axes=plt.subplots(2,2,figsize=(12,7),constrained_layout=True)
        for col,cid in enumerate((0,1)):
            axes[0,col].plot(t,data[f'case{cid}_voltage_mV'],color='black',lw=.8,label='Observed analog')
            axes[1,col].plot(t,data[f'case{cid}_truth_components_mV'][1],color='black',lw=.8,label='True weak ancestry')
            for method,color in [('old','tab:blue'),('corrected','tab:orange')]:
                components=data[f'case{cid}_{method}_components_mV']
                axes[0,col].plot(t,components.sum(0)+background,color=color,lw=1,label=method+' fitted')
                axes[1,col].plot(t,components[1],color=color,lw=1,label=method+' fitted weak')
            axes[0,col].set(title=f'({chr(97+col)}) First selected pair: {[8,40][col]} ns',ylabel='Total voltage / mV')
            axes[1,col].set(title=f'({chr(99+col)}) Weak component: diagnostic truth overlay',ylabel='Component / mV')
            for axis in axes[:,col]:axis.set(xlim=(0,2048),xlabel='Aligned time / ns');axis.legend(fontsize=7)
        fig.suptitle('Figure 0007-139-1: Full observed window; repeat 0; fixed first hash-selected pair')
        p=ART/'figures/fig_0007_139_1_waveforms.png';fig.savefig(p);plt.close(fig);paths['139']=p
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for axis,gap in zip(axes,(8.,40.)):
        for kind,marker,color in [('individual','o','tab:blue'),('mean','x','tab:orange')]:
            rows={m:[next(r for r in details if r['method']==m and r['kind']==kind and r['block']==b and
                  r['condition']=='high_weak' and r['gap_ns']==gap) for b in summary['blocks']] for m in ('old','corrected')}
            axis.scatter([r['weak_error_percent'] for r in rows['old']],[r['weak_error_percent'] for r in rows['corrected']],
                         marker=marker,color=color,label=kind,s=35)
        all_errors=[r['weak_error_percent'] for r in details if r['gap_ns']==gap and r['condition']=='high_weak']
        limits=[min(-10,min(all_errors)),max(10,max(all_errors))]
        axis.plot(limits,limits,'k--',lw=.8);axis.axhline(0,color='grey',lw=.5);axis.axvline(0,color='grey',lw=.5)
        axis.set(xlabel='Old weak E error / %',ylabel='Corrected weak E error / %',title=f'H-W {gap:g} ns; same 8 mother pairs')
        axis.legend(fontsize=8)
    axes[0].set_title('(a) H-W 8 ns; same 8 mother pairs');axes[1].set_title('(b) H-W 40 ns; same 8 mother pairs')
    fig.suptitle('Figure 0007-140-1: Mother-pair variability and overlapping mean observations')
    p=ART/'figures/fig_0007_140_1_scatter.png';fig.savefig(p);plt.close(fig);paths['140']=p
    return paths


def main():
    study.verify();assert read(EXP/'replay.json')['passed']
    summary=read(EXP/'result-summary.json');details=read(EXP/'fit-details.json');gate=read(EXP/'gate.json')
    protocol=read(EXP/'protocol.json');numeric=read(EXP/'numerical.json')
    paths=figures(summary,details)
    csv_table('energy_groups',summary['groups']);csv_table('short_pooled',summary['short_pooled']);csv_table('fit_details',details)
    mothers=read(EXP/'input-mothers.json');cases=read(EXP/'cases.json')
    catalog=[{k:r[k] for k in ('id','root','macro','snapshot_macro','seed','sipm_seed','readout_seed',
                              'source_energy_MeV','edep_MeV','root_sha256','log')} for r in mothers]
    csv_table('reused_root_catalog',catalog);csv_table('waveform_cases',cases)
    individual=[r for r in summary['groups'] if r['kind']=='individual'];means=[r for r in summary['groups'] if r['kind']=='mean']
    energy_table=table(['顺序/间隔ns','整对±10%成功：旧→修正','弱E偏差/%','弱E RMSE/%','弱误差SD/%','边界/未收敛：旧→修正'],[
        [f"{r['condition']}/{r['gap_ns']:g}",f"{r['old']['strict_pair_count']}/8→{r['corrected']['strict_pair_count']}/8",
         f"{r['old']['weak_bias_percent']:+.3f}→{r['corrected']['weak_bias_percent']:+.3f}",
         f"{r['old']['weak_rmse_percent']:.3f}→{r['corrected']['weak_rmse_percent']:.3f}",
         f"{r['old']['weak_error_sd_percent']:.3f}→{r['corrected']['weak_error_sd_percent']:.3f}",
         f"{r['old']['boundary']}/{r['old']['nonconverged']}→{r['corrected']['boundary']}/{r['corrected']['nonconverged']}"] for r in individual])
    strong_table=table(['顺序/间隔ns','强E偏差/%：旧→修正','强E RMSE/%','强误差SD/%'],[
        [f"{r['condition']}/{r['gap_ns']:g}",f"{r['old']['strong_bias_percent']:+.3f}→{r['corrected']['strong_bias_percent']:+.3f}",
         f"{r['old']['strong_rmse_percent']:.3f}→{r['corrected']['strong_rmse_percent']:.3f}",
         f"{r['old']['strong_error_sd_percent']:.3f}→{r['corrected']['strong_error_sd_percent']:.3f}"] for r in individual])
    mean_table=table(['强→弱间隔ns','均值弱E偏差/%：旧→修正','均值弱E RMSE/%','均值强E偏差/%','均值整对±10%计数'],[
        [int(r['gap_ns']),f"{r['old']['weak_bias_percent']:+.3f}→{r['corrected']['weak_bias_percent']:+.3f}",
         f"{r['old']['weak_rmse_percent']:.3f}→{r['corrected']['weak_rmse_percent']:.3f}",
         f"{r['old']['strong_bias_percent']:+.3f}→{r['corrected']['strong_bias_percent']:+.3f}",
         f"{r['old']['strict_pair_count']}/8→{r['corrected']['strict_pair_count']}/8"] for r in means])
    uncertainty_table=table(['8/40ns合并','弱E RMSE比及95%区间','整对成功增益/百分点及95%区间'],[
        [r['kind'],f"{r['comparison']['weak_rmse_ratio']:.4f} [{r['comparison']['weak_rmse_ratio_ci95'][0]:.4f}, {r['comparison']['weak_rmse_ratio_ci95'][1]:.4f}]",
         f"{r['comparison']['strict_gain_pp']:+.2f} [{r['comparison']['strict_gain_ci95'][0]:+.2f}, {r['comparison']['strict_gain_ci95'][1]:+.2f}]"] for r in summary['short_pooled']])
    primary=summary['short_pooled'][0]
    key=(f"在L97保存的电压上完成40条单次、16条八次平均双脉冲的旧/修正模型比较，共112次已知K=2/时刻的能量拟合。"
         f"8/40 ns单次整对E±10%成功为{primary['old']['strict_pair_count']}/16→{primary['corrected']['strict_pair_count']}/16，"
         f"弱E RMSE为{primary['old']['weak_rmse_percent']:.2f}%→{primary['corrected']['weak_rmse_percent']:.2f}%；"
         f"冻结推进门{'通过' if gate['passed'] else '未通过'}。")
    next_text=('下一步使用预留母对作冻结扩展，并另留新生成母事件确认；扩展前保持模型、评价口径和默认流程。'
               '当前只有8对已作响应检查的母事件，不能据此认定盲重建、高率或实物收益。' if gate['passed'] else
               '按冻结规则停止扩样和默认替换。下一步在开发资料中分解强/弱孤立残差、共享恢复残差与噪声的作用，'
               '解释本轮失败条件；任何模型或权重修订另建快照，并另留确认样本，避免用这8对母反复选参。')
    failure_text='；'.join(gate['failures']) or '无'
    counts=('正式112次双脉冲＋32次孤立拟合；排除在正式评价外的耗时先导6次，合成反演4次，重放2次双脉冲＋4次孤立。'
            '总计160次已完成并保存的能量拟合；新ROOT、随机数字化和训练均为0。实际输入为48个旧电压文件，追溯20个既有ROOT。'
            '首次正式进程收到SIGTERM，29条完整结果保留；从缺失83条任务继续，未改模型或计划。'
            '中断原因及未完成的在途调用数无法从现有日志确定，未保存的计算不计入完成数。')
    validation=(f"数值导数最大相对差{max(r['gradient_relative_error'] for r in numeric['gradient_rows']):.3g}，"
                f"合成反演最大相对能量差{max(r['max_relative_error'] for r in numeric['inverse_rows']):.3g}。"
                '32次孤立控制与L97相同输入的结果逐字段一致（计时除外）；6次拟合重放同样一致。'
                '两项母对统计/抽样测试通过，已验证空目录重跑输入完整且拒绝覆盖旧目录。')
    gates=('推进门沿用L89短间隔要求：单次8/40 ns合并严格整对成功增加≥10个百分点，弱RMSE比≤0.95；'
           '各短间隔弱RMSE比≤1.1，八次均值的修正弱偏差绝对值≤10%。另要求修正强均值偏差绝对值≤5%、双脉冲未收敛0。'
           '这些是项目决策门，不是文献性能保证；全部门须同时满足，不按条件择优。')
    evidence=read(EXP/'literature-evidence.json')
    literature_text='\n\n'.join(
        f"{'; '.join(r['authors'])}，*{r['title']}*，[DOI](https://doi.org/{r['doi']})。据已有精读笔记：{r['source_positions']}。"
        f"笔记位置 `{r['note']}`；原PDF位置 `{r['pdf']}`。本轮未重新读PDF，DOI与页码据笔记。"
        +('TES非线性/非平稳噪声模型的条件均值与协方差原则，不是PCA或GAGG双脉冲性能证明；本轮固定共同GLS工作度量，不称完整参数相关似然。'
           if 'Fixsen' in r['note'] else
           'LaBr3(Ce)/R6233 PMT、500 MHz、14 bit；同模板高斯噪声模拟不能验证SiPM共享微胞，能峰效率不能替代逐事件能量成功率。')
        +'原状态AI草稿、待人工审核；没有转移Busch第658项的限定确认。' for r in evidence['entries'])
    sections=[('136',f'''## 136. 已知K/t双脉冲能量先导：设计与证据（2026-11-02）

{key}

L97孤立/共享响应门通过后固定本设计；模型与物理参数不变，源1.625/8 MeV、实际逐母E_dep分别记录。当前G4 11.3.2/ROOT6.34.10输入已用于L97前向检查，所以本轮是同批电压的能量反演扩展，不称全新独立确认。没有调用旧环境、生成新ROOT或重新数字化。

耗时先导仅用L97排除母0/6的8/40/1200 ns、相位0、第一重复，各两模型，共6次；实际每拟合约60–131 s。先导能量效果未用于选择方案。按固定字符串SHA256排序，选择母对block {summary['blocks']}；只用索引，不看能量、残差或成功率。8对父母互不重叠，相位0；五条件各用重复0（40条单次），仅8/40 ns再用同母全部8重复均值（16条）。多条件共用母，对应的56条波形不是56对独立母。

另选四能区各2个孤立母，用重复0和八次均值，32次拟合作输入/实现控制。每条都与L97已保存结果逐字段对照；这8条孤立案例不是新独立能标确认。

{literature_text}
'''),('137',f'''## 137. 单次双脉冲沉积能量结果

{energy_table}

{strong_table}

每行8对母、每对一次完整含噪电压。严格成功要求两份E_dep相对误差均≤10%且所选解收敛；边界与未收敛都保留在分母中。误差SD描述这批相对能量误差散布，不是单能峰FWHM；条件成功率不等于盲检出、类型鉴别或完整重建效率。

图 0007-137-1：(a)(b)弱能量偏差与RMSE；(c)(d)强能量偏差与RMSE。蓝色旧空间模型、橙色冻结均值修正；五种条件完整展示。
'''),('138',f'''## 138. 配对不确定度与推进门

{uncertainty_table}

10000次整母对配对重采样，种子741900000；同一母对的全部条件、旧/修正以及单次/均值一起移动。只有8个独立母对，16个短间隔观测并非16个独立样本。区间逐终点、条件于冻结模型和已观察重复；没有传播旧能标/空间先验/噪声模型及形状修正训练误差，不宣称同时覆盖或总体等效性。推进门使用预先规定的点估计，区间另列，不事后改门。

{gates}

失败项：{failure_text}。

图 0007-138-1：(a)各条件单次整对能量成功率；(b)8/40 ns合并弱E RMSE比及整母对95%区间。八次均值只作偏差诊断。
'''),('139',f'''## 139. 同输入、同权重和预算的完整波形比较

两模型均只输入2048 ns完整含噪电压与明确给定的K=2/t；时刻是首有效面光子对齐128 ns及已知间隔，不是盲触发估计。电压导出的共同旧模型NNLS初值、两个0.1/0.9电荷分配起点和GLS权重在同一观测的两方法间逐项相同。固定logE边界，每起点max_nfev=80，ftol/xtol=1e-8、gtol=1e-5；选择最低代价候选，保留未收敛候选。有限搜索预算不证明全局最优。

权重沿用冻结标定，是工作度量，没有重训或在拟合过程中随候选E更新，不声称已校准共享器件的完整概率似然。推断工作进程的文件审计禁止读取ROOT、真实E、母标签、光子/祖先/预测数组；真值只在评分、合成数值控制和本节诊断图中使用。

图 0007-139-1：(a)(b)固定哈希排序第一对母、重复0、8/40 ns的完整观测和拟合；(c)(d)拟合弱分量与真实祖先信号对照，真值未进入拟合。示例没有按成功或失败挑选。
'''),('140',f'''## 140. 八次平均波形与母间变化

{mean_table}

均值是同母8次含噪数字化平均，仍包含母事件的光学实现变化及有限重复噪声，并非无噪声条件期望；其中包含单次评价使用的重复0，两组不能当作独立重复证据。均值沿用原工作度量的构造，没有另行标定八次平均波形的协方差。200/1200 ns和反向40 ns未做均值能量拟合，不补造相应数字。相位0之外、未知K/t、自由形状参数和不同器件条件均未在本轮拟合。

图 0007-140-1：(a)(b)8/40 ns逐母弱E误差，横轴旧、纵轴修正；圆点为单次、叉号为八次均值，同一母在两种波形下相关。对角线表示能量误差不变。
'''),('141',f'''## 141. 复核、计数及下一步

{counts}

{validation}

{next_text}

盲漏检、误拆分/误合并、能谱选择效率、背景和真实器件性能未测，不写为0。采用记录L98与研究计划第39节同步。结果位置pileup/experiments/alpha_current_pair_energy_20261102；电压副本/拟合/图在jupyter/artifacts同名目录，耗时先导和旧笔记备份在jupyter/run_records同名目录。20个输入ROOT沿用原RunInfo与实际宏，Run_Catalog明确标记复用，没有增加新运行。旧Notebook全部计算和图像保留。
''')]
    append_notebooks(sections,paths,key,next_text,catalog,counts)
    documents(key,next_text,energy_table,strong_table,mean_table,uncertainty_table,gates,failure_text,counts,validation,evidence)
    outputs=[*paths.values(),EXP/'README.md',ROOT/'README.md',ROOT/'docs/literature-use.md',
             *sorted((EXP/'tables').glob('*.csv')),*[ROOT/'jupyter'/n for n in
             ('codex_0007_SiPM_Energy_Calibration.ipynb','Research_Plan_20260921.ipynb','Run_Catalog.ipynb')]]
    save(EXP/'report-output-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in outputs})
    print('Wrote4figures,sections136–141,plan39,L98 and20reused-ROOT mappings',flush=True)


def append_notebooks(sections,paths,key,next_text,catalog,counts):
    names=('codex_0007_SiPM_Energy_Calibration.ipynb','Research_Plan_20260921.ipynb','Run_Catalog.ipynb')
    notebooks={name:read(ROOT/'jupyter'/name) for name in names}
    for name in names:assert notebooks[name]==read(REC/(name+'.before')),name
    assert [len(notebooks[name]['cells']) for name in names]==[220,107,68]
    cells=[]
    for number,body in sections:
        cells.append(nbformat.v4.new_markdown_cell(body,id=f'current-pair-energy-{number}'))
        if number in paths:
            source=("from pathlib import Path\nfrom IPython.display import display, Image\n"
                    "root = Path.cwd()\nif root.name == 'jupyter':\n    root = root.parent\n"
                    f"display(Image(filename=str(root / '{paths[number].relative_to(ROOT).as_posix()}')))\n")
            cells.append(nbformat.v4.new_code_cell(source,id=f'current-pair-energy-figure-{number}'))
    new=nbformat.v4.new_notebook(cells=cells)
    NotebookClient(new,timeout=120,kernel_name='python3',resources={'metadata':{'path':str(ROOT)}}).execute()
    notebooks[names[0]]['cells'].extend(new['cells'])
    body=('## 39. 当前环境已知K/t双脉冲能量先导（2026-11-02）\n\n'+key+'\n\n'+next_text+'\n\n'+counts+
          '\n\n固定8对母、相位0、完整五种条件；8/40 ns追加同母八次均值诊断。输入已作L97前向检查，不能称全新确认。'
          '拟合只读电压与已知K/t；真值在评分及数值控制使用。四幅固定编号图和完整门/区间见codex_0007第136–141节，L98同步。')
    notebooks[names[1]]['cells'].append(nbformat.v4.new_markdown_cell(body,id='current-pair-energy-plan-20261102'))
    mapping=table(['原母ID','复用ROOT','实际宏副本','G4/雪崩/读出种子','E_source/E_dep MeV'],[
        [r['id'],r['root'],r['snapshot_macro'],f"{r['seed']}/{r['sipm_seed']}/{r['readout_seed']}",
         f"{r['source_energy_MeV']:.6f}/{r['edep_MeV']:.9f}"] for r in catalog])
    body=('## 当前环境双脉冲能量先导：复用20个ROOT，无新运行（2026-11-02）\n\n'+key+'\n\n'
          '实验alpha_current_pair_energy_20261102复用L97的20个原生ROOT；原来单α20 μs窗、5 μs中央正入射，'
          '源能量1.625/2.5/5.486/8 MeV，三组种子与实际宏逐条列于下表，不把当前宏推作历史条件。'
          'Geant4 11.3.2/ROOT6.34.10原环境、真实RunInfo与源宏不变；本轮未调用G4。\n\n'+mapping+'\n\n'
          '48个已有完整含噪电压文件逐字节复制并冻结：40双脉冲条件、8孤立控制。双脉冲只取相位0，'
          '五条件均用重复0，8/40 ns另用全部8重复均值。离线输入仍为4096 ns数字化的前2048 ns、1024点，'
          '2 ns采样/0.5 ns孔径。原双脉冲随机种子740500000+母对block×1000+重复ID，孤立740300000+母ID×100+重复ID，读出加10000000；'
          '本轮只是读取这些波形，没有再次消耗仿真随机流。\n\n'+counts+'\n\n'
          '代码基于8231f62及新冻结驱动；旧/修正L93模型、age.so和相关源码哈希保持。推断最多8进程、BLAS单线程，'
          '已完成29条后收到SIGTERM，原因未确定；完整结果哈希保留，从余下83条恢复，模型/预算/抽样不变。'
          '新增source/实际宏/RunInfo副本与输入哈希在实验目录；拟合、图及电压副本在jupyter/artifacts同名目录；'
          '耗时先导、恢复日志和Notebook备份在jupyter/run_records同名目录。\n\n'+next_text)
    notebooks[names[2]]['cells'].append(nbformat.v4.new_markdown_cell(body,id='current-pair-energy-runs-20261102'))
    for name,nb in notebooks.items():
        nbformat.validate(nbformat.from_dict(nb));(ROOT/'jupyter'/name).write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n')


def documents(key,next_text,energy_table,strong_table,mean_table,uncertainty_table,gates,failure_text,counts,validation,evidence):
    literature_table=table(['来源与位置（据已有笔记）','本轮采用范围与审核边界'],[
        [f"{'; '.join(r['authors'])}, *{r['title']}*, [DOI](https://doi.org/{r['doi']})；"
         f"[笔记]({r['note']})、[原PDF]({r['pdf']})；{r['source_positions']}",
         ('TES条件均值/完整协方差的局部最小二乘原则；非PCA、非本GAGG实测保证。本轮固定共同GLS工作度量，未声称完整参数相关似然。'
          if 'Fixsen' in r['note'] else
          'LaBr3(Ce)/R6233 PMT、500 MHz/14 bit、同模板高斯噪声模拟；借鉴模板条件与失配检查，峰效率不当逐事件能量成功。')+
         'AI草稿待人工审核；未重读PDF或联网核DOI，未转移Busch第658项限定确认。'] for r in evidence['entries']])
    (EXP/'README.md').write_text(f'''# 当前环境固定K/t双脉冲能量先导

{key}

{energy_table}

{strong_table}

{mean_table}

{uncertainty_table}

{gates}

失败项：{failure_text}。

输入为L97已检查的完整含噪电压，模型保持L93旧空间12×12与面积保持均值修正。按固定SHA256排序选择8对母，phase0，五种顺序/间隔；每条件使用重复0，仅8/40 ns补同母八次平均。56条双脉冲波形仍只有8对独立母，单次与均值还共享重复0；均值不是无噪声理论期望，仍沿用原工作度量，没有另标定均值协方差。此轮是已暴露响应样本上的能量先导，不是新增母事件的独立确认。

两模型使用同一电压导出初值/固定GLS权重、相同三起点和搜索预算；真实K=2及到达时刻明确给定，真实E不进推断。推断文件审计禁止读取ROOT、母标签、光子和真值数组。严格整对成功要求两能量误差均≤10%且收敛；所有边界和未收敛结果保留。误差SD不是实测单能峰FWHM；盲漏检、误拆分/误合并、能谱选择效率、背景和真实器件性能未测。

10000次整母对配对重采样，种子741900000；条件/方法/单次和均值均随同一母对移动。8母的区间仅作条件描述，不传播原模型/能标/空间先验和均值训练误差，不宣称联合覆盖。推进门是项目设计，沿用L89短间隔标准并增加强均值/收敛控制，没有按本轮输出改门。

{counts}

{validation}

48个电压文件和20个ROOT/实际宏映射分别在tables/waveform_cases.csv、tables/reused_root_catalog.csv。实际源码、配置、原生RunInfo及宏副本保留；全部模型与age.so沿用已验证二进制哈希，物理源码未改，未重新构建G4。当前输入环境为Geant411.3.2/ROOT6.34.10，分析NumPy1.26.4/SciPy1.15.3；版本与运行前未提交状态在protocol记录。

项目根目录只读验证：

```bash
python3 -B pileup/experiments/{study.NAME}/run.py verify
python3 -B pileup/experiments/{study.NAME}/validate.py
```

同环境重跑先用`python3 -B pileup/experiments/{study.NAME}/reproduce.py --target-root /path/to/new-project`创建不存在的新目录。它复制冻结源码、模型、电压及必要上游输入，保留原输出，然后可在新目录按numerical→single→pair→score→replay运行run.py。numerical、single、pair和replay阶段完成后均拒绝重复运行；score会从既有拟合重新生成派生评分文件，宜在新复制目录执行。部分完成的pair可用resume.py继续缺失任务，旧结果哈希保持，最多8进程。长任务宜由终端或持久会话运行；resume.py把状态和日志写入jupyter/run_records对应目录。输入检查/拒绝覆盖已实测，不依赖切换Git。Notebook报告需要已有本地笔记上下文，不能在空复制目录直接生成。

ROOT、NPZ、二进制与Notebook不随Git上传，跨机器需另提供输入与权重并按build记录重建库。图与拟合在jupyter/artifacts/{study.NAME}，日志与旧笔记备份在jupyter/run_records同名目录。报告追加至codex_0007第136–141节、四幅固定图号、研究计划第39节及Run_Catalog；L98保存采用范围。

{next_text}
''')
    literature=ROOT/'docs/literature-use.md';before=(REC/'literature-use.md.before').read_text()
    assert literature.read_text()==before and '## L98 ' not in before
    literature.write_text(before+f'''

<a id="current-pair-energy-20261102"></a>
## L98 — 2026-11-02：当前环境已知K/t双脉冲能量先导

**项目结果。** {key}

{uncertainty_table}

**范围和输入。** 复用L97已经作过前向响应检查的电压，不称全新独立确认。固定哈希选8对母，相位0；强→弱8/40/200/1200 ns及弱→强40 ns均用重复0，8/40 ns另用同母8重复均值。40条单次＋16条均值只有8对独立母，均值包含单次所用重复0。所有波形、L93两模型、共同电压初值/GLS权重、三起点和每起点80评估预算提前冻结。真实K/t明确给定，真实E只作评分和数值合成控制；推断文件审计禁止读取ROOT、真实母标签及光子/祖先数组。

**判定。** {gates} 失败项：{failure_text}。严格整对成功还要求收敛，边界与未收敛保留；完整两分量偏差/RMSE/误差SD及逐母散点见Notebook和表格。8母整对配对区间只是冻结模型下的描述，不把重复或条件当作新增独立母，也不传播先验/模型训练误差。当前能量误差SD不写成实测FWHM；已知K/t条件结果不能用于盲效率、漏检或误拆并。

**本轮文献核验。** 按paper-search只读现有主题导航、两份完整来源笔记及原审核记录。来源身份、作者、DOI和原页位置据既有精读，没有重新读PDF或联网核验。库内默认检索脚本尝试写可选缓存时被只读挂载阻止，随后禁用缓存完成检索；已核实际来源/导航哈希未变。Fixsen检索无命中不作为未审核证明，审核边界以来源状态和原始记录为准。

{literature_table}

**计数与复核。** {counts} {validation} 数值梯度/反演、共同权重与起点、读取隔离、输入哈希、实际宏/RunInfo及原Notebook保留都核对。完整统计不按结果筛样或替换母事件。

**采用与下一步。** {next_text} 盲漏检、误拆分/误合并、能谱选择效率、高率和真实器件未测，不填零；本轮不更换默认方法。

位置：`pileup/experiments/{study.NAME}/`；模块`alpha_current_pair_energy.py`、`alpha_current_pair_energy_study.py`、`alpha_current_pair_energy_report.py`与中断恢复适配器。已有codex_0007第136–141节、图0007-137-1/138-1/139-1/140-1，研究计划第39节与Run_Catalog同步；原计算及图像保留。
''')
    readme=ROOT/'README.md';before=(REC/'README.md.before').read_text();assert readme.read_text()==before
    readme.write_text(before+f'\n\n2026-11-02完成[当前环境已知K/t双脉冲能量先导](pileup/experiments/{study.NAME}/README.md)：'
        +key+' 固定8对母、相位0，复用L97完整电压与原模型；没有新增ROOT或数字化。'+next_text+
        ' 完整结果见本地codex_0007第136–141节及[采用记录L98](docs/literature-use.md#current-pair-energy-20261102)。\n')


if __name__=='__main__':main()
