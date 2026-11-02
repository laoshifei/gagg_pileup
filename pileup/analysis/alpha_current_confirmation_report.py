"""Report the engineering controls and independent current-environment transfer."""
from __future__ import annotations
import csv
import json
import os

os.environ.setdefault('MPLCONFIGDIR','/tmp/gagg_current_confirmation_mpl')
os.environ.setdefault('IPYTHONDIR','/tmp/gagg_current_confirmation_ipython')
os.environ.setdefault('JUPYTER_RUNTIME_DIR','/tmp/gagg_current_confirmation_jupyter')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nbformat
from nbclient import NotebookClient
import numpy as np

from . import alpha_current_confirmation_study as study
from . import alpha_current_response_study as control
ROOT, EXP, ART, REC = study.ROOT, study.EXP, study.ART, study.REC
read, save, sha = study.read, study.save, study.sha


def table(headers, rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(map(str,row))+' |' for row in rows])


def write_csv(experiment, name, rows):
    (experiment/'tables').mkdir(exist_ok=True)
    with (experiment/'tables'/f'{name}.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader()
        writer.writerows({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in row.items()} for row in rows)


def figures(controls, summary, details):
    plt.rcParams.update({'font.size':10,'figure.dpi':140,'axes.grid':True,'grid.alpha':.2})
    paths={};groups=np.arange(4)
    fig,axes=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    for env,style in [('historical','--'),('current','-')]:
        for method,color in [('old','tab:blue'),('corrected','tab:orange')]:
            values=[next(r['bias_percent'] for r in controls['energy'] if r['environment']==env and r['group']==g
                         and r['phase_ns']==0 and r['method']==method and r['kind']=='mean') for g in groups]
            axes[0].plot(groups,values,style+'o',color=color,label=env+' / '+method)
    axes[0].axhline(5,color='grey',ls=':',lw=.8);axes[0].axhline(-5,color='grey',ls=':',lw=.8)
    axes[0].set(xticks=groups,xticklabels=['1.625','2.5','5.486','8'],xlabel='Source energy / MeV',
                ylabel='Mean-waveform fitted E bias / %',title='(a) 4 selected mothers per source / environment')
    axes[0].legend(fontsize=7)
    for g in groups:
        rows=[r for r in controls['energy'] if r['environment']=='current' and r['group']==g and
              r['method']=='corrected' and r['kind']=='mean']
        baseline=next(r['bias_percent'] for r in rows if r['phase_ns']==0)
        axes[1].plot([r['phase_ns'] for r in rows],[r['bias_percent']-baseline for r in rows],'o-',
                     label=['1.625','2.5','5.486','8'][g]+' MeV')
    axes[1].set(xlabel='ADC phase shift / ns',ylabel='Bias relative to phase 0 / percentage points',
                title='(b) Current environment, frozen correction')
    axes[1].legend(fontsize=8)
    fig.suptitle('Figure 0007-124-1: Exposed-sample isolated engineering controls')
    p=ART/'figures/fig_0007_124_1_engineering.png';fig.savefig(p);plt.close(fig);paths['124']=p
    fig,axes=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    x=np.arange(3);rows=summary['projection']
    ci=np.array([r['mse_ratio_ci95'] for r in rows])
    axes[0].vlines(x,ci[:,0],ci[:,1],color='tab:purple',lw=2)
    axes[0].scatter(x,[r['mse_ratio'] for r in rows],color='tab:purple',zorder=3)
    axes[0].axhline(1,color='black',ls='--',lw=.8)
    axes[0].set(xticks=x,xticklabels=['8','40','200 (secondary)'],xlabel='Nominal weak delay / ns',
                ylabel='Corrected / old per-mother projected MSE',title='(a) Paired mother-bootstrap 95% intervals')
    for method,offset,color in [('old',-.08,'tab:blue'),('corrected',.08,'tab:orange')]:
        ci=np.array([r[method+'_mean_ci95'] for r in rows]);means=[r[method+'_mean'] for r in rows]
        axes[1].vlines(x+offset,ci[:,0],ci[:,1],color=color,lw=2)
        axes[1].scatter(x+offset,means,color=color,label=method,zorder=3)
    axes[1].axhline(0,color='black',lw=.8)
    axes[1].set(xticks=x,xticklabels=['8','40','200'],xlabel='Nominal weak delay / ns',
                ylabel='Mean projection / 100 delta ln E',title='(b) Mean of 336 independent mother means')
    axes[1].legend(fontsize=8)
    fig.suptitle('Figure 0007-126-1: Fresh current-environment transfer; fixed models')
    p=ART/'figures/fig_0007_126_1_fresh_projection.png';fig.savefig(p);plt.close(fig);paths['126']=p
    with np.load(ART/'summary.npz') as data:
        t=data['time_ns'];actual=data['actual_mV'].mean(0)
        fig,axes=plt.subplots(2,2,figsize=(12,7),constrained_layout=True)
        axes[0,0].plot(t,actual,color='black',lw=1,label='Fresh 8-repeat ancestor mean')
        all_errors=[r['relative_error_percent'] for r in details if r['repeat']>=0]
        bins=np.linspace(min(all_errors)-.05,max(all_errors)+.05,61)
        for method,color in [('old','tab:blue'),('corrected','tab:orange')]:
            residual=data[method+'_residual_mV'].mean(0)
            axes[0,0].plot(t,actual-residual,color=color,lw=1,label=method+' at true E')
            axes[0,1].plot(t,residual,color=color,lw=1,label=method)
            error=[r['relative_error_percent'] for r in details if r['method']==method and r['repeat']>=0]
            axes[1,0].hist(error,bins=bins,histtype='step',color=color,label=method)
            mean_error=[r['relative_error_percent'] for r in details if r['method']==method and r['repeat']==-1]
            axes[1,1].scatter(data['edep_MeV'],mean_error,s=9,alpha=.5,color=color,label=method)
        axes[0,0].set(xlim=(0,2048),xlabel='Aligned time / ns',ylabel='Signal / mV',title='(a) Complete observed signal window')
        axes[0,0].legend(fontsize=7)
        axes[0,1].set(xlim=(0,2048),xlabel='Aligned time / ns',ylabel='Observed - predicted / mV',title='(b) Mean forward residual at true E')
        axes[0,1].legend(fontsize=8)
        axes[1,0].set(xlabel='Individual fitted E error / %',ylabel='Count',title='(c) 2688 observations in 336 mother clusters')
        axes[1,0].legend(fontsize=8)
        axes[1,1].set(xlabel='True deposited energy / MeV',ylabel='8-repeat mean fitted E error / %',title='(d) Known K=1 / time isolated fits')
        axes[1,1].axhline(0,color='black',lw=.8)
        fig.suptitle('Figure 0007-127-1: Frozen isolated energy response on new strong mothers')
        p=ART/'figures/fig_0007_127_1_energy_waveforms.png';fig.savefig(p);plt.close(fig);paths['127']=p
        fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
        for j,axis in enumerate(axes):
            old=data['old_projection_100logE'][:,j];new=data['corrected_projection_100logE'][:,j]
            axis.scatter(old,new,s=10,alpha=.55)
            limits=[min(old.min(),new.min()),max(old.max(),new.max())]
            axis.plot(limits,limits,'k--',lw=.8);axis.axhline(0,color='grey',lw=.5);axis.axvline(0,color='grey',lw=.5)
            axis.set(xlabel='Old / 100 delta ln E',ylabel='Corrected / 100 delta ln E',
                     title=f'({chr(97+j)}) {[8,40][j]} ns: each point is one fresh mother')
        fig.suptitle('Figure 0007-128-1: Mother variability remains after the fixed mean correction')
        p=ART/'figures/fig_0007_128_1_mother_scatter.png';fig.savefig(p);plt.close(fig);paths['128']=p
    return paths


def main():
    study.verify();control.verify();assert read(EXP/'replay.json')['passed']
    summary=read(EXP/'result-summary.json');controls=read(control.EXP/'result-summary.json')
    gate=read(EXP/'gate.json');details=read(EXP/'fit-details.json')
    paths=figures(controls,summary,details)
    for name in ('energy','projection','waveform'):write_csv(EXP,name,summary[name])
    for name in ('energy','paired_environment_shift','waveform','projection'):write_csv(control.EXP,name,controls[name])
    mothers=read(EXP/'mothers.json')
    inventory=[{k:r[k] for k in ('id','role','root','macro','seed','sipm_seed','readout_seed',
                                 'source_energy_MeV','edep_MeV','root_sha256','log')} for r in mothers]
    write_csv(EXP,'root_catalog',inventory)
    write_csv(control.EXP,'root_catalog',[{k:r[k] for k in ('id','pair','environment','root','macro','seed',
                                                         'source_energy_MeV','edep_MeV','root_sha256')}
                                        for r in read(control.EXP/'mothers.json')])
    def ce(g,method,kind):
        return next(r for r in controls['energy'] if r['environment']=='current' and r['group']==g and
                    r['phase_ns']==0 and r['method']==method and r['kind']==kind)
    engineering_table=table(['源能量/MeV','均值波形E偏差/%：旧→修正','单次E RMSE/%：旧→修正','当前−历史偏差/百分点：修正'],[
        [energy,f"{ce(g,'old','mean')['bias_percent']:+.3f}→{ce(g,'corrected','mean')['bias_percent']:+.3f}",
         f"{ce(g,'old','individual')['rmse_percent']:.3f}→{ce(g,'corrected','individual')['rmse_percent']:.3f}",
         f"{next(r['mean_change'] for r in controls['paired_environment_shift'] if r['group']==g and r['method']=='corrected'):+.3f}"]
         for g,energy in enumerate((1.625,2.5,5.486,8.))])
    projection_table=table(['间隔/ns','平均投影：旧→修正','逐母MSE比及配对95%区间','冻结门'],[
        [int(r['gap_ns']),f"{r['old_mean']:+.2f}→{r['corrected_mean']:+.2f}",
         f"{r['mse_ratio']:.4f} [{r['mse_ratio_ci95'][0]:.4f}, {r['mse_ratio_ci95'][1]:.4f}]",
         ('通过' if r['passed'] else '未通过') if r['primary'] else '次要终点'] for r in summary['projection']])
    energy_table=table(['模型','波形','母数/观测数','偏差/%','RMSE/%','误差SD/%','E±10%','边界/未收敛'],[
        [r['method'],r['kind'],f"336/{r['rows']}",f"{r['bias_percent']:+.4f}",f"{r['rmse_percent']:.4f}",
         f"{r['error_sd_percent']:.4f}",f"{r['strict10_count']}/{r['rows']}",f"{r['boundary']}/{r['nonconverged']}"]
        for r in summary['energy']])
    variance_table=table(['间隔/ns','旧→修正均值平方','旧→修正逐母中心二阶矩'],[
        [int(r['gap_ns']),f"{r['old_mean']**2:.2f}→{r['corrected_mean']**2:.2f}",
         f"{r['old_mse']-r['old_mean']**2:.2f}→{r['corrected_mse']-r['corrected_mean']**2:.2f}"]
         for r in summary['projection']])
    audit=read(control.EXP/'native-audit.json')
    signal_error=max(r['python_signal_max_error_mV'] for r in audit['rows'])
    dt_error=max(r['half_dt_max_error_over_peak'] for r in audit['rows'])
    spans=[np.ptp([r['bias_percent'] for r in controls['energy'] if r['environment']=='current' and
                   r['group']==g and r['method']=='corrected' and r['kind']=='mean']) for g in range(4)]
    key=(f"当前环境四档孤立响应/读出检查通过后，另生成344个新种子单α ROOT，其中8个先导、336个正式强母。"
         f"冻结模型的8/40 ns逐母弱方向投影MSE比为{summary['projection'][0]['mse_ratio']:.4f}/"
         f"{summary['projection'][1]['mse_ratio']:.4f}，独立迁移检验{'通过' if gate['passed'] else '未通过'}全部预设门。")
    if gate['passed']:
        next_text=('下一步在当前环境新增其他源能区、采样相位和共享微胞受控双脉冲的独立响应检查。'
                   '本次强源单事件通过只支持固定模型在本条件下的迁移，仍不足以升级默认重建；'
                   '共享控制通过后才比较同输入、同权重和搜索预算的双脉冲能量拟合。')
    else:
        next_text=('本次新母数据未支持按既定要求推广该均值修正，停止双脉冲能量和高率扩展。'
                   '保留这一独立结果，后续将其作为已暴露的验证记录；若调整条件均值或权重，'
                   '应另建开发设计并保留新的确认样本，不能用本批结果选择幅度、平滑宽度或最佳间隔。'
                   '下一步优先分析冻结平均修正与母间涨落的关系，区分能标通过与弱方向误差改善。')
    failed='无' if gate['passed'] else '；'.join(gate['failures'])
    sections=[('123',f'''## 123. 当前环境的读出与孤立响应控制（2026-11-02）

{key}

先使用L95固定的16对历史/当前同种子光子记录，源能量1.625/2.5/5.486/8 MeV各4对。旧分支曾参与L93训练，当前分支复用了这些种子，因此这里只是已暴露样本的工程迁移检查。L93旧空间12×12响应和面积保持均值修正均原样冻结；没有重训或按结果调整器件参数。

32个原生ROOT的雪崩索引、微胞恢复、电荷守恒、原生采样及Python读出全部通过。Python从同一保存雪崩重算的最大信号差{signal_error:.3g} mV；内部时间步长0.25→0.125 ns的最大差占信号峰值{100*dt_error:.5f}%，低于既有0.2%检查要求。这些是确定性读出复算，不计为新的随机数字化。原生SiPM、离线SiPM及到达光子仍分开记录。

四相位0/0.5/1/1.5 ns，每条件8次独立数字化，共1024次正式离线数字化、2304次已知K=1/时刻的能量拟合。环境/相位间复用739000000+配对ID×100+重复ID种子，读出加10000000；不同光子输入会改变随机数消耗，不保证同一随机雪崩。只用完整含噪电压构造共同初值/GLS权重，真实E_dep仅用于评价。
'''),('124',f'''## 124. 四档能标和相位检查通过

{engineering_table}

上表为相位0；每能区/环境仅4母，32个单次观测不能视为32个独立输运样本。历史→当前偏差变化由相同母对与重复配对求得；完整4000次整母配对重采样的描述性95%区间在表格文件中保留，4母不足以建立总体统计等效性。最大相位偏差跨度{max(spans):.5f}个百分点。

预先沿用工程门：当前环境各能区/相位旧与修正均值E偏差≤5%，修正/旧单次E RMSE≤1.1，修正相位跨度≤1个百分点，未收敛0。全部通过。强母弱方向投影仍保留作小样本诊断，不用仅4个强母决定统计收益。通过只允许进入新母事件的固定模型检验，不代表均值修正已被独立确认。

图 0007-124-1：(a)四档相位0均值波形E偏差，虚线为历史环境、实线为当前环境，颜色区分旧/修正；(b)当前修正模型相位偏差相对相位0的变化。
'''),('125',f'''## 125. 冻结后新增336个正式强母

另建alpha_current_confirmation_20261102快照，原39份G4源码/配置、L93全部能量/均值/空间/噪声模型及投影算子逐项保留。运行环境为Geant4 11.3.2、ROOT 6.34.10、g++11.4.0；模型原开发数据来自11.4.2环境，所以此轮是当前环境迁移检验，不能写成旧环境内的重复确认。

先导ID 0–7不进入正式统计；正式ID 8–343共336母，固定8 MeV源能量、中央正入射、无背景。实际E_dep范围为[{summary['edep_range_MeV'][0]:.6f}, {summary['edep_range_MeV'][1]:.6f}] MeV，不等于源能量。G4种子738000000+ID，原生雪崩/读出分别加100000/200000；已与历史manifest中的种子流检查无重叠。全部样本保留，没有按Q、E或拟合结果删选或替换种子。

每ROOT一个20 μs窗，单α在5 μs发射；不是泊松注入率扫描。G4内部串行，最多8独立任务、启动间隔≥1.05 s，秒级文件名统一写g4/data。8个先导的完整读出核对、64次离线数字化和36次拟合只用于结构/成本检查，未用效果决定模型或门限。实际代码基于eaa8563和运行前未提交的新驱动，源码、配置、实际宏、RunInfo及哈希均保存在新快照，不依赖切换Git才能复查。

正式每母8次离线数字化，共2688次；种子738300000+ID×100+重复ID，读出另加10000000。保留全部记录面光子，首个有效面光子对齐128 ns，4096 ns数字化、前2048 ns/1024点观测，2 ns采样、0.5 ns孔径。本轮固定相位0，未新增其他能区、共享双脉冲或高率数据。
'''),('126',f'''## 126. 独立母事件的弱方向投影结果

{projection_table}

每母先对8次祖先信号求平均，再减去该母真实E_dep处的固定模型响应；以旧模型名义强7.5/弱0.25 MeV、8/40/200 ns间隔的logE影响算子投影。单位为100δlnE，MSE单位(100δlnE)²；这是前向残差诊断，不是实际双脉冲能量百分误差。

主要终点是336个逐母投影平方的平均之比，而非先把336母残差平均再平方。10000次按母配对重采样、种子738900000；旧/修正始终使用同一重采样索引。95%区间逐终点给出，条件于已冻结模型及当前环境，没有传播旧能标、空间先验和均值训练不确定度，也不是联合同时置信区间。

8/40 ns各要求MSE比<1、配对95%上界<1，且|修正均值|≤max(10,0.7|旧均值|)；200 ns为次要终点。门和样本量在新ROOT产生前冻结。总体失败项：{failed}。336来自历史开发方差给出的334粗估向上取整，不保证均值精度或检验通过。

图 0007-126-1：(a)逐母投影MSE比及配对母重采样95%区间；(b)旧/修正平均投影及各自95%区间。336为独立新母数量，重复与重采样不增加此数量。
'''),('127',f'''## 127. 孤立能量与完整波形

{energy_table}

每母8条单次含噪电压和1条8次均值波形，各拟合旧与修正两模型，共6048次正式拟合。相同观测、已知K=1/时刻、同电压导出的初值和固定GLS权重、相同边界与搜索预算。拟合进程禁止读取ROOT、mothers.json、真实分量/光子及原始RunInfo，逐条检查共同权重签名。真值只用于本节评分和前向残差图。

单次误差SD描述当前8 MeV源、不同E_dep母的相对误差分布，不是单能实测峰的FWHM。E±10%为给定K/t的条件能量计数，不能当作盲恢复效率；盲漏检、误拆分/误合并、能谱选择效率及高率性能均未测，不填为零。完整每母/每重复误差、边界和未收敛标志全部保留。

图 0007-127-1：(a)336母的完整2048 ns祖先信号均值及各自真实E处的模型均值；(b)相应平均残差；(c)全部2688条单次能量误差分布，保留所有观测；(d)每母8次均值的拟合E误差随真实E_dep变化。拟合使用含噪完整电压，(a–b)只作祖先真值辅助的响应诊断。
'''),('128',f'''## 128. 均值修正与母间涨落的区别

均值形状修正是从历史训练母估计的固定函数，只随候选能量和平移时刻变化；没有利用新母的真实光子或祖先残差逐事件回归。面积保持不保证沿弱能量敏感方向的误差变小，孤立E门与弱方向投影门也检验不同问题。

{variance_table}

以上为结果产生后的描述性分解：MSE=均值平方＋按N归一的中心二阶矩，单位均为(100δlnE)²。主要变化来自均值项，母间残差散布基本保持；这不是新增筛选门，也没有据此再调修正。孤立单次E RMSE由0.6568%变为0.6636%，通过原1.1倍容差，但没有能量精度改善证据。

图 0007-128-1：(a)8 ns、(b)40 ns，逐母旧投影与修正投影散点，虚线y=x。每点是一个正式新母的8次均值，336点均保留；图展示固定修正如何移动已有散布，不将散点相关性解释为弱脉冲可辨识能力。

历史L93/L94样本已参与开发且环境不同，与本轮不是配对效果估计。本轮能判断固定旧模型在本机新强母中的表现，不能仅凭它分离依赖版本变化、历史有限样本和原模型失配的因果贡献。两阶段结果不作为调整修正幅度的依据。
'''),('129',f'''## 129. 验证、完整计数和下一步

{next_text}

两项新增统计测试验证重复次数不扩增母样本量、整母配对重采样及能量失败/边界计数。四档阶段4次先导数字化与正式前两次数组一致，36次先导拟合与正式逐字段一致（仅计时除外）。新母阶段另重放正式ID 8/343各两次数字化，6个数组完全一致；ID8的18次拟合输出除计时外逐项相同。

四档阶段正式1024＋先导4=1028次离线随机数字化，正式2304＋先导36=2340次能量拟合。新母阶段正式2688＋先导64＋重放4=2756次离线随机数字化，正式6048＋先导36＋重放18=6102次拟合。合计3784次离线数字化、8442次拟合；原生344次SiPM数字化另列。32个既有ROOT及8个新先导的各两次确定性读出复算不增加随机数字化数量。

344个新ROOT的宏、三组种子、E_dep、状态、日志和SHA256在tables/root_catalog.csv及Run_Catalog逐条登记。32个复用控制输入也保留本地索引，没有重复计为新ROOT。当前门{'通过' if gate['passed'] else '未通过'}，失败原样保留，默认重建未升级。

本轮未新增外部论文检索、原PDF核查或文献性能主张，沿用L91保存的条件均值方法证据及人工审核边界。完整代码/冻结模型/源宏/运行器在两个新实验目录，原始ROOT在g4/data，图和数组在jupyter/artifacts对应目录，日志和旧笔记备份在jupyter/run_records。L96和研究计划第37节同步，原Notebook全部计算与图像保留。
''')]
    append_notebooks(sections,paths,key,next_text,inventory)
    write_documents(key,next_text,engineering_table,projection_table,energy_table,failed)
    outputs=[*paths.values(),EXP/'README.md',control.EXP/'README.md',ROOT/'README.md',ROOT/'docs/literature-use.md',
             *sorted((EXP/'tables').glob('*.csv')),*sorted((control.EXP/'tables').glob('*.csv')),
             *[ROOT/'jupyter'/name for name in ('codex_0007_SiPM_Energy_Calibration.ipynb','Research_Plan_20260921.ipynb','Run_Catalog.ipynb')]]
    save(EXP/'report-output-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in outputs})
    print('Wrote4figures,sections123–129,plan37,L96 and344-ROOT catalog',flush=True)


def append_notebooks(sections, paths, key, next_text, inventory):
    names=('codex_0007_SiPM_Energy_Calibration.ipynb','Research_Plan_20260921.ipynb','Run_Catalog.ipynb')
    notebooks={}
    for name in names:
        notebooks[name]=read(ROOT/'jupyter'/name)
        assert notebooks[name]==read(control.REC/(name+'.before')),name
    assert [len(notebooks[name]['cells']) for name in names]==[199,105,66]
    cells=[]
    for number,body in sections:
        cells.append(nbformat.v4.new_markdown_cell(body,id=f'current-confirmation-{number}'))
        if number in paths:
            source=("from pathlib import Path\nfrom IPython.display import display, Image\n"
                    "root = Path.cwd()\nif root.name == 'jupyter':\n    root = root.parent\n"
                    f"display(Image(filename=str(root / '{paths[number].relative_to(ROOT).as_posix()}')))\n")
            cells.append(nbformat.v4.new_code_cell(source,id=f'current-confirmation-figure-{number}'))
    new=nbformat.v4.new_notebook(cells=cells)
    NotebookClient(new,timeout=120,kernel_name='python3',resources={'metadata':{'path':str(ROOT)}}).execute()
    notebooks[names[0]]['cells'].extend(new['cells'])
    body='## 37. 当前环境四档控制与新强母固定模型确认（2026-11-02）\n\n'+key+'\n\n'+next_text
    body+='\n\n当前工程控制不是独立验证；新336正式强母与8先导严格分开。投影单位100δlnE，不是实际弱能量误差。已知K/t孤立拟合不能给出盲效率。详见codex_0007第123–129节、四幅固定编号图和L96。'
    notebooks[names[1]]['cells'].append(nbformat.v4.new_markdown_cell(body,id='current-confirmation-plan-20261102'))
    mapping=table(['ID/用途','ROOT','实际宏','G4/雪崩/读出种子','E_dep/MeV'],[
        [f"{r['id']}/{r['role']}",r['root'],r['macro'],f"{r['seed']}/{r['sipm_seed']}/{r['readout_seed']}",
         f"{r['edep_MeV']:.9f}"] for r in inventory])
    body=('## 当前环境四档响应控制与新强α独立母事件（2026-11-02）\n\n'+key+'\n\n'+
          '四档控制复用alpha_transport_environment_20261102中的16个新环境对照及其16个历史原ROOT；32个文件原样保留，未生成新ROOT。'
          '原ROOT/宏/种子/E_dep映射在alpha_current_response_20261102/tables/root_catalog.csv，本地original-runinfo和original-macros保留实际来源。'
          '正式1024/先导4次离线数字化，正式2304/先导36次能量拟合。32个原生雪崩记录各做原时间步长及半步长的确定性读出核对。\n\n'
          '新实验alpha_current_confirmation_20261102：344个ROOT全部完成，无失败或替换；8先导、336正式。每个固定8 MeV单α源在20 μs窗的5 μs中央正入射，'
          '无背景，不对应泊松注入率。三组种子逐条列于下表，原始RunInfo及实际宏随实验保留。G4串行、最多8进程、启动≥1.05 s，ROOT使用秒级文件名。'
          '环境Geant4 11.3.2、ROOT6.34.10、g++11.4.0，NumPy1.26.4/SciPy1.15.3，实际源码基于eaa8563及当时未提交驱动，design-freeze.json保存源码/配置哈希。'
          '二进制复用此前已验证的当前构建，物理源码39文件未变。\n\n'+mapping+'\n\n'
          '新母离线数字化正式2688、先导64、重放4；能量拟合正式6048、先导36、重放18。'
          '数字化4096 ns、观测前2048 ns；每母8重复，首有效光子对齐128 ns、相位0。种子738300000+ID×100+重复ID，读出加10000000。'
          '8先导另做原生读出/电荷检查，不进入正式评价。全部ROOT实际宏、日志、SHA256见tables/root_catalog.csv；'
          '输出数组与图在jupyter/artifacts/alpha_current_confirmation_20261102，日志在jupyter/run_records同名目录。\n\n'+next_text)
    notebooks[names[2]]['cells'].append(nbformat.v4.new_markdown_cell(body,id='current-confirmation-runs-20261102'))
    for name,nb in notebooks.items():
        nbformat.validate(nbformat.from_dict(nb));(ROOT/'jupyter'/name).write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n')


def write_documents(key, next_text, engineering_table, projection_table, energy_table, failed):
    response_readme=f'''# 当前环境四档孤立响应与读出检查

32个既有ROOT（16对历史/当前同种子）全部通过原生雪崩、电荷、Python读出和半时间步长核对。四档/四相位工程能量门通过，固定L93旧模型与均值修正，没有训练或调整参数。

{engineering_table}

每档/环境仅4母；同母8次数字化不是独立输运。旧分支参与过模型训练，当前分支也复用旧种子，因此本阶段不是独立确认或总体等效性检验。源能量与真实E_dep分别记录。已知K=1/时刻，拟合只读观测电压；真值用于评分。

正式1024次离线数字化、2304次拟合；4次数字化先导及36次拟合先导独立记账并与正式对应结果核对。32个原生记录的确定性读出复算不计新随机数字化，无新ROOT。actual模型和分析源码在source，原始ROOT/宏映射在tables/root_catalog.csv，RunInfo在original-runinfo，原模型/输入/设计哈希单独冻结。

本阶段通过后启动了alpha_current_confirmation_20261102的新强母迁移检验；其独立结果不能由本控制通过代替。两阶段合并报告保存在本地codex_0007第123–129节和L96。

项目根目录验证本阶段输入和冻结源码：

```bash
python3 -B pileup/experiments/alpha_current_response_20261102/run.py verify
python3 -B pileup/experiments/alpha_current_confirmation_20261102/validate.py
```

run.py导入本目录实际源码副本。原模型NPZ、ROOT、二进制、数组及Notebook保留在本地，不随Git发布。新统计使用新输出目录；既有stage已完成时拒绝覆盖。图和数组在jupyter/artifacts/alpha_current_response_20261102，日志及旧笔记备份在jupyter/run_records同名目录。
'''
    (control.EXP/'README.md').write_text(response_readme)
    confirmation_readme=f'''# 当前环境独立强母的固定均值模型迁移检验

{key}

{projection_table}

主要终点为逐母投影MSE，而非合并均值的平方。投影单位100δlnE，不是实际双脉冲能量百分误差；200 ns为次要终点。10000次整母配对重采样，区间条件于冻结模型、当前环境和正式样本。冻结门失败项：{failed}。

{energy_table}

已知K=1/时刻，单次2688观测只有336个母；均值波形是同母8次平均。误差SD不是实测峰FWHM，E±10%是条件能量计数，不是盲检出或完整重建效率。没有双脉冲能量拟合或默认替换。

样本：8个先导ID0–7，336正式ID8–343。固定8 MeV源，中央正入射，单α在20 μs窗的5 μs发射。G4种子738000000+ID，原生雪崩/读出加100000/200000；离线种子738300000+ID×100+重复ID，读出再加10000000。所有光子保留，首有效面光子对齐128 ns，相位0；4096 ns数字化、前2048 ns观测。独立新母检验了旧模型向Geant4 11.3.2/ROOT6.34.10当前环境的迁移，不能写成11.4.2旧环境内的再次确认。

全部344个ROOT成功、无替换或筛选。G4串行、至多8任务、秒级文件名、启动≥1.05 s；实际源码、宏、种子、RunInfo、模型/输入/输出哈希保存于source、macros、original-runinfo与各freeze.json。tables/root_catalog.csv逐条映射ROOT及实际宏，ROOT统一保存在g4/data。原生雪崩/读出结果与离线随机数字化分开记账。

本阶段正式2688＋先导64＋重放4=2756次离线数字化；正式6048＋先导36＋重放18=6102次能量拟合。两条正式母的6个数组精确重放，18次拟合除计时外相同。完整图和数组在jupyter/artifacts/alpha_current_confirmation_20261102，日志在jupyter/run_records同名目录，Notebook本地第123–129节，L96记录研究范围。

项目根目录执行：

```bash
python3 -B pileup/experiments/alpha_current_confirmation_20261102/run.py verify
python3 -B pileup/experiments/alpha_current_confirmation_20261102/validate.py
```

validate.py只读核对既有结果，不重跑生产。run.py直接加载本实验的source副本，原生二进制与年龄模型库来自本机验证构建，build.json记录来源和编译命令。新机器可用source/g4独立配置编译，再编译source/pileup/analysis中的两个年龄模型C++文件；外部数据集及环境在transport-environment.json/analysis-environment.json。二进制、模型NPZ、ROOT和Notebook不随Git分发，须另行提供；重算应另建实验/输出，不能改写这份已冻结的记录。

{next_text}
'''
    (EXP/'README.md').write_text(confirmation_readme)
    literature=ROOT/'docs/literature-use.md';before=(control.REC/'literature-use.md.before').read_text()
    assert literature.read_text()==before and '## L96 ' not in before
    literature.write_text(before+f'''

<a id="current-confirmation-20261102"></a>
## L96 — 2026-11-02：当前环境四档工程检查与新强母独立迁移检验

**项目结果。** {key}

{projection_table}

**先决控制。** 复用16对历史/当前同种子ROOT、四源能区各4对，旧分支已有L93训练暴露。32个原生雪崩/电荷及Python读出、半时间步长核对全部通过。四相位×8重复=1024次正式离线数字化，2304次能量拟合；当前各档/相位均值E偏差≤5%、修正/旧单次RMSE≤1.1、相位跨度≤1个百分点、未收敛0，工程门通过。该控制不是独立验证或总体等效性证明。

{engineering_table}

**独立输入与冻结条件。** 在控制通过后固定模型和门限，另生成8个先导、336个正式单α光学母；G4种子738000000+ID，已与历史manifest种子流去重。源能量8 MeV，E_dep逐母记录；中央正入射，20 μs窗内5 μs一个α，无背景和泊松注入率。全部344个ROOT保留，无Q/E筛选、无替换种子；先导ID0–7不进入正式统计。39份物理源码、95/351 ns快慢常数及0.3快比例、29 ns恢复代理、14400微胞、旧空间12×12先验和L93面积保持形状修正均原样。实际环境Geant4 11.3.2/ROOT6.34.10，模型开发数据来自11.4.2，所以是跨环境的冻结模型迁移验证。

**主要终点。** 每母8次祖先信号均值减去真实E处预测，使用固定旧名义强7.5/弱0.25 MeV的logE影响方向。8/40 ns要求逐母投影MSE比<1且配对母bootstrap95%上界<1，并要求|修正均值|≤max(10,0.7|旧均值|)；200 ns次要。10000次按母配对重采样、种子738900000；区间逐终点且条件于冻结模型，不传播旧模型/均值训练误差。投影单位100δlnE、MSE单位(100δlnE)²，比值无量纲；不是实际双脉冲E百分误差。失败项：{failed}。

**孤立能量。** 正式2688单次波形及336条同母8次平均，各比较旧/修正模型；6048次正式拟合均给定K=1/时刻，相同电压初值/GLS权重/搜索预算，拟合进程禁止读取ROOT、真实E和祖先/光子。真值只用于评分。

{energy_table}

误差SD不是实测单能峰FWHM，±10%计数不代表盲效率。盲漏检、误拆分/误合并、能谱选择效率、高率及真实器件性能尚未评价。正式样本只覆盖强源和相位0，不把工程控制4母/能区当作新独立的全能区证明。

**验证和计数。** 两项统计测试、四档成本先导与正式对应数组/拟合一致；新母阶段重放2母各2次，6数组精确相同，另18拟合除计时外一致。两阶段合计3784次离线数字化、8442次能量拟合，344次原生SiPM另列。40个原生记录的确定性读出交叉检查不增加随机数字化数。全部新ROOT的单窗、种子、实际宏、RunInfo和哈希核对，最多8任务、启动≥1.05 s。旧Notebook全部计算和图像保留。

**采用与下一步。** {next_text}

**文献层级。** 沿用L91据既有精读笔记保存的Fixsen条件均值思想及Busch/Zhang适用范围，没有新联网检索、原PDF核查、论文性能主张或人工审核变更。样本量、冻结门和重采样是项目设计；文献不保证本GAGG/SiPM模型或当前环境达到这些门限。

位置：`pileup/experiments/alpha_current_response_20261102/`与`pileup/experiments/alpha_current_confirmation_20261102/`。模块`pileup/analysis/alpha_current_response.py`、两个study及`alpha_current_confirmation_report.py`。本地codex_0007第123–129节、图0007-124-1/126-1/127-1/128-1，研究计划第37节、Run_Catalog与ROOT逐条清单同步。
''')
    readme=ROOT/'README.md';before=(control.REC/'README.md.before').read_text();assert readme.read_text()==before
    readme.write_text(before+f'\n\n2026-11-02完成[当前环境的四档控制与新强母迁移检验](pileup/experiments/{study.NAME}/README.md)：{key} '
        '336正式强母与8先导分开，保持L93模型、物理参数和预设门不变；投影不是实际双脉冲能量误差。'
        f'全部344个ROOT保留，原生/离线读出核对及局部精确重放通过。{next_text} '
        '结果与四幅编号图见本地codex_0007第123–129节及[采用记录L96](docs/literature-use.md#current-confirmation-20261102)。\n')


if __name__=='__main__':main()
