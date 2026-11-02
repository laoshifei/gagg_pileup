"""Report the frozen residual diagnosis while preserving prior Notebook cells."""
from __future__ import annotations
import csv
import json
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/gagg_energy_diagnosis_mpl')
os.environ.setdefault('IPYTHONDIR','/tmp/gagg_energy_diagnosis_ipython')
os.environ.setdefault('JUPYTER_RUNTIME_DIR','/tmp/gagg_energy_diagnosis_jupyter')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nbformat
from nbclient import NotebookClient
import numpy as np
from . import alpha_current_energy_diagnosis_study as study
from .alpha_current_energy_diagnosis import TERM_NAMES

ROOT,EXP,ART,REC=study.ROOT,study.EXP,study.ART,study.REC
read,save,sha=study.read,study.save,study.sha
CONDITIONS=[('high_weak',8.),('high_weak',40.),('high_weak',200.),('high_weak',1200.),('weak_high',40.)]
LABELS=['H-W 8','H-W 40','H-W 200','H-W 1200','W-H 40']
TERM_LABELS=['Strong isolated','Weak isolated','Shared difference','Background/readout']


def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(map(str,row))+' |' for row in rows])


def csv_table(name,rows):
    (EXP/'tables').mkdir(exist_ok=True)
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (EXP/'tables'/f'{name}.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,lineterminator='\n');writer.writeheader()
        writer.writerows({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in rows)


def select(summary,method,kind):
    return [next(r for r in summary['groups'] if (r['method'],r['kind'],r['condition'],r['gap_ns'])==
                 (method,kind,condition,gap)) for condition,gap in CONDITIONS]


def figures(summary,details):
    plt.rcParams.update({'font.size':10,'figure.dpi':140,'axes.grid':True,'grid.alpha':.2})
    paths={};x=np.arange(5)
    fig,axes=plt.subplots(2,2,figsize=(12,7),constrained_layout=True)
    for method,color,offset in [('old','tab:blue',-.12),('corrected','tab:orange',.12)]:
        rows=select(summary,method,'individual')
        for axis,key,ylabel in [(axes[0,0],'median_cosine','Median derivative cosine'),
             (axes[0,1],'median_weak_retained_fraction','Median weak norm retained'),
             (axes[1,0],'median_normalized_condition','Median normalized condition'),
             (axes[1,1],'residual_rms_mV','RMS residual at true E / mV')]:
            axis.bar(x+offset,[r[key] for r in rows],width=.24,color=color,label=method)
            axis.set(xticks=x,xticklabels=LABELS,xlabel='Order / gap (ns)',ylabel=ylabel)
    for i,axis in enumerate(axes.ravel()):axis.set_title(f'({chr(97+i)}) 8 mother pairs; individual metric')
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Figure 0007-143-1: True-energy local geometry in the inherited working metric')
    p=ART/'figures/fig_0007_143_1_geometry.png';fig.savefig(p);plt.close(fig);paths['143']=p

    fig,axes=plt.subplots(1,2,figsize=(13,4.8),sharey=True,constrained_layout=True)
    for j,(axis,method) in enumerate(zip(axes,('old','corrected'))):
        rows=select(summary,method,'mean');values=np.array([r['weak_terms_mean'] for r in rows])
        for k,label in enumerate(TERM_LABELS):axis.bar(x+(k-1.5)*.17,values[:,k],width=.17,label=label,color=f'C{k}')
        axis.scatter(x,[r['weak_projection_mean'] for r in rows],color='black',marker='D',s=25,label='Total',zorder=3)
        axis.axhline(0,color='grey',lw=.6)
        axis.set(xticks=x,xticklabels=LABELS,xlabel='Order / gap (ns)',ylabel=r'Mean weak projection / $100\delta\ln E$',
                 title=f'({chr(97+j)}) {method}; 8-repeat mean metric')
    axes[0].legend(fontsize=7,ncol=2)
    fig.suptitle('Figure 0007-144-1: Exact residual partition; local projections are not energy errors')
    p=ART/'figures/fig_0007_144_1_terms.png';fig.savefig(p);plt.close(fig);paths['144']=p

    matrices=[]
    for gap in (8.,40.):
        rows=[r for r in details if (r['method'],r['kind'],r['condition'],r['gap_ns'])==('corrected','mean','high_weak',gap)]
        matrices.append(np.mean([r['repeat_covariance']['covariance'] for r in rows],axis=0))
    limit=max(abs(m).max() for m in matrices)
    fig,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
    ticks=['Strong iso','Weak iso','Shared diff','Background']
    for j,(axis,matrix,gap) in enumerate(zip(axes,matrices,(8,40))):
        view=axis.imshow(matrix,cmap='coolwarm',vmin=-limit,vmax=limit)
        for row in range(4):
            for col in range(4):axis.text(col,row,f'{matrix[row,col]:.1f}',ha='center',va='center',fontsize=8)
        axis.grid(False);axis.set(xticks=range(4),yticks=range(4),xticklabels=ticks,yticklabels=ticks,
                                  title=f'({chr(97+j)}) H-W {gap} ns; corrected')
        axis.tick_params(axis='x',labelrotation=25)
    fig.colorbar(view,ax=axes,label=r'Within-mother covariance / $(100\delta\ln E)^2$',shrink=.8)
    fig.suptitle('Figure 0007-145-1: Average repeat covariance; cross terms must be retained')
    p=ART/'figures/fig_0007_145_1_covariance.png';fig.savefig(p);plt.close(fig);paths['145']=p

    fig,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
    for j,(axis,gap) in enumerate(zip(axes,(8.,40.))):
        allx=[];ally=[]
        for method,color in [('old','tab:blue'),('corrected','tab:orange')]:
            for kind,marker in [('individual','o'),('mean','x')]:
                rows=[r for r in details if (r['method'],r['kind'],r['condition'],r['gap_ns'])==(method,kind,'high_weak',gap)]
                xv=[r['projection_total'][1] for r in rows];yv=[r['reference_fit']['energy_log_error'][1] for r in rows]
                allx.extend(xv);ally.extend(yv);axis.scatter(xv,yv,color=color,marker=marker,s=28,label=method+' '+kind)
        limits=[min(allx+ally),max(allx+ally)]
        axis.plot(limits,limits,'k--',lw=.8);axis.axhline(0,color='grey',lw=.5);axis.axvline(0,color='grey',lw=.5)
        axis.set(xlabel=r'Local weak projection / $100\delta\ln E$',ylabel=r'Saved fit error / $100\ln(E_{fit}/E_{true})$',
                 title=f'({chr(97+j)}) H-W {gap:g} ns; same 8 mother pairs')
        axis.legend(fontsize=7)
    fig.suptitle('Figure 0007-146-1: Local linearization compared with saved nonlinear energy fits')
    p=ART/'figures/fig_0007_146_1_linearity.png';fig.savefig(p);plt.close(fig);paths['146']=p
    return paths


def main():
    study.verify();assert read(EXP/'replay.json')['passed'] and read(EXP/'reproduction-check.json')['passed']
    summary=read(EXP/'result-summary.json');details=read(EXP/'details.json');paths=figures(summary,details)
    csv_table('groups',summary['groups']);csv_table('details',details)
    mothers=read(EXP/'input-mothers.json');csv_table('reused_roots',mothers)
    geometry_table=table(['条件/ns','模型','导数余弦中位数','弱范数保留/%','归一化条件数中位数','真值处残差RMS/mV'],[
        [f"{r['condition']}/{r['gap_ns']:g}",r['method'],f"{r['median_cosine']:.4f}",
         f"{100*r['median_weak_retained_fraction']:.2f}",f"{r['median_normalized_condition']:.3f}",f"{r['residual_rms_mV']:.3f}"]
        for r in summary['groups'] if r['kind']=='individual'])
    terms_table=table(['条件/ns','模型','强孤立','弱孤立','共享差','背景/读出','总投影'],[
        [f"{r['condition']}/{r['gap_ns']:g}",r['method'],*[f'{v:+.2f}' for v in r['weak_terms_mean']],f"{r['weak_projection_mean']:+.2f}"]
        for r in summary['groups'] if r['kind']=='mean'])
    covariance_table=table(['条件/ns','模型','重复总方差','仅对角和','交叉项和','均值MC标准误中位数'],[
        [f"{r['condition']}/{r['gap_ns']:g}",r['method'],f"{r['mean_repeat_total_variance']:.2f}",
         f"{r['mean_repeat_diagonal_variance']:.2f}",f"{r['mean_repeat_cross_variance']:+.2f}",f"{r['median_mean_mc_se']:.2f}"]
        for r in summary['groups'] if r['kind']=='mean'])
    linear_table=table(['条件/ns','波形','模型','局部弱投影均值','实际弱log误差均值','二者差RMSE','边界/未收敛'],[
        [f"{r['condition']}/{r['gap_ns']:g}",r['kind'],r['method'],f"{r['weak_projection_mean']:+.2f}",
         f"{r['weak_fitted_log_mean']:+.2f}",f"{r['linear_minus_fit_rmse']:.2f}",f"{r['boundary']}/{r['nonconverged']}"]
        for r in summary['groups'] if r['reference_available']])
    m8,m40=select(summary,'corrected','mean')[:2]
    interval_table=table(['间隔/ns','项','均值','整母95%描述区间'],[
        [int(r['gap_ns']),name,f'{value:+.2f}',f'[{ci[0]:+.2f}, {ci[1]:+.2f}]']
        for r in (m8,m40) for name,value,ci in zip(('强孤立','弱孤立','共享差','背景/读出'),r['weak_terms_mean'],r['weak_terms_ci95'])])
    key=(f"复用L98的8对母、40条件完成160条真值辅助局部诊断；112个原拟合点的工作代价均低于真值点。"
         f"修正模型8 ns的弱方向范数保留中位数{100*m8['median_weak_retained_fraction']:.1f}%，40 ns为{100*m40['median_weak_retained_fraction']:.1f}%。"
         f"8 ns八次均值的弱投影由强孤立{m8['weak_terms_mean'][0]:+.2f}、弱孤立{m8['weak_terms_mean'][1]:+.2f}、"
         f"共享差{m8['weak_terms_mean'][2]:+.2f}、背景/读出{m8['weak_terms_mean'][3]:+.2f}组成，单位均为100δlnE。")
    next_text=('下一步先在已暴露的开发母上固定增加数字化重复，检查弱方向投影随重复数的收敛，'
               '区分有限数字化波动与跨母光学/模板差异，再提出模型修改。此8对仅作失败诊断，不用于选择修正参数；'
               '任何后续模型改变另建快照并另留确认母，当前默认模型和L98失败判定保持。')
    counts=('40个既有双脉冲条件×旧/修正×单次/八次均值，共160条局部诊断；320份既有共享观测重复。'
            '112个既有能量拟合只重算代价，另2条件诊断精确重放、1条件在新复制目录重跑。'
            '新能量拟合、ROOT、随机数字化和训练均为0；120个输入数组文件按字节复制，追溯16个既有ROOT。')
    checks=summary['checks']
    validation=(f"残差闭合最大差{checks['max_residual_closure_mV']:.3g} mV，投影闭合差{checks['max_projection_closure']:.3g}；"
                f"数值导数最大相对差{checks['max_gradient_relative_error']:.3g}，与L97均值最大差{checks['max_mean_identity_mV']:.3g} mV，"
                f"112次代价重算最大相对差{checks['max_cost_replay_relative_error']:.3g}。三项代数/协方差测试通过；"
                '重放与新目录重跑结果除计时外逐字段及数组一致，新目录入口拒绝覆盖。')
    interpretation=(f"8 ns均值强孤立项{m8['weak_terms_mean'][0]:+.2f}是本批最大的带符号组均值投影，"
        f"整母描述区间[{m8['weak_terms_ci95'][0][0]:+.2f}, {m8['weak_terms_ci95'][0][1]:+.2f}]很宽；"
        '它混合了单母光学实现、有限雪崩重复及固定模型差异，不能全归为确定模板偏差。'
        f"同一8 ns条件的均值MC标准误中位数仍为{m8['median_mean_mc_se']:.2f}，40 ns为{m40['median_mean_mc_se']:.2f}。"
        '背景/读出项的组均值小，不等于全部随机噪声小：光子来源雪崩的统计涨落已包含在前三项。')
    sections=[('142',f'''## 142. 已暴露失败样本的残差诊断：设计与推导（2026-11-02）

{key}

本轮只诊断L98已暴露的8对母，不使用保留母或做模型选择。五条件、相位0和全部8次数字化重复固定；两种模型、物理参数、能标、共享恢复项、原112个拟合及其共同权重保持。对所有40条件计算单次/均值的局部投影；200/1200 ns及反向40 ns的均值没有实际能量拟合，不能补造能量结果。

记y为观测、b为冻结背景，S/I为共享/孤立的祖先信号，μ/μ_iso为真E_dep处的联合/孤立预测。本项目直接推导恒等式：r=y−b−Σμ = (I_s−μ_iso,s)+(I_w−μ_iso,w)+[Σ(S−I)−Σ(μ−μ_iso)]+(y−ΣS−b)。四项依次为强孤立、弱孤立、共享差、背景/读出。每个重复均闭合，允许负投影和抵消。

在每条观测的固定工作度量中，D=W·T·∂μ/∂lnE；以SVD计算δlnE=D⁺W·T·r。弱方向先消去强导数投影，保留比例=‖D_w−proj(D_s)D_w‖/‖D_w‖。T沿用原FIR变换，W沿用电压导出的共同GLS；没有在真E上重算权重。真E只用于本轮明确的局部诊断，结果不进入盲推断。100δlnE不是100(E_fit/E_true−1)，大位移下不能互换。

文献范围沿用第136节及L98的Fixsen2004条件均值/协方差原则和Zhang2024模板适用条件；原TES与LaBr3/PMT条件、AI草稿与审核边界保持。本轮未重新检索文献、查PDF或核DOI。上述残差分解与SVD计算是项目推导，不冒充论文给出的GAGG性能结论。
'''),('143',f'''## 143. 8 ns条件的局部能量方向

{geometry_table}

表中权重来自重复0；每个条件8对母，两模型权重完全相同。消去强参数后，弱方向的可区分分量减小；25%范数保留不是25%探测效率，也不是经实验标定的误差下界。归一化条件数剔除两列自身幅度比例，未归一化量和原导数保存在结果数组。没有判定严格秩亏或普适8 ns分辨极限。

图 0007-143-1：(a)导数余弦；(b)弱范数保留；(c)归一化条件数；(d)真值处电压残差RMS。前向残差小不自动意味着弱能量误差小。
'''),('144',f'''## 144. 残差沿弱能量方向的来源

以下均为8重复均值、8对母的组均值投影，单位100δlnE；不是能量误差或独立因果效应。

{terms_table}

{interval_table}

{interpretation}

10000次整母对重采样，种子742100000；同一抽样索引跨全部条件、方法与各项使用。只有8个已暴露母对，区间条件于原模型和8次已观察重复，不传播模型训练误差，不宣称同时覆盖或独立确认。

图 0007-144-1：(a)旧模型、(b)修正模型的四项均值投影及总和；带符号展示，不把各项绝对值归一化成因果百分比。
'''),('145',f'''## 145. 有限重复与相关项不能省略

共享S和孤立I来自不同数字化，虽使用相同种子，但事件排序、占用及随机消耗不同；S−I不是逐次无噪声恢复效应。强/弱孤立误差与共享差包含共用量，存在相关或抵消。每母用同一个均值工作度量投影其8次重复，完整协方差Σ满足Var(Σ项)=1ᵀΣ1；不能只加对角方差。

{covariance_table}

方差单位(100δlnE)²；标准误单位100δlnE。表中方差先逐母计算再取均值，MC标准误先逐母取sqrt(1ᵀΣ1/8)再取中位数，所以二者不能直接开方对应。该标准误只描述固定母光学输入和固定度量下的有限数字化重复，不是总体能量分辨率；均值工作协方差未重新标定。重复0与八次均值相关，同母在不同条件也相关。

图 0007-145-1：(a)8 ns、(b)40 ns修正模型的4×4重复协方差逐母平均，含全部非对角项；颜色范围一致。
'''),('146',f'''## 146. 局部投影与原非线性能量拟合的差异

{linear_table}

局部预测和实际误差均以log能量坐标比较，单位100δlnE或100ln(E_fit/E_true)。原能量百分误差仍见第137/140节；二者的组均值不通过指数直接互换。8 ns位移大且部分解触边界，局部近似明显失效；它能定位敏感方向，不能代替真实非线性拟合或精确误差归因。

112个已保存拟合的电压初值、三起点及度量哈希均复核；代价重新计算一致，全部比对应真值E点低。说明本批观测在冻结工作目标下可偏好其他E，不证明全局最优，也未排除优化器的其他局部解。没有新增拟合或把真实E作为拟合起点。

图 0007-146-1：(a)8 ns、(b)40 ns逐母局部弱投影与原拟合log误差；旧/修正与单次/均值完整保留，虚线为相等。
'''),('147',f'''## 147. 复核与下一步：先检查条件均值的收敛

{counts}

{validation}

{next_text}

本轮没有新的盲漏检、误拆分/误合并、能谱效率、背景实验、高率或实物结果，未测项不填零。L99、研究计划第40节及Run_Catalog同步；源码、宏、RunInfo、模型与输入哈希保存在pileup/experiments/alpha_current_energy_diagnosis_20261102。既有Notebook计算和图像保留，当前图号143-1至146-1固定。
''')]
    append_notebooks(sections,paths,key,next_text,counts,mothers)
    documents(key,next_text,counts,validation,interpretation,geometry_table,terms_table,interval_table,covariance_table,linear_table)
    outputs=[*paths.values(),EXP/'README.md',ROOT/'README.md',ROOT/'docs/literature-use.md',*sorted((EXP/'tables').glob('*.csv')),
        *[ROOT/'jupyter'/n for n in ('codex_0007_SiPM_Energy_Calibration.ipynb','Research_Plan_20260921.ipynb','Run_Catalog.ipynb')]]
    save(EXP/'report-output-freeze.json',{str(p.relative_to(ROOT)):sha(p) for p in outputs})
    print('Wrote4figures,sections142–147,plan40,L99 and16reusedROOT mappings',flush=True)


def append_notebooks(sections,paths,key,next_text,counts,mothers):
    names=('codex_0007_SiPM_Energy_Calibration.ipynb','Research_Plan_20260921.ipynb','Run_Catalog.ipynb')
    notebooks={n:read(ROOT/'jupyter'/n) for n in names}
    for name in names:assert notebooks[name]==read(REC/(name+'.before'))
    assert [len(notebooks[n]['cells']) for n in names]==[230,108,69]
    cells=[]
    for number,body in sections:
        cells.append(nbformat.v4.new_markdown_cell(body,id='current-energy-diagnosis-'+number))
        if number in paths:
            code=("from pathlib import Path\nfrom IPython.display import display, Image\nroot = Path.cwd()\n"
                  "if root.name == 'jupyter':\n    root = root.parent\n"
                  f"display(Image(filename=str(root / '{paths[number].relative_to(ROOT).as_posix()}')))\n")
            cells.append(nbformat.v4.new_code_cell(code,id='current-energy-diagnosis-figure-'+number))
    new=nbformat.v4.new_notebook(cells=cells)
    NotebookClient(new,timeout=120,kernel_name='python3',resources={'metadata':{'path':str(ROOT)}}).execute()
    notebooks[names[0]]['cells'].extend(new['cells'])
    body='## 40. 短间隔能量失败的真值辅助诊断（2026-11-02）\n\n'+key+'\n\n'+next_text+'\n\n'+counts+(
        '\n\n8次重复仍不足以把全部残差归为确定模型偏差；8 ns线性外推误差明显。投影单位100δlnE，不能当实际能量百分误差或盲效率。'
        '详见codex_0007第142–147节、四幅固定图与L99。')
    notebooks[names[1]]['cells'].append(nbformat.v4.new_markdown_cell(body,id='current-energy-diagnosis-plan'))
    mapping=table(['母ID','既有ROOT','实际宏副本','G4/雪崩/读出种子','源/沉积能量MeV'],[
        [r['id'],r['root'],r['diagnostic_macro'],f"{r['seed']}/{r['sipm_seed']}/{r['readout_seed']}",
         f"{r['source_energy_MeV']:.6f}/{r['edep_MeV']:.9f}"] for r in mothers])
    body='## 短间隔能量诊断：复用16个ROOT，无新运行（2026-11-02）\n\n'+key+'\n\n'+mapping+'\n\n'+counts+(
        '\n\n原生G4 11.3.2/ROOT6.34.10、20μs单α窗、5μs中央正入射条件沿用原RunInfo。离线共享输入原为4096ns数字化的前2048ns，'
        '1024点、2ns采样/0.5ns孔径、相位0；原随机种子740500000+block×1000+重复，读出加10000000。'
        '本轮只读已有数组、真E及原拟合，不调用G4或随机数字化。原源码、模型、二进制和实际宏均有副本/哈希；'
        '诊断以bb130af和新冻结源码为基线，最多8进程、BLAS单线程。\n\n'+next_text)
    notebooks[names[2]]['cells'].append(nbformat.v4.new_markdown_cell(body,id='current-energy-diagnosis-runs'))
    for name,nb in notebooks.items():
        nbformat.validate(nbformat.from_dict(nb));(ROOT/'jupyter'/name).write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n')


def documents(key,next_text,counts,validation,interpretation,geometry,terms,intervals,covariance,linearity):
    (EXP/'README.md').write_text(f'''# 当前环境短间隔能量失败的真值辅助诊断

{key}

{interpretation}

这是对L98已暴露样本的事后诊断，沿用8对母、五条件、相位0；保留母未被读取，模型与参数未调整。局部量单位100δlnE，实际比较量100ln(E_fit/E_true)，均不能替代能量百分误差或盲重建指标。

## 残差与局部方向

逐重复恒等分解：观测−冻结背景−真E处联合均值 = 强孤立残差 + 弱孤立残差 + 共享减孤立的模型差 + 背景/读出残差。共享与孤立来自不同数字化，共用种子不保证同一随机实现；共享差不是无噪声因果效应。

D为经原FIR变换T及固定GLS变换W作用后的dμ/dlnE，δlnE=D⁺WTr；用SVD求逆，并与消去强方向后的弱投影核对。权重仍由观测电压导出，没有用真E改权重。八次均值也沿用原工作度量，没有另标定均值协方差。

{geometry}

## 八次均值的弱方向分解

下表单位100δlnE；所有条件完整保留。

{terms}

{intervals}

10000次整母对配对重采样，seed742100000；跨条件、模型和各项共用索引。8个已暴露母的区间只是条件描述，没有传播模型训练误差或提供同时覆盖。

## 重复协方差与非线性范围

{covariance}

方差单位(100δlnE)²，MC标准误单位100δlnE。保留全部非对角项；逐母8重复总方差除8只估计有限数字化均值的波动，不包括跨母光学变化或系统模型误差。统计涨落也在前三个祖先信号项中，不只在背景/读出项。

{linearity}

112个原拟合点重算代价均低于对应真值点；原度量哈希、初值和起点逐项一致。有限搜索仍不证明全局最优。8 ns大位移和边界使局部外推明显失效；不能把局部拆分解释为精确非线性误差归因。其他三个条件的均值只有投影，没有新增能量拟合。

## 计数、验证与复现

{counts}

{validation}

项目根目录运行`python3 -B pileup/experiments/{study.NAME}/run.py verify`检查冻结输入，`python3 -B pileup/experiments/{study.NAME}/validate.py`核对完整本地结果。重跑先执行`python3 -B pileup/experiments/{study.NAME}/reproduce.py --target-root /path/to/new-project`，再在新目录按pilot→run→score→replay运行run.py。新目录必须不存在；完整阶段拒绝覆盖。新复制目录的首条件已实测逐字段/数组一致，不依赖切换Git，也不需要原ROOT参与诊断计算。

计算使用NumPy1.26.4/SciPy1.15.3与旧age.so；未改C++或物理模型、未重建G4。120个本地NPZ及模型/二进制不随Git上传，跨机器需另提供；ROOT宏、RunInfo和源码快照随实验保留。完整本地验证还会核对原16个ROOT，Notebook报告需要已有本地笔记上下文。

四图与详细数组在jupyter/artifacts/{study.NAME}，codex_0007第142–147节、研究计划第40节和Run_Catalog均同步。复用L98文献证据边界，没有新PDF/DOI核验；恒等分解和SVD为项目推导。

{next_text}
''')
    literature=ROOT/'docs/literature-use.md';before=(REC/'literature-use.md.before').read_text()
    assert literature.read_text()==before and '## L99 ' not in before
    literature.write_text(before+f'''

<a id="current-energy-diagnosis-20261102"></a>
## L99 — 2026-11-02：短间隔能量失败的残差与局部方向诊断

**项目结果。** {key}

{interpretation}

**证据和推导范围。** 复用L98及codex_0007第136节已经核过的Fixsen2004与Zhang2024来源笔记、条件和审核状态，本轮没有新文献搜索、原PDF复核、DOI验证或共享库修改。TES条件均值/协方差及LaBr3/PMT模板适用条件不构成GAGG共享微胞性能证明，原AI草稿和限定人工审核范围保持。四项残差恒等式、SVD局部投影与完整重复协方差为本项目计算，不冒充外部文献结论。

**条件与限制。** 只读L98已暴露8对母的40条件、相位0和8重复，模型/物理参数不变；没有接触预留确认母。权重由原电压初值导出并逐条与原112拟合核对；真E用于当前明确的事后诊断，未进入原推断。共享/孤立来自不同数字化，差项包含随机波动；须保留四项完整协方差。100δlnE投影不是实际百分能量误差；8 ns线性外推失效明显。112点代价都比真值低说明工作目标可偏向别的能量，不证明全局搜索最优。组均值背景项小不意味着全部噪声小，祖先信号也含统计涨落。

**不确定度。** 10000次整母重采样，seed742100000；相同索引保留方法、条件和各项依赖。8母区间只作条件描述，不传播训练误差或宣称同时覆盖。重复MC标准误只量化固定母/固定度量下的数字化均值波动，不是总能量分辨率，未重新标定均值协方差。

**计数与检查。** {counts} {validation}

**下一步。** {next_text} 盲漏检、误拆分/误合并、能谱效率、高率及实物结果仍未测，未测不填零。

位置：`pileup/experiments/{study.NAME}/`及`pileup/analysis/alpha_current_energy_diagnosis*.py`；本地codex_0007第142–147节、图0007-143-1/144-1/145-1/146-1，研究计划第40节和Run_Catalog。原Notebook全部计算与图像保留。
''')
    p=ROOT/'README.md';before=(REC/'README.md.before').read_text();assert p.read_text()==before
    p.write_text(before+f'\n\n2026-11-02完成[短间隔能量失败诊断](pileup/experiments/{study.NAME}/README.md)：'+key+
        ' 8次重复仍不足以把全部残差视为确定模型偏差，8 ns局部线性外推误差明显。'+next_text+
        ' 结果见本地codex_0007第142–147节与[采用记录L99](docs/literature-use.md#current-energy-diagnosis-20261102)。\n')


if __name__=='__main__':main()
