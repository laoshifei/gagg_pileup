"""Paired optimizer comparison, figures, and executed Notebook delivery."""
from __future__ import annotations
import argparse
import ast
import csv
import hashlib
import io
import json
import os
import re
import numpy as np
from scipy.stats import binomtest
from .alpha_pair_optimizer_polish_study import (
    PROJECT,NAME,EXPERIMENT,OUTPUT,RECORD,NOTEBOOK,METHODS,
    load,save,digest,rows,write_new,now,models,check_freeze)

LABELS=dict(baseline='Original',restart='Wide restart',multistart='Multistart',polished='+ Powell')
COLORS=dict(baseline='#557799',restart='#dd9944',multistart='#55aa77',polished='#aa5577')
MARKER='optimizer-polish-20261028'
LITERATURE='''据共享库既有精读笔记，Haselman等（2012，DOI10.1109/TNS.2012.2207403，PDF3–7图2–11及表I）支持分开评价电荷、定时与保留计数；原条件为LFS-3/MAPD-N孤立实测脉冲构成的合成流，65–1000MHz采样、100kcps–1Mcps，不是本项目GAGG性能保证。Žugec等（2016，DOI10.1016/j.nima.2015.12.054，PDF7–9式13–23、图9–10）提供模板拟合、亚样点相位及质量检查的依据；n_TOF多种探测器、典型100MHz–1GHz采样、形状稳定模板和顺序扣尾，不等同本轮联合高斯优化。

本轮只读核对主题、审核入口、两篇来源笔记及global-review中N29-010所在原始AI审核整理；两篇仍为AI待审，没有新增人工认可，未重读PDF、联网或修改共享库。延长迭代、多起点和Powell组合是项目数值消融设计，不是这两篇的算法复现；不把较低代价解释成真值误差必然更小或全局最优。'''

def freeze_extra():
    check_freeze();assert not (OUTPUT/'test').exists()
    folder=EXPERIMENT/'report-source';folder.mkdir(exist_ok=False)
    p=PROJECT/'pileup/analysis/alpha_pair_optimizer_polish_report.py'
    write_new(folder/p.name,p.read_bytes())
    wrapper='''from pathlib import Path
import os,sys
folder=Path(__file__).resolve().parent
root=str(folder.parents[2])
os.environ['GAGG_POLISH_PROJECT']=root
os.environ['GAGG_UNCERTAINTY_PROJECT']=root
sys.path.insert(0,str(folder/'source'))
import pileup.analysis
pileup.analysis.__path__.insert(0,str(folder/'report-source'))
if __name__=='__main__':
    from pileup.analysis.alpha_pair_optimizer_polish_report import main,check_extra
    check_extra();main()
'''
    write_new(EXPERIMENT/'report.py',wrapper.encode())
    save(EXPERIMENT/'report-freeze.json',dict(created=now(),files={
        'report-source/'+p.name:digest(folder/p.name),'report.py':digest(EXPERIMENT/'report.py')}))
    print('REPORT_FROZEN',flush=True)

def check_extra():
    for name,h in load(EXPERIMENT/'report-freeze.json')['files'].items():assert digest(EXPERIMENT/name)==h,name

def fraction(k,n):
    if not n:return dict(numerator=int(k),denominator=0,value=None,binomial95=None)
    ci=binomtest(k,n).proportion_ci(.95,method='exact')
    return dict(numerator=int(k),denominator=int(n),value=k/n,binomial95=[ci.low,ci.high])

def spread(values):
    if not values:return None
    a=np.percentile(values,[16,50,84]);return dict(median=float(a[1]),sigma68=float((a[2]-a[0])/2),p16=float(a[0]),p84=float(a[2]))

def summarize(group,band,method):
    e=[r['evaluation'] for r in group]
    return dict(band=band,method=method,frames=len(group),strict=fraction(sum(r['strict_success'] for r in group),len(group)),
        output_K2=sum(r['predicted_count']==2 for r in group),target_events=sum(v['target_count'] for v in e),
        matched_events=sum(v['matched_count'] for v in e),recovered_events=sum(v['recovered_count'] for v in e),
        missed=sum(v['missed_count'] for v in e),extra=sum(v['unmatched_predictions'] for v in e),
        split=sum(v['split_flags'] for v in e),merge=sum(v['merge_flags'] for v in e),
        Q_error_matched=spread([x for v in e for x in v['relative_charge_errors']]),
        t_error_matched_ns=spread([x for v in e for x in v['time_errors_ns']]),
        baseline_failures=sum(bool(r['baseline_failure']) for r in group),polish_failures=sum(bool(r['polish_failure']) for r in group))

def analyze():
    check_freeze();scored={}
    for mode in ('development','test'):
        v=load(RECORD/(mode+'-validation.json'));assert v['passed']
        assert digest(OUTPUT/mode/'scores.truth.json')==v['scores_sha256']
        scored[mode]=load(OUTPUT/mode/'scores.truth.json')
    test=scored['test'];bands=list(dict.fromkeys(r['band'] for r in test));summaries=[];paired=[]
    by={(r['index'],r['method']):r for r in test}
    for band in bands:
        for method in METHODS:
            group=[r for r in test if r['band']==band and r['method']==method]
            summaries.append(summarize(group,band,method))
            if method=='baseline':continue
            improved=sum(r['strict_success'] and not by[r['index'],'baseline']['strict_success'] for r in group)
            harmed=sum(not r['strict_success'] and by[r['index'],'baseline']['strict_success'] for r in group)
            selected=[r for r in group if r['predicted_count']==2]
            paired.append(dict(band=band,method=method,frames=len(group),improved=improved,harmed=harmed,
                unchanged=len(group)-improved-harmed,net=improved-harmed,
                exact_discordance_p=float(binomtest(improved,improved+harmed,.5).pvalue) if improved+harmed else 1.,
                selected_K2=len(selected),cost_improved_over005=sum(r['cost_gain']>.05 for r in selected),
                selected_cost_gain=spread([r['cost_gain'] for r in selected])))
    dev=[];original=load(OUTPUT/'development/baseline.json')
    for r in scored['development']:
        if not r['known_witness']:continue
        baseline=original[r['index']]['selected'];witness=original[r['index']]['prior_exterior']
        dev.append(dict(original_index=r['original_index'],method=r['method'],strict=r['strict_success'],
            cost_gain=r['cost_gain'],accuracy_distance=r['accuracy_distance'],Q_errors=r['ordered_Q_errors'],
            t_errors_ns=r['ordered_t_errors_ns'],prior_witness_gain=baseline['cost']-witness['cost']))
    predictions=rows(OUTPUT/'test/predictions.jsonl');searched=[r for r in predictions if r['diagnostic'] and r['diagnostic']['searched']]
    report=dict(created=now(),primary='polished',summary=summaries,paired=paired,development_six=dev,
        test_frames=288,development_frames=15,selected_double_frames=len(searched),
        primary_cost_improvements=sum(r['fits']['baseline']['cost']-r['fits']['polished']['cost']>.05 for r in searched),
        secondary_trial_failures=sum(len(r['diagnostic'].get('failures',[])) for r in searched),
        selected_origin_counts={name:sum(r['fits']['polished'].get('origin')==name for r in searched)
                               for name in sorted({r['fits']['polished'].get('origin','original') for r in searched})},
        actual_inference_seconds=load(RECORD/'test-inference.json')['seconds'],
        count_unchanged_by_design=True,non_K2_exact_copies=True,not_a_global_optimum_proof=True,
        no_new_ROOTs=True,no_training=True,default_changed=False,
        secondary_methods_exploratory=True,exact_tests_not_adjusted_across_conditions=True)
    save(OUTPUT/'report.json',report)
    stream=io.StringIO();writer=csv.writer(stream);writer.writerow(['band','method','frames','K2','strict','missed','extra','split','merge'])
    for r in summaries:writer.writerow([r['band'],r['method'],r['frames'],r['output_K2'],r['strict']['numerator'],r['missed'],r['extra'],r['split'],r['merge']])
    write_new(OUTPUT/'summary.csv',stream.getvalue().encode())
    save(RECORD/'report-validation.json',dict(passed=True,report_sha256=digest(OUTPUT/'report.json'),
        all_288_test_frames_and_15_development_frames_retained=True,methods=4,created=now()))
    print(conclusions(report),flush=True)

def table(report):
    lines=['| 间隔/ns | 原流程 | 宽域重启 | 多起点 | 再加Powell（主方案） | 改善/退化 |',
           '| --- | --- | --- | --- | --- | --- |']
    for band in ('pair40','pair200'):
        row=[]
        for m in METHODS:
            s=next(r for r in report['summary'] if r['band']==band and r['method']==m)
            row.append(f"{s['strict']['numerator']}/{s['frames']}")
        p=next(r for r in report['paired'] if r['band']==band and r['method']=='polished')
        lines.append('| '+band[4:]+' | '+' | '.join(row)+f" | {p['improved']}/{p['harmed']} |")
    return '\n'.join(lines)

def conclusions(report):
    parts=['独立288条确认采用同一波形的配对比较，原事件数门和K判决保持不变；主方案为预先固定的多起点后Powell，另外两种仅作过程消融。']
    for band in ('pair40','pair200'):
        b=next(r for r in report['summary'] if r['band']==band and r['method']=='baseline')
        p=next(r for r in report['summary'] if r['band']==band and r['method']=='polished')
        d=next(r for r in report['paired'] if r['band']==band and r['method']=='polished')
        parts.append(f"{band[4:]}ns：整对严格恢复{b['strict']['numerator']}/{b['frames']}→{p['strict']['numerator']}/{p['frames']}；改善{d['improved']}、退化{d['harmed']}，配对精确检验p={d['exact_discordance_p']:.4g}。时间漏检{b['missed']}→{p['missed']}，额外输出{b['extra']}→{p['extra']}。")
    controls=[r for r in report['summary'] if not r['band'].startswith('pair') and r['method']=='baseline']
    parts.append(f"{report['selected_double_frames']}条被盲选为K2的确认波形中，{report['primary_cost_improvements']}条代价降低超过0.05。控制组原误报双事件{sum(r['output_K2'] for r in controls)}/160，K误报数保持不变是设计约束，不能作为新筛选器有效的证据。")
    parts.append('较低代价与逐事件准确度须分别解释；有限搜索未找到更优解不证明物理不可辨识。本轮不替换默认算法，也不外推为原生SiPM、α/γ、free-τ或高计数率性能。')
    return '\n\n'.join(parts)

def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.2})
    report=load(OUTPUT/'report.json');test=load(OUTPUT/'test/scores.truth.json')
    folder=OUTPUT/'figures';folder.mkdir(exist_ok=False);made=[]
    def finish(fig,number,title):
        fig.suptitle('Figure '+number+': '+title,fontsize=12)
        fig.tight_layout(rect=(0,0,1,.95));stream=io.BytesIO()
        fig.savefig(stream,format='png',dpi=150);plt.close(fig)
        write_new(folder/(number+'.png'),stream.getvalue());made.append(number)
    dev=report['development_six'];indices=list(dict.fromkeys(r['original_index'] for r in dev));x=np.arange(len(indices))
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    for j,m in enumerate(METHODS[1:]):
        rr=[next(r for r in dev if r['original_index']==i and r['method']==m) for i in indices]
        axes[0].bar(x+(j-1)*.24,[r['cost_gain'] for r in rr],.23,label=LABELS[m],color=COLORS[m])
    refs=[next(r for r in dev if r['original_index']==i) for i in indices]
    axes[0].scatter(x,[r['prior_witness_gain'] for r in refs],marker='_',s=180,color='black',label='Prior exterior witness')
    for m in ('baseline','polished'):
        rr=[next(r for r in dev if r['original_index']==i and r['method']==m) for i in indices]
        axes[1].plot(x,[r['accuracy_distance'] for r in rr],'o-',label=LABELS[m],color=COLORS[m])
    axes[1].axhline(1,color='black',ls='--');axes[1].set_yscale('symlog',linthresh=1.)
    for ax in axes:ax.set_xticks(x,[str(i) for i in indices]);ax.set_xlabel('Exposed old frame index');ax.legend(fontsize=8)
    axes[0].set(title='(a) Objective improvement',ylabel='Original cost - candidate cost')
    axes[1].set(title='(b) Actual parameter accuracy',ylabel='Largest normalized Q/t error (pass <= 1)')
    finish(fig,'0004-133-1','Six exposed optimization failures: cost and truth error are separate')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    for ax,band,panel in zip(axes,('pair40','pair200'),('a','b')):
        rr=[next(r for r in report['summary'] if r['band']==band and r['method']==m) for m in METHODS]
        for j,r in enumerate(rr):
            f=r['strict'];ax.bar(j,f['value'],color=COLORS[r['method']])
            ax.errorbar(j,f['value'],yerr=[[f['value']-f['binomial95'][0]],[f['binomial95'][1]-f['value']]],color='black',capsize=3)
            ax.text(j,min(1.07,f['binomial95'][1]+.025),f"{f['numerator']}/{f['denominator']}",ha='center')
        ax.set_xticks(range(4),[LABELS[m] for m in METHODS]);ax.set_ylim(0,1.14)
        ax.set(title=f'({panel}) '+band[4:]+' ns',ylabel='Strict pair recovery / all true pairs')
    finish(fig,'0004-135-1','Fresh paired confirmation; exact binomial sampling intervals')
    by={(r['index'],r['method']):r for r in test}
    fig,axes=plt.subplots(1,2,figsize=(11,5))
    for ax,band,panel in zip(axes,('pair40','pair200'),('a','b')):
        rr=[r for r in test if r['band']==band and r['method']=='polished' and r['accuracy_distance'] is not None]
        before=[by[r['index'],'baseline']['accuracy_distance'] for r in rr];after=[r['accuracy_distance'] for r in rr]
        ax.scatter(before,after,c=[COLORS['polished'] if r['cost_gain']>.05 else '.65' for r in rr],alpha=.8)
        maximum=max([1.5]+before+after)*1.12
        ax.plot([0,maximum],[0,maximum],'k:',lw=1);ax.axvline(1,color='gray',ls='--');ax.axhline(1,color='gray',ls='--')
        ax.set_xscale('symlog',linthresh=1.);ax.set_yscale('symlog',linthresh=1.)
        ax.set(xlim=(-.05,maximum),ylim=(-.05,maximum),title=f'({panel}) '+band[4:]+' ns',
               xlabel='Original normalized Q/t error',ylabel='Polished normalized Q/t error')
    finish(fig,'0004-135-2','Selected pairs: lower cost can improve or worsen parameter recovery')
    response,noise=models();pred={r['index']:r for r in rows(OUTPUT/'test/predictions.jsonl')}
    with np.load(OUTPUT/'test/observations.npz') as f:waves=f['voltage_mV']
    cases=load(OUTPUT/'test/labels.cases.json');grid=response.sample_grid
    for j,index in enumerate(load(EXPERIMENT/'protocol.json')['example_indices'],1):
        case=cases[index];p=pred[index];truth=response.components(case['times_ns'],case['charge_pe'])
        fits={m:response.components(p['fits'][m]['times_ns'],p['fits'][m]['charge_pe']) for m in ('baseline','polished')}
        offset=response.baseline+noise.background_mV
        fig,axes=plt.subplots(3,1,figsize=(11,8),sharex=True)
        axes[0].plot(grid,waves[index],color='.7',lw=.6,label='Observed')
        axes[0].plot(grid,truth.sum(0)+offset,'k--',lw=1,label='Truth sum (display only)')
        for m in ('baseline','polished'):
            axes[0].plot(grid,fits[m].sum(0)+offset,color=COLORS[m],lw=1.1,label=f'{LABELS[m]} K={len(fits[m])}')
        for k in (0,1):
            axes[k+1].plot(grid,truth[k],'k--',lw=1,label=f'True component {k+1}')
            for m in ('baseline','polished'):
                if len(fits[m])>k:
                    f=p['fits'][m];axes[k+1].plot(grid,fits[m][k],color=COLORS[m],lw=1.1,
                        label=f"{LABELS[m]} Q={f['charge_pe'][k]:.0f} pe, t={f['times_ns'][k]:.2f} ns")
            if k==1 and len(fits['baseline'])<2:
                axes[k+1].text(.55,.55,'No second event selected',transform=axes[k+1].transAxes)
        for ax,title in zip(axes,('(a) Complete waveform','(b) First output and strong truth','(c) Second output and weak truth')):
            ax.set_title(title,fontsize=10,loc='left');ax.set_ylabel('mV');ax.set_xlim(0,float(grid[-1]));ax.legend(fontsize=8)
        axes[-1].set_xlabel('Time in complete record / ns')
        finish(fig,f'0004-136-{j}',f'Prespecified fresh frame {index}, gap {case["times_ns"][1]-case["times_ns"][0]:g} ns')
    save(RECORD/'figures-validation.json',dict(passed=True,figures=made,created=now(),
        files={p.name:digest(p) for p in folder.glob('*.png')}))
    print('FIGURES',len(made),flush=True)

def replace_verified(path,payload,expected):
    temporary=path.with_name(path.name+'.optimizer-writing');write_new(temporary,payload)
    assert digest(path)==expected,'Concurrent edit; prepared file retained: '+str(path)
    os.replace(temporary,path);assert path.read_bytes()==payload

def notebook():
    import nbformat
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    check_freeze();report=load(OUTPUT/'report.json')
    assert load(RECORD/'report-validation.json')['passed'] and load(RECORD/'figures-validation.json')['passed']
    old=(RECORD/'notebook.before.ipynb').read_bytes();assert NOTEBOOK.read_bytes()==old
    book=nbformat.reads(old.decode(),as_version=4);assert len(book.cells)==627;cells=[]
    def md(text):cells.append(nbformat.v4.new_markdown_cell(text))
    def code(text):cells.append(nbformat.v4.new_code_cell(text))
    def fig(number,caption):
        md(f'<a id="codex-figure-{number}"></a>\n**图{number}：{caption}**')
        code(f'display(Image(filename=str(polish_output / "figures/{number}.png")))')
    md('''<a id="optimizer-polish-20261028"></a>
## 133. 原优化器“已收敛”，是否仍能找到更好的拟合？

上一轮79条盲选双事件中有6条存在代价降低超过0.05的容差外替代点。它们的原优化器都报告收敛；该标志只是局部停止条件，不能证明全局最优。本轮检查数值搜索不足是否实际影响逐事件恢复。

先复查15条已公开旧波形，其中6条为上述问题案例，另外包含固定选择的双脉冲、空窗、孤立脉冲和误拆案例。它们不是随机性能样本，旧6条是按代价问题挑出的诊断案例。随后在源码、参数和主方案冻结后另生成288条独立噪声波形，不用新测试结果调参。

旧真值与此前的替代解只在评分和展示中读取，不作为重拟合起点。推断worker实际尝试并阻断truth/cases/ROOT文件访问；只输入观测、冻结模板/噪声/门及观测产生的原候选。
''')
    code(f'''from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Image
polish_root=Path.cwd().resolve()
while not (polish_root / "pileup/experiments/{NAME}").is_dir():
    if polish_root.parent==polish_root: raise RuntimeError("请从项目根目录或jupyter目录运行")
    polish_root=polish_root.parent
polish_output=polish_root / "jupyter/artifacts/{NAME}"
polish_report=json.loads((polish_output / "report.json").read_text())
display(pd.DataFrame(polish_report["development_six"]))
''')
    fig('0004-133-1','六条已暴露问题案例：(a)代价降低，黑色横线为上一轮找到的替代解；(b)相对真值的最大归一化Q/t误差，≤1才通过参数容差。两者不能互相替代。')
    md('''## 134. 保持盲判事件数的优化消融

完整流程仍先用原waveform_candidates和原生控制标定的0.025门判断未知K。仅对原流程选出的K=2重拟合；K=0/1输出逐项保留，因此本轮无法挽回已判成单事件的漏检。重拟合后不另加Q≥500pe筛选，避免悄悄改变K；Q仍限制50–50000pe。不能把“不增加误报双事件数”当成新选门的效果。

比较四个阶段：原流程；从原点在更宽域重启L-BFGS-B；六个分散观测起点；再以最佳结果做Powell方向搜索（预先指定的主方案）。第一阶段同时改变搜索域/坐标并增加迭代，不把其作用全部归因于延长迭代。全过程保留此前的较优候选，实际代价必须不升。

参数域沿用上一轮的观测定义：原K1时刻前96至后512ns，仅为包含原K2扩展，再截至完整记录。两事件时刻有序，电荷取对数坐标。起点来自原拟合、各时刻±1/2/4/8/16/32ns，以及观测K1锚点与4/8/16/32/64/128/256/512ns试探间隔、强弱两种顺序和七种电荷份额。未输入真40/200ns或真电荷。L-BFGS-B最多160迭代，Powell最多100迭代/2400次调用；有限预算不能证明全局最优。

目标函数仍是包含logdet的候选相关高斯代价，物理响应、协方差、固定95/351ns双衰减均不改变；没有加入真值先验或网络训练。本轮是优化消融，不宣称复现文献中的新物理模型。
''')
    md('''## 135. 新随机样本的逐事件配对结果

确认集含40/200ns先强后弱双脉冲各64条，Q=[18000,800]pe、t1=200.3ns；另空窗及800/3000/8000/18800pe孤立脉冲各32条，共288条。噪声种子为SeedSequence[202610281501,3,index,997]，与旧432条波形逐行哈希不重复。都是条件高斯模型内独立抽样，尚不是新原生SiPM或光学母事件验证。

完整1024×2ns float64电压、足够动态范围；FIR首三个残差为零。时间一对一匹配≤10ns，电荷相对误差以真实Q为分母且≤20%，整对严格成功还要求没有额外输出。Q为等效pe，没有标定成沉积能量。失败、漏检、错配和误拆均留在分母。

每种方法在同一条波形上比较；“改善/退化”分别是原失败变成功、原成功变失败。主方案在测试前固定，其他阶段为解释优化过程的探索性消融。逐固定参数组的精确二项95%区间不包含物理模型误差；配对精确检验基于不一致结果，未作多条件检验校正，不据单个小p值宣布普适改进。
''')
    md(table(report))
    fig('0004-135-1','(a)40ns、(b)200ns。四个优化阶段的整对严格恢复，均以全部64条真实双脉冲作分母；误差条仅表示当前固定模型的重复抽样不确定性。')
    fig('0004-135-2','(a)40ns、(b)200ns中被原流程选为双事件的样本。横纵轴为优化前后最大归一化Q/t误差；紫色表示代价降低>0.05，灰色为其余。左下角通过，右下角改善，左上角退化；漏检仍在上一张图分母内。')
    code('''display(pd.DataFrame(polish_report["summary"]))
display(pd.DataFrame(polish_report["paired"]))
''')
    md('''匹配后的Q误差中位数和sigma68=(p84−p16)/2、时间误差，须结合漏检和额外输出阅读；sigma68不是能谱FWHM。事件数保持不变，但拟合时刻改变后，时间错配造成的漏检/额外输出仍可能改变。误合并/误拆分沿用既有邻近匹配定义，不冒充所有失败的物理原因。''')
    md('''## 136. 完整波形及分量展示

预先指定新测试序号0和64，分别来自40ns和200ns组，不根据成功与否选图。展示完整0–2046ns，强、弱分量分轴便于观察。K=1时第一输出代表原单脉冲拟合，不强行把它解释成真实强分量；第二输出按实际缺失保留。真值只用于显示。
''')
    for j,index in enumerate((0,64),1):
        fig(f'0004-136-{j}',f'新测试序号{index}：(a)观测及拟合总波形；(b)第一输出与真实强分量；(c)第二输出与真实弱分量。原流程及优化结果均保持实际事件数。')
    md('## 137. 本轮结论与适用范围\n\n'+conclusions(report)+'\n\n'+LITERATURE+f'''

可复现源码、模型、原事件数门、设计与运行器保存在pileup/experiments/{NAME}/；数据、预测、完整评分、CSV及PNG保存在jupyter/artifacts/{NAME}/；测试、环境、SHA256及旧笔记备份保存在jupyter/run_records/{NAME}/。旧627个单元和输出完整保留。没有生成ROOT、修改共享库或替换默认拟合。

下一阶段应区分选门漏检与被选中后的参数误差，再扩展幅度比、间隔、采样相位和窗口。更低的目标函数、局部停止标志或本轮两个间隔，都不能推出普适物理分辨极限。
''')
    shell=InteractiveShell.instance();shell.user_ns.clear();shell.user_ns['__name__']='__main__';executed=0
    for i,cell in enumerate(cells):
        cell.id=f'optimizer-polish-{i:03d}';cell.metadata['gagg_section_marker']=MARKER
        if cell.cell_type=='markdown':
            assert not re.search(r'\]\([^)]*\.(?:md|ipynb|pdf)(?:[#?][^)]*)?\)',cell.source,re.I)
        else:
            executed+=1
            with capture_output(display=True) as capture:result=shell.run_cell(cell.source,store_history=False)
            assert result.error_before_exec is None and result.error_in_exec is None
            output=[]
            if capture.stdout:output.append(nbformat.v4.new_output('stream',name='stdout',text=capture.stdout))
            if capture.stderr:output.append(nbformat.v4.new_output('stream',name='stderr',text=capture.stderr))
            output.extend(nbformat.v4.new_output('display_data',data=o.data,metadata=o.metadata) for o in capture.outputs)
            cell.outputs=output;cell.execution_count=executed
    book.cells.extend(cells);nbformat.validate(book);payload=nbformat.writes(book).encode()
    assert json.loads(payload)['cells'][:627]==json.loads(old)['cells']
    replace_verified(NOTEBOOK,payload,hashlib.sha256(old).hexdigest())
    save(RECORD/'notebook-validation.json',dict(passed=True,previous_cells=627,total_cells=len(book.cells),
        added_cells=len(cells),executed_cells=executed,old_cells_outputs_preserved=True,notebook_sha256=digest(NOTEBOOK)))
    print('NOTEBOOK_UPDATED',len(book.cells),executed,flush=True)

def documents():
    report=load(OUTPUT/'report.json');text=conclusions(report)
    for p,b in ((PROJECT/'README.md',RECORD/'README.before.md'),(PROJECT/'docs/literature-use.md',RECORD/'literature.before.md')):
        assert digest(p)==digest(b)
    readme=(RECORD/'README.before.md').read_text()+f'\n\n### 2026-10-28：保持事件数判决的优化器独立消融\n\n{text}\n\n{table(report)}\n\n先检查15条已暴露开发波形，再冻结方案生成288条新确认波形；比较原流程、宽域重启、多起点和Powell。完整记录见[实验说明](pileup/experiments/{NAME}/README.md)；本地0004笔记第133–137节新增5幅固定编号图，原627个单元与输出保留。全部结果仍限固定双τ的纯α条件高斯模型。\n'
    literature=(RECORD/'literature.before.md').read_text();assert not re.search(r'^## L70\b',literature,re.M)
    literature+='\n\n## L70 — 2026-10-28：原点重启、多起点与Powell的独立优化消融\n\n'+LITERATURE+'''

| 采用依据 | 本轮核查层级与边界 |
| --- | --- |
| Haselman等，FPGA-Based Pulse Pile-Up Correction With Energy and Timing Recovery，TNS59(5),1823–1830(2012)；[DOI](https://doi.org/10.1109/TNS.2012.2207403)；[既有笔记](../../papers/sources/2012_Haselman_TNS_FPGA-Based-Pulse-Pile-Up-Correction-With-Energy-and-Timing-Recovery.md)；[原PDF](../../papers/raw/2012_Haselman_TNS_FPGA-Based-Pulse-Pile-Up-Correction-With-Energy-and-Timing-Recovery.pdf) | 仅读既有完整来源笔记；PDF3–7图2–11/表I页码沿用笔记。global-review的N29-010所属段明确新增判断仍AI待审；约30%额外保留和能量/定时变化有原合成流条件，不作当前优化器效果依据。 |
| Žugec等，Pulse processing routines for neutron time-of-flight data，NIMA812,134–144(2016)；[DOI](https://doi.org/10.1016/j.nima.2015.12.054)；[既有笔记](../../papers/sources/2016_Zugec_NIMA_Pulse-processing-routines-for-neutron-time-of-flight-data.md)；[原PDF](../../papers/raw/2016_Zugec_NIMA_Pulse-processing-routines-for-neutron-time-of-flight-data.pdf) | 仅读既有完整来源笔记，定位PDF7–9式13–23/图9–10；形状质量和重构曲线不代替逐事件效率。AI待审，旧555意见不扩展为本篇人工认可。 |

**项目设计**：不改变物理响应、候选相关协方差和完整Gaussian quadratic+logdet。原未知K选门固定，仅原K2参数重拟合；K0/K1逐项保留，不在优化后重新按Q≥500筛选。宽域重启同时改变参数化/边界和迭代预算，不能把效果全归因于迭代次数。观测起点、六起点L-BFGS-B及Powell配置在新测试生成前冻结；此前替代解与真值不输入优化器。

**实际结果**：

'''+text+f'''

{table(report)}

**评价与限制**：全体64条/间隔的严格整对恢复、电荷误差中位数/sigma68、定时误差、漏检、额外输出及误拆/误合并同时保留。比较同一波形的改善与退化，不拿不同随机批次的恢复率直接评算法。主方案预定为polished；其他阶段探索性，逐组精确区间/配对检验不提供模型误差或普适保证。事件数不变不能证明新误报控制有效；未选中双事件的漏检无法被本轮后处理挽回。

**复现**：纯α固定95/351ns、完整1024×2ns电压、充足动态范围，Q为等效pe而非沉积能量；FIR首三残差零，沿用条件Gaussian模型。15条旧开发样本明确暴露；新288条使用SeedSequence[202610281501,3,index,997]且与旧432条哈希不重复。源码/模型/协议在pileup/experiments/{NAME}/，本地结果jupyter/artifacts/{NAME}/、记录jupyter/run_records/{NAME}/，Notebook为jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb第133–137节。未修改共享库、原实验或默认流程，未commit/push。

**待验证**：本次有限搜索不是全局优化或物理不可辨识证据。后续需分别检查选门漏检和参数估计，再扩展幅比/间隔/采样相位/窗口；原生独立母事件、α/γ、free-τ和连续高率仍需独立验证。若提出全局性或统计覆盖保证，需另行在共享库整理对应一手依据。
'''
    for p,content,b in ((PROJECT/'README.md',readme,RECORD/'README.before.md'),
                        (PROJECT/'docs/literature-use.md',literature,RECORD/'literature.before.md')):
        replace_verified(p,content.encode(),digest(b))
    experiment=f'''# Count-preserving optimizer ablation

{text}

{table(report)}

The original waveform-only candidate generator and native count gate remain unchanged. K=0/1 outputs are exact copies; only selected K=2 parameters are polished. The primary configuration was frozen before fresh test generation: observed starts, six L-BFGS-B starts and bounded Powell, retaining earlier lower-cost fits. Wide restart changes the parameterization/domain as well as iteration budget. No post-polish Q>=500 cut, true count, truth parameters or prior exterior witness enters the optimizer. Improved objective is not guaranteed improved event recovery, and finite searches do not prove a global optimum.

Fifteen previously exposed development frames include six cases selected for suboptimal original objective. They are not an unbiased performance sample. Independent confirmation uses 288 fresh model draws and same-waveform paired comparisons. See protocol.json and test-design.truth.json for fixed design and random seeds. Full per-event and denominator-aware scoring is saved locally. No default algorithm promotion, training or new ROOT files.

Pure alpha, fixed95/351ns, full1024x2ns float64 voltage, ample dynamic range; conditional Gaussian model with first three FIR residuals fixed at zero. Q is equivalent photoelectrons, not deposited energy. No native SiPM, mixed alpha/gamma, free-tau or continuous-rate claim. Models and arrays are ignored local assets, required separately for replay. Source snapshots contain actual dependencies; no Git checkout is required.

From project root using the project environment and one BLAS thread:

```bash
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py tests
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py development-inputs
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py infer-development --workers 4
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py score-development
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py generate
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py infer-test --workers 4
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/run.py score-test
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/report.py analyze
bash tools/run_gagg_python.sh pileup/experiments/{NAME}/report.py figures
```

Existing outputs are protected. Replay in a separate project copy with original local development assets and models. Primary and reporting source freezes are separate manifests. Notebook delivery also requires the recorded before-version of the local Notebook, README and literature-use document; old 627 cells and outputs are preserved. Full commands/environment/hashes and validations are in jupyter/run_records/{NAME}/. No commit or push was performed.
'''
    write_new(EXPERIMENT/'README.md',experiment.encode())
    save(RECORD/'documentation-validation.json',dict(passed=True,created=now(),files={str(p.relative_to(PROJECT)):digest(p)
        for p in (PROJECT/'README.md',PROJECT/'docs/literature-use.md',EXPERIMENT/'README.md')}))
    print('DOCUMENTS_UPDATED',flush=True)

def validate():
    import nbformat
    check_freeze();check_extra()
    names=('tests','development-inputs','development-inference','development-validation','test-generation',
           'test-inference','test-validation','report-validation','figures-validation','notebook-validation','documentation-validation')
    checks={n:load(RECORD/(n+'.json')) for n in names};assert all(v['passed'] for v in checks.values())
    assert checks['tests']['tests']==11 and checks['development-inference']['frames']==15 and checks['test-inference']['frames']==288
    for mode,n in (('development',15),('test',288)):
        assert digest(OUTPUT/mode/'predictions.jsonl')==checks[mode+'-inference']['sha256']
        rr=rows(OUTPUT/mode/'predictions.jsonl');assert len(rr)==n and {r['index'] for r in rr}==set(range(n))
        assert all(r['truth_inputs']==[] for r in rr)
        assert digest(OUTPUT/mode/'scores.truth.json')==checks[mode+'-validation']['scores_sha256']
        generation=checks['development-inputs' if mode=='development' else 'test-generation']
        for name,h in generation['files'].items():assert digest(OUTPUT/mode/name)==h
    assert load(EXPERIMENT/'freeze.json')['created']<checks['test-generation']['created']
    assert load(EXPERIMENT/'report-freeze.json')['created']<checks['test-generation']['created']
    assert digest(OUTPUT/'report.json')==checks['report-validation']['report_sha256']
    book=nbformat.read(NOTEBOOK,as_version=4);before=nbformat.read(RECORD/'notebook.before.ipynb',as_version=4)
    nbformat.validate(book);assert len(before.cells)==627 and book.cells[:627]==before.cells
    assert len({c.id for c in book.cells})==len(book.cells)
    codes=[c for c in book.cells[627:] if c.cell_type=='code']
    assert all(c.execution_count is not None and not any(o.output_type=='error' for o in c.outputs) for c in codes)
    pngs=sum('image/png' in o.get('data',{}) for c in codes for o in c.outputs);assert pngs==5
    assert digest(NOTEBOOK)==checks['notebook-validation']['notebook_sha256']
    for name,h in checks['figures-validation']['files'].items():assert digest(OUTPUT/'figures'/name)==h
    for name,h in checks['documentation-validation']['files'].items():assert digest(PROJECT/name)==h
    for name,h in load(RECORD/'environment.json')['original_inputs'].items():assert digest(PROJECT/name)==h
    for p in (EXPERIMENT/'source').rglob('*.py'):ast.parse(p.read_text())
    ast.parse((EXPERIMENT/'report-source/alpha_pair_optimizer_polish_report.py').read_text())
    result=dict(passed=True,created=now(),tests=11,development_frames=15,confirmation_frames=288,
        figures=pngs,previous_notebook_cells=627,total_notebook_cells=len(book.cells),executed_new_cells=len(codes),
        original_notebook_cells_outputs_preserved=True,original_count_decisions_preserved=True,
        predictions_before_truth_scoring=True,source_and_report_frozen_before_fresh_generation=True,
        all_frames_retained=True,notebook_sha256=digest(NOTEBOOK),report_sha256=digest(OUTPUT/'report.json'),
        no_default_change=True,no_new_ROOTs=True,no_training=True,no_commit_or_push=True)
    save(RECORD/'final-validation.json',result);print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['freeze-extra','analyze','figures','notebook','documents','validate'])
    args=parser.parse_args();globals()[args.command.replace('-','_')]()

if __name__=='__main__':main()
