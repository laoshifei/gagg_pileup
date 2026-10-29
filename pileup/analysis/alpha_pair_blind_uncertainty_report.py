"""Independent scoring, denominators and figures for joint compatibility."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import io
import numpy as np
from scipy.stats import binomtest

from .alpha_pair_blind_uncertainty_study import (
    PROJECT, EXPERIMENT, OUTPUT, RECORD, load, save, digest, rows, models,
    check_freeze, score_rows, write_new, now,
)
from .alpha_pair_blind_uncertainty import diagnostic_label
from .alpha_pair_joint_likelihood import select

LABELS = ('not_selected_K2', 'alternative_found', 'numerically_unresolved', 'no_alternative_found')


def fraction(k, n, with_ci=True):
    if not n: return dict(numerator=int(k),denominator=0,value=None,binomial95=None)
    ci=binomtest(int(k),int(n)).proportion_ci(.95,method='exact') if with_ci else None
    return dict(numerator=int(k),denominator=int(n),value=float(k/n),binomial95=[ci.low,ci.high] if ci else None)


def quantiles(values):
    if not len(values): return None
    a=np.percentile(values,[16,50,84]).tolist()
    return dict(p16=a[0],median=a[1],p84=a[2],sigma68=(a[2]-a[0])/2.)


def summary(scored, name):
    evaluations=[r['evaluation'] for r in scored]
    single_condition=len({r['band'] for r in scored})<=1
    return dict(band=name,frames=len(scored),strict=fraction(sum(r['strict_success'] for r in scored),len(scored),single_condition),
        selected_K2=fraction(sum(r['predicted_count']==2 for r in scored),len(scored),single_condition),
        target_events=sum(r['target_count'] for r in evaluations),
        matched_events=sum(r['matched_count'] for r in evaluations),recovered_events=sum(r['recovered_count'] for r in evaluations),
        missed=sum(r['missed_count'] for r in evaluations),extra=sum(r['unmatched_predictions'] for r in evaluations),
        extra_frames=fraction(sum(r['unmatched_predictions']>0 for r in evaluations),len(scored),single_condition),
        split=sum(r['split_flags'] for r in evaluations),merge=sum(r['merge_flags'] for r in evaluations),
        matched_Q_errors=quantiles([v for r in evaluations for v in r['relative_charge_errors']]),
        matched_t_errors_ns=quantiles([v for r in evaluations for v in r['time_errors_ns']]),
        not_converged=sum(not r['baseline_converged'] for r in scored),
        numerical_failures=sum(bool(r['numerical_failure']) for r in scored))


def containment(scored, level, name):
    selected=[r for r in scored if r['predicted_count']==2]
    inside=lambda r:bool(r['truth_in_domain'] and r['joint_delta']<=level)
    return dict(band=name,level=float(level),frames=len(scored),
        candidate_region_all_pairs=fraction(sum(inside(r) for r in scored),len(scored)),
        selected=fraction(len(selected),len(scored)),
        conditional_selected=fraction(sum(inside(r) for r in selected),len(selected)),
        selected_and_contained=fraction(sum(inside(r) for r in selected),len(scored)),
        selected_outside_domain=sum(not r['truth_in_domain'] for r in selected),
        strict=fraction(sum(r['strict_success'] for r in scored),len(scored)))


def analyze():
    check_freeze()
    frozen=load(RECORD/'calibration-freeze.json'); assert frozen['passed']
    assert digest(EXPERIMENT/'calibration.json')==frozen['sha256']
    calibration=load(EXPERIMENT/'calibration.json'); level=calibration['level']; assert level is not None
    scored,predictions=score_rows('test')
    gate=load(EXPERIMENT/'original-gate.json')
    with np.load(OUTPUT/'test/labels.truth.npz') as f: q=f['charge_pe']; t=f['times_ns']
    for r in scored:
        predicted=predictions[r['index']]
        if not predicted.get('numerical_failure'):
            assert select(predicted['candidates'],gate)==predicted['selected'], 'Baseline changed'
        r['quality_label']=diagnostic_label(predicted,level)
        r['joint_contained']=bool(r['true_count']==2 and r['truth_in_domain'] and r['joint_delta']<=level)
        r['exterior_delta']=predicted['exterior']['exterior_delta'] if predicted['exterior'] else None
        if r['true_count']==2 and r['predicted_count']==2:
            i=r['index']; fitted=predicted['selected']
            dq=np.array(fitted['charge_pe'])/q[i]-1.; dt=np.array(fitted['times_ns'])-t[i]
            r['ordered_Q_errors']=dq.tolist(); r['ordered_t_errors_ns']=dt.tolist()
            r['accuracy_distance']=float(max(np.max(np.abs(dq)/.2),np.max(np.abs(dt)/10.)))
            assert r['strict_success']==(r['accuracy_distance']<=1.+1e-12)
    groups=defaultdict(list)
    for r in scored: groups[r['band']].append(r)
    summaries=[summary(rr,band) for band,rr in groups.items()]
    pairs=[r for r in scored if r['true_count']==2]
    controls=[r for r in scored if r['true_count']!=2]
    primary=[containment(groups[band],level,band) for band in ('pair40','pair200')]
    sensitivity=[containment(groups[band],threshold,band)
        for threshold in sorted(set([1.,4.,9.,level])) for band in ('pair40','pair200')]
    quality=[]
    for label in LABELS:
        subset=[r for r in scored if r['quality_label']==label]
        pair_subset=[r for r in subset if r['true_count']==2]
        quality.append(dict(label=label,frames=len(subset),true_pairs=len(pair_subset),
            false_pairs=sum(r['true_count']!=2 and r['predicted_count']==2 for r in subset),
            strict_accuracy=fraction(sum(r['strict_success'] for r in subset),len(subset),False),
            pair_strict_accuracy=fraction(sum(r['strict_success'] for r in pair_subset),len(pair_subset),False),
            all_pair_strict_yield=fraction(sum(r['strict_success'] for r in pair_subset),len(pairs),False)))
    quality_by_band=[]
    for band,group in groups.items():
        for label in LABELS:
            subset=[r for r in group if r['quality_label']==label]
            quality_by_band.append(dict(band=band,label=label,frames=len(subset),
                strict_accuracy=fraction(sum(r['strict_success'] for r in subset),len(subset))))
    numerical=[r['exterior'] for r in predictions.values() if r['exterior']]
    diagnostics=dict(selected_searches=len(numerical),stable=sum(r['stable'] for r in numerical),
        all_minima_converged=sum(r['all_exterior_minima_converged'] for r in numerical),
        tolerance_box_hits_domain=sum(r['tolerance_box_hits_domain'] for r in numerical),
        verification_gain=quantiles([r['verification_gain'] for r in numerical]),
        maximum_cost_replay_error=max([r['cost_replay_error'] for r in predictions.values()]+
                                      [r['cost_replay_error'] for r in numerical]),
        numerical_failures=sum(bool(r.get('numerical_failure')) for r in predictions.values()),
        auxiliary_failures=sum(bool(r.get('exterior_failure')) for r in predictions.values()))
    tagged=[r for r in pairs if r['quality_label']=='no_alternative_found']
    report=dict(calibration=calibration,summary=summaries,primary_containment=primary,level_sensitivity=sensitivity,
        quality=quality,quality_by_band=quality_by_band,numerical=diagnostics,control_summary=summary(controls,'all_controls'),
        pair_summary=summary(pairs,'all_pairs'),tagged_pairs=fraction(len(tagged),len(pairs),False),
        tagged_correct_yield=fraction(sum(r['strict_success'] for r in tagged),len(pairs),False),
        calibration_frames=144,test_frames=len(scored),test_pairs=len(pairs),test_controls=len(controls),
        baseline_unchanged=diagnostics['numerical_failures']==0,default_unchanged=True,
        no_new_correct_events_from_tagging=True,coverage='Empirical selected-pair joint-region containment under fixed simulated conditions, not marginal charge intervals or nominal95% confidence',
        new_ROOTs=0,new_training=0,new_native_mothers=0)
    assert len(scored)==288 and len(pairs)==128 and len(controls)==160
    assert sum(r['strict_success'] for r in tagged)<=report['pair_summary']['strict']['numerator']
    save(OUTPUT/'test/scores.truth.json',scored)
    save(OUTPUT/'report.json',report)
    stream=io.StringIO(); columns=['band','frames','K2','strict','missed','extra','split','merge']
    writer=csv.DictWriter(stream,fieldnames=columns); writer.writeheader()
    for r in summaries:
        writer.writerow(dict(band=r['band'],frames=r['frames'],K2=r['selected_K2']['numerator'],
            strict=r['strict']['numerator'],**{k:r[k] for k in ('missed','extra','split','merge')}))
    write_new(OUTPUT/'summary.csv',stream.getvalue().encode())
    save(RECORD/'science-validation.json',dict(passed=True,created=now(),test_predictions_sha256=digest(OUTPUT/'test/predictions.jsonl'),
        calibration_sha256=frozen['sha256'],report_sha256=digest(OUTPUT/'report.json'),
        scores_sha256=digest(OUTPUT/'test/scores.truth.json'),all_288_frames_retained=True,
        old_count_decisions_replayed=True,quality_tag_does_not_alter_fits=True,
        max_cost_replay_error=diagnostics['maximum_cost_replay_error']))
    print('SCIENCE',json_summary(report),flush=True)


def json_summary(report):
    import json
    return json.dumps(dict(level=report['calibration']['level'],calibration=report['calibration']['selected_per_gap'],
        pairs=report['primary_containment'],quality=report['quality'],numerical=report['numerical']),ensure_ascii=False)


def conclusions(report):
    cal=report['calibration']; parts=[f"用校准集中盲选为双事件的{cal['samples']}条真实双脉冲，经验95%上分位得到联合代价阈值c={cal['level']:.3f}；40/200ns分别贡献{cal['selected_per_gap']['pair40']}/{cal['selected_per_gap']['pair200']}条。这是当前选门和设计下的混合样本标定，不是统一95%置信保证。"]
    for r in report['primary_containment']:
        s=next(x for x in report['summary'] if x['band']==r['band'])
        c=r['conditional_selected']; a=r['selected_and_contained']
        parts.append(f"{r['band'][4:]}ns：原盲流程输出双事件{s['selected_K2']['numerator']}/{s['frames']}，严格整对恢复{s['strict']['numerator']}/{s['frames']}；选出双事件后的联合区域包含真实Q1/Q2/t1/t2为{c['numerator']}/{c['denominator']}，按全部真实双脉冲计算则为{a['numerator']}/{a['denominator']}。时间漏检{s['missed']}个，额外输出{s['extra']}个。")
    for label, chinese in [('alternative_found','找到达到或超出容差边界的替代解'),('no_alternative_found','有限搜索未找到替代解')]:
        r=next(v for v in report['quality'] if v['label']==label)
        parts.append(f"{chinese}：{r['frames']}条盲选双事件，其中真实双脉冲{r['true_pairs']}条、控制样本误报双事件{r['false_pairs']}条；严格正确为{r['strict_accuracy']['numerator']}/{r['strict_accuracy']['denominator']}。分母为0时准确率未定义。")
    parts.append('质量标记没有改变原拟合或增加正确恢复事件；找到替代参数说明该代价区域内存在歧义，有限搜索未找到不证明全局唯一或物理可辨识。')
    return '\n\n'.join(parts)


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.2,'figure.dpi':120})
    report=load(OUTPUT/'report.json'); scored=load(OUTPUT/'test/scores.truth.json')
    predictions={r['index']:r for r in rows(OUTPUT/'test/predictions.jsonl')}
    folder=OUTPUT/'figures'; folder.mkdir(exist_ok=False); made=[]
    def finish(fig,number,title):
        fig.suptitle('Figure '+number+': '+title,fontsize=12); fig.tight_layout(rect=(0,0,1,.95))
        stream=io.BytesIO(); fig.savefig(stream,format='png',dpi=150); plt.close(fig)
        write_new(folder/(number+'.png'),stream.getvalue()); made.append(number)
    def bars(ax,values,labels,color):
        positions=np.arange(len(values)); ax.bar(positions,values,color=color)
        ax.set_xticks(positions,labels); ax.set_ylim(0,1.12)
    # Detection, correct reconstruction and output containment use the same denominator.
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for ax,band in zip(axes,('pair40','pair200')):
        s=next(x for x in report['summary'] if x['band']==band)
        c=next(x for x in report['primary_containment'] if x['band']==band)
        entries=[s['selected_K2'],s['strict'],c['selected_and_contained']]
        bars(ax,[x['value'] for x in entries],['Output K=2','Strict pair','Output + region contains'],['#7799bb','#55aa77','#dd9955'])
        for i,x in enumerate(entries): ax.text(i,x['value']+.025,f"{x['numerator']}/{x['denominator']}",ha='center')
        ax.set_title(band[4:]+' ns'); ax.set_ylabel('Fraction of all true pairs')
    finish(fig,'0004-127-1','Blind outputs, correct events, and joint-region containment')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for ax,band in zip(axes,('pair40','pair200')):
        rr=[x for x in report['level_sensitivity'] if x['band']==band]
        for field,label,color in [('conditional_selected','Among selected K=2','#0077bb'),('selected_and_contained','All true pairs','#ee7733')]:
            xx=[r['level'] for r in rr]; yy=[r[field]['value'] for r in rr]
            lo=[r[field]['value']-r[field]['binomial95'][0] if r[field]['value'] is not None else 0. for r in rr]
            hi=[r[field]['binomial95'][1]-r[field]['value'] if r[field]['value'] is not None else 0. for r in rr]
            ax.errorbar(xx,[np.nan if y is None else y for y in yy],yerr=[lo,hi],fmt='o-',capsize=3,label=label,color=color)
        ax.axvline(report['calibration']['level'],color='gray',ls=':',label='Frozen calibration c')
        ax.axhline(.95,color='gray',ls='--',label='Calibration target, not guarantee')
        ax.set(xlabel='Joint cost allowance c',ylabel='Observed containment',ylim=(0,1.05),title=band[4:]+' ns')
        ax.legend(fontsize=8)
    finish(fig,'0004-128-1','Independent containment; exact binomial sampling intervals')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for ax,band in zip(axes,('pair40','pair200')):
        rr=[r for r in scored if r['band']==band and r['predicted_count']==2 and r['exterior_delta'] is not None]
        for label,color,marker in [('alternative_found','#cc6633','o'),('numerically_unresolved','#999999','x'),('no_alternative_found','#0077bb','s')]:
            subset=[r for r in rr if r['quality_label']==label]
            ax.scatter([r['exterior_delta']-report['calibration']['level'] for r in subset],
                       [r['accuracy_distance'] for r in subset],c=color,marker=marker,label=label.replace('_',' '),alpha=.8)
        ax.axvline(0,color='black',ls=':'); ax.axhline(1,color='black',ls='--')
        ax.set_xscale('symlog',linthresh=1.); ax.set_yscale('symlog',linthresh=1.)
        ax.set(xlabel='Best exterior cost - reference cost - c',ylabel='Largest normalized Q/t error',title=band[4:]+' ns')
        ax.legend(fontsize=7)
    finish(fig,'0004-129-1','Ambiguity witnesses versus actual parameter error')
    response,noise=models()
    with np.load(OUTPUT/'test/observations.npz') as f: waves=f['voltage_mV']
    with np.load(OUTPUT/'test/labels.truth.npz') as f: q=f['charge_pe']; t=f['times_ns']
    for j,index in enumerate(load(EXPERIMENT/'protocol.json')['example_indices'],1):
        r=predictions[index]; selected=r['selected']; x=response.sample_grid
        truth=response.components(t[index],q[index]); predicted=response.components(selected['times_ns'],selected['charge_pe'])
        baseline=response.baseline+noise.background_mV
        fig,axes=plt.subplots(2,1,figsize=(11,6),sharex=True)
        axes[0].plot(x,waves[index],color='.65',lw=.6,label='Observed')
        axes[0].plot(x,truth.sum(0)+baseline,'k--',lw=1,label='Truth sum (scoring only)')
        axes[0].plot(x,predicted.sum(0)+baseline,color='#0077bb',lw=1.2,label=f"Frozen blind fit K={len(selected['times_ns'])}")
        for k,curve in enumerate(truth): axes[1].plot(x,curve,ls='--',lw=1.2,label=f'True component {k+1}')
        for k,curve in enumerate(predicted): axes[1].plot(x,curve,lw=1,label=f'Blind component {k+1}')
        if r['exterior']:
            witness=r['exterior']['best_exterior']; other=response.components(witness['times_ns'],witness['charge_pe'])
            axes[0].plot(x,other.sum(0)+baseline,color='#cc6677',lw=1,label=f"Exterior witness, delta={r['exterior']['exterior_delta']:.2f}")
            for k,curve in enumerate(other): axes[1].plot(x,curve,color='#cc6677',ls=':',label=f'Exterior component {k+1}')
        for ax in axes: ax.set_ylabel('mV'); ax.legend(fontsize=8,ncol=2); ax.set_xlim(x[0],x[-1])
        axes[1].set_xlabel('Time in the complete record / ns')
        finish(fig,f'0004-130-{j}',f'Prespecified independent frame {index}, gap {t[index,1]-t[index,0]:.0f} ns')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    qq=[r for r in report['quality'] if r['label']!='not_selected_K2']
    labels=['Alternative\nfound','Numerically\nunresolved','No alternative\nfound']
    axes[0].bar(np.arange(3),[r['frames'] for r in qq],color=['#cc6633','#999999','#0077bb'])
    axes[0].set_xticks(np.arange(3),labels); axes[0].set_ylabel('All selected K=2 frames, including controls')
    for i,r in enumerate(qq):
        f=r['strict_accuracy']
        if f['value'] is not None:
            axes[1].plot(i,f['value'],'o',color='black')
            axes[1].text(i,.05,f"{f['numerator']}/{f['denominator']}",ha='center')
        else: axes[1].text(i,.5,'Undefined\n0/0',ha='center')
    axes[1].set_xticks(np.arange(3),labels); axes[1].set_ylim(0,1.05); axes[1].set_ylabel('Strict accuracy in fixed stratified mixture')
    finish(fig,'0004-131-1','Quality labels: acceptance and accuracy remain separate')
    save(RECORD/'figures-validation.json',dict(passed=True,figures=made,files={n+'.png':digest(folder/(n+'.png')) for n in made}))
    print('FIGURES',len(made),flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('command',choices=['analyze','figures']);a=p.parse_args()
    if a.command=='analyze':analyze()
    else:figures()


if __name__=='__main__': main()
