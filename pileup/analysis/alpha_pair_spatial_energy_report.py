"""Paired energy endpoints for a fixed known-K/t historical-mother pilot."""
from __future__ import annotations

import argparse
import csv
import numpy as np

from . import alpha_pair_spatial_energy_study as study
from .alpha_pair_causal_controls import read, save, sha


def table(name, rows):
    folder = study.EXP/'tables';folder.mkdir(exist_ok=True)
    with (folder/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def load_rows(kind):
    study.verify(True)
    cases={c['id']:c for c in read(study.EXP/'cases.json')}
    tasks=[t for t in read(study.EXP/'fit-plan.json') if t['kind']==kind]
    folder=study.ART/('single' if kind=='single' else 'pairs')
    rows=[];by_key={}
    for task in tasks:
        row=read(folder/f"{task['id']:04d}.json")
        case=cases[row['case']]
        assert row['read_inputs']==[task['waveform']]
        assert not row['true_energy_input']
        e=np.asarray(row['energy_MeV']);truth=np.asarray(case['energy_dep_MeV'])
        assert np.isfinite(e).all() and np.isfinite(row['cost'])
        error=100*(e/truth-1)
        row.update(condition=case['condition'],gap_ns=case['gap_ns'],block=case['block'],
                   true_energy_MeV=truth.tolist(),relative_error_percent=error.tolist(),
                   strict=bool(row['converged'] and np.all(abs(error)<=10.)))
        if kind=='pair':
            weak=0 if case['condition']=='weak_high' else 1
            row.update(weak_index=weak,weak_bias_percent=float(error[weak]),
                       high_bias_percent=float(error[1-weak]))
        rows.append(row);by_key[(row['case'],row['repeat'],row['method'])]=row
    for case,repeat,method in by_key:
        if method!='uniform':continue
        a=by_key[(case,repeat,'uniform')];b=by_key[(case,repeat,'spatial12')]
        assert a['metric_signature']==b['metric_signature']
        assert a['initial_MeV']==b['initial_MeV'] and a['start_energy_MeV']==b['start_energy_MeV']
        if kind=='single':
            assert a['energy_MeV']==b['energy_MeV'] and a['cost']==b['cost'] and a['converged']==b['converged']
    return rows,by_key


def single_report():
    rows,_=load_rows('single')
    groups=[]
    for condition in ('weak','low','mid','high'):
        avg=[r for r in rows if r['condition']==condition and r['method']=='spatial12' and r['repeat']==-1]
        noisy=[r for r in rows if r['condition']==condition and r['method']=='spatial12' and r['repeat']>=0]
        a=np.array([r['relative_error_percent'][0] for r in avg])
        e=np.array([r['relative_error_percent'][0] for r in noisy])
        groups.append(dict(condition=condition,mothers=4,individual_n=len(noisy),mean_bias_percent=float(a.mean()),
            individual_bias_percent=float(e.mean()),individual_sd_percent=float(e.std(ddof=1)),
            individual_rmse_percent=float(np.sqrt(np.mean(e*e))),strict=sum(r['strict'] for r in noisy),
            boundary=sum(r['boundary'] for r in noisy),nonconverged=sum(not r['converged'] for r in noisy)))
    limit=read(study.EXP/'protocol.json')['gates']['single_group_mean_bias_percent']
    passed=all(abs(g['mean_bias_percent'])<=limit for g in groups) and all(r['converged'] for r in rows)
    legacy=[]
    cases=read(study.EXP/'cases.json')
    for case in cases:
        if case['kind']!='single':continue
        old=read(study.PREVART/'single/single'/f"{case['source_id']:04d}.json")
        row=next(r for r in rows if r['case']==case['id'] and r['repeat']==-1 and r['method']=='spatial12')
        legacy.append(abs(row['energy_MeV'][0]/old['fitted_energy_MeV']-1))
    assert max(legacy)<1e-10
    report=dict(passed=bool(passed),groups=groups,method_outputs_exact=True,
        mean_fit_legacy_max_relative=float(max(legacy)),fit_rows=len(rows),individual_waveforms=128,
        mean_waveforms=16,known_K_t=True,blind_detection_measured=False)
    save(study.EXP/'single-gate.json',report);save(study.ART/'single-report.json',report)
    save(study.ART/'single-rows.json',rows);table('single-summary.csv',groups)
    print('Single energy controls',passed,groups,flush=True)


def paired_intervals(errors, success, seed):
    """Shape[4 mothers, repetitions, 2 paired methods]; no row independence."""
    rng=np.random.default_rng(seed)
    mothers=rng.integers(0,4,(10000,4))
    repeats=rng.integers(0,errors.shape[1],(10000,4,errors.shape[1]))
    sampled=errors[mothers[:,:,None],repeats]
    ok=success[mothers[:,:,None],repeats]
    bias=sampled.mean((1,2));rmse=np.sqrt((sampled*sampled).mean((1,2)))
    gain=100*(ok[:,:, :,1]-ok[:,:, :,0]).mean((1,2))
    ci=lambda v:np.percentile(v,[2.5,97.5]).tolist()
    return dict(weak_bias_ci95=[ci(bias[:,i]) for i in range(2)],
                weak_rmse_ci95=[ci(rmse[:,i]) for i in range(2)],
                strict_gain_ci95_pp=ci(gain))


def validate_statistics():
    errors=np.empty((4,2,2));errors[:,:,0]=20.;errors[:,:,1]=-5.
    success=np.zeros_like(errors);success[:,:,1]=1.
    result=paired_intervals(errors,success,716900001)
    assert np.allclose(result['weak_bias_ci95'],[[20,20],[-5,-5]])
    assert np.allclose(result['weak_rmse_ci95'],[[20,20],[5,5]])
    assert np.allclose(result['strict_gain_ci95_pp'],[100,100])
    save(study.REC/'statistics-check.json',dict(passed=True,controls=result))


def pair_report():
    rows,index=load_rows('pair')
    groups=[]
    conditions=[('high_weak',8.),('high_weak',40.),('high_weak',200.),('high_weak',1200.),('weak_high',40.)]
    for domain in ('individual','mean'):
        for gid,(condition,gap) in enumerate(conditions):
            selected=[r for r in rows if r['condition']==condition and r['gap_ns']==gap and
                      ((r['repeat']>=0) if domain=='individual' else (r['repeat']==-1))]
            case_ids=sorted({r['case'] for r in selected})
            repeats=[0,1] if domain=='individual' else [-1]
            errors=np.array([[[index[(c,j,m)]['weak_bias_percent'] for m in study.METHODS]
                              for j in repeats] for c in case_ids])
            success=np.array([[[float(index[(c,j,m)]['strict']) for m in study.METHODS]
                               for j in repeats] for c in case_ids])
            ci=paired_intervals(errors,success,716910000+gid+(0 if domain=='individual' else 100))
            methods={}
            for i,method in enumerate(study.METHODS):
                method_rows=[r for r in selected if r['method']==method]
                weak=np.array([r['weak_bias_percent'] for r in method_rows]);high=np.array([r['high_bias_percent'] for r in method_rows])
                methods[method]=dict(n=len(method_rows),strict=sum(r['strict'] for r in method_rows),
                    weak_bias_percent=float(weak.mean()),weak_sd_percent=float(weak.std(ddof=1)),
                    weak_rmse_percent=float(np.sqrt(np.mean(weak**2))),weak_bias_ci95=ci['weak_bias_ci95'][i],
                    weak_rmse_ci95=ci['weak_rmse_ci95'][i],high_bias_percent=float(high.mean()),
                    high_sd_percent=float(high.std(ddof=1)),high_rmse_percent=float(np.sqrt(np.mean(high**2))),
                    nonconverged=sum(not r['converged'] for r in method_rows),boundary=sum(r['boundary'] for r in method_rows))
            improved=int(np.sum((success[:,:,1]==1)&(success[:,:,0]==0)))
            worsened=int(np.sum((success[:,:,0]==1)&(success[:,:,1]==0)))
            groups.append(dict(domain=domain,condition=condition,gap_ns=gap,mothers=4,repeats=len(repeats),methods=methods,
                improved=improved,worsened=worsened,strict_gain_pp=float(100*(success[:,:,1]-success[:,:,0]).mean()),
                strict_gain_ci95_pp=ci['strict_gain_ci95_pp']))
    gates=read(study.EXP/'protocol.json')['gates']
    short=[r for r in rows if r['condition']=='high_weak' and r['gap_ns'] in (8.,40.) and r['repeat']>=0]
    metrics={}
    for method in study.METHODS:
        s=[r for r in short if r['method']==method]
        metrics[method]=dict(n=len(s),strict=sum(r['strict'] for r in s),
                             weak_rmse_percent=float(np.sqrt(np.mean([r['weak_bias_percent']**2 for r in s]))))
    gain=100*(metrics['spatial12']['strict']-metrics['uniform']['strict'])/metrics['uniform']['n']
    ratio=metrics['spatial12']['weak_rmse_percent']/metrics['uniform']['weak_rmse_percent']
    failures=[]
    if gain<gates['next_stage_short_strict_gain_pp']:failures.append('short strict gain <10 pp')
    if ratio>gates['next_stage_short_rmse_ratio']:failures.append('short pooled RMSE ratio >0.95')
    for g in groups:
        if g['condition']!='high_weak' or g['gap_ns'] not in (8.,40.):continue
        if g['domain']=='individual':
            if g['methods']['spatial12']['weak_rmse_percent']/g['methods']['uniform']['weak_rmse_percent']>gates['next_stage_each_short_rmse_ratio']:
                failures.append(f"{g['gap_ns']:g} ns individual RMSE ratio >1.1")
        elif abs(g['methods']['spatial12']['weak_bias_percent'])>gates['next_stage_short_mean_abs_bias_percent']:
            failures.append(f"{g['gap_ns']:g} ns mean weak bias >10%")
    pilot=read(study.REC/'pilot.json')
    pilot_checks=[]
    for p in pilot['rows']:
        row=next(r for r in rows if r['id']==p['id'])
        assert p['energy_MeV']==row['energy_MeV'] and p['cost']==row['cost']
        pilot_checks.append(p['id'])
    report=dict(groups=groups,short_individual=metrics,short_strict_gain_pp=gain,short_weak_rmse_ratio=ratio,
        recommend_independent_confirmation=not failures,decision_failures=failures,
        fit_rows=len(rows),individual_waveforms=40,mean_waveforms=20,pilot_repeated_fits_exact=pilot_checks,
        common_weights_starts_exact=True,known_K_t=True,new_ROOT=0,default_changed=False,
        blind_detection_metrics=None,scope='Four historically exposed mother pairs;2 repeats/condition for individual fits;8-repeat means. No independent confirmation or physical resolution limit.')
    save(study.ART/'pair-rows.json',rows);save(study.ART/'pair-report.json',report)
    save(study.EXP/'decision.json',report)
    flat=[dict(domain=g['domain'],condition=g['condition'],gap_ns=g['gap_ns'],method=m,**{k:v for k,v in stats.items() if not isinstance(v,list)})
          for g in groups for m,stats in g['methods'].items()]
    table('pair-summary.csv',flat)
    detailed=[dict(id=r['id'],case=r['case'],condition=r['condition'],gap_ns=r['gap_ns'],repeat=r['repeat'],method=r['method'],
        E1_true_MeV=r['true_energy_MeV'][0],E2_true_MeV=r['true_energy_MeV'][1],E1_fit_MeV=r['energy_MeV'][0],E2_fit_MeV=r['energy_MeV'][1],
        weak_bias_percent=r['weak_bias_percent'],high_bias_percent=r['high_bias_percent'],strict=r['strict'],boundary=r['boundary'],
        converged=r['converged'],cost=r['cost'],evaluations=r['evaluations'],seconds=r['seconds']) for r in rows]
    table('pair-fits.csv',detailed)
    print('Independent-confirmation recommendation:',not failures,failures,flush=True)
    for g in groups:
        print(g['domain'],g['condition'],g['gap_ns'],g['methods'],flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['validate','single','pair','plots'])
    args=parser.parse_args()
    {'validate':validate_statistics,'single':single_report,'pair':pair_report,'plots':plots}[args.mode]()


def plots():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    study.verify(True);study.initialize()
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    figures=[]
    def finish(fig,stem):
        for extension in ('png','svg'):
            path=study.ART/(stem+'.'+extension)
            fig.savefig(path,dpi=180,bbox_inches='tight')
            figures.append(dict(path=str(path.relative_to(study.ROOT)),sha256=sha(path)))
        plt.close(fig)
    single=read(study.ART/'single-report.json')
    single_rows=read(study.ART/'single-rows.json')
    labels=['weak','low','mid','high'];x=np.arange(4)
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for i,condition in enumerate(labels):
        group=single['groups'][i]
        values=[r['relative_error_percent'][0] for r in single_rows if r['condition']==condition and
                r['method']=='spatial12' and r['repeat']>=0]
        axes[0].scatter(i+np.linspace(-.15,.15,len(values)),values,s=12,alpha=.55,color='#355C9B')
        axes[0].plot(i,group['individual_bias_percent'],'D',color='black',ms=6)
    axes[0].set(ylabel='Individual E error (%)',title='(a) 4 mothers x 8 repeats per group')
    axes[1].bar(x-.16,[g['mean_bias_percent'] for g in single['groups']],.32,label='8-repeat mean fit',color='#238B67')
    axes[1].bar(x+.16,[g['individual_bias_percent'] for g in single['groups']],.32,label='Individual-fit mean',color='#355C9B')
    axes[1].set(ylabel='Group mean E error (%)',title='(b) Uniform and spatial outputs identical')
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.axhline(0,color='.5',lw=.8);ax.set_xticks(x,labels);ax.grid(axis='y',alpha=.15)
    for y in (-5,5):axes[1].axhline(y,color='.6',ls='--',lw=.8)
    fig.suptitle('Fig. 0007-88-1: Isolated energy controls with known event times')
    finish(fig,'0007-88-1_single_energy')

    report=read(study.ART/'pair-report.json')
    labels=['H-W 8','H-W 40','H-W 200','H-W 1200','W-H 40'];x=np.arange(5)
    colors={'uniform':'#355C9B','spatial12':'#238B67'}
    names={'uniform':'Uniform','spatial12':'Spatial 12 x 12'}
    groups=[g for g in report['groups'] if g['domain']=='individual']
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for k,method in enumerate(study.METHODS):
        xx=x+(k-.5)*.22
        for ax,field in zip(axes.ravel()[:3],['weak_bias_percent','weak_sd_percent','weak_rmse_percent']):
            ax.plot(xx,[g['methods'][method][field] for g in groups],'o-',color=colors[method],label=names[method])
        stats=[g['methods'][method] for g in groups]
        values=[100*s['strict']/s['n'] for s in stats]
        axes[1,1].bar(x+(k-.5)*.34,values,.32,color=colors[method],label=names[method])
        for pos,val,s in zip(x+(k-.5)*.34,values,stats):
            axes[1,1].text(pos,val+2,f"{s['strict']}/{s['n']}",ha='center',fontsize=8)
    for ax,title,ylabel in zip(axes.ravel(),
        ['(a) Weak energy bias','(b) Weak relative-error spread','(c) Weak energy RMSE','(d) Both energies within 10%, converged'],
        ['Mean E error (%)','SD of E error (%)','RMSE (%)','Strict pair success (%)']):
        ax.set(title=title,ylabel=ylabel);ax.set_xticks(x,labels,rotation=18);ax.grid(axis='y',alpha=.15)
    axes[0,0].axhline(0,color='.5',lw=.8);axes[1,1].set_ylim(0,112)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Fig. 0007-89-1: Individual noisy pairs; four mothers x two repeats per condition')
    finish(fig,'0007-89-1_noisy_pair_energy')

    groups=[g for g in report['groups'] if g['domain']=='mean']
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
    for k,method in enumerate(study.METHODS):
        for ax,field in zip(axes,['weak_bias_percent','high_bias_percent']):
            ax.plot(x+(k-.5)*.1,[g['methods'][method][field] for g in groups],'o-',color=colors[method],label=names[method])
    for ax,title in zip(axes,['(a) Weak component','(b) High component']):
        ax.set(title=title,ylabel='Mean E error (%)');ax.set_xticks(x,labels,rotation=18)
        ax.axhline(0,color='.5',lw=.8);ax.grid(axis='y',alpha=.15)
        for y in (-10,10):ax.axhline(y,color='.6',lw=.8,ls='--')
    axes[0].legend(fontsize=8)
    fig.suptitle('Fig. 0007-90-1: Fits to eight-repeat mean voltages; four mothers per condition')
    finish(fig,'0007-90-1_mean_pair_energy')

    rows=read(study.ART/'pair-rows.json');cases={c['id']:c for c in read(study.EXP/'cases.json')}
    fig,axes=plt.subplots(2,3,figsize=(13,7),layout='constrained')
    time=study.MODELS['uniform'].readout.adc_centers
    for col,cid in enumerate((16,17,20)):
        case=cases[cid];weak=0 if case['condition']=='weak_high' else 1
        with np.load(study.ART/'waveforms'/f'{cid:04d}.npz') as f:
            observed=f['voltage_mV'].mean(0)
        with np.load(study.ART/'diagnostics'/f'{cid:04d}.npz') as f:
            truth=f['event_signal_mV'].mean(0)
        axes[0,col].plot(time,observed-study.MODELS['uniform'].background,color='.2',lw=1,label='Physical mean')
        axes[1,col].plot(time,truth[weak],color='.2',lw=1,label='Truth ancestry mean')
        for method in study.METHODS:
            row=next(r for r in rows if r['case']==cid and r['repeat']==-1 and r['method']==method)
            curve=study.MODELS[method].components(case['times_ns'],row['energy_MeV'])
            axes[0,col].plot(time,curve.sum(0),ls='--',color=colors[method],lw=1.3,label=names[method])
            axes[1,col].plot(time,curve[weak],ls='--',color=colors[method],lw=1.3,label=names[method])
        label=('H-W' if weak else 'W-H')+f" {case['gap_ns']:g} ns"
        axes[0,col].set_title(f'({chr(97+col)}) {label}: total')
        axes[1,col].set_title(f'({chr(100+col)}) {label}: weak component')
        for ax in axes[:,col]:
            ax.set(xlabel='Time (ns)',ylabel='Voltage (mV)',xlim=(0,2048));ax.grid(alpha=.15)
    axes[0,0].legend(fontsize=8);axes[1,0].legend(fontsize=8)
    fig.suptitle('Fig. 0007-91-1: First fixed mother pair; complete observed window, fitted E and known t')
    finish(fig,'0007-91-1_energy_waveforms')
    save(study.REC/'figures.json',figures)


if __name__=='__main__':main()
