"""Staged scoring of spatial response controls, with no pair energy inference."""
from __future__ import annotations

import argparse
import csv

import numpy as np

from . import alpha_pair_spatial_study as study
from .alpha_pair_causal_controls import read, save, sha


def table(name,rows):
    folder=study.EXP/'tables';folder.mkdir(exist_ok=True)
    with (folder/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def single_report():
    count=study.verify(True);study.initialize()
    truth={r['parent']:r for r in read(study.EXP/'evaluation-mothers.json')}
    rows=[];raw=[]
    for task in read(study.EXP/'single-plan.json'):
        result=read(study.ART/'single/single'/f"{task['id']:04d}.json")
        with np.load(study.ART/'single/single'/f"{task['id']:04d}.npz") as f:
            assert all(np.isfinite(f[k]).all() for k in f.files)
            event=f['event_signal_mV'].mean(0)
        energy=truth[task['parent']]['energy_dep_MeV']
        model=study.MODELS[study.PRIMARY]
        curve=model.components([task['time_ns']],[energy])[0]
        before=model.readout.adc_centers+.25<=task['time_ns']
        row=dict(**task,energy_dep_MeV=energy,fitted_energy_MeV=result['fitted_energy_MeV'],
                 energy_bias_percent=100*(result['fitted_energy_MeV']/energy-1),converged=result['converged'],
                 anchor_area_error_percent=float(100*(curve.sum()/event.sum()-1)),
                 pre_event_mV=float(abs(curve[before]).max()))
        rows.append(row)
        if task['phase_ns']==0:
            item=dict(parent=task['parent'],condition=task['condition'],energy_dep_MeV=energy)
            with np.load(study.ART/'predictions'/f"single-{task['id']:04d}.npz") as predictions:
                for name in study.MODELS:
                    item[name+'_raw_area_error_percent']=float(100*(predictions[name+'_raw'].sum()/event.sum()-1))
            raw.append(item)
    groups=[];failures=[];gates=read(study.EXP/'protocol.json')['gates']
    for condition in ('weak','low','mid','high'):
        selected=[r for r in rows if r['condition']==condition]
        phases=[float(np.mean([r['energy_bias_percent'] for r in selected if r['phase_ns']==phase])) for phase in study.PHASES]
        selected_raw=[r for r in raw if r['condition']==condition]
        group=dict(condition=condition,mothers=4,energy_bias_percent=float(np.mean(phases)),
             phase_bias_percent=phases,phase_span_pp=float(np.ptp(phases)),
             raw_area_errors={name:float(np.mean([r[name+'_raw_area_error_percent'] for r in selected_raw])) for name in study.MODELS})
        groups.append(group)
        if abs(group['energy_bias_percent'])>gates['single_group_mean_energy_percent']:failures.append(condition+': isolated energy anchor')
        if group['phase_span_pp']>gates['single_phase_span_pp']:failures.append(condition+': ADC phase')
        if abs(group['raw_area_errors'][study.PRIMARY])>gates['raw_single_group_area_error_percent']:failures.append(condition+': raw spatial response')
    if max(r['pre_event_mV'] for r in rows)>gates['pre_event_mV']:failures.append('pre-event support')
    if not all(r['converged'] for r in rows):failures.append('single optimizer nonconvergence')
    with np.load(study.ART/'pilot/single/0000.npz') as pilot,np.load(study.ART/'single/single/0000.npz') as full:
        pilot_exact=all(np.array_equal(pilot[k],full[k][:2]) for k in pilot.files)
    assert pilot_exact
    summary=dict(passed=not failures,failures=failures,groups=groups,parents=16,phases=4,digitizations=512,
        source_input_model_hashes=count,pilot_arrays_exact=pilot_exact,
        limitation='Isolated anchored voltage invariant by construction; raw truth-time diagnostic additionally exercises spatial operator. Four exposed mothers per group, not individual resolution or fresh confirmation.')
    save(study.EXP/'single-gate.json',summary);save(study.ART/'single-report.json',summary)
    save(study.ART/'single-rows.json',rows);save(study.ART/'single-raw-rows.json',raw)
    table('single-controls.csv',rows);table('single-raw-controls.csv',raw)
    print('single gate',summary['passed'],failures,flush=True)
    for g in groups:print(g['condition'],'E bias',g['energy_bias_percent'],'raw',g['raw_area_errors'],flush=True)


def paired_intervals(areas,predictions,seed):
    rng=np.random.default_rng(seed)
    mother=rng.integers(0,4,(10000,4))
    repeat=rng.integers(0,8,(10000,4,8))
    samples=areas[mother[:,:,None],repeat].mean(axis=2)
    observed=(100*(samples[:,:,0]/samples[:,:,1]-1)).mean(axis=1)
    interval=lambda x:np.percentile(x,[2.5,97.5]).tolist()
    return interval(observed),{name:interval(np.asarray(values)[mother].mean(axis=1)-observed) for name,values in predictions.items()}


def pair_report():
    count=study.verify(True);study.initialize()
    predicted={r['id']:r for r in read(study.EXP/'predictions.json')['rows'] if r['kind']=='pair'}
    rows=[];display={};areas={}
    for task in read(study.EXP/'pair-plan.json'):
        cid=task['id'];weak=task['parents'].index(task['block'])
        with np.load(study.ART/'pair/pair'/f'{cid:04d}.npz') as f:
            arrays={k:f[k].copy() for k in f.files}
            assert all(np.isfinite(v).all() for v in arrays.values())
        a=np.column_stack([arrays['on_shared'][:,weak].sum(-1),arrays['on_isolated'][:,weak].sum(-1)])
        areas[cid]=a
        actual=float(100*(a[:,0].mean()/a[:,1].mean()-1))
        off=float(100*(arrays['off_shared'][:,weak].sum(-1).mean()/arrays['off_isolated'][:,weak].sum(-1).mean()-1))
        p=predicted[cid]['predictions']
        with np.load(study.ART/'predictions'/f'pair-{cid:04d}.npz') as f:
            minimum=min(float(f[name+'_joint'].min()) for name in study.MODELS)
        row=dict(**task,actual_percent=actual,recovery_off_percent=off,min_model_mV=minimum,
            spatial_grid_change_pp=p['spatial12']-p['spatial6'],time_grid_change_pp=p['spatial12']-p['spatial12_fine'])
        for name,value in p.items():row[name+'_prediction_percent']=value;row[name+'_error_pp']=value-actual
        rows.append(row)
        if task['block']==5:
            display[f'case{cid}_shared']=arrays['on_shared'].mean(0)
            display[f'case{cid}_isolated']=arrays['on_isolated'].mean(0)
    groups=[];failures=[];gates=read(study.EXP/'protocol.json')['gates']
    for groupid,(condition,gap) in enumerate((('high_weak',8.),('high_weak',40.),('high_weak',200.),('high_weak',1200.),('weak_high',40.))):
        selected=[r for r in rows if (r['condition'],r['gap_ns'])==(condition,gap)]
        values={name:[r[name+'_prediction_percent'] for r in selected] for name in study.MODELS}
        ci,error_ci=paired_intervals(np.stack([areas[r['id']] for r in selected]),values,715900000+groupid)
        actual=float(np.mean([r['actual_percent'] for r in selected]))
        group=dict(condition=condition,gap_ns=gap,mothers=4,repeats=8,actual_percent=actual,actual_ci95=ci,
            recovery_off_percent=float(np.mean([r['recovery_off_percent'] for r in selected])),
            models={name:dict(prediction_percent=float(np.mean(v)),error_pp=float(np.mean(v)-actual),error_ci95=error_ci[name]) for name,v in values.items()},
            spatial_grid_max_abs_pp=max(abs(r['spatial_grid_change_pp']) for r in selected),
            time_grid_max_abs_pp=max(abs(r['time_grid_change_pp']) for r in selected))
        groups.append(group)
        tag=condition+'_'+str(int(gap))
        if abs(group['models'][study.PRIMARY]['error_pp'])>gates['shared_weak_area_error_pp']:failures.append(tag+': spatial shared response')
        if abs(group['recovery_off_percent'])>gates['recovery_off_abs_pp']:failures.append(tag+': recovery-off')
        if group['spatial_grid_max_abs_pp']>gates['spatial_grid_change_pp']:failures.append(tag+': spatial grid')
        if group['time_grid_max_abs_pp']>gates['time_grid_change_pp']:failures.append(tag+': time grid')
    if min(r['min_model_mV'] for r in rows)<-1e-10:failures.append('negative model component')
    with np.load(study.ART/'pilot/pair/0001.npz') as pilot,np.load(study.ART/'pair/pair/0001.npz') as full:
        exact=all(np.array_equal(pilot[k],full[k][:2]) for k in pilot.files)
    assert exact
    summary=dict(passed=not failures,failures=failures,groups=groups,primary=study.PRIMARY,
        source_input_model_hashes=count,pilot_arrays_exact=exact,production_digitizations=960,
        default_changed=False,pair_energy_inference=False,new_ROOT=0,
        uncertainty='10000 paired mother/digitizer nested resamples; frozen128-mother prior uncertainty not propagated; four historically exposed evaluation mother pairs only.')
    save(study.EXP/'shared-gate.json',summary);save(study.ART/'pair-report.json',summary)
    save(study.ART/'pair-rows.json',rows);table('pair-controls.csv',rows)
    np.savez_compressed(study.ART/'pair-display.npz',**display)
    print('shared gate',summary['passed'],failures,flush=True)
    for g in groups:print(g['condition'],g['gap_ns'],'actual',g['actual_percent'],'models',g['models'],flush=True)


def validate_statistics():
    areas=np.full((4,8,2),100.)
    areas[:,:,0]=85.
    ci,errors=paired_intervals(areas,{'uniform':np.full(4,-8.),'spatial':np.full(4,-15.)},715910001)
    assert np.allclose(ci,[-15.,-15.],atol=1e-12)
    assert np.allclose(errors['uniform'],[7.,7.],atol=1e-12)
    assert np.allclose(errors['spatial'],[0.,0.],atol=1e-12)
    save(study.REC/'resampling-check.json',dict(passed=True,constant_response_ci=ci,known_model_error_ci=errors))
    print('Paired response statistics controls passed',flush=True)


def plots():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    study.verify(True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    figures=[]
    def finish(fig,stem):
        for extension in ('png','svg'):
            p=study.ART/(stem+'.'+extension)
            fig.savefig(p,dpi=180,bbox_inches='tight')
            figures.append(dict(path=str(p.relative_to(study.ROOT)),sha256=sha(p)))
        plt.close(fig)
    with np.load(study.EXP/'spatial-profile.npz') as f:
        maps=[f['cell_probability']/.0025,f['profile_6'],f['profile_12']/.25]
    fig,axes=plt.subplots(1,3,figsize=(12,4.2),layout='constrained')
    vmax=max(p.max() for p in maps)
    for ax,p,label in zip(axes,maps,('(a) Development mean','(b) 6 x 6 D4 prior','(c) 12 x 12 D4 prior')):
        im=ax.imshow(p,origin='lower',extent=(-3,3,-3,3),vmin=0,vmax=vmax,cmap='magma')
        ax.set(title=label,xlabel='x (mm)',ylabel='y (mm)')
    fig.colorbar(im,ax=axes,label='PDE-weighted density (mm$^{-2}$)',shrink=.8)
    fig.suptitle('Fig. 0007-82-1: Frozen illumination from 128 disjoint development mothers')
    finish(fig,'0007-82-1_spatial_prior')

    single=read(study.ART/'single-report.json')
    single_rows=read(study.ART/'single-rows.json')
    raw_rows=read(study.ART/'single-raw-rows.json')
    labels=['weak','low','mid','high'];x=np.arange(4)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    axes[0].bar(x,[g['energy_bias_percent'] for g in single['groups']],color='#5479AD',alpha=.6,label='Group mean')
    for i,condition in enumerate(labels):
        parents=sorted({r['parent'] for r in single_rows if r['condition']==condition})
        vals=[np.mean([r['energy_bias_percent'] for r in single_rows if r['parent']==parent]) for parent in parents]
        axes[0].scatter(i+np.linspace(-.13,.13,4),vals,s=20,color='#203D6A',label='Individual mother' if i==0 else None)
    for k,(name,color,label) in enumerate((('uniform','#355C9B','Uniform raw'),('spatial12','#238B67','Spatial raw'))):
        xx=x+(k-.5)*.22
        axes[1].plot(xx,[g['raw_area_errors'][name] for g in single['groups']],'o',color=color,label=label,ms=8)
        for i,condition in enumerate(labels):
            vals=[r[name+'_raw_area_error_percent'] for r in raw_rows if r['condition']==condition]
            axes[1].scatter(xx[i]+np.linspace(-.04,.04,4),vals,s=13,color=color,alpha=.5)
    for ax in axes:
        ax.axhline(0,color='.5',lw=.8)
        for y in (-5,5):ax.axhline(y,color='.6',ls='--',lw=.8)
        ax.set_xticks(x,labels);ax.grid(axis='y',alpha=.15);ax.legend(fontsize=8)
    axes[0].set(title='(a) Isolated mean-voltage energy fit',ylabel='Energy bias (%)')
    axes[1].set(title='(b) Raw truth-time forward check',ylabel='Predicted / physical area - 1 (%)')
    fig.suptitle('Fig. 0007-83-1: Separate anchored energy and raw spatial-response controls')
    finish(fig,'0007-83-1_single_controls')
    if not (study.ART/'pair-report.json').exists():
        save(study.REC/'figures.json',figures)
        return

    pair=read(study.ART/'pair-report.json');groups=pair['groups'];x=np.arange(5)
    labels=['H-W 8','H-W 40','H-W 200','H-W 1200','W-H 40']
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    observed=np.array([g['actual_percent'] for g in groups]);ci=np.array([g['actual_ci95'] for g in groups])
    axes[0].vlines(x,ci[:,0],ci[:,1],color='#303030',lw=2)
    axes[0].plot(x,observed,'o',color='#303030',label='Physical digitizer')
    for index,(name,color,label) in enumerate((('uniform','#355C9B','Uniform'),('spatial12','#238B67','Spatial 12 x 12'))):
        axes[0].plot(x,[g['models'][name]['prediction_percent'] for g in groups],'--',color=color,label=label)
        xx=x+(index-.5)*.16
        error=np.array([g['models'][name]['error_pp'] for g in groups])
        interval=np.array([g['models'][name]['error_ci95'] for g in groups])
        axes[1].vlines(xx,interval[:,0],interval[:,1],color=color,lw=2)
        axes[1].plot(xx,error,'o',color=color,label=label)
    for ax in axes:
        ax.axhline(0,color='.5',lw=.8)
        ax.set_xticks(x,labels,rotation=20);ax.grid(axis='y',alpha=.15)
        ax.set_xlabel('Order and gap (ns)')
    axes[0].set(title='(a) Weak-component area change',ylabel='Shared / isolated - 1 (%)')
    axes[1].set(title='(b) Prediction error; paired descriptive intervals',ylabel='Error (percentage points)')
    axes[0].legend(loc='lower center',bbox_to_anchor=(.58,.01),fontsize=8)
    axes[1].legend(fontsize=8)
    for value in (-5,5):axes[1].axhline(value,color='.6',ls=':',lw=.8)
    fig.suptitle('Fig. 0007-84-1: Shared response on four held-out historical mother pairs')
    finish(fig,'0007-84-1_shared_response')

    with np.load(study.ART/'pair-display.npz') as f:display={k:f[k].copy() for k in f.files}
    fig,axes=plt.subplots(2,3,figsize=(13,6.8),layout='constrained');t=np.arange(1024)*2+.25
    for column,(cid,weak,label) in enumerate(((0,1,'H-W 8 ns'),(1,1,'H-W 40 ns'),(4,0,'W-H 40 ns'))):
        shared=display[f'case{cid}_shared'];isolated=display[f'case{cid}_isolated']
        axes[0,column].plot(t,shared.sum(0),color='#303030',label='Physical mean',lw=1.2)
        axes[1,column].plot(t,shared[weak]-isolated[weak],color='#303030',label='Physical mean',lw=1.2)
        with np.load(study.ART/'predictions'/f'pair-{cid:04d}.npz') as f:
            for name,color,label_model in (('uniform','#355C9B','Uniform'),('spatial12','#238B67','Spatial12')):
                joint=f[name+'_joint'];solo=f[name+'_solo']
                axes[0,column].plot(t,joint.sum(0),'--',color=color,label=label_model,lw=1.1)
                axes[1,column].plot(t,joint[weak]-solo[weak],'--',color=color,label=label_model,lw=1.1)
        for row in range(2):
            axes[row,column].set(xlim=(0,2048),xlabel='Time (ns)',ylabel='Voltage (mV)')
            axes[row,column].grid(alpha=.15)
        axes[0,column].set_title(f'({chr(97+column)}) {label}: total')
        axes[1,column].set_title(f'({chr(100+column)}) {label}: weak difference')
        axes[1,column].axhline(0,color='.6',lw=.7)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Fig. 0007-85-1: Complete observed window; first fixed evaluation mother pair')
    finish(fig,'0007-85-1_waveforms')
    save(study.REC/'figures.json',figures)
    print('Saved four numbered spatial-control figures',flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=('single','pair','plots','validate'))
    args=parser.parse_args()
    if args.mode=='validate':validate_statistics()
    elif args.mode=='single':single_report()
    elif args.mode=='pair':pair_report()
    else:plots()


if __name__=='__main__':main()
