"""Score native-forward diagnostics; label access is confined to reporting."""
import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from . import alpha_pair_recovery_forward_study as s

METHODS=('new_mean_new_noise','recovery128',*s.METHODS)
NAMES=['Existing voltage GLS','RL128 + recovery','Native OLS','Native GLS','Native calibrated GLS']
COLORS=['#20252b','#bd9cbb','#4799bc','#d7a242','#d66543']
KEYS=['weak','low','mid','high','high_weak_8','high_weak_40','high_weak_200','weak_high_40']
LABELS=['Single weak','Single low','Single mid','Single high','High→weak 8 ns','High→weak 40 ns','High→weak 200 ns','Weak→high 40 ns']


def key(r):return r['condition'] if r['kind']=='single' else r['condition']+'_'+str(int(r['gap_ns']))
def csv_write(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def save(fig,n,name):
    for ext in ('png','svg'):fig.savefig(s.ART/f'0007-{n}-1_{name}.{ext}',dpi=160,bbox_inches='tight')
    plt.close(fig)


def main():
    labels=s.read(s.ART/'labels.json');rows=[];components=[];results={};oldfits={};seed_diffs=[]
    for r in labels:
        i=r['id'];new=s.read(s.ART/'cases'/f'{i:04d}.json');results[i]=new
        old=s.read(s.ROOT/'jupyter/artifacts/alpha_pair_energy_calibrated_20261030/cases'/f'{i:04d}.json')['fits']['new_mean_new_noise']
        rl=s.read(s.PREVART/'cases'/f'{i:04d}.json')['fits']['recovery128'];oldfits[i]=old
        assert new['worker_truth_read_guard'] and not new['true_E_Q_input']
        seed_diffs.append(float(np.max(np.abs(np.array(new['seed']['energy_MeV'])/rl['energy_MeV']-1))))
        fits={'new_mean_new_noise':old,'recovery128':rl,**new['fits']}
        for method,fit in fits.items():
            errors=np.array(fit['energy_MeV'])/r['energy_MeV']-1
            rows.append(dict(id=i,condition=key(r),method=method,success5=bool(np.all(abs(errors)<=.05)),
                             success10=bool(np.all(abs(errors)<=.1)),success20=bool(np.all(abs(errors)<=.2)),
                             boundary=fit['boundary'],converged=fit['converged']))
            for j,error in enumerate(errors):
                components.append(dict(id=i,condition=key(r),method=method,component=j,weak=j==int(np.argmin(r['energy_MeV'])),
                    true_energy_MeV=r['energy_MeV'][j],fitted_energy_MeV=fit['energy_MeV'][j],relative_error=float(error)))
    outcomes={(r['id'],r['method']):r['success10'] for r in rows};summary=[]
    for k in KEYS:
        for method in METHODS:
            selected=[r for r in rows if r['condition']==k and r['method']==method]
            err=np.array([r['relative_error'] for r in components if r['condition']==k and r['method']==method and r['weak']])*100
            ids=[r['id'] for r in selected]
            summary.append(dict(condition=k,method=method,frames=len(selected),
                success5=sum(r['success5'] for r in selected),success10=sum(r['success10'] for r in selected),success20=sum(r['success20'] for r in selected),
                bias_pct=float(err.mean()),median_pct=float(np.median(err)),sigma68_pct=float((np.percentile(err,84)-np.percentile(err,16))/2),
                rmse_pct=float(np.sqrt(np.mean(err**2))),boundary=sum(r['boundary'] for r in selected),nonconverged=sum(not r['converged'] for r in selected),
                gained_vs_gls=sum(outcomes[i,method] and not outcomes[i,METHODS[0]] for i in ids),
                lost_vs_gls=sum(not outcomes[i,method] and outcomes[i,METHODS[0]] for i in ids)))
    lookup={(r['condition'],r['method']):r for r in summary}
    totals={m:sum(lookup[k,m]['success10'] for k in KEYS[4:]) for m in METHODS}
    report=dict(scope='Pure alpha, known TRUE K/t, same exposed64 diagnostic frames; not new validation',frames=len(labels),
        new_fits=len(s.METHODS)*len(labels),reference_fits=2*len(labels),component_scores=len(components),summary=summary,
        pair_totals=totals,primary='native_gls_calibrated',seed_max_relative_difference_from_python_reference=max(seed_diffs),
        new_nonconverged=sum(not r['converged'] for r in rows if r['method'] in s.METHODS),
        default_changed=False,expanded_to480=False)
    s.write(s.ART/'report.json',report);table=s.EXP/'tables';table.mkdir(exist_ok=True)
    csv_write(table/'summary.csv',summary);csv_write(table/'frames.csv',rows);csv_write(table/'components.csv',components)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(14,5),constrained_layout=True)
    for panel,(ax,keys) in enumerate(zip(axes,(KEYS[:4],KEYS[4:]))):
        for j,m in enumerate(METHODS):
            x=np.arange(4)+(j-2)*.15;values=[lookup[k,m]['success10'] for k in keys]
            ax.bar(x,values,width=.145,color=COLORS[j],label=NAMES[j])
            for xx,v in zip(x,values):ax.text(xx,v+.08,str(v),ha='center',fontsize=8)
        ax.set(xticks=np.arange(4),xticklabels=[LABELS[KEYS.index(k)] for k in keys],ylim=(0,9),ylabel='All E within ±10% / 8',
               title=f'({chr(97+panel)}) '+('Isolated controls' if panel==0 else 'Known-K/t pairs'))
        ax.tick_params(axis='x',rotation=18);ax.grid(axis='y',alpha=.2)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='outside lower center',ncol=5,fontsize=9)
    fig.suptitle('Figure 0007-64-1: Fitting space, weights and independent mean calibration')
    save(fig,64,'success')
    fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for panel,(ax,k) in enumerate(zip(axes.flat,KEYS[4:])):
        values=[[r['relative_error']*100 for r in components if r['condition']==k and r['method']==m and r['weak']] for m in METHODS]
        boxes=ax.boxplot(values,patch_artist=True,tick_labels=['Old GLS','RL+recovery','Native OLS','Native GLS','+mean cal.'])
        for patch,color in zip(boxes['boxes'],COLORS):patch.set_facecolor(color);patch.set_alpha(.5)
        for j,v in enumerate(values):ax.scatter(np.full(len(v),j+1),v,s=10,color=COLORS[j])
        ax.axhspan(-10,10,color='green',alpha=.1);ax.axhline(0,color='grey',lw=.7)
        ax.set(title=f'({chr(97+panel)}) '+LABELS[KEYS.index(k)],ylabel='Weak E error [%]');ax.tick_params(axis='x',rotation=15)
    fig.suptitle('Figure 0007-65-1: Weak-component errors; all cases and independent y scales')
    save(fig,65,'weak_errors')
    s.initialize()
    with np.load(s.ART/'inputs/observations.npz') as f:voltage=f['voltage_mV'];ids=f['ids']
    time_ns=np.arange(1024)*2.+.25
    fig,axes=plt.subplots(2,2,figsize=(14,7),constrained_layout=True)
    for col,identifier in enumerate((40,120)):
        idx=list(ids).index(identifier);label=labels[idx];old=oldfits[identifier]
        old_curve=s.LEGACY.components(old['times_ns'],old['energy_MeV'])[0].sum(0)
        axes[0,col].plot(time_ns,voltage[idx]-s.BACKGROUND,color='k',lw=.8,label='Observed - background')
        axes[0,col].plot(time_ns,old_curve,color='#555555',lw=.9,label='Existing voltage GLS')
        for method,color in [('native_gls_centered',COLORS[3]),('native_gls_calibrated',COLORS[4])]:
            fit=results[identifier]['fits'][method];curve=s.MODEL.signal(fit['times_ns'],fit['energy_MeV'],corrected=method.endswith('calibrated'))
            axes[0,col].plot(time_ns,curve,color=color,lw=.9,label=NAMES[METHODS.index(method)])
            axes[1,col].plot(time_ns,voltage[idx]-s.BACKGROUND-curve,color=color,lw=.8,label=NAMES[METHODS.index(method)])
        axes[0,col].set(title=f'({chr(97+col)}) id{identifier}: {key(label)}',ylabel='Signal [mV]',xlim=(0,2048));axes[0,col].legend(fontsize=8)
        axes[1,col].set(title=f'({chr(99+col)}) Full-window residual',ylabel='Residual [mV]',xlabel='Time from window start [ns]',xlim=(0,2048));axes[1,col].legend(fontsize=8)
    fig.suptitle('Figure 0007-66-1: Full native voltage and forward-model residuals')
    save(fig,66,'full_waveforms')
    fig,axes=plt.subplots(2,2,figsize=(12,7),constrained_layout=True)
    correction=s.MODEL.correction
    for j,ax in enumerate(axes.flat):
        ax.plot(correction['relative_ns'],1000*correction['profiles'][j],color=COLORS[4])
        ax.axhline(0,color='grey',lw=.7)
        ax.set(title=f'({chr(97+j)}) '+LABELS[j]+f", mean E={correction['energy_knots'][j]:.3f} MeV",
               xlabel='Time relative to event [ns]',ylabel='Voltage correction per 1000 q0 [mV]',xlim=(-64,1720))
    fig.suptitle('Figure 0007-67-1: Mean residual calibration from isolated development data only')
    save(fig,67,'mean_correction')
    paths=list(table.glob('*.csv'))+[s.ART/'report.json']+list(s.ART.glob('0007-*.png'))+list(s.ART.glob('0007-*.svg'))
    s.write(s.EXP/'report-freeze.json',s.hashes(paths))
    print(json.dumps({k:v for k,v in report.items() if k!='summary'},indent=2),flush=True)
    for k in KEYS:print(k,[(m,lookup[k,m]['success10'],round(lookup[k,m]['bias_pct'],2),round(lookup[k,m]['rmse_pct'],2),lookup[k,m]['boundary']) for m in METHODS],flush=True)


if __name__=='__main__':main()
