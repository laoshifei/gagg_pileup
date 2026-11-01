"""Post-fit truth-point causality diagnostic; never imported by inference."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from . import alpha_pair_recovery_forward_study as s
from .alpha_pair_energy_calibrated import learned_metric


def main():
    s.initialize();labels=s.read(s.ART/'labels.json');rows=[]
    with np.load(s.ART/'inputs/observations.npz') as f:voltage=f['voltage_mV']
    methods=['native_gls_centered','native_gls_calibrated','existing_gls_mean']
    for i,r in enumerate(labels):
        if r['kind']!='single':continue
        fit=s.read(s.ART/'cases'/f"{r['id']:04d}.json");t=np.asarray(r['times_ns']);e=np.asarray(r['energy_MeV'])
        metric=learned_metric(s.LEGACY,t,fit['seed']['energy_MeV'],s.NOISE)
        before=s.LEGACY.sample_grid<t[0];before_fir=s.LEGACY.sample_grid[3:]<t[0]
        for method in methods:
            model=s.LEGACY.components(t,e)[0].sum(0) if method=='existing_gls_mean' else s.MODEL.signal(t,e,corrected=method=='native_gls_calibrated')
            residual=metric.apply(s.LEGACY.transform(model-(voltage[i]-s.BACKGROUND)))
            model_only=metric.apply(s.LEGACY.transform(model))
            rows.append(dict(id=r['id'],condition=r['condition'],method=method,
                pre_event_max_model_mV=float(np.max(abs(model[before]))),
                pre_event_model_white_norm2=float(np.sum(model_only[before_fir]**2)),
                pre_event_residual_chi2_fraction=float(np.sum(residual[before_fir]**2)/np.sum(residual**2))))
    summary=[]
    for condition in ['weak','low','mid','high']:
        for method in methods:
            values=[x for x in rows if x['condition']==condition and x['method']==method]
            summary.append(dict(condition=condition,method=method,frames=len(values),
                pre_event_max_model_mV_mean=float(np.mean([x['pre_event_max_model_mV'] for x in values])),
                pre_event_residual_chi2_fraction_mean=float(np.mean([x['pre_event_residual_chi2_fraction'] for x in values]))))
    s.write(s.REC/'pre_event_diagnostic.json',dict(scope='Truth E/t used ONLY for post-fit physical-admissibility diagnostic; not inference or calibration',rows=rows,summary=summary))
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
    # Post-fit diagnostic illustration: first isolated high-energy record.
    idx=next(i for i,r in enumerate(labels) if r['condition']=='high');r=labels[idx]
    relative=s.LEGACY.sample_grid-r['times_ns'][0];mask=(relative>=-24)&(relative<0)
    axes[0].plot(relative[mask],(voltage[idx]-s.BACKGROUND)[mask],color='grey',marker='.',label='Observed - mean background')
    names=['Uncalibrated native mean','Calibrated native mean','Existing GLS mean'];colors=['#d7a242','#d66543','#20252b']
    for method,name,color in zip(methods,names,colors):
        model=s.LEGACY.components(r['times_ns'],r['energy_MeV'])[0].sum(0) if method=='existing_gls_mean' else s.MODEL.signal(r['times_ns'],r['energy_MeV'],corrected=method=='native_gls_calibrated')
        axes[0].plot(relative[mask],model[mask],marker='.',label=name,color=color)
    axes[0].axhline(0,color='black',lw=.6);axes[0].set(title=f"(a) Isolated high id{r['id']}, at TRUE E/t",xlabel='Time BEFORE event [ns]',ylabel='Mean event signal [mV]',xlim=(-24,0))
    axes[0].legend(fontsize=8)
    for j,method in enumerate(methods):
        values=[next(x['pre_event_residual_chi2_fraction_mean'] for x in summary if x['condition']==c and x['method']==method)*100 for c in ('weak','low','mid','high')]
        axes[1].bar(np.arange(4)+(j-1)*.24,values,width=.23,color=colors[j],label=names[j])
    axes[1].set(title='(b) Before-event share of whitened residual squared',xticks=np.arange(4),xticklabels=['Weak','Low','Mid','High'],ylabel='Mean per-frame share [%]',ylim=(0,100));axes[1].legend(fontsize=8)
    fig.suptitle('Figure 0007-67-2: Physical support failure, not a pile-up resolution limit')
    for extension in ('png','svg'):fig.savefig(s.ART/f'0007-67-2_pre_event.{extension}',dpi=160,bbox_inches='tight')
    plt.close(fig)
    s.write(s.EXP/'diagnostic-freeze.json',s.hashes([s.REC/'pre_event_diagnostic.json',s.ART/'0007-67-2_pre_event.png',s.ART/'0007-67-2_pre_event.svg']))
    print('CAUSALITY_DIAGNOSTIC_RECORDED',flush=True)


if __name__=='__main__':main()
