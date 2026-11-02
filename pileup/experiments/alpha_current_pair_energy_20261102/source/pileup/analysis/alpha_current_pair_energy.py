"""Frozen, cost-limited energy comparison on previously checked waveforms."""
import hashlib
import numpy as np

CONDITIONS=(('high_weak',8.),('high_weak',40.),('high_weak',200.),('high_weak',1200.),('weak_high',40.))


def selected_blocks():
    return sorted(range(32),key=lambda b:hashlib.sha256(f'current-pair-energy-20261102:{b}'.encode()).hexdigest())[:8]


def summarize(errors,success,boundary,converged):
    x=np.asarray(errors,float)
    if x.ndim!=2 or x.shape[1]!=2 or not np.isfinite(x).all():raise ValueError('Expected finite weak/strong errors')
    return dict(rows=len(x),weak_bias_percent=float(x[:,0].mean()),strong_bias_percent=float(x[:,1].mean()),
        weak_rmse_percent=float(np.sqrt(np.mean(x[:,0]**2))),strong_rmse_percent=float(np.sqrt(np.mean(x[:,1]**2))),
        weak_error_sd_percent=float(x[:,0].std(ddof=1)),strong_error_sd_percent=float(x[:,1].std(ddof=1)),
        strict_pair_count=int(np.sum(success)),strict_pair_percent=float(100*np.mean(success)),
        boundary=int(np.sum(boundary)),nonconverged=int(np.sum(~np.asarray(converged,bool))))


def paired_statistics(old,new,old_success,new_success,index):
    """Each resampled row is an optical mother pair, carrying every condition."""
    a,b=np.asarray(old,float),np.asarray(new,float)
    u,v=np.asarray(old_success,float),np.asarray(new_success,float)
    idx=np.asarray(index)
    if a.shape!=b.shape or a.ndim!=2 or u.shape!=a.shape or v.shape!=a.shape:
        raise ValueError('Mother-by-condition arrays must agree')
    if idx.ndim!=2 or idx.shape[1]!=len(a):raise ValueError('Resample whole mother pairs')
    denominator=np.mean(a*a)
    if denominator<=0:raise ValueError('Positive reference MSE required')
    ratios=np.sqrt(np.mean(b[idx]**2,axis=(1,2))/np.mean(a[idx]**2,axis=(1,2)))
    gains=100*np.mean(v[idx]-u[idx],axis=(1,2))
    return dict(weak_rmse_ratio=float(np.sqrt(np.mean(b*b)/denominator)),
        weak_rmse_ratio_ci95=np.percentile(ratios,[2.5,97.5]).tolist(),
        strict_gain_pp=float(100*np.mean(v-u)),strict_gain_ci95=np.percentile(gains,[2.5,97.5]).tolist())
