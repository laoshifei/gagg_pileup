"""Plot the exact plus-one, tie-aware empirical rank used by the gate."""
from __future__ import annotations

import shutil

import numpy as np

from . import alpha_pair_null_study as study


def rank_curve(values, right):
    values=np.sort(np.asarray(values,float))
    assert len(values) and right>values[-1]
    x=np.r_[np.unique(values),right]
    rank=(1.+len(values)-np.searchsorted(values,x,side='left'))/(len(values)+1.)
    return x,rank


def draw():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    x,p=rank_curve([0.,0.,.1],.11)
    np.testing.assert_array_equal(x,[0.,.1,.11])
    np.testing.assert_allclose(p,[1.,.5,.25],atol=0.,rtol=0.)
    calibration=study.load(study.EXPERIMENT/'model-null-calibration.json')
    folder=study.OUTPUT/'figures';folder.mkdir(exist_ok=True)
    path=folder/'0004-49-1.png'
    if path.exists():
        history=study.OUTPUT/'figure-history';history.mkdir(exist_ok=True)
        original=history/'0004-49-1.before-exact-rank.png'
        if not original.exists():shutil.copyfile(path,original)
    fig,axes=plt.subplots(1,2,figsize=(13,4.6))
    for ax,(mode,cal),letter in zip(axes,calibration.items(),'ab'):
        right=max(max(v) for v in cal['scores'])*1.06+1e-12
        for label,values in zip(('Weak-Q bin','Low-Q bin','Mid-Q bin','High-Q bin'),cal['scores']):
            x,p=rank_curve(values,right)
            line,=ax.step(x,p,where='pre',label=f'{label}, n={len(values)}')
            if len(np.unique(values))==1:ax.plot(x[0],p[0],'.',color=line.get_color())
            rank=int(np.floor(cal['alpha']*(len(values)+1)))
            if rank:ax.axvline(sorted(values)[-rank],color=line.get_color(),lw=.8,ls=':',alpha=.6)
        ax.axhline(.025,color='black',ls='--',lw=1.)
        ax.set(title=f'({letter}) {mode}',xlabel='Post-search evidence score',
               ylabel='Exact plus-one upper-tail rank',yscale='log',ylim=(.003,1.1))
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Figure 0004-49-1 | Isolated-event calibration; dashed line: rank 0.025')
    fig.tight_layout();fig.savefig(path,dpi=150);plt.close(fig)
    study.save(study.RECORD/'calibration-figure-review.json',dict(passed=True,
        figure_sha256=study.digest(path),plus_one=True,ties_use_greater_or_equal=True,
        original_plot_preserved=True,inference_and_selection_unchanged=True,
        explanation='The initial display omitted plus one and did not collapse ties; the inference rank was already correct'))
    print('EXACT_RANK_FIGURE_WRITTEN',flush=True)


if __name__=='__main__':draw()
