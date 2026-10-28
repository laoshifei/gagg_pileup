"""Separate recovery stages and compare errors with local Gaussian references."""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict

import numpy as np

from .alpha_pair_information import invert_information, standardized_weak_error
from .alpha_pair_information_study import PROJECT, EXPERIMENT, OUTPUT, RECORD, INPUT, load, save, rows, digest, check_freeze

CONDITIONS = ('mid_mid', 'high_weak', 'weak_high')
GAPS = (8., 40., 200.)


def quantiles(values):
    return np.percentile(values, [16, 50, 84]).tolist() if len(values) else None


def csv_table(path, records):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def local_summary(row):
    f = row['information']
    weak = f['weak_index']
    strong = 1-weak
    F = np.asarray(f['total'])
    inverse = np.asarray(f['inverse'])
    corr = np.asarray(f['correlation'])
    # Exploratory decomposition of the predeclared known-strong-Q/t comparison.
    q_known_indices = [0, 1, weak+2]
    t_known_indices = [weak, 2, 3]
    q_known = invert_information(F[np.ix_(q_known_indices, q_known_indices)])
    t_known = invert_information(F[np.ix_(t_known_indices, t_known_indices)])
    return dict(task_id=row['task_id'], role=row['role'], native_index=row['native_index'],
        condition=row['condition'], gap_ns=row['gap_ns'], weak_charge_pe=row['charge_pe'][weak],
        strong_charge_pe=row['charge_pe'][strong], phase_ns=row.get('phase_ns'),
        weak_sigma_time_ns=f['sigma_time_ns'][weak], weak_sigma_log_charge=f['sigma_log_charge'][weak],
        fixed_times_sigma_log_charge=float(np.sqrt(f['fixed_times_charge'][weak][weak])),
        known_strong_sigma_time_ns=float(10*np.sqrt(f['known_strong'][0][0])),
        known_strong_sigma_log_charge=float(np.sqrt(f['known_strong'][1][1])),
        known_strong_Q_sigma_log_charge=float(np.sqrt(q_known[-1, -1])),
        known_strong_t_sigma_log_charge=float(np.sqrt(t_known[weak+1, weak+1])),
        correlation_weakQ_strongQ=float(corr[weak+2, strong+2]),
        correlation_weakQ_weakt=float(corr[weak+2, weak]),
        condition_number=f['condition_number'], derivative_step=f['step'],
        derivative_converged=f['derivative_converged'],
        mean_trace_fraction=float(np.trace(f['mean'])/np.trace(F)),
        charge_nuisance_variance_ratio=float(inverse[weak+2, weak+2]/f['known_strong'][1][1]))


def analyze():
    check_freeze()
    for name in ('numerical-tests', 'pilot-validation', 'information-run'):
        assert load(RECORD/(name+'.json'))['passed']
    information = sorted(rows(OUTPUT/'information.truth.jsonl'), key=lambda r:r['task_id'])
    tasks = load(EXPERIMENT/'tasks.truth.json')
    assert [r['task_id'] for r in information] == list(range(444))
    for task, r in zip(tasks, information):
        assert all(r[k] == v for k, v in task.items())
    assert digest(OUTPUT/'information.truth.jsonl') == load(RECORD/'information-run.json')['output_sha256']
    available = {r['native_index']:r for r in information if r['role']=='native_parameters' and r['status']=='ok'}
    unresolved = [dict(task_id=r['task_id'], error=r.get('error')) for r in information if r['status']!='ok']
    unconverged = [r['task_id'] for r in information if r['status']=='ok' and not r['information']['derivative_converged']]
    values = [local_summary(r) for r in information if r['status']=='ok']
    info_summary = []
    fields = [k for k in values[0] if k.startswith(('weak_sigma', 'fixed_times_', 'known_strong', 'correlation_', 'charge_nuisance'))]
    for condition in CONDITIONS:
        for gap in GAPS:
            subset = [r for r in values if r['role']=='native_parameters' and (r['condition'], r['gap_ns'])==(condition, gap)]
            info_summary.append(dict(condition=condition, gap_ns=gap, frames=len(subset),
                **{k:quantiles([r[k] for r in subset]) for k in fields}))
    save(OUTPUT/'local-information-summary.json', dict(rows=values, summary=info_summary,
        unresolved=unresolved, derivative_unconverged=unconverged,
        quantile_order=[16,50,84], warning='Known K and true parameters; local information scales, not calibrated finite-fit uncertainties',
        exploratory='Splitting known strong Q/t into Q-only and t-only conditioning is an exploratory derived-matrix comparison'))
    csv_table(OUTPUT/'local-information.csv', values)

    source_scores = rows(INPUT/'scores.jsonl')
    primary = [r for r in source_scores if r['method'] in ('native/native_gate','matched/native_gate')]
    assert len(primary) == 2280
    stage_rows = []
    for r in primary:
        if r['kind'] != 'pair':
            continue
        q = np.asarray(r['true_q']); t = np.asarray(r['true_t'])
        assert np.all(q>=500.) and len(q)==2 and np.all(np.diff(t)>0.)
        fit = r['fit']; K2 = len(fit['times_ns'])==2
        time10 = timestrict = success = False
        if K2:
            order = np.argsort(fit['times_ns'])
            dt = np.asarray(fit['times_ns'])[order]-t
            dq = np.asarray(fit['charge_pe'])[order]/q-1.
            time10 = bool(np.all(abs(dt)<=10.))
            timestrict = bool(np.all(abs(dt)<=min(10.,r['gap_ns']/4.)))
            success = timestrict and bool(np.all(abs(dq)<=.20)) and r['extra']==0
        assert success == r['success'], (r['domain'],r['index'])
        assert int(time10)*2 == (2 if r['matched']==2 else 0)
        stage_rows.append(dict(domain=r['domain'], index=r['index'], native_index=r['native_index'],
            condition=r['condition'], gap_ns=r['gap_ns'], block=r['block'], draw=r['draw'],
            K2_selected=K2, both_times10=time10, both_times_strict=timestrict,
            strict_success=success, missed=r['missed'], extra=r['extra'], split=r['split'], merge=r['merge']))
    assert len(stage_rows) == 1080
    stages = []
    for domain in ('native','matched'):
        for condition in CONDITIONS:
            for gap in GAPS:
                rr = [r for r in stage_rows if (r['domain'],r['condition'],r['gap_ns'])==(domain,condition,gap)]
                stage = dict(domain=domain,condition=condition,gap_ns=gap,frames=len(rr))
                for key in ('K2_selected','both_times10','both_times_strict','strict_success','missed','extra','split','merge'):
                    stage[key]=sum(r[key] for r in rr)
                assert stage['strict_success']<=stage['both_times_strict']<=stage['both_times10']<=stage['K2_selected']
                stages.append(stage)
    save(OUTPUT/'stage-summary.json', dict(rows=stage_rows,summary=stages,
        warning='K=2 output alone is not correct second-event detection; all denominators include every pair, no selection conditioning'))
    csv_table(OUTPUT/'stage-summary.csv',stages)

    diagnostic = load(INPUT/'diagnostic-summary.json')
    errors = []; excluded = []
    for r in diagnostic['rows']:
        if r['method'] == 'blind_native_gate':
            continue
        if r['native_index'] not in available:
            excluded.append(dict(domain=r['domain'],index=r['index'],method=r['method'],reason='unresolved information'))
            continue
        point = available[r['native_index']]
        normalization = standardized_weak_error(r['fit'],point['times_ns'],point['charge_pe'],point['information'],
                                               fixed_times=r['method']=='fixed_times')
        error = dict(domain=r['domain'],index=r['index'],native_index=r['native_index'],method=r['method'],
            condition=r['condition'],gap_ns=r['gap_ns'],block=r['block'],strict_success=r['success'],
            converged=r['fit'].get('converged',True),**normalization)
        if r['method'] == 'fixed_times':
            error['inside_reference95'] = abs(error['standardized_log_charge_error'])<=1.959963984540054
            assert abs(error['time_error_ns']) < 1e-9
        else:
            error['inside_reference95'] = error['mahalanobis_squared']<=-2*np.log(.05)
        errors.append(error)
    assert len(errors)+len(excluded)==2160
    grouped = defaultdict(list)
    for r in errors:
        grouped[r['domain'],r['condition'],r['gap_ns'],r['method']].append(r)
    error_summary = []
    for (domain,condition,gap,method), rr in sorted(grouped.items()):
        error_summary.append(dict(domain=domain,condition=condition,gap_ns=gap,method=method,frames=len(rr),
            strict_success=sum(r['strict_success'] for r in rr), inside_reference95=sum(r['inside_reference95'] for r in rr),
            not_converged=sum(not r['converged'] for r in rr),
            weak_charge_fractional_error=quantiles([r['charge_fractional_error'] for r in rr]),
            weak_time_error_ns=quantiles([r['time_error_ns'] for r in rr]),
            log_charge_rmse=float(np.sqrt(np.mean([r['log_charge_error']**2 for r in rr]))),
            time_rmse_ns=float(np.sqrt(np.mean([r['time_error_ns']**2 for r in rr]))),
            normalized_squared_error=quantiles([r['standardized_log_charge_error']**2 if method=='fixed_times' else r['mahalanobis_squared'] for r in rr])))
    save(OUTPUT/'error-summary.json',dict(rows=errors,summary=error_summary,excluded=excluded,
        warning='Truth-conditioned local Gaussian reference regions, not operational confidence intervals or independent empirical coverage estimates'))
    csv_table(OUTPUT/'error-summary.csv',error_summary)
    phase = [r for r in values if r['role']=='phase_check']
    same_phase = next(r for r in phase if r['phase_ns']==.37)
    nominal = next(r for r in values if r['role']=='nominal_grid' and r['condition']=='high_weak' and r['gap_ns']==40. and r['weak_charge_pe']==700.)
    assert all(abs(same_phase[k]-nominal[k])<1e-12 for k in fields)
    good = [r['information'] for r in information if r['status']=='ok']
    for f in good:
        F=np.asarray(f['total']); V=np.asarray(f['inverse']); w=f['weak_index']
        np.testing.assert_allclose(F@V,np.eye(4),atol=1e-8)
        assert np.linalg.eigvalsh(np.asarray(f['weak_joint'])-f['known_strong']).min()>-1e-9
        assert f['fixed_times_charge'][w][w]<=V[w+2,w+2]+1e-10
        np.testing.assert_allclose(F,np.asarray(f['mean'])+f['covariance'],rtol=1e-12)
    report = dict(stages=stages,information=info_summary,errors=error_summary,phase_check=phase,
        primary_cached_summary=[r for r in load(INPUT/'summary.json') if r['method'] in ('native/native_gate','matched/native_gate')],
        tasks=len(information),native_points=360,nominal_points=80,phase_points=4,
        unresolved=unresolved,derivative_unconverged=unconverged,new_blind_fits=0,new_roots=0,
        warning='Diagnostic only: true parameters and K, local regular Gaussian information; no universal separation or unbiasedness claim')
    save(OUTPUT/'report.json',report)
    save(RECORD/'science-validation.json',dict(passed=True,tasks=444,stage_scores=1080,cached_primary_scores=2280,
        normalized_diagnostic_scores=len(errors),excluded=excluded,unresolved=unresolved,derivative_unconverged=unconverged,
        strict_scores_match_cache=True,inputs_unchanged=True,freeze_unchanged=True,
        maximum_whitening_error=max(f['whitening_relative_error'] for f in good),
        maximum_derivative_matrix_change=max(f['derivative_checks'][-1]['relative_information_change'] for f in good),
        maximum_derivative_sigma_change=max(f['derivative_checks'][-1]['maximum_sigma_fractional_change'] for f in good),
        maximum_condition_number=max(f['condition_number'] for f in good),
        phase_duplicate_matches=True,known_parameters_reduce_variance=True,
        output_hashes={name:digest(OUTPUT/name) for name in ('report.json','local-information-summary.json','stage-summary.json','error-summary.json','information.truth.jsonl')},
        source_sha256=digest(__file__)))
    print('INFORMATION_ANALYSIS_COMPLETE',len(errors),'errors',len(unresolved),'unresolved',len(unconverged),'unstable',flush=True)


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.special import erf
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.2,'figure.dpi':120})
    assert load(RECORD/'science-validation.json')['passed']
    report = load(OUTPUT/'report.json')
    local = load(OUTPUT/'local-information-summary.json')
    errors = load(OUTPUT/'error-summary.json')['rows']
    destination=OUTPUT/'figures';destination.mkdir(exist_ok=True)
    hashes={}

    def finish(fig,number,title):
        fig.suptitle('Figure '+number+' | '+title,fontsize=13)
        fig.tight_layout(rect=(0,0,1,.95))
        path=destination/(number+'.png');fig.savefig(path,dpi=160);plt.close(fig)
        hashes[number]=digest(path)

    colors=['#888888','#0077bb','#ee7733','#009988']
    fig,axes=plt.subplots(2,3,figsize=(13.5,7.8),sharey=True)
    stage_keys=['K2_selected','both_times10','both_times_strict','strict_success']
    names=['Select K=2','Both t within 10 ns','Both t within strict limit','Both t + Q correct']
    for row,domain in enumerate(('native','matched')):
        for col,condition in enumerate(CONDITIONS):
            ax=axes[row,col]
            rr=[next(r for r in report['stages'] if (r['domain'],r['condition'],r['gap_ns'])==(domain,condition,g)) for g in GAPS]
            for j,(key,label,color) in enumerate(zip(stage_keys,names,colors)):
                ax.bar(np.arange(3)+(j-1.5)*.2,[r[key]/r['frames'] for r in rr],width=.19,label=label,color=color)
            ax.set(xticks=np.arange(3),xticklabels=['8','40','200'],ylim=(0,1.05),xlabel='True separation (ns)',
                   title=f'({chr(97+row*3+col)}) {domain}: {condition}')
            if col==0:ax.set_ylabel('Fraction of all pairs')
    axes[0,0].legend(fontsize=8,loc='lower right')
    finish(fig,'0004-77-1','Two outputs, correct times, and correct charges are distinct stages')

    fig,axes=plt.subplots(2,3,figsize=(13.5,7.8))
    for row,condition in enumerate(('high_weak','weak_high')):
        rr=[next(r for r in report['information'] if (r['condition'],r['gap_ns'])==(condition,g)) for g in GAPS]
        specifications=[([('weak_sigma_log_charge','All four unknown','#0077bb'),('fixed_times_sigma_log_charge','Both times fixed','#ee7733'),('known_strong_sigma_log_charge','Strong Q,t fixed','#009988')],'Local scale of ln Q(weak)'),
            ([('weak_sigma_time_ns','All four unknown','#0077bb'),('known_strong_sigma_time_ns','Strong Q,t fixed','#009988')],'Local time scale (ns)'),
            ([('correlation_weakQ_strongQ','Weak Q vs strong Q','#aa3377'),('correlation_weakQ_weakt','Weak Q vs weak t','#228833')],'Local parameter correlation')]
        for col,(series,ylabel) in enumerate(specifications):
            ax=axes[row,col]
            for key,label,color in series:
                a=np.array([r[key] for r in rr])
                ax.plot(GAPS,a[:,1],'-o',label=label,color=color)
                ax.fill_between(GAPS,a[:,0],a[:,2],color=color,alpha=.12)
            ax.set(xscale='log',xticks=GAPS,xticklabels=['8','40','200'],xlabel='True separation (ns)',ylabel=ylabel,
                   title=f'({chr(97+row*3+col)}) {condition}')
            if col<2:ax.set_yscale('log')
            else:ax.set_ylim(-1.05,1.05)
            ax.legend(fontsize=8)
    finish(fig,'0004-78-1','Known K and true parameters; median and 16-84% spread of local scales')

    fig,axes=plt.subplots(1,2,figsize=(12.5,4.8),sharey=True)
    for col,condition in enumerate(('high_weak','weak_high')):
        ax=axes[col]
        for q,color in zip((500.,700.,1500.,3000.),('#aa3377','#ee7733','#0077bb','#009988')):
            rr=sorted([r for r in local['rows'] if r['role']=='nominal_grid' and r['condition']==condition and r['weak_charge_pe']==q],key=lambda r:r['gap_ns'])
            ax.plot([r['gap_ns'] for r in rr],[r['weak_sigma_log_charge'] for r in rr],'-o',color=color,label=f'Weak Q = {q:.0f} pe')
        ax.axhline(np.log(1.2),color='gray',ls=':',label='ln(1.2): tolerance scale only')
        ax.set(xscale='log',yscale='log',xlabel='True separation (ns)',ylabel='Local scale of ln Q(weak)',
               title=f'({chr(97+col)}) {condition}; strong Q = 18000 pe')
        ax.legend(fontsize=8)
    finish(fig,'0004-79-1','Model design grid, fixed phase; not measured recovery or a separation boundary')

    fig,axes=plt.subplots(2,3,figsize=(13.5,7.8),sharey=True)
    for row,condition in enumerate(('high_weak','weak_high')):
        for col,gap in enumerate(GAPS):
            ax=axes[row,col]
            for domain,color in [('native','#0077bb'),('matched','#ee7733')]:
                for method,style in [('known_K','--'),('true_start','-')]:
                    rr=[r for r in errors if (r['domain'],r['condition'],r['gap_ns'],r['method'])==(domain,condition,gap,method)]
                    x=np.sort([r['mahalanobis_squared'] for r in rr])
                    ax.step(x,np.arange(1,len(x)+1)/len(x),where='post',color=color,ls=style,label=domain+' / '+method)
            lo,hi=ax.get_xlim();x=np.geomspace(.001,max(20.,hi),400)
            ax.plot(x,1-np.exp(-x/2),'k:',label='Local Gaussian chi-square(2)')
            ax.axvline(-2*np.log(.05),color='gray',lw=.7)
            ax.set(xscale='log',ylim=(0,1.02),xlabel='Squared normalized weak (t, ln Q) error',ylabel='Empirical cumulative fraction',
                   title=f'({chr(97+row*3+col)}) {condition}; {gap:.0f} ns')
    axes[0,0].legend(fontsize=7,loc='lower right')
    finish(fig,'0004-80-1','Truth-conditioned error reference; finite fitted errors need not follow chi-square')

    fig,axes=plt.subplots(2,3,figsize=(13.5,7.8),sharey=True)
    for row,condition in enumerate(('high_weak','weak_high')):
        for col,gap in enumerate(GAPS):
            ax=axes[row,col]
            for domain,color in [('native','#0077bb'),('matched','#ee7733')]:
                rr=[r for r in errors if (r['domain'],r['condition'],r['gap_ns'],r['method'])==(domain,condition,gap,'fixed_times')]
                x=np.sort([abs(r['standardized_log_charge_error']) for r in rr])
                ax.step(x,np.arange(1,len(x)+1)/len(x),where='post',color=color,label=domain)
            lo,hi=ax.get_xlim();x=np.geomspace(.001,max(5.,hi),400)
            ax.plot(x,erf(x/np.sqrt(2.)),'k:',label='Local Gaussian |N(0,1)|')
            ax.axvline(1.959963984540054,color='gray',lw=.7)
            ax.set(xscale='log',ylim=(0,1.02),xlabel='Absolute normalized weak ln Q error',ylabel='Empirical cumulative fraction',
                   title=f'({chr(97+row*3+col)}) {condition}; {gap:.0f} ns')
    axes[0,0].legend(fontsize=8,loc='lower right')
    finish(fig,'0004-80-2','Both true times fixed; local Gaussian reference is not a calibrated confidence interval')
    assert list(hashes)==load(EXPERIMENT/'protocol.json')['figures']
    save(RECORD/'figure-delivery.json',dict(passed=True,figures=hashes,source_sha256=digest(__file__)))
    print('INFORMATION_FIGURES_COMPLETE',len(hashes),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('analyze','figures'))
    args=parser.parse_args()
    if args.command=='analyze':analyze()
    else:figures()


if __name__=='__main__':
    main()
