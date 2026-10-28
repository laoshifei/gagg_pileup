"""Mean-only estimating equations versus the complete joint Gaussian score.

IRGLS freezes C within each solve and updates it between solves. This is not
minimization of r.T C(theta)^-1 r with logdet simply removed. No truth inputs.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from .alpha_pair_model_sampler import GaussianFactor

METHODS=('baseline','wide_nll','frozen_gls','iterated_gls')
DEFAULTS=dict(starts_per_metric=3,maxiter=150,max_outer=12,damping=.5,
              charge_log_tolerance=.002,time_tolerance_ns=.02)


def encode(times,charges,bounds):
    t=np.asarray(times,float);q=np.asarray(charges,float)
    if len(t)!=len(q) or len(t) not in (1,2):raise ValueError('One or two components required')
    order=np.argsort(t);t=np.clip(t[order],*bounds);q=np.clip(q[order],50.,50000.)
    if len(t)==1:return np.r_[t/10.,np.log(q)]
    gap=(t[1]-t[0])/max(bounds[1]-t[0],1e-12)
    return np.r_[t[0]/10.,100.*gap,np.log(q)]


def decode(z,count,bounds):
    z=np.asarray(z,float)
    if count==1:return np.array([10*z[0]]),np.exp(z[1:])
    if count!=2:raise ValueError('One or two components required')
    t1=10*z[0];return np.array([t1,t1+(bounds[1]-t1)*z[1]/100.]),np.exp(z[2:])


def parameter_bounds(count,bounds):
    if count==1:return [(bounds[0]/10.,bounds[1]/10.),(np.log(50.),np.log(50000.))]
    return [(bounds[0]/10.,bounds[1]/10.),(0.,100.)]+[(np.log(50.),np.log(50000.))]*2


def local_bounds(response,one,base):
    anchor=one['times_ns'][0]
    return [max(0.,min(anchor-96.,min(base['times_ns']))),
            min(float(response.sample_grid[-1]),max(anchor+512.,max(base['times_ns'])))]


def optimize(objective,initial,bounds,maxiter):
    lower,upper=np.asarray(bounds).T;z0=np.clip(initial,lower,upper);initial_cost=float(objective(z0))
    fit=minimize(objective,z0,method='L-BFGS-B',bounds=bounds,
        options=dict(maxiter=maxiter,ftol=1e-10,gtol=1e-5,maxls=30,eps=1e-5))
    use=bool(np.isfinite(fit.fun) and fit.fun<=initial_cost)
    z=fit.x if use else z0;value=float(objective(z))
    if not np.isfinite(value):raise ValueError('Non-finite retained objective')
    return dict(z=z,value=value,initial_cost=initial_cost,converged=bool(fit.success),
        retained_initial=not use,evaluations=int(fit.nfev),iterations=int(fit.nit),message=str(fit.message))


class Objectives:
    def __init__(self,response,noise,voltage,count,bounds):
        self.response=response;self.noise=noise;self.count=count;self.bounds=bounds
        self.observed=np.asarray(voltage,float)-response.baseline-noise.background_mV
        self.zy=noise.transform(self.observed)

    def curves(self,z):
        t,q=decode(z,self.count,self.bounds)
        m=self.response.components(t,q)
        return t,q,m,self.noise.transform(m)

    def nll(self,z):
        _,q,_,zm=self.curves(z)
        return float(self.noise.score(self.zy-zm.sum(0),zm,q,transformed=True))

    def metric(self,z):
        _,q,m,_=self.curves(z)
        return GaussianFactor.from_noise(self.noise,m,q)

    def fixed_objective(self,factor):
        def objective(z):
            _,_,_,zm=self.curves(z)
            w=factor.whiten(self.zy-zm.sum(0))
            return float(w@w)
        return objective

    def result(self,z,**metadata):
        t,q,m,_=self.curves(z)
        return dict(times_ns=t.tolist(),charge_pe=q.tolist(),
            **self.noise.score(self.observed-m.sum(0),m,q,details=True),**metadata)


def observed_bank(one,base,count,bounds):
    seeds=[encode(base['times_ns'],base['charge_pe'],bounds)]
    anchor=float(one['times_ns'][0]);total=float(sum(one['charge_pe']))
    if count==1:
        for dt in (-8.,0.,8.):
            for scale in (.9,1.,1.1):seeds.append(encode([anchor+dt],[total*scale],bounds))
    else:
        for direction in (-1,1):
            for gap in (4.,8.,16.,32.,64.,128.,256.,512.):
                t=np.sort([anchor,anchor+direction*gap])
                if t[0]<bounds[0] or t[1]>bounds[1]:continue
                for share in (.01,.03,.1,.3,.5,.7,.9,.97,.99):
                    seeds.append(encode(t,[total*share,total*(1-share)],bounds))
    return list({tuple(z):z for z in seeds}.values())


def diverse_indices(bank,values,count,bounds,number):
    scaled=[]
    for z in bank:
        t,q=decode(z,count,bounds);scaled.append(np.r_[t/16.,np.log(q)/.35])
    chosen=[]
    for i in np.argsort(values,kind='stable'):
        if not chosen or min(np.linalg.norm(scaled[i]-scaled[j]) for j in chosen)>=1.:
            chosen.append(int(i))
        if len(chosen)==number:break
    return chosen


def damped_step(objective,current,proposal,damping):
    old=float(objective(current))
    for fraction in (damping,damping/2.,damping/4.,damping/8.):
        trial=current+fraction*(proposal-current);value=float(objective(trial))
        if np.isfinite(value) and value<=old+1e-8:
            return trial,dict(accepted=True,fraction=fraction,before=old,after=value)
    return current.copy(),dict(accepted=False,fraction=0.,before=old,after=old)


def fit_estimators(response,noise,voltage,one,base,settings=None):
    cfg=dict(DEFAULTS if settings is None else settings);count=len(base['charge_pe'])
    bounds=local_bounds(response,one,base);limits=parameter_bounds(count,bounds)
    objective=Objectives(response,noise,voltage,count,bounds)
    z0=encode(base['times_ns'],base['charge_pe'],bounds);factor=objective.metric(z0)
    frozen=objective.fixed_objective(factor);bank=observed_bank(one,base,count,bounds)
    nll_values=[objective.nll(z) for z in bank];gls_values=[frozen(z) for z in bank]
    ids=list(dict.fromkeys([0]+diverse_indices(bank,nll_values,count,bounds,cfg['starts_per_metric'])+
                           diverse_indices(bank,gls_values,count,bounds,cfg['starts_per_metric'])))
    seeds=[bank[i] for i in ids]
    trials={name:[optimize(function,z,limits,cfg['maxiter']) for z in seeds]
            for name,function in (('wide_nll',objective.nll),('frozen_gls',frozen))}
    best={name:min(rr,key=lambda r:r['value']) for name,rr in trials.items()}
    fits={'baseline':objective.result(z0,converged=base.get('converged',True))}
    for name in best:
        fits[name]=objective.result(best[name]['z'],converged=best[name]['converged'],
            objective_value=best[name]['value'],objective_kind='full_nll' if name=='wide_nll' else 'frozen_quadratic')
    assert fits['wide_nll']['cost']<=fits['baseline']['cost']+1e-7
    current=best['frozen_gls']['z'];history=[];converged=False
    for iteration in range(cfg['max_outer']):
        fixed=objective.fixed_objective(objective.metric(current))
        starts=list({tuple(z):z for z in (current,z0,best['wide_nll']['z'])}.values())
        rr=[optimize(fixed,z,limits,cfg['maxiter']) for z in starts]
        proposal=min(rr,key=lambda r:r['value']);old_t,old_q=decode(current,count,bounds)
        new_t,new_q=decode(proposal['z'],count,bounds)
        qchange=float(np.max(np.abs(np.log(new_q/old_q))));tchange=float(np.max(np.abs(new_t-old_t)))
        updated,step=damped_step(fixed,current,proposal['z'],cfg['damping'])
        converged=bool(step['accepted'] and proposal['converged'] and qchange<=cfg['charge_log_tolerance'] and tchange<=cfg['time_tolerance_ns'])
        history.append(dict(iteration=iteration+1,charge_log_change=qchange,time_change_ns=tchange,
            fixed_point_converged=converged,step=step,proposal=proposal,trials=rr,updated=updated))
        current=updated
        if converged or not step['accepted']:break
    fits['iterated_gls']=objective.result(current,converged=converged,outer_iterations=len(history),
        objective_kind='mean_score_fixed_point',stop_reason='fixed_point' if converged else 'no_descent' if not history[-1]['step']['accepted'] else 'outer_budget')
    replay=max(abs(objective.nll(encode(f['times_ns'],f['charge_pe'],bounds))-f['cost'])/max(1.,abs(f['cost'])) for f in fits.values())
    assert replay<1e-9
    return dict(fits=fits,bounds_ns=bounds,starts=[z.tolist() for z in seeds],seed_indices=ids,
        seed_bank_size=len(bank),trials=trials,iteration_history=history,cost_replay_error=replay,truth_inputs=[])


def oracle_diagnostics(response,noise,voltage,one,base,primary,true_times,true_charges,settings=None):
    """Explicitly truth-assisted scoring diagnostics; never a primary input."""
    cfg=dict(DEFAULTS if settings is None else settings);count=len(true_charges)
    bounds=primary['bounds_ns'];limits=parameter_bounds(count,bounds)
    objective=Objectives(response,noise,voltage,count,bounds)
    truth_z=encode(true_times,true_charges,bounds)
    if not np.allclose(decode(truth_z,count,bounds)[0],true_times,atol=1e-8):
        return dict(outside_time_domain=True,truth_assisted=True)
    starts=[encode(f['times_ns'],f['charge_pe'],bounds) for f in primary['fits'].values()]
    starts.extend(np.asarray(z) for z in primary['starts'])
    starts=list({tuple(z):z for z in starts}.values())
    nt=count
    def constrained(kind):
        # Known times use only observed charges; known charges use only observed
        # times. Truth starts for all parameters belong to the separate branch.
        ids=list(range(nt,2*nt)) if kind=='known_times' else list(range(nt))
        def expand(x):
            z=truth_z.copy();z[ids]=x;return z
        unique=list({tuple(z[ids]):z[ids] for z in starts}.values())
        ranked=sorted(unique,key=lambda x:objective.nll(expand(x)))[:6]
        rr=[optimize(lambda x:objective.nll(expand(x)),z,[limits[i] for i in ids],cfg['maxiter']*2) for z in ranked]
        best=min(rr,key=lambda r:r['value'])
        return objective.result(expand(best['z']),converged=best['converged'],truth_assisted=kind),rr
    fits={};trials={}
    for kind in ('known_times','known_charges'):fits[kind],trials[kind]=constrained(kind)
    r=optimize(objective.nll,truth_z,limits,cfg['maxiter']*2)
    fits['truth_initialized_free']=objective.result(r['z'],converged=r['converged'],truth_assisted='all_parameters_initialization_only')
    fits['true_parameters']=objective.result(truth_z,converged=True,truth_assisted='true_parameters_no_fit')
    return dict(fits=fits,trials=trials,truth_initialized_trial=r,outside_time_domain=False,truth_assisted=True)
