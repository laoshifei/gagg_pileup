"""Energy-conditioned mean correction and separately fitted residual noise.

Calibration accepts deposited-energy labels and actual ancestry components.
Inference accepts only voltage and explicitly supplied count/times. This is a
local pure-alpha model, not a microscopic SiPM response or a blind detector.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from types import SimpleNamespace

import numpy as np
from scipy.linalg import cholesky_banded, solve_banded
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize

from .alpha_pair_recovery_energy import shape_features
from .alpha_pair_energy_covariance import covariance_band, candidate_factor
from .alpha_pair_energy_noise import FixedMetricObjective, fit_known_times
from .alpha_pair_energy_oracle import EnergyObjective, VoltageSeeds

METHODS = ("old_mean_old_noise", "new_mean_old_noise", "old_mean_new_noise", "new_mean_new_noise")


def taper(grid, end):
    onset = np.clip(np.asarray(grid)/8., 0., 1.)
    return onset*onset*(3.-2.*onset)*np.clip((end-np.asarray(grid))/120., 0., 1.)


class CalibratedEnergyResponse:
    """A causal residual mean calibrated on E, never on realized pulse charge."""
    def __init__(self, base, grid, coefficients):
        self.base = base
        self.grid = np.asarray(grid, float)
        self.coefficients = np.asarray(coefficients, float)
        if self.coefficients.shape != (len(self.grid), 9):
            raise ValueError("Expected nine smooth mean features")
        self._taper = taper(self.grid, self.grid[-1])

    def __getattr__(self, name):
        return getattr(self.base, name)

    def components(self, times, energies):
        t, e = np.asarray(times, float), np.asarray(energies, float)
        z, neighbor, q0 = self.base.drivers(t, e)
        full, q = self.base.from_drivers(z, neighbor)
        changed = full.copy()
        for j in range(len(e)):
            other = 1-j if len(e)==2 else None
            x = shape_features(self.grid, q0[j], q0[other] if other is not None else 0.,
                               t[other]-t[j] if other is not None else 0.)
            correction = np.sum(x*self.coefficients, axis=1)*self._taper
            values = q0[j]*np.interp(self.base.full_grid-t[j], self.grid, correction, left=0., right=0.)
            changed[j] = np.maximum(full[j]+values, 0.)
        # Preserve the calibrated energy label; do not renormalize to realized Q.
        # The correction has finite support; its voltage area updates model Q.
        q = q+(changed-full).sum(-1)*self.dt/self.area_per_pe
        return changed[:, :self.observed_samples], q

    def save(self, path):
        np.savez_compressed(path, grid=self.grid, coefficients=self.coefficients)

    @classmethod
    def load(cls, base, path):
        with np.load(path) as f:
            return cls(base, f["grid"], f["coefficients"])


def calibrate_mean(base, records, *, ridge=.03):
    """Balanced strata, normalized by frozen q0(E); no realized-Q regressor."""
    grid = np.arange(-16., 1720.01, 2.)
    window = taper(grid, grid[-1])
    counts = Counter((r["condition"], r["gap_ns"], j) for r in records for j in range(len(r["energy_MeV"])))
    gram = np.zeros((len(grid), 9, 9)); rhs = np.zeros((len(grid), 9)); total = np.zeros(len(grid))
    component_rows = 0
    for record in records:
        t, e = np.asarray(record["times_ns"]), np.asarray(record["energy_MeV"])
        predicted, _ = base.components(t, e); q0 = base.q_isolated(e)
        for j in range(len(e)):
            other = 1-j if len(e)==2 else None
            x = shape_features(grid, q0[j], q0[other] if other is not None else 0.,
                               t[other]-t[j] if other is not None else 0.)*window[:,None]
            residual = np.asarray(record["components_mV"])[j]-predicted[j]
            target = np.interp(grid+t[j], base.sample_grid, residual, left=0., right=0.)/q0[j]
            valid = (grid+t[j]>=base.sample_grid[0]) & (grid+t[j]<=base.sample_grid[-1])
            weight = valid.astype(float)/counts[record["condition"], record["gap_ns"], j]
            gram += np.einsum("si,sj,s->sij", x, x, weight)
            rhs += x*(target*weight)[:,None]; total += weight
            component_rows += 1
    gram /= np.maximum(total,1e-12)[:,None,None]
    rhs /= np.maximum(total,1e-12)[:,None]
    regularizer = np.eye(9)*ridge; regularizer[0,0] *= .01
    coefficients = np.linalg.solve(gram+regularizer[None],rhs[...,None])[...,0]
    coefficients = gaussian_filter1d(coefficients, 2., axis=0)
    return CalibratedEnergyResponse(base, grid, coefficients), dict(
        frames=len(records), component_rows=component_rows, strata=len(counts),
        ridge=ridge, smoothing_sigma_ns=4., coordinate="Frozen q0(E_dep), not realized Q",
        waveform_interpolation="Used only for smooth mean regression; unobserved tails excluded",
        coefficient_max_abs=float(abs(coefficients).max()))


class BandedModes:
    """Exact whitening for positive C = L L.T + U U.T, with rank <= 2."""
    def __init__(self, factor, modes=None):
        self.factor = np.asarray(factor, float)
        n = self.factor.shape[1]
        self.modes = np.empty((n,0)) if modes is None else np.asarray(modes,float)
        if self.modes.ndim!=2 or self.modes.shape[0]!=n:
            raise ValueError("Mode/factor dimensions differ")
        self.A = solve_banded((len(self.factor)-1,0),self.factor,self.modes,check_finite=True)
        if self.A.shape[1]:
            values, vectors = np.linalg.eigh(self.A.T@self.A)
            values = np.maximum(values,0.); root = np.sqrt(1.+values)
            self.shrink = (vectors/(root*(1.+root)))@vectors.T
            self.logdet = float(2*np.log(self.factor[0]).sum()+np.log1p(values).sum())
        else:
            self.shrink = np.empty((0,0))
            self.logdet = float(2*np.log(self.factor[0]).sum())

    def apply(self, values):
        w = solve_banded((len(self.factor)-1,0),self.factor,np.asarray(values,float),check_finite=False)
        return w-self.A@(self.shrink@(self.A.T@w))

    def marginal_variance(self):
        n = self.factor.shape[1]; result = np.zeros(n)
        for lag,row in enumerate(self.factor):result[lag:] += row[:n-lag]**2
        return result+np.sum(self.modes*self.modes,axis=1)


def noise_proxy(response, parameters):
    return SimpleNamespace(observed_samples=response.observed_samples, transform=response.transform,
        noise_lags=np.asarray(response.noise_lags)*parameters["background_scale"],
        noise_parameters=np.asarray(parameters["shot_parameters"]),
        lag1=parameters["lag1"], cross=parameters["cross"])


def learned_metric(response, times, energies, parameters):
    curves, _ = response.components(times, energies)
    q0 = response.q_isolated(energies)
    band = covariance_band(noise_proxy(response, parameters), curves, q0)
    factor = cholesky_banded(band, lower=True, check_finite=True)
    variance = np.interp(np.log(q0),np.log(parameters["mode_charge_knots_pe"]),
                         parameters["fractional_mode_variance"])
    modes = response.transform(curves).T*np.sqrt(np.maximum(variance,0.))[None,:]
    return BandedModes(factor,modes)


def calibrate_noise(response, records, empty_voltage):
    """Fit noise on separate optical mothers, on the un-interpolated FIR grid.

    Short-range parameters describe residual second moments and can include
    remaining mean error. Positive amplitude modes use within-stratum centered
    integrated residual variance. Neither is an exact microscopic covariance.
    """
    empty = response.transform(np.asarray(empty_voltage)-response.baseline-response.background)
    centered = empty-empty.mean()
    background_scale = float(np.mean(centered**2)/response.noise_lags[0])
    groups = Counter((r["condition"],r["gap_ns"],j) for r in records for j in range(len(r["energy_MeV"])))
    xrows=[]; targets=[]; weights=[]; component_data=[]; paired_data=[]
    for record in records:
        t,e = record["times_ns"],record["energy_MeV"]
        mean,_ = response.components(t,e); q0=response.q_isolated(e)
        zmean=response.transform(mean); residual=response.transform(np.asarray(record["components_mV"])-mean)
        activity=np.maximum(zmean,0.)
        data=[]
        for j in range(len(e)):
            other=activity[1-j] if len(e)==2 else np.zeros_like(activity[j])
            mask=activity[j]>.02
            x=np.column_stack((np.ones(mask.sum()),np.full(mask.sum(),np.log(q0[j]/3000.)),other[mask]/10.))
            target=residual[j,mask]**2/activity[j,mask]
            weight=np.full(mask.sum(),1./(groups[record["condition"],record["gap_ns"],j]*max(mask.sum(),1)))
            xrows.append(x);targets.append(target);weights.append(weight)
            d=dict(condition=record["condition"],gap_ns=record["gap_ns"],component=j,
                   q0=float(q0[j]), activity=activity[j],neighbor=other,residual=residual[j],
                   integrated_mean=float(zmean[j].sum()),energy_MeV=float(e[j]))
            data.append(d);component_data.append(d)
        if len(data)==2:paired_data.append(data)
    x=np.concatenate(xrows);target=np.concatenate(targets);weight=np.concatenate(weights);weight/=weight.sum()
    def objective(theta):
        raw=x@theta; eta=np.clip(raw,-8.,4.)
        ratio=target*np.exp(-eta)
        value=float(np.sum(weight*(eta+ratio)))
        gradient=x.T@(weight*(1.-ratio)*((raw>-8.)&(raw<4.)))
        return value,gradient
    initial=np.asarray(response.noise_parameters,float)
    fit=minimize(objective,initial,jac=True,method="L-BFGS-B",bounds=[(-8.,2.),(-2.,2.),(-.5,.5)],
                 options=dict(maxiter=120,ftol=1e-11,gtol=1e-7))
    theta=fit.x if fit.fun<=objective(initial)[0] else initial
    def variance(d):
        a,b,c=theta
        return d["activity"]*np.exp(np.clip(a+b*np.log(d["q0"]/3000.)+c*d["neighbor"]/10.,-8.,4.))
    numerator=denoma=denomb=0.
    for d in component_data:
        v=variance(d);mask=(v[1:]>1e-7)&(v[:-1]>1e-7)
        u=d["residual"]/np.sqrt(np.maximum(v,1e-12));w=1./max(mask.sum(),1)
        numerator += w*np.sum(u[1:][mask]*u[:-1][mask])
        denoma += w*np.sum(u[1:][mask]**2);denomb += w*np.sum(u[:-1][mask]**2)
    raw_lag=numerator/np.sqrt(denoma*denomb); lag=float(np.clip(raw_lag,-.45,.45))
    numerator=denoma=denomb=0.
    for pair in paired_data:
        va,vb=map(variance,pair);mask=(va>1e-7)&(vb>1e-7)
        ua=pair[0]["residual"]/np.sqrt(np.maximum(va,1e-12));ub=pair[1]["residual"]/np.sqrt(np.maximum(vb,1e-12))
        w=1./max(mask.sum(),1)
        numerator += w*np.sum(ua[mask]*ub[mask]);denoma += w*np.sum(ua[mask]**2);denomb += w*np.sum(ub[mask]**2)
    raw_cross=numerator/np.sqrt(denoma*denomb) if denoma*denomb>0. else 0.
    cross=float(np.clip(raw_cross,-.45,.45))
    # Four frozen q0 bins identify the four deposited-energy clusters.
    bins=(1500.,5000.,12500.); by_stratum=defaultdict(list)
    for d in component_data:
        v=variance(d); mean=d["integrated_mean"]
        if mean<=0.:raise ValueError("Nonpositive integrated model innovation")
        d["fractional_total_residual"]=float(d["residual"].sum()/mean)
        d["band_fractional_total_variance"]=float((v.sum()+2*lag*np.sqrt(v[1:]*v[:-1]).sum())/mean**2)
        by_stratum[np.searchsorted(bins,d["q0"]),d["condition"],d["gap_ns"],d["component"]].append(d)
    estimates=defaultdict(list)
    for key,rows in by_stratum.items():
        value=float(np.var([d["fractional_total_residual"] for d in rows],ddof=1)
                    -np.mean([d["band_fractional_total_variance"] for d in rows]))
        estimates[int(key[0])].append(value)
    raw_modes=[float(np.mean(estimates[i])) for i in range(4)]
    mode_variance=np.maximum(raw_modes,0.)
    knots=[float(np.median([d["q0"] for d in component_data if np.searchsorted(bins,d["q0"])==i])) for i in range(4)]
    parameters=dict(background_scale=background_scale,shot_parameters=theta.tolist(),lag1=lag,cross=cross,
        mode_charge_knots_pe=knots,fractional_mode_variance=mode_variance.tolist(),
        conditioning="q0(E_dep), never realized avalanche charge",mean_model="Frozen separately fitted energy mean")
    report=dict(frames=len(records),components=len(component_data),empty_frames=len(empty_voltage),
        active_native_FIR_samples=len(target),optimizer_converged=bool(fit.success),optimizer_message=str(fit.message),
        loss_before=objective(initial)[0],loss_after=objective(theta)[0],raw_lag1=float(raw_lag),raw_cross=float(raw_cross),
        raw_fractional_mode_variance=raw_modes,mode_floor_count=int(np.count_nonzero(np.asarray(raw_modes)<0)),
        residual_moment_caveat="Band scale uses uncentered residual second moments; residual mean mismatch may remain. Amplitude mode is centered within condition/gap/component strata.")
    return parameters,report


class MetricObjective(EnergyObjective):
    def __init__(self,response,voltage,metric):
        self.r=response;self.metric=metric
        self.signal=np.asarray(voltage,float)-response.baseline-response.background
        if self.signal.shape!=(1024,) or not np.isfinite(self.signal).all():
            raise ValueError("Expected finite1024-sample voltage")
        self.target=self.whiten(response.transform(self.signal))

    def whiten(self,values):
        return self.metric.apply(values)


def compare_known_times(base,calibrated,noise_parameters,voltage,known_times):
    t=np.asarray(known_times,float)
    if t.ndim!=1 or len(t) not in (1,2) or not np.isfinite(t).all() or np.any(np.diff(t)<0):
        raise ValueError("Expected ordered one/two explicitly supplied times")
    y=np.asarray(voltage,float)-base.baseline-base.background
    crossed=np.flatnonzero(y>max(.5,.03*y.max()))
    onset=max(0.,float(base.sample_grid[crossed[0]])-3.) if len(crossed) else float(base.sample_grid[np.argmax(y)])
    initial=VoltageSeeds(base).refine(y,[onset])
    seed=fit_known_times(base,voltage,t,base.fixed_metric(initial),initial)
    e0=np.asarray(seed["energy_MeV"])
    factor=candidate_factor(base,t,e0)
    old_metric=BandedModes(factor)
    new_metric=learned_metric(calibrated,t,e0,noise_parameters)
    linear=FixedMetricObjective(base,voltage,factor).linear_energy_start(t)
    starts=[e0,linear]
    if len(t)==2:
        total=max(sum(initial["charge_pe"]),100.)
        starts.extend(base.initial_energy(total*np.array([s,1.-s])) for s in (.03,.15,.5,.85,.97))
    rows={}
    for name,model,metric in ((METHODS[0],base,old_metric),(METHODS[1],calibrated,old_metric),
                               (METHODS[2],base,new_metric),(METHODS[3],calibrated,new_metric)):
        objective=MetricObjective(model,voltage,metric)
        fits=[objective.fit(t,e,free_times=False) for e in starts]
        rows[name]=dict(min(fits,key=lambda f:f["cost"]),starts=len(starts),
                       total_evaluations=sum(f["evaluations"] for f in fits),covariance_energy_MeV=e0.tolist())
    return dict(fits=rows,seed=seed,initial=initial,
                scope="Known K/t; fitted energy only; common fixed metrics and common starts within2x2 ablation")
