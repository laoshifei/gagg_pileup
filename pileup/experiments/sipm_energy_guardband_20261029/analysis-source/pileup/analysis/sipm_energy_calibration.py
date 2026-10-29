"""Monotone alpha charge-to-deposited-energy calibration and blind scoring.

Calibration may use independent labels. ``invert_charge`` accepts charge and a
frozen calibration only. Energy outside the observed calibration range is not
clipped or silently extrapolated. No second Birks correction is applied.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares, linear_sum_assignment


METHODS = ("linear_energy", "yield_scale", "yield_saturation", "empirical_pchip")


def fit_calibration(energy, charge, groups, table_energy, table_light, method):
    e, q, g = np.asarray(energy, float), np.asarray(charge, float), np.asarray(groups)
    if method not in METHODS or len(e) < 8 or not (len(e) == len(q) == len(g)):
        raise ValueError("Invalid calibration input")
    if not np.all(np.isfinite(e) & np.isfinite(q) & (e > 0) & (q > 0)):
        raise ValueError("Calibration must retain valid positive deposited energy/charge")
    te, tl = np.asarray(table_energy), np.asarray(table_light)
    if np.any(np.diff(te) <= 0) or np.any(np.diff(tl) <= 0):
        raise ValueError("Cumulative yield must be strictly increasing")
    if e.min() < te.min() or e.max() > te.max():
        raise ValueError("Energy outside frozen material table")
    levels = np.unique(g)
    if len(levels) < 4:
        raise ValueError("Require four independent source-energy groups")
    # Equal group weight and a fixed group mean scale, never 1/observed Q.
    weights = np.zeros_like(e)
    for level in levels:
        mask = g == level
        weights[mask] = 1 / (np.sqrt(mask.sum()) * q[mask].mean())
    grid = np.linspace(float(e.min()), float(e.max()), 10001)
    light = np.interp(e, te, tl)
    grid_light = np.interp(grid, te, tl)
    parameters = None
    if method in ("linear_energy", "yield_scale"):
        x = e if method == "linear_energy" else light
        scale = float(np.sum(weights**2 * x * q) / np.sum(weights**2 * x*x))
        curve = scale * (grid if method == "linear_energy" else grid_light)
        parameters = dict(scale=scale)
    elif method == "yield_saturation":
        norm = float(light.max())
        x = light / norm
        initial = np.sum(weights**2*x*q) / np.sum(weights**2*x*x)
        fit = least_squares(lambda v: weights * (v[0]*x/(1+v[1]*x)-q),
                            [initial, .01], bounds=([1e-12, 0.], [np.inf, 20.]),
                            xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=2000)
        if not fit.success:
            raise RuntimeError("Calibration optimizer did not converge")
        a, b = fit.x
        gx = grid_light/norm
        curve = a*gx/(1+b*gx)
        parameters = dict(scale_pe=float(a), saturation=float(b), light_norm=norm,
                          bound_touched=bool(b > 19.999))
    else:
        knots = sorted((float(e[g == l].mean()), float(q[g == l].mean())) for l in levels)
        ke, kq = np.array(knots).T
        if np.any(np.diff(kq) <= 0):
            raise ValueError("Empirical group means are not monotone")
        curve = PchipInterpolator(ke, kq, extrapolate=True)(grid)
        parameters = dict(knots_energy_MeV=ke.tolist(), knots_charge_pe=kq.tolist(),
                          edge_rule="Endpoint cubic evaluated only within observed calibration E range")
    if not np.all(np.isfinite(curve)) or np.any(np.diff(curve) <= 0) or curve[0] <= 0:
        raise ValueError("Non-monotone/invalid forward calibration curve")
    return dict(method=method, parameters=parameters, energy_MeV=grid.tolist(),
                charge_pe=curve.tolist(), calibration_energy_range_MeV=[float(e.min()), float(e.max())],
                training_events=len(e), groups=levels.tolist(), inference_inputs=["charge_pe"],
                output="deposited_energy_MeV", out_of_range="NaN; preserve denominator/status")


def invert_charge(model, charge):
    q = np.asarray(charge, float)
    return np.interp(q, model["charge_pe"], model["energy_MeV"], left=np.nan, right=np.nan)


def error_summary(truth, estimate):
    e, p = np.asarray(truth, float), np.asarray(estimate, float)
    if len(e) != len(p) or not np.all(np.isfinite(e) & (e > 0)):
        raise ValueError("Positive true deposited energies required")
    valid = np.isfinite(p)
    r = (p[valid]-e[valid])/e[valid]
    delta = p[valid]-e[valid]
    n = len(e)
    return dict(n=n, valid=int(valid.sum()), out_of_domain=int((~valid).sum()),
                median_relative=float(np.median(r)) if len(r) else None,
                sigma68_relative=float((np.quantile(r,.84)-np.quantile(r,.16))/2) if len(r) else None,
                rmse_MeV=float(np.sqrt(np.mean(delta**2))) if len(r) else None,
                within_5=int(np.sum(np.abs(r)<=.05)), within_10=int(np.sum(np.abs(r)<=.10)),
                within_20=int(np.sum(np.abs(r)<=.20)),
                selection_loss=float(np.mean(np.where(valid, ((np.where(valid,p,e)-e)/e)**2, 1.))) if n else None)


def selection_score(model, energy, charge, groups):
    e, q, g = np.asarray(energy), np.asarray(charge), np.asarray(groups)
    pred = invert_charge(model, q)
    summaries = {str(l): error_summary(e[g == l], pred[g == l]) for l in np.unique(g)}
    return float(np.mean([r["selection_loss"] for r in summaries.values()])), summaries


def time_matches(true_t, estimated_t, tolerance=10.):
    """Maximum-cardinality, then minimum-distance, one-to-one matching."""
    t, p = np.asarray(true_t, float), np.asarray(estimated_t, float)
    if not len(t) or not len(p):
        return []
    distance = np.abs(t[:,None]-p[None,:])
    forbidden = (min(len(t),len(p))+1)*(tolerance+1)
    rows, cols = linear_sum_assignment(np.where(distance<=tolerance, distance, forbidden))
    return [(int(i), int(j)) for i,j in zip(rows,cols) if distance[i,j]<=tolerance]


def score_events(true_t, true_q, true_e, estimated_t, estimated_e, gap_ns=None):
    t,q,e,p,ep = map(lambda x: np.asarray(x,float), (true_t,true_q,true_e,estimated_t,estimated_e))
    target = (q >= 500.) & (e > 0)
    t,e = t[target],e[target]
    matched = time_matches(t,p)
    strict_t = min(10.,gap_ns/4) if gap_ns is not None and gap_ns>0 else 10.
    errors=[]; recovered={v:0 for v in (5,10,20)}
    for i,j in matched:
        relative=float((ep[j]-e[i])/e[i]) if np.isfinite(ep[j]) else None
        dt=float(p[j]-t[i])
        errors.append(dict(target_index=int(i), prediction_index=j, true_energy_MeV=float(e[i]),
                           estimated_energy_MeV=float(ep[j]) if np.isfinite(ep[j]) else None,
                           relative_energy_error=relative, time_error_ns=dt))
    for i,j in time_matches(t,p,strict_t):
        relative=float((ep[j]-e[i])/e[i]) if np.isfinite(ep[j]) else None
        for percent in recovered:
            recovered[percent] += bool(relative is not None and abs(relative)<=percent/100)
    return dict(target=int(target.sum()), matched=len(matched), missed=int(target.sum())-len(matched),
                extra=len(p)-len(matched), output_count=len(p),
                energy_out_of_domain=int(np.sum(~np.isfinite(ep))),
                strict_5=recovered[5], strict_10=recovered[10], strict_20=recovered[20],
                frame_success_10=bool(len(p)==target.sum()==recovered[10]), matches=errors)
