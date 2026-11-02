"""Fixed-mean energy diagnostics with explicitly separated truth interfaces.

Candidate search accepts voltage only. Known count is applied by selecting an
already frozen candidate; known times are accepted only by oracle_diagnostics.
All candidates in a frame share one observation-derived weighting matrix.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import solve_banded
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares, minimize, nnls
from scipy.signal import find_peaks


class VoltageSeeds:
    """Existing PairFitter variable-projection starts using frozen kernel arrays."""
    def __init__(self, response):
        self.r = response
        self.slope = np.gradient(response.base_kernel, response.base_grid)

    def columns(self, times):
        return np.interp(self.r.sample_grid[:, None]-np.asarray(times), self.r.base_grid,
                         self.r.base_kernel, left=0., right=0.)

    def refine(self, signal, times):
        times = np.asarray(times, float)
        def objective(t):
            h = self.columns(t); q = nnls(h, signal, maxiter=100)[0]; residual = h@q-signal
            slope = np.interp(self.r.sample_grid[:, None]-t, self.r.base_grid, self.slope, left=0., right=0.)
            return float(residual@residual), -2*np.sum(slope*q[None]*residual[:, None], axis=0)
        fit = minimize(objective, times, method='L-BFGS-B', jac=True,
            bounds=[(0., float(self.r.sample_grid[-1]))]*len(times),
            options=dict(maxiter=100, ftol=1e-12, gtol=1e-6, maxls=30))
        t = fit.x if np.isfinite(fit.fun) and fit.fun <= objective(times)[0] else times
        q = nnls(self.columns(t), signal, maxiter=100)[0]; order = np.argsort(t)
        return dict(times_ns=t[order].tolist(), charge_pe=q[order].tolist(),
                    cost=objective(t)[0], converged=bool(fit.success))

    def candidates(self, voltage):
        y = np.asarray(voltage)-self.r.baseline-self.r.background
        if y.max() <= 1.: return []
        crossed = np.flatnonzero(y > max(.5, .03*y.max()))
        onset = max(0., float(self.r.sample_grid[crossed[0]])-3.)
        one = self.refine(y, [onset])
        derivative = gaussian_filter1d(np.gradient(y, 2.), 1.)
        peaks, _ = find_peaks(derivative, distance=2, prominence=max(.05, .1*derivative.max()))
        peaks = peaks[np.argsort(derivative[peaks])[-4:]]
        starts = [[one['times_ns'][0], one['times_ns'][0]+g] for g in (4., 12., 40., 100., 200.)]
        starts += [[max(0., float(self.r.sample_grid[i])-4.), max(0., float(self.r.sample_grid[j])-4.)]
                   for k, i in enumerate(peaks) for j in peaks[k+1:]]
        fits = [self.refine(y, np.clip(s, 0., self.r.sample_grid[-1])) for s in starts]
        return [one, min(fits, key=lambda f: f['cost'])]


class EnergyObjective:
    def __init__(self, response, voltage, weight_initial):
        self.r = response
        self.signal = np.asarray(voltage, float)-response.baseline-response.background
        if self.signal.shape != (1024,) or not np.isfinite(self.signal).all():
            raise ValueError('Expected one finite 1024-sample voltage')
        self.factor = response.fixed_metric(weight_initial)
        self.target = self.whiten(response.transform(self.signal))

    def whiten(self, values):
        return solve_banded((len(self.factor)-1, 0), self.factor, values, check_finite=False)

    def residual(self, times, energies):
        curves, _ = self.r.components(times, energies)
        return self.whiten(self.r.transform(curves.sum(0)))-self.target

    def describe(self, times, energies):
        t, e = np.asarray(times, float), np.asarray(energies, float)
        _, q = self.r.components(t, e); order = np.argsort(t)
        residual = self.residual(t, e)
        return dict(times_ns=t[order].tolist(), energy_MeV=e[order].tolist(), charge_pe=q[order].tolist(),
                    cost=float(residual@residual))

    def fit(self, times, energies, *, free_times, max_nfev=120):
        times, energies = np.asarray(times, float), np.asarray(energies, float)
        n = len(times)
        if not n: return dict(self.describe([], []), converged=True, evaluations=0, boundary=False, initial_cost=float(self.target@self.target))
        lo, hi = np.log(self.r.energy_bounds)
        if free_times:
            lower = np.r_[np.zeros(n), np.full(n, lo)]
            upper = np.r_[np.full(n, self.r.sample_grid[-1]/10.), np.full(n, hi)]
            x0 = np.r_[times/10., np.log(energies)]
            unpack = lambda x: (x[:n]*10., np.exp(x[n:]))
        else:
            lower, upper = np.full(n, lo), np.full(n, hi)
            x0 = np.log(energies)
            unpack = lambda x: (times, np.exp(x))
        x0 = np.clip(x0, lower+1e-9, upper-1e-9)
        def residual(x): return self.residual(*unpack(x))
        r0 = residual(x0); cost0 = float(r0@r0)
        fit = least_squares(residual, x0, bounds=(lower, upper), max_nfev=max_nfev,
            ftol=1e-9, xtol=1e-9, gtol=1e-5, diff_step=2e-5)
        use_fit = np.isfinite(fit.fun).all() and float(fit.fun@fit.fun) <= cost0
        x = fit.x if use_fit else x0
        return dict(self.describe(*unpack(x)), initial_cost=cost0, converged=bool(fit.success and use_fit),
            evaluations=int(fit.nfev), boundary=bool(np.any(x-lower < 1e-5) or np.any(upper-x < 1e-5)),
            optimality=float(fit.optimality), message=str(fit.message))

    def linear_energy_start(self, times):
        # Only known/observed times and voltage are used; no true Q/E.
        columns = self.r.legacy_components(times, np.full(len(times), 3000.)).T/3000.
        transformed = np.column_stack([self.whiten(self.r.transform(columns[:, j])) for j in range(len(times))])
        q = nnls(transformed, self.target, maxiter=100)[0]
        return self.r.initial_energy(q)


def blind_candidates(response, voltage):
    """Voltage-only finite search at K=0/1/2; no count or true timing input."""
    raw = VoltageSeeds(response).candidates(voltage)
    initial = raw[0] if raw else dict(times_ns=[], charge_pe=[])
    objective = EnergyObjective(response, voltage, initial)
    zero = objective.fit([], [], free_times=False)
    if not raw:
        return dict(weight_initial=initial, candidates=[zero], triggered=False, search_count=0)
    one = objective.fit(raw[0]['times_ns'], response.initial_energy(raw[0]['charge_pe']), free_times=True)
    starts = [(raw[1]['times_ns'], response.initial_energy(raw[1]['charge_pe']))]
    origin, total = float(one['times_ns'][0]), float(sum(raw[0]['charge_pe']))
    for gap in [4., 8., 12., 20., 40., 80., 160., 320.]:
        for sign in [-1, 1]:
            times = np.sort([origin, origin+sign*gap])
            if times[0] < 0 or times[-1] > response.sample_grid[-1]: continue
            options = [response.initial_energy(total*np.array([share, 1-share])) for share in [.03, .15, .5, .85, .97]]
            energy = min(options, key=lambda e: objective.describe(times, e)['cost'])
            starts.append((times, energy))
    ranked = sorted(starts, key=lambda s: objective.describe(*s)['cost'])
    # Preserve the raw two-pulse start even when coarse mean scores rank it low.
    shortlist = [starts[0]]+ranked[:4]
    screened = [objective.fit(t, e, free_times=True, max_nfev=24) for t, e in shortlist]
    polished = [objective.fit(s['times_ns'], s['energy_MeV'], free_times=True, max_nfev=120)
                for s in sorted(screened, key=lambda s: s['cost'])[:2]]
    two = min(screened+polished, key=lambda s: s['cost'])
    return dict(weight_initial=initial, candidates=[zero, one, two], triggered=True,
                search_count=len(shortlist)+len(polished), two_start_min_cost=float(min(s['initial_cost'] for s in screened)))


def oracle_diagnostics(response, voltage, blind, known_times):
    """Truth-assisted controls; true energy/charge is never accepted.

    fixed_time: true count and fixed true times. free_from_true_time: additional
    truth initialization, with full observed time bounds. known_count: selected
    from previously frozen voltage-only candidates by true count only.
    """
    t = np.sort(np.asarray(known_times, float)); n = len(t)
    objective = EnergyObjective(response, voltage, blind['weight_initial'])
    if n == 0:
        fixed = objective.fit([], [], free_times=False)
    else:
        starts = [objective.linear_energy_start(t)]
        total = max(sum(blind['weight_initial']['charge_pe']), 100.)
        if n == 2:
            starts += [response.initial_energy(total*np.array([s, 1-s])) for s in [.03, .15, .5, .85, .97]]
        fits = [objective.fit(t, e, free_times=False) for e in starts]
        fixed = min(fits, key=lambda f: f['cost'])
    free = objective.fit(fixed['times_ns'], fixed['energy_MeV'], free_times=True)
    assert free['cost'] <= fixed['cost']+1e-7
    selected = blind['candidates'][n] if len(blind['candidates']) > n else blind['candidates'][0]
    # Do not replace a blind candidate with a truth-initialized solution.
    return dict(fixed_true_time=fixed, true_count_blind_start=selected, free_from_true_time=free,
                supplied_true_count=n, missing_blind_candidate=bool(len(blind['candidates']) <= n))


def evidence(candidates):
    if len(candidates) < 3 or any(q < 500. for q in candidates[2]['charge_pe']): return 0.
    return max(0., candidates[1]['cost']-candidates[2]['cost'])


def calibrate_gate(single_candidates):
    edges = [1500., 5000., 12500.]; scores = [[] for _ in range(4)]
    for candidates in single_candidates:
        charge = sum(candidates[1]['charge_pe']) if len(candidates) > 1 else 0.
        scores[int(np.searchsorted(edges, charge, side='right'))].append(evidence(candidates))
    assert all(len(s) >= 39 for s in scores)
    return dict(edges_pe=edges, scores=scores, alpha=.025,
                scope='Development single-event empirical rank; no finite-sample calibration claim under fitted response/mother reuse')


def select(blind, gate):
    candidates = blind['candidates']
    if len(candidates) < 2 or candidates[1]['charge_pe'][0] < 500. or candidates[1]['cost'] >= candidates[0]['cost']:
        return dict(candidates[0], selected_k=0, rank=None, delta_cost=0.)
    score = evidence(candidates)
    group = int(np.searchsorted(gate['edges_pe'], sum(candidates[1]['charge_pe']), side='right'))
    reference = np.asarray(gate['scores'][group]); rank = float((1+np.sum(reference >= score))/(len(reference)+1))
    k = 2 if len(candidates) == 3 and score > 0 and rank <= gate['alpha'] else 1
    return dict(candidates[k], selected_k=k, rank=rank, delta_cost=score)
