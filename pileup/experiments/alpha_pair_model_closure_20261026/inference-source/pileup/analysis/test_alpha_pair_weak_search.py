"""Physical invariants and a separately soluble nuisance-profile check."""
import unittest
import numpy as np

from .alpha_pair_weak_search import split_proposals, SHARES
from .alpha_pair_profile_diagnostics import profile_weak


class GaussianResponse:
    baseline = 0.
    sample_grid = np.arange(0., 128., 2.)

    def components(self, times, charges):
        return np.asarray([q * np.exp(-.5 * ((self.sample_grid - t) / 6.) ** 2) / 1000.
                           for t, q in zip(times, charges)])


class WhiteNoise:
    background_mV = 0.
    def transform(self, value): return np.asarray(value)
    def score(self, residual, curves, charges, **kwargs): return float(np.asarray(residual) @ np.asarray(residual))


class WeakSearchTests(unittest.TestCase):
    def test_proposals_conserve_charge_at_optimizer_boundary_and_both_orders(self):
        for total in (100., 600., 20000., 50000.):
            bank = split_proposals(1000., total, 2046., SHARES['weak_extra'])
            self.assertEqual({p['direction'] for p in bank}, {-1, 1})
            for p in bank:
                self.assertLess(p['times_ns'][0], p['times_ns'][1])
                self.assertGreaterEqual(min(p['charge_pe']), 50.)
                self.assertAlmostEqual(sum(p['charge_pe']), total)
        self.assertEqual(split_proposals(100., 99., 2046., SHARES['weak_extra']), [])
        bank = split_proposals(1000., 20000., 2046., SHARES['weak_extra'])
        self.assertTrue(any(np.allclose(p['charge_pe'], [200., 19800.]) for p in bank))
        self.assertTrue(any(np.allclose(p['charge_pe'], [19800., 200.]) for p in bank))

    def test_profile_recovers_separated_nuisance_without_moving_fixed_weak_parameters(self):
        response = GaussianResponse(); noise = WhiteNoise()
        y = response.components([28., 90.], [9000., 650.]).sum(0)
        fit = profile_weak(response, noise, y, 90., 650., [(24., 7000.), (32., 11000.)], (12., 44.))
        self.assertEqual(fit['times_ns'][1], 90.)
        self.assertEqual(fit['charge_pe'][1], 650.)
        self.assertAlmostEqual(fit['times_ns'][0], 28., places=3)
        self.assertAlmostEqual(fit['charge_pe'][0], 9000., delta=.2)
        self.assertLess(fit['cost'], 1e-6)
        # At fixed true time the nuisance amplitude has an independent analytic LS solution.
        unit = response.components([28.], [1.])[0]
        residual = y - response.components([90.], [650.])[0]
        analytic_q = (unit @ residual) / (unit @ unit)
        self.assertAlmostEqual(fit['charge_pe'][0], analytic_q, delta=.2)

    def test_profile_recomputes_candidate_covariance_and_logdet(self):
        class ChargeNoise(WhiteNoise):
            def score(self, residual, curves, charges, **kwargs):
                variance = 1. + np.asarray(curves).sum(0) ** 2 * .1
                return float(np.sum(np.asarray(residual) ** 2 / variance + np.log(variance)))
        response = GaussianResponse(); noise = ChargeNoise()
        y = response.components([28., 90.], [9000., 650.]).sum(0)
        fit = profile_weak(response, noise, y, 90., 650., [(28., 9000.)], (27.999, 28.001))
        grid = np.linspace(8000., 9500., 1501)
        costs = [noise.score(y - response.components([28., 90.], [q, 650.]).sum(0),
                            response.components([28., 90.], [q, 650.]), [q, 650.]) for q in grid]
        self.assertAlmostEqual(fit['charge_pe'][0], grid[np.argmin(costs)], delta=2.)
        self.assertLess(fit['charge_pe'][0], 9000.)


if __name__ == '__main__': unittest.main()
