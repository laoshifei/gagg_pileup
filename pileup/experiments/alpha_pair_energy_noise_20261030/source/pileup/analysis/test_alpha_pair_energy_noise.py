"""Numerical controls for diagnostic covariance and no-truth-energy fitting."""
import inspect
import unittest

import numpy as np

from .alpha_pair_energy_noise import (diagonal_factor, fit_known_times,
    replace_conditional_mean, gaussian_waveforms)
from .alpha_pair_model_sampler import lower_band_product
from .test_alpha_pair_recovery_energy import fixture


class EnergyNoiseTests(unittest.TestCase):
    def test_diagonal_preserves_marginal_variance(self):
        band = np.array([[2., 3., 4., 5.], [.4, .6, .8, 0.], [.1, .2, 0., 0.]])
        dense = lower_band_product(band, np.eye(4)).T
        np.testing.assert_allclose(diagonal_factor(band)[0]**2, np.diag(dense@dense.T))

    def test_draws_match_whitened_covariance_and_prefix(self):
        r = fixture(); t = [128., 168.]; e = [.35, 7.]
        curves, q = r.components(t, e); mean = curves.sum(0)+r.baseline+r.background
        factor = r.fixed_metric(dict(times_ns=t, charge_pe=q))
        draws = gaussian_waveforms(r, mean, factor, 2026103001, 2000)
        from scipy.linalg import solve_banded
        white = solve_banded((len(factor)-1, 0), factor, r.transform(draws-mean).T).T
        self.assertLess(abs(white.mean()), .004)
        self.assertLess(abs(white.var()-1), .006)
        self.assertLess(abs(np.mean(white[:, 1:]*white[:, :-1])), .006)
        np.testing.assert_allclose(draws[:, :len(r.ar)-1], mean[None, :len(r.ar)-1]+np.zeros((2000, len(r.ar)-1)))

    def test_reference_is_independent_and_preserves_differences(self):
        rng = np.random.default_rng(7)
        ref, obs = rng.normal(size=(10, 20)), rng.normal(size=(12, 20))
        new = replace_conditional_mean(obs, ref, np.ones(20))
        np.testing.assert_allclose(new-new[0], obs-obs[0], atol=1e-14)
        np.testing.assert_allclose(new.mean(0)-1., obs.mean(0)-ref.mean(0))
        self.assertGreater(np.linalg.norm(new.mean(0)-1.), .1)

    def test_noiseless_inverse_both_metrics_and_truth_api(self):
        r = fixture(); t = [128., 168.]; e = [.35, 7.]
        curves, q = r.components(t, e); y = curves.sum(0)+r.baseline+r.background
        initial = dict(times_ns=[130.], charge_pe=[float(q.sum())])
        factor = r.fixed_metric(initial)
        for f in [factor, diagonal_factor(factor)]:
            result = fit_known_times(r, y, t, f, initial)
            np.testing.assert_allclose(result['energy_MeV'], e, rtol=2e-5)
            np.testing.assert_array_equal(result['times_ns'], t)
        self.assertEqual(list(inspect.signature(fit_known_times).parameters),
                         ['response', 'voltage', 'known_times', 'factor', 'weight_initial'])


if __name__ == '__main__':
    unittest.main()
