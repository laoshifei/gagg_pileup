"""Numerical controls for candidate covariance and normalized likelihood."""
import inspect
import unittest

import numpy as np
from scipy.linalg import cholesky_banded

from .alpha_pair_energy_covariance import (
    candidate_factor, compare_known_times, covariance_band, gaussian_terms,
)
from .alpha_pair_energy_noise import diagonal_factor
from .alpha_pair_model_sampler import lower_band_product
from .test_alpha_pair_recovery_energy import fixture


class CandidateCovarianceTests(unittest.TestCase):
    def setUp(self):
        self.r = fixture()

    def test_legacy_equivalence_and_component_permutation(self):
        t, e = [128., 168.], [.3, 7.]
        curves, q = self.r.components(t, e)
        legacy = self.r.legacy_components(t, q)
        actual = cholesky_banded(covariance_band(self.r, legacy, q), lower=True)
        expected = self.r.fixed_metric(dict(times_ns=t, charge_pe=q))
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-14)
        np.testing.assert_allclose(covariance_band(self.r, curves, q),
                                   covariance_band(self.r, curves[::-1], q[::-1]), atol=1e-14)

    def test_weak_front_variance_and_causal_support(self):
        one = candidate_factor(self.r, [168.], [7.])
        both = candidate_factor(self.r, [128., 168.], [.3, 7.])
        early = candidate_factor(self.r, [128.], [.3])
        grid = self.r.sample_grid[len(self.r.ar)-1:]
        v1, v2, ve = [diagonal_factor(f)[0]**2 for f in (one, both, early)]
        pre = grid < 128.
        middle = (grid >= 132.) & (grid <= 164.)
        np.testing.assert_allclose(v1[pre], v2[pre], rtol=0, atol=1e-13)
        np.testing.assert_allclose(v2[middle], ve[middle], rtol=0, atol=1e-12)
        self.assertGreater(float(v2[middle].mean()), float(v1[middle].mean()))

    def test_normalized_terms_match_dense_gaussian(self):
        curves, q = self.r.components([128., 168.], [.3, 7.])
        band = covariance_band(self.r, curves, q)[:, 64:96].copy()
        factor = cholesky_banded(band, lower=True)
        lower = lower_band_product(factor, np.eye(32)).T
        dense = lower @ lower.T
        residual = np.random.default_rng(723).normal(size=32)
        quadratic, logdet = gaussian_terms(residual, factor)
        self.assertAlmostEqual(float(quadratic), float(residual @ np.linalg.solve(dense, residual)), places=9)
        self.assertAlmostEqual(logdet, float(np.linalg.slogdet(dense)[1]), places=10)

    def test_expected_gaussian_score_needs_logdet(self):
        t, energy = [128., 168.], np.array([.35, 6.])
        region = slice(61, 157)
        c, q = self.r.components(t, energy)
        mean0 = self.r.transform(c.sum(0))[region]
        factor0 = cholesky_banded(covariance_band(self.r, c, q)[:, region], lower=True)
        n = factor0.shape[1]
        # Symmetric sigma points have exactly zero mean and covariance C0.
        basis = lower_band_product(factor0, np.eye(n))*np.sqrt(n)
        observations = mean0 + np.concatenate([basis, -basis])
        h = 1e-4
        for j in range(2):
            means, quadratics = [], []
            for sign in (-1, 1):
                proposed = energy.copy(); proposed[j] *= np.exp(sign*h)
                curves, charge = self.r.components(t, proposed)
                f = cholesky_banded(covariance_band(self.r, curves, charge)[:, region], lower=True)
                a, b = gaussian_terms(observations-self.r.transform(curves.sum(0))[region], f)
                means.append(float(np.mean(a)+b)); quadratics.append(float(np.mean(a)))
            score = (means[1]-means[0])/(2*h)
            self.assertLess(abs(score), 1e-2)
            self.assertGreater(abs((quadratics[1]-quadratics[0])/(2*h)), .01)

    def test_noiseless_fixed_weight_inverse_and_truth_free_energy_api(self):
        t, e = [128., 168.], [.4, 6.]
        voltage = self.r.components(t, e)[0].sum(0)
        initial = dict(times_ns=[166.], charge_pe=[12000.])
        fits = compare_known_times(self.r, voltage, t, self.r.fixed_metric(initial), initial)
        for name in ("reference_gls", "all_legacy_once", "energy_once", "energy_fgls"):
            np.testing.assert_allclose(fits[name]["energy_MeV"], e, rtol=2e-4, atol=1e-6)
            self.assertEqual(fits[name]["times_ns"], t)
        self.assertAlmostEqual(fits["energy_nll"]["cost"],
                               fits["energy_nll"]["quadratic"]+fits["energy_nll"]["logdet"], places=10)
        params = set(inspect.signature(compare_known_times).parameters)
        self.assertEqual(params, {"response", "voltage", "known_times", "reference_factor", "reference_initial"})


if __name__ == "__main__":
    unittest.main()
