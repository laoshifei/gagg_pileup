"""Dense linear-algebra and distribution checks for the matched-model sampler."""
import unittest

import numpy as np

from .alpha_pair_joint_likelihood import JointNoise
from .alpha_pair_model_sampler import GaussianFactor, lift_innovations


class ModelSamplerTests(unittest.TestCase):
    def make(self, count, n=19):
        rng = np.random.default_rng(2602601)
        noise = JointNoise(np.array([1.]), np.array([.3, -.08, .01, .0]),
                           np.array([-.7, .08, -.1]), .39, -.2, .95)
        curves = rng.uniform(.01, 3., (count, n)); charges = np.array([800., 18000.])[:count]
        band, V = noise.covariance(curves, charges)
        dense = sum(np.diag(band[k, :n-k], -k) + (np.diag(band[k, :n-k], k) if k else 0)
                    for k in range(len(band))) - V @ V.T
        return noise, curves, charges, GaussianFactor.from_covariance(band, V), dense

    def test_full_square_root_and_inverse_against_dense_covariance(self):
        for count in (0, 1, 2):
            _, _, _, factor, covariance = self.make(count)
            square_root = factor.color(np.eye(factor.size)).T
            np.testing.assert_allclose(square_root @ square_root.T, covariance, atol=2e-12)
            np.testing.assert_allclose(factor.whiten(square_root.T), np.eye(factor.size), atol=2e-12)
            self.assertAlmostEqual(factor.logdet, np.linalg.slogdet(covariance)[1], places=11)

    def test_quadratic_matches_dense_solve_and_frozen_likelihood(self):
        noise, curves, charges, factor, covariance = self.make(2)
        r = np.random.default_rng(2602602).normal(size=(4, factor.size))
        quadratic = np.sum(factor.whiten(r)**2, axis=1)
        np.testing.assert_allclose(quadratic, np.sum(r * np.linalg.solve(covariance, r.T).T, axis=1), atol=1e-10)
        for rr, value in zip(r, quadratic):
            self.assertAlmostEqual(noise.score(rr, curves, charges), value + factor.logdet, places=9)

    def test_monte_carlo_covariance(self):
        _, _, _, factor, covariance = self.make(2, n=9)
        draws = factor.draw(np.random.default_rng(2602603), 50000)
        normalized = np.cov(draws, rowvar=False) / np.sqrt(np.outer(np.diag(covariance), np.diag(covariance)))
        expected = covariance / np.sqrt(np.outer(np.diag(covariance), np.diag(covariance)))
        self.assertLess(np.max(np.abs(normalized - expected)), .025)
        self.assertLess(np.max(np.abs(draws.mean(0)) / np.sqrt(np.diag(covariance))), .02)

    def test_fir_right_inverse_and_supplied_boundary(self):
        noise = JointNoise(np.array([1., -1.6405963003788442, .8010774506327718, -.12245642825298192]),
                           np.array([1., 0.]), np.zeros(3), 0., 0., 0.)
        z = np.random.default_rng(2602604).normal(size=1021)
        for prefix in (None, np.array([.1, -.2, .3])):
            y = lift_innovations(z, noise.coefficients, prefix)
            np.testing.assert_allclose(noise.transform(y), z, atol=2e-13)
            np.testing.assert_allclose(y[:3], np.zeros(3) if prefix is None else prefix, atol=1e-15)

    def test_invalid_downdate_is_rejected(self):
        with self.assertRaises(np.linalg.LinAlgError):
            GaussianFactor.from_covariance(np.array([[1., 1., 1.]]), np.ones((3, 1)))


if __name__ == '__main__':
    unittest.main()
