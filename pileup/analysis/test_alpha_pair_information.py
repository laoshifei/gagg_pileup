"""Independent algebra and Monte Carlo checks for Gaussian information."""
import unittest

import numpy as np

from .alpha_pair_information import gaussian_information, invert_information, dense_covariance


class InformationTests(unittest.TestCase):
    def test_scalar_mean_and_log_variance(self):
        # N(mu, exp(eta)): information diag(exp(-eta), 1/2).
        variance = 3.
        r = gaussian_information([[1.], [0.]], [[[0.]], [[variance]]], [[1./np.sqrt(variance)]])
        np.testing.assert_allclose(r['total'], np.diag([1./variance, .5]), atol=1e-14)

    def test_dense_precision_identity_and_coordinate_change(self):
        rng = np.random.default_rng(73101)
        a = rng.normal(size=(7, 7))
        C = a@a.T + np.eye(7)
        inv = np.linalg.inv(C)
        W = np.linalg.inv(np.linalg.cholesky(C))
        dm = rng.normal(size=(3, 7))
        dc = rng.normal(size=(3, 7, 7)); dc = (dc + dc.swapaxes(1, 2))/2.
        actual = gaussian_information(dm, dc, W)['total']
        expected = np.array([[dm[i]@inv@dm[j] + .5*np.trace(inv@dc[i]@inv@dc[j]) for j in range(3)] for i in range(3)])
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
        scaling = np.array([10., .3, 2.])
        changed = gaussian_information(dm*scaling[:, None], dc*scaling[:, None, None], W)['total']
        np.testing.assert_allclose(changed, scaling[:, None]*expected*scaling[None, :], rtol=1e-12)

    def test_score_covariance_from_finite_likelihood_differences(self):
        rng = np.random.default_rng(73102)
        n, p, count = 5, 3, 140000
        a = rng.normal(size=(n, n)); C = a@a.T + 2.*np.eye(n)
        dm = rng.normal(size=(p, n)); dc = rng.normal(size=(p, n, n)); dc = (dc+dc.swapaxes(1, 2))/6.
        y = rng.multivariate_normal(np.zeros(n), C, size=count)
        h = 1e-4
        scores = []
        for i in range(p):
            costs = []
            for sign in (-1., 1.):
                covariance = C + sign*h*dc[i]; residual = y - sign*h*dm[i]
                costs.append(.5*(np.einsum('ni,ij,nj->n', residual, np.linalg.inv(covariance), residual)
                                 + np.linalg.slogdet(covariance)[1]))
            scores.append(-(costs[1]-costs[0])/(2.*h))
        scores = np.asarray(scores).T
        F = gaussian_information(dm, dc, np.linalg.inv(np.linalg.cholesky(C)))['total']
        standardized = scores @ np.linalg.inv(np.linalg.cholesky(F)).T
        np.testing.assert_allclose(np.cov(standardized, rowvar=False), np.eye(p), atol=.015)
        self.assertLess(float(np.max(abs(standardized.mean(0)))), .015)

    def test_nuisance_schur_complement(self):
        F = np.array([[4., 1., .7], [1., 3., .2], [.7, .2, 2.]])
        full = invert_information(F)
        schur = F[1:, 1:] - np.outer(F[1:, 0], F[0, 1:])/F[0, 0]
        np.testing.assert_allclose(full[1:, 1:], np.linalg.inv(schur))
        self.assertTrue(np.all(np.linalg.eigvalsh(full[1:, 1:]-np.linalg.inv(F[1:, 1:])) > -1e-13))

    def test_fixed_covariance_linear_gls(self):
        rng = np.random.default_rng(73103)
        design = rng.normal(size=(8, 3)); C = np.diag(np.linspace(1., 3., 8))
        W = np.diag(1./np.sqrt(np.diag(C)))
        F = gaussian_information(design.T, np.zeros((3, 8, 8)), W)['total']
        expected = np.linalg.inv(design.T@np.linalg.solve(C, design))
        np.testing.assert_allclose(invert_information(F), expected)
        self.assertRaises(np.linalg.LinAlgError, invert_information, np.ones((3, 3)))

    def test_symmetric_band_and_downdate(self):
        band = np.array([[5., 6., 7., 8.], [.2, .3, .4, 0.], [.1, .15, 0., 0.]])
        V = np.array([[.1], [.2], [.3], [.4]])
        expected = np.array([[5., .2, .1, 0.], [.2, 6., .3, .15], [.1, .3, 7., .4], [0., .15, .4, 8.]]) - V@V.T
        np.testing.assert_allclose(dense_covariance(band, V), expected)


if __name__ == '__main__':
    unittest.main()
