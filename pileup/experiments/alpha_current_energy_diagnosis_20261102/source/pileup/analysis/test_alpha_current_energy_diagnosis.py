import unittest
import numpy as np
from .alpha_current_energy_diagnosis import residual_terms, local_operator, repeat_covariance


class ResidualDiagnostics(unittest.TestCase):
    def test_exact_partition_with_shared_changes(self):
        rng = np.random.default_rng(17)
        shared, isolated = rng.normal(size=(2, 8, 2, 31))
        predicted, solo = rng.normal(size=(2, 2, 31))
        y = shared.sum(1)+rng.normal(size=(8, 31))+3.
        for weak in (0, 1):
            pieces, total = residual_terms(y, shared, isolated, predicted, solo, 3., weak)
            np.testing.assert_allclose(pieces.sum(1), total, atol=1e-13)

    def test_local_inverse_and_nuisance_rejection(self):
        x = np.linspace(0., 1., 41)
        d = np.column_stack((1+x, 1+x+.01*np.sin(5*x)))
        operator, geometry = local_operator(d)
        np.testing.assert_allclose(operator@d, np.eye(2), atol=1e-11)
        np.testing.assert_allclose(operator@(d@np.array([.03, -.08])), [.03, -.08], atol=1e-11)
        self.assertGreater(geometry['cosine'], .99)
        with self.assertRaises(ValueError):local_operator(np.column_stack((x, x)))

    def test_correlated_cancellation_not_independent_variances(self):
        x = np.arange(8.)
        terms = np.column_stack((x, -x, np.zeros(8), np.ones(8)))
        result = repeat_covariance(terms)
        self.assertEqual(result['total_variance'], 0.)
        self.assertGreater(result['diagonal_variance'], 0.)
        self.assertAlmostEqual(result['cross_variance'], -result['diagonal_variance'])


if __name__ == '__main__':unittest.main()
