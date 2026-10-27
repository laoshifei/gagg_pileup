"""Numerical covariance checks independent of the production pulse fixtures."""
import unittest

import numpy as np

from .alpha_pair_empirical_noise import LowRankPrecision, calibrate


class EmpiricalNoiseTests(unittest.TestCase):
    def test_whitening_agrees_with_dense_inverse(self):
        rng=np.random.default_rng(267001)
        diagonal=rng.uniform(.2,2.,31);modes=rng.normal(size=(31,5))
        metric=LowRankPrecision(diagonal,modes)
        x=rng.normal(size=31);y=rng.normal(size=31)
        covariance=np.diag(diagonal)+modes@modes.T
        self.assertAlmostEqual(float(metric.apply(x)@metric.apply(y)),
                               float(x@np.linalg.solve(covariance,y)),places=10)
        np.testing.assert_allclose(metric.apply(np.column_stack([x,y]))[:,0],metric.apply(x),atol=1e-13)

    def test_diagonal_limit_and_fixed_pilot_normalization(self):
        diagonal=np.arange(1.,18.);x=np.linspace(1.,3.,17)
        metric=LowRankPrecision(diagonal,np.empty((17,0)),pilot=x)
        self.assertAlmostEqual(float(metric.apply(x)@metric.apply(x)),1e6,places=7)
        np.testing.assert_allclose(metric.apply(2.*x),2.*metric.apply(x),rtol=1e-14)

    def test_covariance_calibration_does_not_apply_mean_as_a_signal(self):
        rng=np.random.default_rng(267002);grid=np.arange(19.)
        residuals=rng.normal(size=(80,19))+5.
        calibration=calibrate(grid,residuals,np.repeat([700.,2700.,9000.,17000.],20),rank=4)
        self.assertEqual(calibration.counts.tolist(),[20]*4)
        self.assertGreater(float(calibration.mean_residual.mean()),4.)
        self.assertLess(float(calibration.variance.mean()),2.)
        self.assertTrue(np.all(calibration.variance>0.))


if __name__=='__main__':
    unittest.main()
