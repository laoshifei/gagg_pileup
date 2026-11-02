import unittest
from types import SimpleNamespace

import numpy as np

from .alpha_single_residual_sources import (variance_components, opposite_half_targets,
    prefix_features, select_ridge, ridge_predict)


class ResidualSourceTests(unittest.TestCase):
    def test_random_effects_and_finite_repeat_subtraction(self):
        rng = np.random.default_rng(735001)
        mother = rng.normal(0, 3, (5000, 1, 2))
        values = 2+mother+rng.normal(0, 4, (5000, 16, 2))
        result = variance_components(values)
        np.testing.assert_allclose(result['within'], 16, rtol=.02)
        np.testing.assert_allclose(result['between'], 9, rtol=.05)
        self.assertLess(result['closure_error'], 1e-10)
        # Equal mother means with independent-looking within fluctuations may
        # give a negative estimate. It must not be silently clipped to zero.
        centered = np.tile([[-1., 1.]], (8, 1))
        self.assertLess(variance_components(centered)['between'], 0)

    def test_opposite_half_excludes_input_half(self):
        x = np.arange(32.).reshape(2, 16)
        target = opposite_half_targets(x)
        changed = x.copy(); changed[0, :8] += 100
        np.testing.assert_array_equal(opposite_half_targets(changed)[0, :8], target[0, :8])
        self.assertNotEqual(opposite_half_targets(changed)[0, 8], target[0, 8])

    def test_prefix_has_no_future_samples(self):
        class Model:
            background = 10.
            energy_bounds = (.1, 10.)
            readout = SimpleNamespace(adc_centers=np.arange(1024)*2+.25, aperture=.5)
            def isolated(self, times, energies):
                t = self.readout.adc_centers-times[0]
                shape = np.where(t >= 0, np.maximum(t, 0)*np.exp(-np.maximum(t, 0)/20), 0)
                return np.array([shape*energies[0]]), np.array(energies)
        model = Model()
        y = model.isolated([128.], [3.])[0][0]+10
        a, meta = prefix_features(model, y, 8.)
        changed = y.copy(); changed[model.readout.adc_centers >= 136] = np.nan
        b, _ = prefix_features(model, changed, 8.)
        np.testing.assert_array_equal(a, b)
        self.assertLessEqual(meta['aperture_end_ns'], 136.)
        self.assertAlmostEqual(meta['energy_MeV'], 3., places=5)

    def test_grouped_regression_generalizes_to_new_mothers(self):
        rng = np.random.default_rng(735002)
        x = rng.normal(size=(32, 8, 2))
        y = 1+2*x[..., 0]-.7*x[..., 1]
        fit, _ = select_ridge(x[:16], y[:16])
        error = np.mean((ridge_predict(fit, x[16:])-y[16:])**2)
        self.assertLess(error, .01)


if __name__ == '__main__':
    unittest.main()
