"""Noiseless and conditioning checks for the truth-assisted diagnostics."""
from pathlib import Path
import os
import unittest

import numpy as np

from .alpha_pair_conditional import ConditionalResponse
from .alpha_pair_weak_diagnostics import fixed_time_charge_fit, model_pilot_metric, oracle_component_refit


class WeakDiagnosticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(os.environ.get('GAGG_WEAK_PROJECT', Path(__file__).resolve().parents[2]))
        cls.response = ConditionalResponse.load(root/'pileup/experiments/alpha_pair_null_gate_20261023/response.npz')

    def test_known_times_recover_asymmetric_charge_without_true_charge_start(self):
        response = self.response
        t = np.array([151.3, 191.3]); q = np.array([16500., 650.])
        y = response.components(t, q).sum(0)+response.baseline
        metric = model_pilot_metric(response, y)
        result = fixed_time_charge_fit(response, y, t, q.sum(), metric=metric)
        np.testing.assert_allclose(result['times_ns'], t, atol=0., rtol=0.)
        np.testing.assert_allclose(result['charge_pe'], q, rtol=1e-5)
        self.assertLess(result['cost'], 1e-7)

    def test_oracle_subtraction_keeps_neighbor_conditioning(self):
        response = self.response
        t = np.array([153.7, 353.7]); q = np.array([16000., 700.])
        components = response.components(t, q)
        y = components.sum(0)+response.baseline
        result = oracle_component_refit(response, y, t, q, components,
                                        metric=model_pilot_metric(response, y))
        np.testing.assert_allclose(result['times_ns'], t, atol=1e-5)
        np.testing.assert_allclose(result['charge_pe'], q, rtol=1e-6)
        self.assertIn('other_component', result['diagnostic_truth'])

    def test_empty_fixed_count_is_explicitly_truth_assisted(self):
        response = self.response
        y = np.full(1024, response.baseline)+np.sin(np.arange(1024))*.1
        result = fixed_time_charge_fit(response, y, [], 0.)
        self.assertEqual(result['times_ns'], [])
        self.assertEqual(result['diagnostic_truth'], ['K', 'times'])
        self.assertGreater(result['cost'], 0.)


if __name__ == '__main__':
    unittest.main()
