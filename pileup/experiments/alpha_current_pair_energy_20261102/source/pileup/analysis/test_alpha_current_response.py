import unittest
import numpy as np
from .alpha_current_response import paired_mother_interval, energy_metrics


class PairedControls(unittest.TestCase):
    def test_repeats_do_not_inflate_mother_count(self):
        old = np.zeros((4, 1))
        new = np.array([[1.], [2.], [-1.], [4.]])
        first = paired_mother_interval(old, new, seed=17)
        repeated = paired_mother_interval(np.repeat(old, 8, axis=1), np.repeat(new, 8, axis=1), seed=17)
        self.assertEqual(first['ci95'], repeated['ci95'])
        self.assertEqual(repeated['mothers'], 4)
        self.assertEqual(repeated['mean_change'], 1.5)
        with self.assertRaises(ValueError):
            paired_mother_interval(old, new[:3], seed=17)

    def test_energy_metric_preserves_boundary_and_failures(self):
        result = energy_metrics([-20., 0., 20.], [True, False, False], [True, False, True])
        self.assertEqual(result['strict10_count'], 1)
        self.assertEqual(result['boundary'], 1)
        self.assertEqual(result['nonconverged'], 1)
        self.assertAlmostEqual(result['rmse_percent'], np.sqrt(800/3))


if __name__ == '__main__':
    unittest.main()
