"""Analytic and accounting checks for the joint compatibility diagnostic."""
import unittest
import numpy as np

from .alpha_pair_blind_uncertainty import (
    exterior_boxes, box_encode, box_decode, exterior_search, empirical_level,
    in_domain, diagnostic_label,
)


class JointCompatibilityTests(unittest.TestCase):
    def setUp(self):
        self.candidate = dict(times_ns=[100., 300.], charge_pe=[10000., 800.])
        self.domain = dict(time=[0., 800.], charge=[50., 50000.])

    def test_true_denominator_charge_tolerance(self):
        boxes = exterior_boxes(self.candidate, self.domain)
        low, high = boxes[6:]
        self.assertAlmostEqual(np.exp(low['limits'][3][1]), 800/1.2)
        self.assertAlmostEqual(np.exp(high['limits'][3][0]), 800/.8)
        self.assertAlmostEqual(800/np.exp(low['limits'][3][1])-1., .2)
        self.assertAlmostEqual(800/np.exp(high['limits'][3][0])-1., -.2)

    def test_box_transform_respects_order_and_exterior(self):
        for box in exterior_boxes(self.candidate, self.domain):
            for t in ([400., 100.], [-1., 900.], [100., 100.]):
                z = box_encode(t, [10., 1e6], box['limits'])
                tt, qq = box_decode(z, box['limits'])
                self.assertLessEqual(tt[0], tt[1]+1e-10)
                x = np.r_[tt, np.log(qq)]
                for value, (lo, hi) in zip(x, box['limits']):
                    self.assertGreaterEqual(value+1e-9, lo)
                    self.assertLessEqual(value-1e-9, hi)

    def test_separable_quadratic_has_known_exterior_minimum(self):
        reference = np.r_[self.candidate['times_ns'], np.log(self.candidate['charge_pe'])]
        scales = np.array([2., 2., .05, .4])
        def objective(t, q): return float(np.sum(((np.r_[t, np.log(q)]-reference)/scales)**2))
        settings = dict(starts=1, maxiter=100, verify_starts=2, verify_maxiter=180,
                        stability_cost=.05, q_tolerance=.2, t_tolerance_ns=10.)
        one = dict(times_ns=[100.], charge_pe=[10800.])
        result = exterior_search(objective, one, self.candidate, self.domain, settings)
        expected = (np.log(1.2)/.4)**2
        self.assertAlmostEqual(result['exterior_delta'], expected, places=6)
        self.assertLess(result['cost_replay_error'], 1e-12)
        self.assertTrue(result['stable'])
        self.assertTrue(in_domain(result['best_exterior']['times_ns'], result['best_exterior']['charge_pe'], self.domain))

    def test_empty_exterior_at_parameter_boundary(self):
        candidate = dict(times_ns=[0., 800.], charge_pe=[50., 50000.])
        boxes = exterior_boxes(candidate, self.domain)
        self.assertFalse(boxes[0]['possible'])
        self.assertFalse(boxes[3]['possible'])
        self.assertFalse(boxes[4]['possible'])
        self.assertFalse(boxes[7]['possible'])

    def test_calibration_does_not_drop_failures_or_claim_small_n(self):
        result = empirical_level(range(20))
        self.assertEqual(result['level'], 18.)
        self.assertIsNone(empirical_level(range(19))['level'])
        self.assertIsNone(empirical_level([0.]*18+[np.inf, np.inf])['level'])
        self.assertEqual(empirical_level([-1.]*20)['level'], 0.)

    def test_witness_remains_valid_even_if_search_is_unstable(self):
        row = dict(selected=dict(times_ns=[1., 2.], converged=True), exterior=dict(
            exterior_delta=3., stable=False, all_exterior_minima_converged=False,
            tolerance_box_hits_domain=True))
        self.assertEqual(diagnostic_label(row, 4.), 'alternative_found')
        self.assertEqual(diagnostic_label(row, 2.), 'numerically_unresolved')
        row['selected']['times_ns'] = [1.]
        self.assertEqual(diagnostic_label(row, 4.), 'not_selected_K2')


if __name__ == '__main__': unittest.main()
