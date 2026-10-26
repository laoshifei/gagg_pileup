"""Boundary checks for calibrated false-birth decisions."""
import unittest

from .alpha_pair_null_gate import calibrate, evidence, select


def candidates(gain=.2, q=2700., second=700.):
    return [dict(times_ns=[], charge_pe=[], cost=1000.),
            dict(times_ns=[100.], charge_pe=[q], cost=100.),
            dict(times_ns=[100., 108.], charge_pe=[q-second, second], cost=100.*(1.-gain))]


class NullGateTests(unittest.TestCase):
    def test_small_calibration_cannot_claim_tiny_tail(self):
        cal = calibrate([candidates(.01)]*38, edges=[], mode='relative_residual')
        fit = select(candidates(.9), cal)
        self.assertEqual(len(fit['times_ns']), 1)
        self.assertAlmostEqual(fit['null_rank_p'], 1./39.)

    def test_rank_boundary_with_plus_one(self):
        cal = calibrate([candidates(.01)]*39, edges=[], mode='relative_residual')
        fit = select(candidates(.9), cal)
        self.assertEqual(len(fit['times_ns']), 2)
        self.assertAlmostEqual(fit['null_rank_p'], .025)

    def test_ties_are_conservative(self):
        cal = calibrate([candidates(.2)]*80, edges=[], mode='relative_residual')
        fit = select(candidates(.2), cal)
        self.assertEqual(fit['null_rank_p'], 1.)
        self.assertEqual(len(fit['times_ns']), 1)

    def test_subthreshold_proposal_is_not_evidence(self):
        self.assertEqual(evidence(candidates(.99, second=499.)), 0.)
        cal = calibrate([candidates(.01)]*80, edges=[], mode='relative_residual')
        fit = select(candidates(.99, second=499.), cal)
        self.assertEqual(len(fit['times_ns']), 1)

    def test_bin_is_observed_single_charge(self):
        calibration = [candidates(.1, q=1200., second=500.)]*80 + [candidates(.9)]*80
        cal = calibrate(calibration, edges=[1500.], mode='relative_residual')
        a = select(candidates(.5, q=1200., second=500.), cal)
        b = select(candidates(.5), cal)
        self.assertTrue(a['birth_accepted'])
        self.assertFalse(b['birth_accepted'])

    def test_empty_calibration_bin_fails_closed(self):
        with self.assertRaises(ValueError):
            calibrate([candidates()]*80, edges=[1500.], mode='relative_residual')

    def test_calibrated_birth_can_recover_below_old_penalty(self):
        cal = calibrate([candidates(0.)]*80, edges=[], mode='relative_residual', require_birth_penalty=False)
        # Gain 0.05 / zero cost1000 =5e-5, below the old1e-3 penalty.
        fit = select(candidates(.0005), cal)
        self.assertTrue(fit['birth_accepted'])
        self.assertEqual(len(fit['times_ns']), 2)


if __name__ == '__main__':
    unittest.main()
