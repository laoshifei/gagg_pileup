"""Numerical controls for finite-window unfolding and recovery inversion."""
import json
from pathlib import Path
import unittest

import numpy as np

from .alpha_pair_deconvolution import AgeRecovery, EnergyProfiles, FiniteConvolution, spe_kernel


class DeconvolutionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = next(p for p in Path(__file__).resolve().parents if (p/'g4').is_dir())
        cls.config = json.loads((root/'pileup/experiments/alpha_pair_energy_calibrated_20261030/sensor-config.json').read_text())

    def test_finite_adjoint_and_charge_units(self):
        kernel = spe_kernel(self.config, 256)
        operator = FiniteConvolution(kernel, 256)
        rng = np.random.default_rng(31); x, y = rng.normal(size=(2, 256))
        self.assertAlmostEqual(float(operator.forward(x)@y), float(x@operator.adjoint(y)), places=11)
        # Average-phase voltage area must approach R*full_cell_charge.
        area = kernel.sum()*2.
        expected = 50*1.7e6*1.602176634e-19*1e12
        self.assertLess(abs(area/expected-1), .001)
        self.assertGreater(operator.exposure[0], 100*operator.exposure[-1])

    def test_rl_descent_with_background_and_window_edge(self):
        h = np.exp(-np.arange(32)/4.)*.3
        op = FiniteConvolution(h, 64); x = np.zeros(64); x[[12, 55]] = [30, 20]
        y = op.forward(x)+2.
        outputs = op.unfold(y, 2., (32, 128, 512))
        losses = []
        for value in outputs.values():
            self.assertTrue(np.all(value >= 0))
            mu = op.forward(value)+2.
            losses.append(float(np.sum(mu-y+y*np.log(y/mu))))
        self.assertTrue(np.all(np.diff(losses) < 0))
        self.assertLess(losses[-1], losses[0]*.01)
        self.assertLess(abs(outputs[512].sum()/x.sum()-1), .01)

    def test_recovery_causality_low_rate_and_zero(self):
        m = AgeRecovery(self.config['sensor'])
        np.testing.assert_array_equal(m.forward(np.zeros(256)), np.zeros(256))
        x = np.zeros(256); x[40:70] = .001
        y = m.forward(x)
        self.assertEqual(float(y[:40].sum()), 0.)
        self.assertAlmostEqual(float(y.sum()/x.sum()), 1/(1-.03), places=5)
        later = x.copy(); later[120:] = 100.
        np.testing.assert_array_equal(y[:120], m.forward(later)[:120])

    def test_nonlinear_inverse_and_grid_convergence(self):
        m = AgeRecovery(self.config['sensor'])
        def pulse(dt):
            t = (np.arange(int(1200/dt))+.5)*dt
            return dt*(140*np.exp(-t/95)+15*np.exp(-t/351))
        x = pulse(2.); y = m.forward(x)
        inferred, diag = m.invert(y, tolerance=1e-7)
        self.assertTrue(diag['converged'])
        self.assertLess(float(np.abs(inferred-x).sum()/x.sum()), 2e-7)
        fine = AgeRecovery(self.config['sensor'], 1.).forward(pulse(1.)).reshape(-1, 2).sum(1)
        self.assertLess(float(np.abs(fine-y).sum()/fine.sum()), .01)

    def test_known_time_energy_inverse(self):
        relative = np.arange(-64, 1722, 2.)
        shape = np.exp(-np.maximum(relative, 0)/200.)*(relative >= 0)/100.
        knots = np.array([.2, 1.5, 5., 7.5]); q = np.array([500., 3000., 10000., 17000.])
        model = EnergyProfiles(relative, knots, np.stack([shape]*4), knots, q, np.full(1024, .01), np.ones(1024))
        for times, energies in [([180.], [4.]), ([180., 220.], [7., .3]), ([180., 188.], [.3, 7.])]:
            x = model.background+model.components(times, energies).sum(0)
            fit = model.fit(x, times)
            self.assertTrue(fit['converged'])
            np.testing.assert_allclose(fit['energy_MeV'], energies, rtol=1e-5)


if __name__ == '__main__':
    unittest.main()
