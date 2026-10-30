"""Identifiability controls and separation of truth-assisted/voltage-only APIs."""
import inspect
import unittest

import numpy as np
from scipy.linalg import cho_solve_banded

from .alpha_pair_energy_oracle import EnergyObjective, blind_candidates, oracle_diagnostics, calibrate_gate, select
from .test_alpha_pair_recovery_energy import fixture


class OracleTests(unittest.TestCase):
    def setUp(self):
        self.response = fixture()
        self.times = np.array([128.37, 168.37]); self.energy = np.array([.35, 7.])
        comp, q = self.response.components(self.times, self.energy)
        self.voltage = comp.sum(0)
        self.initial = dict(times_ns=self.times.tolist(), charge_pe=q.tolist())

    def test_whitening_matches_banded_quadratic(self):
        obj = EnergyObjective(self.response, self.voltage, self.initial)
        z = np.random.default_rng(2026102901).normal(size=1021)
        expected = z@cho_solve_banded((obj.factor, True), z)
        self.assertAlmostEqual(float(obj.whiten(z)@obj.whiten(z))/expected, 1., places=12)

    def test_known_times_recover_noiseless_energy(self):
        for gap in [8., 40., 200.]:
            t = np.array([128.37, 128.37+gap]); y = self.response.components(t, self.energy)[0].sum(0)
            obj = EnergyObjective(self.response, y, self.initial)
            f = obj.fit(t, self.energy*[1.3, .7], free_times=False)
            np.testing.assert_array_equal(f['times_ns'], t)
            np.testing.assert_allclose(f['energy_MeV'], self.energy, rtol=1e-5)

    def test_free_time_does_not_increase_cost(self):
        obj = EnergyObjective(self.response, self.voltage, self.initial)
        fixed = obj.fit(self.times+[.4, -.6], self.energy, free_times=False)
        free = obj.fit(fixed['times_ns'], fixed['energy_MeV'], free_times=True)
        self.assertLessEqual(free['cost'], fixed['cost'])
        np.testing.assert_allclose(free['times_ns'], self.times, atol=.01)
        np.testing.assert_allclose(free['energy_MeV'], self.energy, rtol=1e-3)

    def test_known_count_does_not_replace_blind_parameters(self):
        obj = EnergyObjective(self.response, self.voltage, self.initial)
        zero = obj.fit([], [], free_times=False)
        fake = dict(times_ns=[120., 180.], energy_MeV=[.5, 6.], charge_pe=[800., 12000.], cost=1e6)
        blind = dict(weight_initial=self.initial, candidates=[zero, fake, fake])
        result = oracle_diagnostics(self.response, self.voltage, blind, self.times)
        self.assertIs(result['true_count_blind_start'], fake)
        self.assertLess(result['fixed_true_time']['cost'], fake['cost'])
        self.assertNotIn('true_energy', inspect.signature(oracle_diagnostics).parameters)
        self.assertEqual(list(inspect.signature(blind_candidates).parameters), ['response', 'voltage'])

    def test_empty_voltage_and_gate_ties(self):
        blind = blind_candidates(self.response, np.zeros(1024))
        self.assertEqual(len(blind['candidates']), 1)
        single = []
        for q in [700., 2500., 8000., 16000.]:
            for _ in range(40):
                single.append([dict(cost=100., charge_pe=[]), dict(cost=10., charge_pe=[q]),
                               dict(cost=5., charge_pe=[q/2, q/2])])
        gate = calibrate_gate(single)
        self.assertEqual(select(blind, gate)['selected_k'], 0)
        self.assertEqual(select(dict(candidates=single[-1]), gate)['selected_k'], 1)
        result = oracle_diagnostics(self.response, np.zeros(1024), blind, [])
        self.assertEqual(result['fixed_true_time']['energy_MeV'], [])


if __name__ == '__main__': unittest.main()
