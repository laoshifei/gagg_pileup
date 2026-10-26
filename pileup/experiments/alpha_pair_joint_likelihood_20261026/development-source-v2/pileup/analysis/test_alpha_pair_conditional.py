"""Scientific invariants of conditional response fitting and acceptance."""
import json
import os
from pathlib import Path
import unittest

import numpy as np

from .alpha_pair_conditional import ConditionalResponse,features,selective_refit,select_count


class ConditionalTests(unittest.TestCase):
    def test_later_neighbor_is_causal(self):
        grid=np.arange(-16.,100.,2.)
        x=features(grid,6000.,14000.,40.)
        np.testing.assert_array_equal(x[grid<=40.,6:],0.)
        np.testing.assert_array_equal(x[:,3:6],0.)
        self.assertGreater(np.linalg.norm(x[grid>40.,6:]),0.)
        past=features(grid,6000.,14000.,-40.)
        np.testing.assert_array_equal(past[:,6:],0.)
        self.assertGreater(np.linalg.norm(past[:,3:6]),0.)

    def test_empty_and_same_count_gate(self):
        empty=dict(times_ns=[],charge_pe=[])
        self.assertFalse(selective_refit(empty,empty,gap_threshold_ns=16.,max_log_charge_change=.25)['accepted_refit'])
        raw=dict(times_ns=[100.,108.],charge_pe=[6000.,7000.])
        refit=dict(times_ns=[100.,108.],charge_pe=[3000.,10000.],converged=True)
        gated=selective_refit(raw,refit,gap_threshold_ns=16.,max_log_charge_change=.25)
        self.assertFalse(gated['accepted_refit']);self.assertEqual(gated['charge_pe'],raw['charge_pe'])
        with self.assertRaises(ValueError):
            selective_refit(raw,empty,gap_threshold_ns=0.,max_log_charge_change=1.)
        self.assertEqual(select_count([dict(empty,cost=100.),dict(times_ns=[100.],charge_pe=[499.],cost=0.)],0.)['times_ns'],[])

    def test_blind_refit_recovers_off_grid_pulses_and_fixed_metric(self):
        root=Path(os.environ.get('GAGG_CONDITIONAL_PROJECT',Path(__file__).resolve().parents[2]))
        settings=json.loads((root/'pileup/experiments/alpha_pair_resolution_20261022/settings.json').read_text())
        grid=np.arange(-16.,1720.01,2.)
        model=ConditionalResponse(settings,grid,np.zeros((len(grid),9)),
            dict(electronic_sigma_mV=.1,dark_variance_mV2=.01,shot_coefficient_mV=1.,rho=.9))
        model.coefficients[:]=1.
        np.testing.assert_array_equal(model.unit(np.arange(-10.,0.),8000.),0.)
        model.coefficients[:]=0.
        times=np.array([121.31,163.77]);q=np.array([3800.,8200.])
        voltage=model.baseline+model.components(times,q).sum(0)
        metric=model.metric(voltage);before=metric.apply(voltage).copy();logdet=metric.logdet
        fit=model.refine(voltage,dict(times_ns=[123.,161.],charge_pe=[5000.,7000.]),metric=metric)
        np.testing.assert_allclose(fit['times_ns'],times,atol=.05)
        np.testing.assert_allclose(fit['charge_pe'],q,rtol=.005)
        np.testing.assert_array_equal(metric.apply(voltage),before)
        self.assertEqual(metric.logdet,logdet)


if __name__=='__main__':unittest.main()
