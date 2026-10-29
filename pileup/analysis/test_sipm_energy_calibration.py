"""Numerical/denominator checks for energy calibration, independent of G4 data."""
import unittest
import numpy as np
from .sipm_energy_calibration import fit_calibration,invert_charge,error_summary,score_events,time_matches


class EnergyCalibrationTest(unittest.TestCase):
    def test_curved_response_roundtrip(self):
        te=np.linspace(0,10,2001);tl=1000*te**1.3
        e=np.concatenate([np.linspace(v-.1,v+.1,20) for v in (1.,2.,4.,7.)])
        light=np.interp(e,te,tl);q=.3*light/(1+.00002*light)
        model=fit_calibration(e,q,np.repeat([1,2,4,7],20),te,tl,'yield_saturation')
        probe=np.linspace(e.min()+.01,e.max()-.01,57);l=np.interp(probe,te,tl)
        np.testing.assert_allclose(invert_charge(model,.3*l/(1+.00002*l)),probe,atol=3e-5)

    def test_no_endpoint_clipping(self):
        m={'charge_pe':[10.,20.],'energy_MeV':[1.,2.]}
        p=invert_charge(m,[9.,10.,15.,20.,21.,np.nan])
        self.assertTrue(np.isnan(p[[0,4,5]]).all());np.testing.assert_equal(p[1:4],[1.,1.5,2.])

    def test_rejections_keep_denominator(self):
        s=error_summary([1.,2.,3.],[1.,np.nan,3.])
        self.assertEqual((s['n'],s['valid'],s['within_10']),(3,2,2))
        self.assertAlmostEqual(s['selection_loss'],1/3)

    def test_one_to_one_prevents_double_counting(self):
        s=score_events([0.],[1000.],[1.],[0.,1.],[1.,1.])
        self.assertEqual((s['matched'],s['extra'],s['strict_10']),(1,1,1))
        self.assertFalse(s['frame_success_10'])

    def test_close_pair_timing_requirement(self):
        s=score_events([0.,8.],[1000.,1000.],[1.,1.],[0.,11.],[1.,1.],8.)
        self.assertEqual(s['matched'],2);self.assertEqual(s['strict_10'],1)

    def test_maximum_cardinality(self):
        self.assertEqual(len(time_matches([0.,9.],[-9.,1.],10.)),2)

    def test_threshold_does_not_hide_false_outputs(self):
        s=score_events([0.,100.],[1000.,499.],[1.,.4],[0.,100.],[1.,.4],100.)
        self.assertEqual((s['target'],s['extra']),(1,1))

    def test_nonmonotone_yield_rejected(self):
        with self.assertRaises(ValueError):
            fit_calibration(np.ones(8),np.ones(8),np.repeat([1,2,3,4],2),[0,1,2],[0,2,1],'yield_scale')


if __name__=='__main__':unittest.main()
