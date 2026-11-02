import unittest
import numpy as np

from .alpha_mean_shape_correction import MeanShapeCorrection,calibrate_profiles,temporal_window


class MeanShapeTests(unittest.TestCase):
    def make_correction(self):
        grid=np.arange(-16.,1602.,2.)+.25
        base=np.exp(-np.maximum(grid,0)/100)*temporal_window(grid)
        base*=np.sin(grid/30)
        weight=temporal_window(grid)
        base-=base.sum()*weight/weight.sum()
        return MeanShapeCorrection(grid,np.log([500.,2000.,10000.]),np.array([base,base*.5,base*.2]))

    def test_area_at_fractional_phases_and_interpolated_energy(self):
        correction=self.make_correction()
        grid=np.arange(3072)*2.+.25
        for q in (300.,800.,5000.,18000.):
            for t in (128.,128.17,128.5,129.5,1328.):
                value=correction.voltage(grid,t,q)
                self.assertLess(abs(value.sum()*2),1e-8)
                self.assertEqual(float(np.max(abs(value[grid+.25<=t]),initial=0)),0.)

    def test_zero_correction_and_endpoint_clamp(self):
        correction=self.make_correction();correction.scale=0.
        np.testing.assert_array_equal(correction.voltage(np.arange(1000),128.,700.),np.zeros(1000))

    def test_calibration_keeps_means_separate_from_realized_charge(self):
        grid=np.arange(1024)*2.+.25
        s=grid-128.
        shape=np.maximum(s,0)*np.exp(-np.maximum(s,0)/100)*(s>=0)/100
        q=np.array([500.,510.,2000.,2010.]);groups=np.array([0,0,1,1])
        predicted=q[:,None]*shape[None]*.001
        means=predicted+q[:,None]*(.0001*np.sin(s/35)*temporal_window(s))[None]
        correction,records=calibrate_profiles(grid,means,predicted,q,groups)
        self.assertEqual([r['mothers'] for r in records],[2,2])
        np.testing.assert_allclose(correction.profiles.sum(1),0,atol=1e-12)
        self.assertTrue(np.isfinite(correction.profiles).all())

    def test_invalid_nonconserving_profile_rejected(self):
        c=self.make_correction();p=c.profiles.copy();p[0,20]+=1
        with self.assertRaises(ValueError):MeanShapeCorrection(c.grid,c.log_charge_knots,p)


if __name__=='__main__':unittest.main()
