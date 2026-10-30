"""Independent algebra, causal and inverse controls for energy calibration."""
import inspect
import unittest

import numpy as np
from scipy.linalg import cholesky_banded

from .alpha_pair_energy_calibrated import (BandedModes, CalibratedEnergyResponse,
    calibrate_mean, compare_known_times, learned_metric)
from .test_alpha_pair_recovery_energy import fixture


def model_fixture():
    base=fixture();grid=np.arange(-16.,1720.01,2.)
    coefficients=np.zeros((len(grid),9))
    coefficients[:,0]=.00012*np.exp(-((grid-100.)/80.)**2)
    coefficients[:,3]=-.00005*np.exp(-np.maximum(grid,0.)/100.)
    coefficients[:,6]=.00002*np.exp(-np.maximum(grid,0.)/100.)
    return base,CalibratedEnergyResponse(base,grid,coefficients)


def noise_fixture():
    return dict(background_scale=1.1,shot_parameters=[-2.,.1,.01],lag1=.4,cross=-.05,
        mode_charge_knots_pe=[600.,3000.,8000.,16000.],fractional_mode_variance=[.0004,.0002,.0001,.0001])


class EnergyCalibrationTests(unittest.TestCase):
    def test_banded_plus_modes_matches_dense(self):
        n=36;band=np.zeros((3,n));band[0]=2.;band[1,:-1]=.25;band[2,:-2]=-.10
        dense=np.diag(band[0])+np.diag(band[1,:-1],-1)+np.diag(band[1,:-1],1)+np.diag(band[2,:-2],-2)+np.diag(band[2,:-2],2)
        rng=np.random.default_rng(823);u=rng.normal(size=(n,2))*.2;r=rng.normal(size=(n,3))
        metric=BandedModes(cholesky_banded(band,lower=True),u)
        covariance=dense+u@u.T;w=metric.apply(r)
        np.testing.assert_allclose(w.T@w,r.T@np.linalg.solve(covariance,r),rtol=1e-12,atol=1e-12)
        np.testing.assert_allclose(metric.marginal_variance(),np.diag(covariance),rtol=1e-12)
        self.assertAlmostEqual(metric.logdet,float(np.linalg.slogdet(covariance)[1]),places=11)

    def test_future_cannot_change_early_mean_and_permutation(self):
        base,model=model_fixture()
        single,_=model.components([128.],[.3]);pair,q=model.components([128.,328.],[.3,8.])
        mask=base.sample_grid<328.
        np.testing.assert_allclose(pair[0,mask],single[0,mask],rtol=0,atol=1e-12)
        reverse,qr=model.components([328.,128.],[8.,.3])
        np.testing.assert_allclose(reverse,pair[::-1],rtol=0,atol=1e-12)
        np.testing.assert_allclose(qr,q[::-1],rtol=0,atol=1e-10)

    def test_mean_regression_predicts_heldout_energy(self):
        base,truth=model_fixture();records=[]
        for e in (.3,1.,2.,4.,6.,8.):
            t=[128.];curves,_=truth.components(t,[e])
            records.append(dict(condition="single",gap_ns=None,times_ns=t,energy_MeV=[e],components_mV=curves))
        calibrated,_=calibrate_mean(base,records,ridge=.001)
        expected=truth.components([141.7],[3.2])[0]
        before=np.linalg.norm(base.components([141.7],[3.2])[0]-expected)
        after=np.linalg.norm(calibrated.components([141.7],[3.2])[0]-expected)
        self.assertLess(after,before*.25)

    def test_energy_covariance_is_positive_and_causal(self):
        base,model=model_fixture();parameters=noise_fixture()
        one=learned_metric(model,[128.],[.3],parameters)
        pair=learned_metric(model,[128.,328.],[.3,8.],parameters)
        grid=base.sample_grid[len(base.ar)-1:];before=grid<328.
        np.testing.assert_allclose(pair.marginal_variance()[before],one.marginal_variance()[before],rtol=1e-11,atol=1e-12)
        self.assertTrue(np.all(pair.marginal_variance()>0))
        self.assertTrue(np.isfinite(pair.logdet))

    def test_noiseless_inverse_with_fixed_known_times(self):
        base,model=model_fixture();t=[128.,168.];e=[.4,6.]
        voltage=model.components(t,e)[0].sum(0)
        result=compare_known_times(base,model,noise_fixture(),voltage,t)
        for name in ("new_mean_old_noise","new_mean_new_noise"):
            np.testing.assert_allclose(result["fits"][name]["energy_MeV"],e,rtol=3e-4,atol=1e-6)
            self.assertEqual(result["fits"][name]["times_ns"],t)
        self.assertEqual(set(inspect.signature(compare_known_times).parameters),
                         {"base","calibrated","noise_parameters","voltage","known_times"})


if __name__=="__main__":unittest.main()
