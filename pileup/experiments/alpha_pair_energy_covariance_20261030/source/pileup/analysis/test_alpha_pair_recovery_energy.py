"""Physical limits and noiseless inverse consistency for the new response."""
import unittest

import numpy as np

from .alpha_pair_recovery_energy import RecoveryEnergyResponse, coupled_charge_inverse, refine_energy


def fixture():
    grid=np.arange(-24.,6000.01,.25)
    s=np.maximum(grid,0.)
    shape=.6*(np.exp(-s/95.)-np.exp(-s/2.))/93.+.4*(np.exp(-s/351.)-np.exp(-s/2.))/349.
    shape[grid<0]=0.
    gain=13.6;kernel=shape*gain
    area=float(np.trapezoid(np.interp(np.arange(-16.,6000.001),grid,kernel),np.arange(-16.,6000.001)))
    ar=np.poly([np.exp(-1.),np.exp(-.1),np.exp(-1.)])
    noise_lags=np.correlate(ar,ar,'full')[len(ar)-1:]*.01
    noise_lags[0]+=.001
    correction_grid=np.arange(-16.,1720.01,2.)
    arrays=dict(base_grid=grid,base_kernel=kernel,base_area=area,correction_grid=correction_grid,
        mean_coefficients=np.zeros((len(correction_grid),9)),gain=gain,baseline_mV=0.,
        background_mV=0.,ar=ar,noise_lags=noise_lags,noise_parameters=[-2.,0.,0.],lag1=.1,cross=0.)
    energy=dict(energy_MeV=[.1,9.],charge_pe=[200.,18000.])
    return RecoveryEnergyResponse(arrays,energy,dict(strength=2.,recovery_ns=29.,cells=14400))


class RecoveryEnergyTests(unittest.TestCase):
    def setUp(self):self.model=fixture()

    def test_isolated_limit(self):
        v,q=self.model.components([128.],[2.],strength=5.)
        v0,q0=self.model.components([128.],[2.],strength=0.)
        np.testing.assert_array_equal(v,v0);np.testing.assert_array_equal(q,q0)

    def test_later_neighbor_cannot_change_earlier_samples(self):
        alone,_=self.model.components([128.],[.3])
        pair,_=self.model.components([128.,328.],[.3,8.])
        keep=self.model.sample_grid<328.
        np.testing.assert_allclose(pair[0,keep],alone[0,keep],rtol=0,atol=1e-12)
        self.assertGreater(float(np.max(np.abs(pair[0,~keep]-alone[0,~keep]))),.01)

    def test_component_permutation(self):
        a,qa=self.model.components([128.,168.],[8.,.3])
        b,qb=self.model.components([168.,128.],[.3,8.])
        np.testing.assert_allclose(a,b[::-1],rtol=0,atol=1e-12)
        np.testing.assert_allclose(qa,qb[::-1],rtol=0,atol=1e-12)

    def test_shared_response_is_not_an_amplitude_rescaling(self):
        single,q0=self.model.components([128.],[.3])
        pair,q=self.model.components([128.,168.],[.3,8.])
        ratio=pair[0,single[0]>.05]/single[0,single[0]>.05]
        self.assertGreater(float(np.ptp(ratio)),.02)
        self.assertLess(q[0],q0[0])

    def test_coupled_inverse_at_known_time(self):
        t=[128.,168.];e=[8.,.3]
        _,q=self.model.components(t,e)
        fit=coupled_charge_inverse(self.model,t,q)
        np.testing.assert_allclose(fit['energy_MeV'],e,rtol=2e-6,atol=1e-7)

    def test_noiseless_waveform_inverse(self):
        t=np.array([128.,208.]);e=np.array([.4,6.])
        components,q=self.model.components(t,e)
        fit=refine_energy(self.model,components.sum(0),dict(times_ns=t+[.4,-.4],charge_pe=q),self.model.strength,maxiter=90)
        np.testing.assert_allclose(fit['times_ns'],t,atol=.02,rtol=0)
        np.testing.assert_allclose(fit['energy_MeV'],e,atol=0,rtol=2e-3)
        self.assertLessEqual(fit['cost'],fit['initial_cost'])

    def test_empty_and_invalid_domain(self):
        c,q=self.model.components([],[])
        self.assertEqual(c.shape,(0,1024));self.assertEqual(q.shape,(0,))
        with self.assertRaises(ValueError):self.model.components([128.],[10.])
        with self.assertRaises(ValueError):self.model.components([-1.],[1.])


if __name__=='__main__':unittest.main()
