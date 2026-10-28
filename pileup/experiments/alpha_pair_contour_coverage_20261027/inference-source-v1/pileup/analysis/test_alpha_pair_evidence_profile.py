import unittest
import numpy as np
from .alpha_pair_evidence_profile import encode,decode,fit_at_charge,compatibility_segments,rank_grid,choose_endpoint


class ToyResponse:
    baseline=10.
    def components(self,t,q):
        return np.array([[q[0]/1000.,0.,t[0]/10.,0.],[q[1]/1000.,q[1]/10000.,0.,t[1]/10.]])


class ToyNoise:
    background_mV=.03
    def transform(self,y):return np.asarray(y,float)
    def score(self,r,m,q,**kwargs):return float(np.dot(r,r))


class ProfileTests(unittest.TestCase):
    def test_coordinates(self):
        for t in ([0.,0.],[12.,40.],[90.,100.],[100.,100.]):
            got,q=decode(encode(t,900.,(0.,100.)),750.,(0.,100.))
            np.testing.assert_allclose(got,t,atol=1e-10);self.assertAlmostEqual(q[0],900.)

    def test_profile_nuisance(self):
        response=ToyResponse();noise=ToyNoise();seed=dict(times_ns=[20.,60.],charge_pe=[10000.,1000.])
        y=response.components(seed['times_ns'],seed['charge_pe']).sum(0)+10.03
        fit=fit_at_charge(response,noise,y,2000.,[seed],(0.,100.),maxiter=100,starts=1)
        self.assertAlmostEqual(fit['charge_pe'][0],9000.,delta=.5)
        self.assertAlmostEqual(fit['cost'],.01,delta=1e-5)
        np.testing.assert_allclose(fit['times_ns'],[20.,60.],atol=.01)

    def test_float32_offset(self):
        response=ToyResponse();noise=ToyNoise();seed=dict(times_ns=[20.,60.],charge_pe=[10000.,1000.])
        y=(response.components(seed['times_ns'],seed['charge_pe']).sum(0)+10.03).astype(np.float32)
        a=fit_at_charge(response,noise,y,1000.,[seed],(0.,100.),starts=1)
        b=fit_at_charge(response,noise,y.astype(float),1000.,[seed],(0.,100.),starts=1)
        self.assertEqual(a['cost'],b['cost']);self.assertEqual(a['times_ns'],b['times_ns'])

    def test_disconnected_support(self):
        out=compatibility_segments([1.,2.,4.,8.,16.],[0.,4.,0.,4.,0.],0.,1.)
        self.assertEqual(len(out['segments']),3)
        self.assertTrue(out['lower_boundary'] and out['upper_boundary'])
        self.assertAlmostEqual(out['segments'][0][1],2**.25)

    def test_rank_ties_and_endpoints(self):
        grid=rank_grid({'a':{'scores':[[0.]*40,[0.]*80]}})
        self.assertIn(1/41,grid);self.assertIn(1/81,grid);self.assertIn(.025,grid)
        end=dict(eligible=True,rank=2/81)
        self.assertEqual(choose_endpoint(end,2/81),'high')
        self.assertEqual(choose_endpoint(end,np.nextafter(2/81,0.)),'low')
        self.assertEqual(choose_endpoint(dict(end,eligible=False),1.),'low')

    def test_invalid_domain(self):
        with self.assertRaises(ValueError):fit_at_charge(None,None,None,0.,[],[0.,100.])
        with self.assertRaises(ValueError):compatibility_segments([2.,1.],[0.,1.],0.,1.)


if __name__=='__main__':unittest.main()
