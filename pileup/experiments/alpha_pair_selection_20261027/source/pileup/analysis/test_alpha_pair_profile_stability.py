import unittest
import numpy as np

from .alpha_pair_profile_stability import (
    midpoint_grid, contour_midpoints, compare_regions, refine_stage, diverse_starts,
)
from .test_alpha_pair_evidence_profile import ToyResponse, ToyNoise


class StabilityTests(unittest.TestCase):
    def test_log_grid(self):
        np.testing.assert_allclose(midpoint_grid([50., 200., 3200.]), [100., 800.])
        for q in ([1., 1.], [2., 1.], [-1., 2.]):
            with self.assertRaises(ValueError): midpoint_grid(q)

    def test_contour_selection(self):
        curve = dict(q_grid_pe=[1., 4., 16., 64., 256.], reference_cost=0.,
                     fits=[dict(cost=x) for x in (20., 0., 20., 30., 40.)])
        np.testing.assert_allclose(contour_midpoints(curve), [2., 8.])

    def test_linear_log_contour_convergence(self):
        def curve(q):
            costs=np.abs(np.log(np.array(q)/4.))
            return dict(q_grid_pe=q, fits=[dict(cost=float(x)) for x in costs], reference_cost=0.)
        a=curve([1.,4.,16.]);b=curve([1.,2.,4.,8.,16.])
        result=compare_regions(a,b,.5)
        self.assertTrue(result['stable']);self.assertAlmostEqual(result['endpoint_log_change'],0.)
        b['reference_cost']=-.1
        self.assertFalse(compare_regions(a,b,.5)['stable'])

    def test_topology_change_rejected(self):
        a=dict(q_grid_pe=[1.,2.,4.,8.,16.],reference_cost=0.,fits=[dict(cost=x) for x in (0.,4.,0.,4.,0.)])
        b=dict(a,fits=[dict(cost=0.) for _ in range(5)])
        self.assertFalse(compare_regions(a,b,1.)['stable'])

    def test_diverse_seeds_include_different_times(self):
        response=ToyResponse();noise=ToyNoise()
        seed=dict(times_ns=[20.,60.],charge_pe=[10000.,1000.])
        y=response.components(seed['times_ns'],seed['charge_pe']).sum(0)+10.03
        close=dict(times_ns=[20.001,60.001],charge_pe=[10000.,1000.])
        far=dict(times_ns=[30.,80.],charge_pe=[10000.,1000.])
        result=diverse_starts(response,noise,y,1000.,[seed,close,far],(0.,100.),2)
        self.assertEqual(len(result),2)
        self.assertGreater(abs(result[1]['times_ns'][0]-result[0]['times_ns'][0]),8.)

    def test_correlated_profile_and_old_point_retention(self):
        response=ToyResponse();noise=ToyNoise();seed=dict(times_ns=[20.,60.],charge_pe=[10000.,1000.])
        y=(response.components(seed['times_ns'],seed['charge_pe']).sum(0)+10.03).astype(np.float32)
        old=dict(bounds_ns=[0.,100.],q_grid_pe=[1000.,3000.],reference_cost=0.,
                 fits=[dict(seed,cost=0.),dict(times_ns=[20.,60.],charge_pe=[8000.,3000.],cost=.04)])
        one=dict(times_ns=[20.],charge_pe=[11000.])
        row=dict(sets=dict(baseline=[{},one,seed],mlp=[{},one,seed]))
        result=refine_stage(response,noise,y,row,old,[2000.],refit_existing=False,
                            starts=2,maxiter=100,label='test')
        self.assertEqual(result['q_grid_pe'],[1000.,2000.,3000.])
        self.assertLessEqual(result['fits'][0]['cost'],old['fits'][0]['cost'])
        self.assertAlmostEqual(result['fits'][1]['cost'],.01,delta=2e-5)
        self.assertAlmostEqual(result['fits'][1]['charge_pe'][0],9000.,delta=1.)
        self.assertTrue(all(f['cost']<=s+1e-9 for f,s in zip(result['fits'],result['slice_cost'])))


if __name__=='__main__':unittest.main()
