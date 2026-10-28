import unittest
import json
import numpy as np
from .alpha_pair_profile_contours import crossing_cells, refinement_targets, quality, converge_contours
from .alpha_pair_evidence_profile import compatibility_segments


def quadratic(q):
    q=sorted(set(map(float,q)))
    return dict(q_grid_pe=q,reference_cost=0.,fits=[dict(cost=float(np.log(v/1000.)**2)) for v in q])


class ContourTests(unittest.TestCase):
    def test_quality_can_be_saved_as_json(self):
        result=quality(quadratic([50.,1000.,50000.]))
        self.assertEqual(json.loads(json.dumps(result)),result)

    def test_brackets_keep_disjoint_components_and_domain_edges(self):
        c=dict(q_grid_pe=[50.,100.,200.,400.,800.],reference_cost=0.,fits=[dict(cost=x) for x in (0.,2.,0.,2.,0.)])
        cells=crossing_cells(c,[1.]);self.assertEqual(len(cells),4)
        r=compatibility_segments(c['q_grid_pe'],[f['cost'] for f in c['fits']],0.,1.)
        self.assertEqual(len(r['segments']),3);self.assertTrue(r['lower_boundary'] and r['upper_boundary'])
        self.assertFalse(quality(c,[1.])['bracket_pass'])

    def test_no_false_precision_from_stationary_interpolated_endpoint(self):
        c=dict(q_grid_pe=[50.,1000.,50000.],reference_cost=0.,fits=[dict(cost=x) for x in (10.,0.,10.)])
        self.assertFalse(quality(c,[1.,4.,9.])['bracket_pass'])
        self.assertEqual(len(refinement_targets(c)),2)

    def test_analytic_curved_profile_endpoints(self):
        c=quadratic([50.,1000.,50000.])
        for _ in range(12):
            q=refinement_targets(c,[1.,4.],.01)
            if not len(q):break
            c=quadratic(c['q_grid_pe']+q.tolist())
        self.assertTrue(quality(c,[1.,4.],.01)['bracket_pass'])
        for level in (1.,4.):
            r=compatibility_segments(c['q_grid_pe'],[f['cost'] for f in c['fits']],0.,level)
            np.testing.assert_allclose(r['segments'][0],1000.*np.exp([-np.sqrt(level),np.sqrt(level)]),rtol=2e-5)

    def test_budget_exhaustion_does_not_claim_convergence(self):
        def exact(*args,**kwargs):
            previous,targets=args[4:6]
            return quadratic(previous['q_grid_pe']+list(targets))
        settings=dict(levels=[1.,4.],bracket_fraction=.00001,reference_tolerance=.01,endpoint_fraction=.01,
            max_rounds=1,local={},grid={},verify={})
        r=converge_contours(None,None,None,None,quadratic([50.,1000.,50000.]),settings,refiner=exact)
        self.assertFalse(r['stable']);self.assertEqual(r['verification_passes'],0)


if __name__=='__main__':unittest.main()
