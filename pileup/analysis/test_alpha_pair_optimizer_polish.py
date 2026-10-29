import unittest
import numpy as np

from .alpha_pair_optimizer_polish import DEFAULTS,polish_objective,seed_bank,waveform_polish


class PolishingTests(unittest.TestCase):
    def setUp(self):
        self.one=dict(times_ns=[200.],charge_pe=[18800.])
        self.start=dict(times_ns=[199.,420.],charge_pe=[17000.,1400.],converged=True)
        self.domain=dict(time=[100.,712.],charge=[50.,50000.])
        self.settings=dict(DEFAULTS,starts=2,maxiter=100,powell_maxfev=800)

    def test_known_quadratic_optimum_and_monotonic_stages(self):
        def loss(t,q):
            return np.sum((np.asarray(t)-[200.3,400.3])**2)/100.+10*np.sum((np.log(q)-np.log([18000.,800.]))**2)
        r=polish_objective(loss,self.one,self.start,self.domain,self.settings)
        costs=[r['original_cost']]+[x['cost'] for x in r['stages'].values()]
        self.assertTrue(all(b<=a+1e-8 for a,b in zip(costs,costs[1:])))
        self.assertLess(r['stages']['polished']['cost'],1e-6)
        self.assertEqual(self.start['times_ns'],[199.,420.])
        self.assertEqual(r['truth_inputs'],[])

    def test_single_and_empty_decisions_are_preserved_without_model_access(self):
        for s in (dict(times_ns=[],charge_pe=[]),self.one):
            r=waveform_polish(None,None,None,None,s)
            for x in r['stages'].values(): self.assertEqual(x,s)
            self.assertFalse(r['searched'])

    def test_bad_reference_is_not_silently_optimized(self):
        with self.assertRaises(ValueError):
            polish_objective(lambda t,q:np.nan,self.one,self.start,self.domain,self.settings)

    def test_seed_time_order_keeps_charge_identity(self):
        start=dict(times_ns=[200.,205.],charge_pe=[18000.,800.])
        seeds=seed_bank(self.one,start,self.domain)
        swapped=[s for s in seeds if s['times_ns']==[205.,232.]]
        self.assertEqual(swapped[0]['charge_pe'],[800.,18000.])
        self.assertTrue(all(s['times_ns'][0]<=s['times_ns'][1] for s in seeds))

    def test_flat_objective_retains_original(self):
        r=polish_objective(lambda t,q:3.,self.one,self.start,self.domain,
                           dict(self.settings,starts=1,maxiter=2,powell_maxfev=20))
        for s in r['stages'].values():
            np.testing.assert_allclose(s['times_ns'],self.start['times_ns'])
            np.testing.assert_allclose(s['charge_pe'],self.start['charge_pe'])

    def test_boundary_roundoff_does_not_leave_physical_domain(self):
        def loss(t,q):
            return -sum(np.log(q))+sum(np.asarray(t)/1000.)
        r=polish_objective(loss,self.one,self.start,self.domain,self.settings)
        for s in r['stages'].values():
            self.assertTrue(all(50.<=q<=50000. for q in s['charge_pe']))
            self.assertTrue(100.<=s['times_ns'][0]<=s['times_ns'][1]<=712.)


if __name__=='__main__': unittest.main()
