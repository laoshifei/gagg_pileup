"""Physical-domain, objective and identifiable-toy checks for charge coordinates."""
import unittest

import numpy as np

from .alpha_pair_charge_allocation import (encode_charges,decode_charges,refine_allocation,
    profile_starts,extend_allocation,forbid_labels,VARIANTS)
from .alpha_pair_joint_likelihood import refine


class ToyResponse:
    baseline=0.
    sample_grid=np.linspace(0.,2048.,4)

    def components(self,t,q):
        # An identifiable linear mean: accurately measured total charge and
        # a less precise partition measurement, plus separately measured times.
        return np.array([[t[0],q[0],0.,q[0]/100.],[0.,q[1],t[1],0.]])


class ToyNoise:
    background_mV=0.

    def transform(self,x):return np.asarray(x)*np.array([1.,.01,1.,1.])

    def score(self,residual,curves,q,*,transformed=False,details=False):
        r=np.asarray(residual) if transformed else self.transform(residual)
        quadratic=float(r@r)
        return dict(cost=quadratic,quadratic=quadratic,logdet=0.) if details else quadratic


class AllocationTests(unittest.TestCase):
    def test_round_trip_entire_physical_domain(self):
        rng=np.random.default_rng(2026102701)
        values=np.vstack((rng.uniform(50.,50000.,(500,2)),[[50,50],[50000,50000],[50,50000],[50000,50],[18000,700],[700,18000]]))
        for q in values:np.testing.assert_allclose(decode_charges(encode_charges(q)),q,rtol=2e-11,atol=1e-8)

    def test_total_and_original_bounds(self):
        for total in (100.,101.,700.,18700.,50050.,75000.,99999.,100000.):
            for fraction in np.linspace(0.,1.,31):
                q=decode_charges([100.*np.log(total),10.*fraction])
                self.assertGreaterEqual(q.min(),50.-1e-8)
                self.assertLessEqual(q.max(),50000.+1e-8)
                self.assertAlmostEqual(q.sum(),total,places=7)

    def test_identifiable_correlated_charge_problem(self):
        response,noise=ToyResponse(),ToyNoise()
        y=response.components([100.,140.],[18000.,700.]).sum(0)
        initial=dict(times_ns=[99.,142.],charge_pe=[15000.,2000.])
        new=refine_allocation(response,noise,y,initial,maxiter=150)
        old=refine(response,noise,y,initial,maxiter=150)
        self.assertLess(new['cost'],1e-5)
        np.testing.assert_allclose(new['times_ns'],[100.,140.],atol=.003)
        np.testing.assert_allclose(new['charge_pe'],[18000.,700.],atol=.3)
        self.assertLess(abs(old['cost']-new['cost']),1e-4)

    def test_likelihood_invariant_to_representing_charge(self):
        response,noise=ToyResponse(),ToyNoise()
        rng=np.random.default_rng(2026102702);y=rng.normal(size=4)
        for q in ([18000.,700.],[50.,50000.],[40000.,30000.]):
            direct=response.components([100.,140.],q)
            changed=response.components([100.,140.],decode_charges(encode_charges(q)))
            np.testing.assert_allclose(direct,changed,rtol=1e-12,atol=1e-9)
            self.assertAlmostEqual(noise.score(y-direct.sum(0),direct,q),noise.score(y-changed.sum(0),changed,q),places=7)

    def test_profile_includes_unprofiled_value(self):
        response,noise=ToyResponse(),ToyNoise()
        y=response.components([100.,140.],[18000.,700.]).sum(0)
        bank=profile_starts(response,noise,y,dict(times_ns=[100.],charge_pe=[18700.]))
        self.assertEqual(len(bank['starts']),8)
        self.assertTrue(all(r['cost']<=r['unprofiled_cost']+1e-9 for r in bank['profiles']))
        for r in bank['profiles']:
            self.assertTrue(np.all(np.diff(r['times_ns'])>0.))
            self.assertGreaterEqual(min(r['charge_pe']),50.-1e-8)
            self.assertLessEqual(max(r['charge_pe']),50000.+1e-8)

    def test_no_signal_candidate_preserved(self):
        old=[dict(times_ns=[],charge_pe=[],cost=0.)]
        result=extend_allocation(None,None,None,old)
        self.assertEqual(set(result['variants']),set(VARIANTS))
        for variant in result['variants'].values():self.assertEqual(variant['candidates'],old)

    def test_truth_access_guard_and_invalid_input(self):
        for p in ('labels.truth.npz','fit.truth.jsonl','labels.cases.json','event.root','data/native-records.json'):
            with self.assertRaises(PermissionError):forbid_labels('open',(p,))
        forbid_labels('open',('observations.npz',))
        with self.assertRaises(ValueError):encode_charges([0.,700.])


if __name__=='__main__':unittest.main()
