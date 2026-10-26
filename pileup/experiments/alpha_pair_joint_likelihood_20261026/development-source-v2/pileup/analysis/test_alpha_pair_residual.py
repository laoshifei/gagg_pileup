"""Behavioral checks for covariance regularization and missing-event search."""
from pathlib import Path
import json
import os
import unittest

import numpy as np

from .alpha_pair_conditional import ConditionalResponse
from .alpha_pair_residual import PrecisionBlend,choose,fit_candidates,residual_starts


class ResidualTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root=Path(os.environ.get('GAGG_RESIDUAL_PROJECT',Path(__file__).resolve().parents[2]))
        cls.response=ConditionalResponse.load(root/'pileup/experiments/alpha_pair_conditional_20261022/response.npz')

    def test_precision_is_fixed_linear_positive_and_matrix_consistent(self):
        response=self.response;rng=np.random.default_rng(3217)
        y=response.components([153.7],[2700.]).sum(0)+response.baseline
        metric=PrecisionBlend(response,y,.2,rho=.9)
        a=rng.normal(size=1024);b=rng.normal(size=1024)
        np.testing.assert_allclose(metric.apply(a+2*b),metric.apply(a)+2*metric.apply(b),rtol=1e-12,atol=1e-10)
        np.testing.assert_allclose(metric.apply(np.column_stack([a,b]))[:,0],metric.apply(a),rtol=1e-12,atol=1e-10)
        self.assertGreater(float(metric.apply(a)@metric.apply(a)),0.)
        r=metric.apply(y-response.baseline)
        self.assertAlmostEqual(float(r@r),1e6,places=7)

    def test_residual_search_recovers_unreported_later_pulse(self):
        response=self.response
        expected_t=np.array([151.3,331.7]);expected_q=np.array([15000.,2400.])
        y=response.components(expected_t,expected_q).sum(0)+response.baseline
        one=dict(times_ns=[151.3],charge_pe=[17400.])
        seeds=dict(candidates=[dict(times_ns=[],charge_pe=[]),one,
            dict(times_ns=[151.3,159.3],charge_pe=[9000.,8400.])],ml_1=one,ml_2=one)
        found=fit_candidates(response,y,seeds,fraction=0.,rho=.9)
        result=choose(found['expanded'],.001)
        self.assertEqual(len(result['times_ns']),2)
        np.testing.assert_allclose(result['times_ns'],expected_t,atol=.5)
        np.testing.assert_allclose(result['charge_pe'],expected_q,rtol=.03)
        self.assertTrue(any(abs(s['times_ns'][-1]-expected_t[-1])<8. for s in found['birth_starts']))

    def test_correct_single_and_empty_are_not_split(self):
        response=self.response
        for q in (0.,2700.,16000.):
            y=response.baseline+response.components([151.3],[q]).sum(0)
            raw=response.base.candidates(y)
            initial=dict(times_ns=[151.3],charge_pe=[q]) if q else dict(times_ns=[],charge_pe=[])
            seeds=dict(candidates=raw,ml_1=initial,ml_2=initial)
            result=choose(fit_candidates(response,y,seeds,fraction=0.,rho=.9)['expanded'],.001)
            self.assertEqual(len(result['times_ns']),int(q>0.))


if __name__=='__main__':unittest.main()
