"""Analytic controls for coordinates, fixed weights and mean-score semantics."""
import unittest
import numpy as np

from .alpha_pair_mean_score import encode,decode,parameter_bounds,optimize,damped_step
from .alpha_pair_model_sampler import GaussianFactor


class MeanScoreTests(unittest.TestCase):
    def test_ordered_roundtrip_and_bounds(self):
        for t,q in (([210.,150.],[800.,18000.]),([200.,200.],[50.,50000.]),([700.,700.],[700.,700.])):
            z=encode(t,q,[100.,700.]);tt,qq=decode(z,2,[100.,700.]);order=np.argsort(t)
            np.testing.assert_allclose(tt,np.asarray(t)[order],atol=1e-9)
            np.testing.assert_allclose(qq,np.asarray(q)[order],rtol=1e-12)
            a,b=np.asarray(parameter_bounds(2,[100.,700.])).T
            self.assertTrue(np.all(z>=a) and np.all(z<=b))

    def test_single_roundtrip(self):
        z=encode([200.3],[800.],[0.,2046.]);t,q=decode(z,1,[0.,2046.])
        np.testing.assert_allclose(t,[200.3]);np.testing.assert_allclose(q,[800.])
        with self.assertRaises(ValueError):encode([],[],[0.,2046.])

    def test_frozen_quadratic_matches_dense(self):
        n=12;band=np.zeros((2,n));band[0]=2.+np.arange(n)/20;band[1,:-1]=.3
        V=np.arange(n)[:,None]/200.;C=np.diag(band[0])+np.diag(band[1,:-1],1)+np.diag(band[1,:-1],-1)-V@V.T
        metric=GaussianFactor.from_covariance(band,V);r=np.sin(np.arange(n))
        self.assertAlmostEqual(float(metric.whiten(r)@metric.whiten(r)),float(r@np.linalg.solve(C,r)),places=11)

    def test_constant_covariance_same_minimum(self):
        y=np.array([1.,2.,4.,5.]);weights=np.array([1.,2.,3.,4.]);truth=float(weights@y/weights.sum())
        f=lambda z:float(np.sum(weights*(y-z[0])**2))
        a=optimize(f,[1.],[(.1,10.)],100);b=optimize(lambda z:f(z)-18.3,[1.],[(.1,10.)],100)
        self.assertAlmostEqual(a['z'][0],truth,places=4);self.assertAlmostEqual(b['z'][0],truth,places=4)

    def test_heteroscedastic_three_distinct_estimators(self):
        y=np.array([1.,2.,2.,5.]);m=y.mean();s=np.mean(y*y)
        full=lambda z:float(np.sum((y-z[0])**2/z[0]**2)+len(y)*np.log(z[0]**2))
        frozen=lambda z:float(np.sum((y-z[0])**2)/4.)
        omit=lambda z:float(np.sum((y-z[0])**2/z[0]**2))
        expected=[(-m+np.sqrt(m*m+4*s))/2,m,s/m]
        for f,answer in zip((full,frozen,omit),expected):
            result=optimize(f,[2.],[(.1,10.)],200)
            self.assertAlmostEqual(result['z'][0],answer,places=4)
        self.assertGreater(min(np.diff(expected)),.1)

    def test_damping_cannot_fake_success(self):
        f=lambda z:float((z[0]-1.)**2)
        z,meta=damped_step(f,np.array([3.]),np.array([1.]),.5)
        self.assertTrue(meta['accepted']);self.assertEqual(z[0],2.)
        z,meta=damped_step(f,np.array([1.]),np.array([10.]),.5)
        self.assertFalse(meta['accepted']);self.assertEqual(z[0],1.)


if __name__=='__main__':unittest.main()
