"""Independent dense checks of normalized, conditioned joint likelihood."""
import unittest

import numpy as np
from scipy.linalg import toeplitz

from .alpha_pair_joint_likelihood import JointNoise,calibrate_gate,select


class JointLikelihoodTests(unittest.TestCase):
    def test_dense_joint_covariance_and_normalization(self):
        rng=np.random.default_rng(1026);n=19
        for count in (0,1,2):
            for cross in (-.2,.1):
                for kappa in (0.,.95):
                    noise=JointNoise(np.array([1.]),np.array([.3,-.08,.01,.0]),np.array([-.7,.08,-.1]),.39,cross,kappa)
                    curves=rng.uniform(.01,3.,(count,n));q=np.array([800.,18000.])[:count];r=rng.normal(size=n)
                    band,V=noise.covariance(curves,q,transformed=True)
                    matrix=sum(np.diag(band[k,:n-k],-k)+(np.diag(band[k,:n-k],k) if k else 0) for k in range(4))-V@V.T
                    dense=toeplitz(np.r_[noise.noise_lags,np.zeros(n-4)])
                    if count:
                        activity=curves;neighbor=curves[::-1] if count==2 else np.zeros_like(curves)
                        v=np.exp(-.7+.08*np.log(q[:,None]/3000.)-.1*neighbor/10.)*activity
                        R=toeplitz(np.r_[1.,.39,np.zeros(n-2)])
                        blocks=[[np.sqrt(v[i,:,None]*v[j,None,:])*R*(1. if i==j else cross)
                                 for j in range(count)] for i in range(count)]
                        B=np.block(blocks);P=np.tile(np.eye(n),(1,count));L=np.kron(np.eye(count),np.ones((1,n)))
                        U=P@B@L.T;S=L@B@L.T
                        dense+=P@B@P.T-kappa*U@np.linalg.solve(S,U.T)
                    np.testing.assert_allclose(matrix,dense,atol=2e-12)
                    self.assertGreater(np.linalg.eigvalsh(matrix)[0],0.)
                    score=noise.score(r,curves,q,transformed=True,details=True)
                    expected=float(r@np.linalg.solve(dense,r))+np.linalg.slogdet(dense)[1]
                    self.assertAlmostEqual(score['cost'],expected,places=9)
                    self.assertAlmostEqual(score['cost'],score['quadratic']+score['logdet'],places=12)

    def test_logdet_penalizes_covariance_inflation(self):
        small=JointNoise(np.array([1.]),np.array([1.,0.]),np.array([0.,0.,0.]),0.,0.,0.)
        large=JointNoise(np.array([1.]),np.array([100.,0.]),np.array([0.,0.,0.]),0.,0.,0.)
        r=np.zeros(12);curves=np.empty((0,12))
        self.assertLess(small.score(r,curves,[],transformed=True),large.score(r,curves,[],transformed=True))

    def test_rank_works_with_negative_nll_and_ties(self):
        def candidate(gain):
            return [dict(times_ns=[],charge_pe=[],cost=10.),dict(times_ns=[100.],charge_pe=[1000.],cost=-100.),
                    dict(times_ns=[100.,140.],charge_pe=[700.,700.],cost=-100.-gain)]
        gate=calibrate_gate([candidate(10.) for _ in range(40)],edges=(),alpha=.025)
        self.assertFalse(select(candidate(10.),gate)['birth_accepted'])
        self.assertTrue(select(candidate(11.),gate)['birth_accepted'])
        bad=candidate(100.);bad[2]['charge_pe'][1]=499.
        self.assertFalse(select(bad,gate)['birth_accepted'])


if __name__=='__main__':unittest.main()
