import unittest
import numpy as np

from .alpha_mean_shape_correction import MeanShapeCorrection
from .alpha_mean_shape_stability import (mother_folds,bootstrap_mothers,many_voltage,
    cached_positivity_scale,fit_update,projected_variance)


class StabilityControls(unittest.TestCase):
    def test_group_identity_and_bootstrap_isolation(self):
        groups=np.repeat(np.arange(4),32);ids=[f'root{i}' for i in range(128)]
        folds=mother_folds(ids,groups,737000)
        for group in range(4):np.testing.assert_array_equal(np.bincount(folds[groups==group]),[8]*4)
        for seed in range(20):
            sample=bootstrap_mothers(groups,folds,seed)
            for fold in range(4):
                self.assertFalse(set(sample[folds[sample]==fold])&set(sample[folds[sample]!=fold]))
        with self.assertRaises(ValueError):mother_folds(['same']*128,groups,0)

    def test_vectorized_voltage_and_positivity(self):
        grid=np.arange(-16.,1602.,2.)+.25
        profiles=np.zeros((4,len(grid)));profiles[:,100]=[.1,.2,.3,.4];profiles[:,130]=[-.1,-.2,-.3,-.4]
        c=MeanShapeCorrection(grid,np.log([100.,500.,2000.,20000.]),profiles)
        time=np.arange(1024)*2.+.25;charges=np.array([20.,100.,800.,20000.,40000.]);phases=np.array([0.,.5,1.,1.5])
        for phase in phases:
            expected=np.array([c.voltage(time,128.+phase,q) for q in charges])
            np.testing.assert_allclose(many_voltage(c,charges,time,128.+phase),expected,rtol=1e-14,atol=1e-12)
        old=np.ones((4,len(charges),len(time)))*100
        cache=dict(q0=charges,time_ns=time,phases_ns=phases,old_mV=old)
        expected=min(1.,100./(40000*.4))*.99
        self.assertAlmostEqual(cached_positivity_scale(c,cache),expected,places=14)

    def test_heldout_signal_is_never_calibrated(self):
        t=np.arange(1024)*2.+.25;groups=np.repeat(np.arange(4),8)
        q=np.repeat([100.,500.,2000.,20000.],8)
        shape=np.exp(-np.maximum(t-128.,0)/300.)*(t>128.)
        predicted=q[:,None]*shape*.02
        data=dict(time_ns=t,groups=groups,q0=q,old_predictions_mV=predicted,
                  means_mV=predicted*(1+.001*np.sin(t/40.)))
        train=np.flatnonzero(np.arange(32)%4!=0)
        a=fit_update(data,train)
        data['means_mV'][np.arange(32)%4==0]=np.nan
        b=fit_update(data,train)
        np.testing.assert_array_equal(a.profiles,b.profiles)
        np.testing.assert_array_equal(a.log_charge_knots,b.log_charge_knots)

    def test_repeat_correction_keeps_negative_between(self):
        repeats=np.tile(np.array([-1.,1.])[None,:,None],(4,1,3))
        operator=np.eye(3)*.01
        rows=projected_variance(repeats,np.zeros((4,3)),operator)
        for row in rows:
            self.assertAlmostEqual(row['within_variance'],2.)
            self.assertAlmostEqual(row['between_variance_unclipped'],-1.)


if __name__=='__main__':unittest.main()
