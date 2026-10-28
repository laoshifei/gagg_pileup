"""Numerical checks for risks, gradients, decoding and model persistence."""
import tempfile
import unittest
from pathlib import Path

import numpy as np

from .alpha_pair_risk_learning import (Adam, RiskNet, OUTPUT_SIZE, TIME_GRID,
    calibrated_gate, charge_risk, decode, objective, select_candidate)


class RiskLearningTests(unittest.TestCase):
    def test_scalar_decisions(self):
        truth=np.array([500.,2000.])
        answers={'absolute':1250.,'logarithmic':1000.,'relative':np.sum(1/truth)/np.sum(1/truth**2)}
        for kind,value in answers.items():
            loss,gradient=charge_risk(np.full(2,value),truth,kind)
            self.assertLess(abs(gradient.mean()),1e-12)
            for offset in (-.1,.1):
                other,_=charge_risk(np.full(2,value*(1+offset)),truth,kind)
                self.assertLess(loss.mean(),other.mean())

    def test_output_gradient(self):
        rng=np.random.default_rng(3);z=rng.normal(0.,.3,(3,OUTPUT_SIZE))
        base=np.array([7000.,15000.,50.]);anchor=np.array([200.,210.,120.])
        truth=np.array([[1000.,6000.],[12000.,3000.],[0.,0.]])
        times=np.array([[195.,240.],[210.,380.],[0.,0.]])
        count=np.array([2,2,0])
        for kind in ('absolute','logarithmic','relative'):
            value,gradient,_=objective(z,base,truth,times,anchor,count,kind)
            self.assertTrue(np.isfinite(value))
            for row,col in ((0,0),(0,1),(0,2),(0,3+190),(1,3+len(TIME_GRID)+270),(2,0),(2,2)):
                step=1e-5;a=z.copy();b=z.copy();a[row,col]+=step;b[row,col]-=step
                numeric=(objective(a,base,truth,times,anchor,count,kind)[0]-objective(b,base,truth,times,anchor,count,kind)[0])/(2*step)
                self.assertAlmostEqual(numeric,gradient[row,col],places=7)

    def test_network_gradient(self):
        model=RiskNet(np.zeros(5),np.ones(5),4,hidden=(4,3))
        x=np.random.default_rng(5).normal(size=(3,5)).astype(np.float32)
        output,saved=model.forward(x,cache=True)
        gradient=np.random.default_rng(6).normal(size=output.shape).astype(np.float32)/100.
        backward=model.backward(gradient,saved)
        for name,index in (('W0',(2,1)),('b0',(2,)),('W1',(1,1)),('W2',(0,3)),('b2',(4,))):
            old=model.parameters[name][index].copy();step=.002
            model.parameters[name][index]=old+step;a=float(np.sum(model.forward(x)*gradient))
            model.parameters[name][index]=old-step;b=float(np.sum(model.forward(x)*gradient))
            model.parameters[name][index]=old
            self.assertAlmostEqual((a-b)/(2*step),backward[name][index],places=5)

    def test_decode_modes_and_charge_conservation(self):
        z=np.zeros((1,OUTPUT_SIZE));z[:,3:]=-20.
        z[0,3+230]=10.;z[0,3+len(TIME_GRID)+180]=10.;z[0,2]=1.
        fit=decode(z,[10000.],[200.])[0]
        self.assertLess(fit['times_ns'][0],fit['times_ns'][1])
        self.assertAlmostEqual(sum(fit['charge_pe']),10000.,places=8)
        self.assertLess(fit['charge_pe'][0],fit['charge_pe'][1])
        self.assertEqual(fit['truth_inputs'],[])

    def test_null_gate_ties_and_invalid_candidate(self):
        candidate=lambda score,q=700.:dict(times_ns=[100.,140.],charge_pe=[5000.,q],network_score=score)
        gate=calibrated_gate([candidate(0.)]*40,[5000.]*40,edges=())
        base=[dict(cost=100.,times_ns=[],charge_pe=[]),dict(cost=50.,times_ns=[100.],charge_pe=[6000.])]
        self.assertFalse(select_candidate(base,candidate(0.),gate)['birth_accepted'])
        self.assertTrue(select_candidate(base,candidate(1.),gate)['birth_accepted'])
        self.assertFalse(select_candidate(base,candidate(100.,499.),gate)['birth_accepted'])
        self.assertFalse(select_candidate(base[:1],candidate(100.),gate)['birth_accepted'])

    def test_persistence_and_adam(self):
        model=RiskNet(np.zeros(5),np.ones(5),7,hidden=(4,3));x=np.ones((3,5),np.float32)
        before=model.forward(x)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'model.npz';model.save(path)
            np.testing.assert_array_equal(RiskNet.load(path).forward(x),before)
        output,saved=model.forward(x,cache=True);gradient=np.ones_like(output)/output.size
        gradients=model.backward(gradient,saved);Adam(model.parameters).update(model.parameters,gradients)
        self.assertLess(model.forward(x).mean(),before.mean())


if __name__=='__main__':unittest.main()
