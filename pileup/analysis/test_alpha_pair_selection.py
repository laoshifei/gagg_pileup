"""Guard the scientific denominators and tied-score decision edge cases."""
import unittest
import numpy as np

from .alpha_pair_selection_study import inclusion_partition, design
from .alpha_pair_joint_likelihood import select, calibrate_gate, evidence


def candidates(gain=10., q2=800., total=18800.):
    return [dict(cost=-10.,charge_pe=[],times_ns=[]),
            dict(cost=-100.,charge_pe=[total],times_ns=[200.]),
            dict(cost=-100.-gain,charge_pe=[total-q2,q2],times_ns=[200.,240.])]


class SelectionTest(unittest.TestCase):
    def test_conditional_and_joint_denominators(self):
        r=inclusion_partition([True,False,True,True],[True,True,False,False])
        self.assertEqual(r['all']['fraction'],.75)
        self.assertEqual(r['selected']['fraction'],.5)
        self.assertEqual(r['selected_and_included']['fraction'],.25)
        self.assertEqual(r['unselected']['fraction'],1.)

    def test_empty_selection_undefined(self):
        r=inclusion_partition([True,False],[False,False])
        self.assertIsNone(r['selected']['fraction'])
        self.assertIsNone(r['selected']['binomial95'])
        self.assertEqual(r['selected_and_included']['fraction'],0.)
        with self.assertRaises(ValueError):inclusion_partition([True],[True,False])

    def test_negative_costs_and_ties(self):
        c=candidates();self.assertEqual(evidence(c),10.)
        gate=dict(edges_pe=[],scores=[[0.]*38+[10.]],alpha=.025)
        self.assertFalse(select(c,gate)['birth_accepted'])
        self.assertEqual(select(c,gate)['null_rank_p'],.05)
        self.assertTrue(select(candidates(10.01),gate)['birth_accepted'])

    def test_presence_and_weak_charge_guards(self):
        gate=dict(edges_pe=[],scores=[[0.]*63],alpha=.1)
        self.assertEqual(len(select(candidates(q2=499.),gate)['charge_pe']),1)
        c=candidates();c[1]['cost']=0.
        self.assertEqual(len(select(c,gate)['charge_pe']),0)
        self.assertEqual(len(select(candidates(gain=0.),gate)['charge_pe']),1)

    def test_fitted_bin_and_alpha_monotonicity(self):
        c=candidates();gate=dict(edges_pe=[1500.,5000.,12500.],scores=[[0.]*64,[0.]*64,[0.]*64,list(range(64))])
        choices=[select(c,dict(gate,alpha=a)) for a in (.025,.05,.1,.9)]
        self.assertTrue(all(r['calibration_bin']==3 for r in choices))
        self.assertEqual([r['birth_accepted'] for r in choices],sorted(r['birth_accepted'] for r in choices))
        with self.assertRaises(ValueError):calibrate_gate([c])

    def test_disjoint_frozen_design(self):
        a,b=design('calibration'),design('test')
        self.assertEqual((len(a),len(b)),(256,448))
        self.assertFalse({tuple(r['seed_words']) for r in a}&{tuple(r['seed_words']) for r in b})
        self.assertEqual(len({tuple(r['seed_words']) for r in a+b}),704)
        self.assertTrue(all(r['count']==1 for r in a))
        self.assertEqual(sum(r['count']==2 for r in b),128)


if __name__=='__main__':unittest.main()
