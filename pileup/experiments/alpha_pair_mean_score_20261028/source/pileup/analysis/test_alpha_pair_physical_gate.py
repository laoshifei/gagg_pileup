"""Controls for search-bank ablation and physical count selection."""
from copy import deepcopy
import unittest

from .alpha_pair_physical_gate import assemble,decide,rejection


def candidate(q,cost,t=None):
    return dict(charge_pe=list(q),times_ns=list(t or [100.+20*i for i in range(len(q))]),cost=float(cost),converged=True)


class PhysicalGateTests(unittest.TestCase):
    def setUp(self):
        self.base=[candidate([],100.),candidate([6000.],50.),candidate([5000.,1000.],20.)]
        self.hybrid=dict(candidate([5100.,900.],10.),penalty=0.,network_score=-1.,
                         trials=[candidate([5050.,950.],15.),candidate([5100.,900.],10.)])
        sets,learned=assemble(self.base,self.hybrid);self.row=dict(sets=sets,learned=learned)
        self.physical=dict(edges_pe=[],alpha=.025,scores=[[5.]*40])
        self.learned=dict(edges_pe=[],alpha=.025,scores=[[0.]*40])

    def test_local_control_and_retained_baseline(self):
        self.assertEqual(self.row['sets']['local'][2]['cost'],15.)
        self.assertEqual(self.row['sets']['mlp'][2]['cost'],10.)
        for group in self.row['sets'].values():self.assertEqual(group[:2],self.base[:2])
        self.assertNotIn('network_score',self.row['sets']['mlp'][2])
        self.row['sets']['local'][0]['cost']=999.
        self.assertEqual(self.base[0]['cost'],100.)

    def test_same_candidate_different_decision(self):
        physical=decide(self.row,'mlp_physical',self.physical)
        learned=decide(self.row,'learned_gate',self.learned)
        self.assertEqual(len(physical['times_ns']),2)
        self.assertEqual(len(learned['times_ns']),1)
        self.assertEqual(rejection(self.row,'learned_gate',learned),'rank_rejected')

    def test_physical_gate_ignores_network_logit(self):
        before=decide(self.row,'mlp_physical',self.physical)
        self.row['learned']['network_score']=1e9
        self.assertEqual(before,decide(self.row,'mlp_physical',self.physical))

    def test_subthreshold_minimum_is_not_silently_replaced(self):
        row=deepcopy(self.row);row['sets']['mlp'][2]['charge_pe']=[5501.,499.]
        fit=decide(row,'mlp_physical',self.physical)
        self.assertEqual(len(fit['times_ns']),1)
        self.assertEqual(rejection(row,'mlp_physical',fit),'subthreshold_Q')

    def test_presence_and_tie_guards(self):
        row=deepcopy(self.row);row['sets']['baseline'][1]['cost']=100.
        row['sets']['mlp'][1]['cost']=100.
        fit=decide(row,'mlp_physical',self.physical)
        self.assertEqual(len(fit['times_ns']),0)
        self.assertEqual(rejection(row,'mlp_physical',fit),'no_presence')
        gate=dict(self.physical,scores=[[40.]*40])
        fit=decide(self.row,'mlp_physical',gate)
        self.assertEqual(len(fit['times_ns']),1)
        self.assertEqual(rejection(self.row,'mlp_physical',fit),'rank_rejected')

    def test_empty_candidate_bank(self):
        sets,learned=assemble(self.base[:1],dict(candidate([300.,300.],10.),network_score=10.))
        row=dict(sets=sets,learned=learned)
        for method in ('baseline','learned_gate','local_physical','mlp_physical','mlp_fixed_gate'):
            self.assertEqual(len(decide(row,method,self.physical)['times_ns']),0)


if __name__=='__main__':unittest.main()
