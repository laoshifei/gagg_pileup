import unittest
import numpy as np
from .alpha_current_shared import make_plans, anchored_pair_predictions


class SharedPlans(unittest.TestCase):
    def test_independent_blocks_and_phase_membership(self):
        native,single,pairs=make_plans()
        self.assertEqual([len(native),len(single),len(pairs)],[136,544,640])
        self.assertEqual(len({r['seed'] for r in native}),136)
        self.assertEqual(sum(r['role']=='formal' for r in single),512)
        self.assertEqual(sum(r['phase_ns']==0 for r in pairs),160)
        blocks={}
        for r in pairs:
            self.assertTrue(all(native[i]['role']=='formal' for i in r['parents']))
            self.assertEqual({native[i]['group'] for i in r['parents']},{0,3})
            self.assertTrue(all(native[i]['block']==r['block'] for i in r['parents']))
            blocks.setdefault(r['block'],set()).update(r['parents'])
            self.assertEqual(r['times_ns'][1]-r['times_ns'][0],r['gap_ns'])
        self.assertEqual(len(set.union(*blocks.values())),64)
        self.assertTrue(all(len(v)==2 for v in blocks.values()))

    def test_reuse_rejects_charge_change(self):
        class Model:
            def __init__(self,delta=0.,charge=1.):self.delta,self.charge=delta,charge
            def isolated(self,t,e):return np.ones((2,3))+self.delta,np.full(2,self.charge)
            def components(self,t,e):return self.isolated(t,e)[0]-.1
        old,new=Model(),Model(.2)
        out=anchored_pair_predictions(old,new,[1,2],[3,4])
        np.testing.assert_allclose(out['corrected_joint'],new.components([1,2],[3,4]),rtol=0,atol=1e-15)
        with self.assertRaises(ValueError):anchored_pair_predictions(old,Model(charge=2.),[1,2],[3,4])


if __name__=='__main__':unittest.main()
