import unittest
import numpy as np
from .alpha_current_pair_energy import selected_blocks,paired_statistics


class PairedEnergyStatistics(unittest.TestCase):
    def test_cluster_cancellation(self):
        old=np.arange(1.,17.).reshape(8,2)
        new=old*.5
        success=np.zeros_like(old)
        index=np.random.default_rng(7).integers(8,size=(100,8))
        result=paired_statistics(old,new,success,success,index)
        self.assertAlmostEqual(result['weak_rmse_ratio'],.5)
        np.testing.assert_allclose(result['weak_rmse_ratio_ci95'],[.5,.5])
        self.assertEqual(result['strict_gain_ci95'],[0.,0.])
        with self.assertRaises(ValueError):paired_statistics(old,new,success,success,index[:,:4])

    def test_cohort_unique_and_stable(self):
        blocks=selected_blocks()
        self.assertEqual(len(blocks),8)
        self.assertEqual(len(set(blocks)),8)
        self.assertTrue(all(0<=b<32 for b in blocks))
        self.assertEqual(blocks,selected_blocks())


if __name__=='__main__':unittest.main()
