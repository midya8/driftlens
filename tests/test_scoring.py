import itertools
import random
import unittest
from driftlens.core import kernel_matrix, statistic, compare
from driftlens.scoring import CachedScorer


class TestScoring(unittest.TestCase):
    def test_all_imbalanced_allocations_match_direct(self):
        rng=random.Random(7)
        rows=[[rng.gauss(0,1) for _ in range(3)] for _ in range(9)]
        k=kernel_matrix(rows[:2],rows[2:],1.)
        scorer=CachedScorer(k,2)
        self.assertTrue(scorer.use_complement)
        for left in itertools.combinations(range(9),2):
            left=list(left);right=[i for i in range(9) if i not in left]
            self.assertAlmostEqual(scorer.score(left,right),statistic(k,left,right),places=13)

    def test_complement_group_swap(self):
        rows=[[i] for i in range(9)];k=kernel_matrix(rows[:7],rows[7:],2.)
        scorer=CachedScorer(k,7)
        left=list(range(7));right=[7,8]
        self.assertAlmostEqual(scorer.score(left,right),statistic(k,left,right),places=13)
        self.assertEqual(scorer.pair_reads,1)

    def test_identical_kernel_cancellation_and_ties(self):
        d={'reference':[[0]]*2,'current':[[0]]*7}
        r=compare(d,exact=True)
        self.assertEqual(r['p_value'],1.)
        self.assertEqual(r['mmd2_unbiased'],0.)
        self.assertEqual(r['scoring_mode'],'complement_sums')

    def test_balanced_keeps_direct(self):
        k=kernel_matrix([[0],[1]],[[2],[3]],1.)
        scorer=CachedScorer(k,2)
        self.assertFalse(scorer.use_complement)
        self.assertEqual(scorer.score([0,1],[2,3]),statistic(k,[0,1],[2,3]))

    def test_exact_tail_matches_direct_at_same_tolerance(self):
        rows=[[i*.3] for i in range(9)];d={'reference':rows[:2],'current':rows[2:]}
        k=kernel_matrix(rows[:2],rows[2:],1.)
        obs=statistic(k,[0,1],list(range(2,9)));count=0
        for chosen in itertools.combinations(range(9),2):
            left=list(chosen);right=[i for i in range(9) if i not in left]
            count+=statistic(k,left,right)>=obs-1e-12
        r=compare(d,exact=True)
        self.assertEqual(r['p_value'],count/36)
        self.assertEqual(r['kernel_pair_reads_per_allocation'],1)
        self.assertEqual(r['scoring_setup_off_diagonal_reads'],72)

    def test_near_constant_and_near_zero_kernel(self):
        rows=[[i*.001] for i in range(10)]
        for bandwidth in [1e6,.000001]:
            k=kernel_matrix(rows[:2],rows[2:],bandwidth);s=CachedScorer(k,2)
            for chosen in itertools.combinations(range(10),2):
                left=list(chosen);right=[i for i in range(10) if i not in left]
                self.assertAlmostEqual(s.score(left,right),statistic(k,left,right),places=12)
