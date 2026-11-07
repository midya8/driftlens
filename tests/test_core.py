import copy
import itertools
import math
import random
import unittest
from unittest.mock import patch
from driftlens.core import compare, kernel_matrix, statistic
from driftlens.demo import run_demo


def data():return {'reference': [[0], [1]], 'current': [[2], [3]]}


class TestMMD(unittest.TestCase):
    def test_kernel_closed_form(self):
        k=kernel_matrix([[0],[1]],[[2],[3]],1)
        self.assertEqual(k[0][0],1)
        self.assertAlmostEqual(k[0][1],math.exp(-.5))
        self.assertEqual(k[0][2],k[2][0])

    def test_unbiased_formula_direct(self):
        k=kernel_matrix([[0],[1]],[[2],[3]],1)
        expected=k[0][1]+k[2][3]-sum(k[i][j] for i in [0,1] for j in [2,3])/2
        self.assertAlmostEqual(statistic(k,[0,1],[2,3]),expected)

    def test_negative_estimates_not_clipped(self):
        r=compare({'reference':[[0],[1]],'current':[[0],[1]]},exact=True)
        self.assertLess(r['mmd2_unbiased'],0)

    def test_exact_tiny_tail_known(self):
        r=compare({'reference':[[0],[.1]],'current':[[3],[3.1]]},exact=True)
        self.assertEqual(r['draws'],6)
        self.assertEqual(r['exceedances'],2)
        self.assertEqual(r['p_value'],1/3)

    def test_exact_p_values_rank_superuniform(self):
        rows=[[0],[.1],[3],[3.1]]
        ps=[]
        for a in itertools.combinations(range(4),2):
            b=[i for i in range(4) if i not in a]
            ps.append(compare({'reference':[rows[i] for i in a], 'current':[rows[i] for i in b]},exact=True)['p_value'])
        for alpha in [.1,1/3,.5,.9]:
            self.assertLessEqual(sum(p<=alpha for p in ps)/6,alpha)

    def test_all_identical_ties_p_one(self):
        for exact in [True,False]:
            self.assertEqual(compare({'reference':[[1],[1]],'current':[[1],[1]]},exact=exact)['p_value'],1)

    def test_seed_and_input_rng_unchanged(self):
        d=data(); saved=copy.deepcopy(d); state=random.getstate()
        self.assertEqual(compare(d,seed=7),compare(d,seed=7))
        self.assertEqual(d,saved);self.assertEqual(random.getstate(),state)

    def test_mc_plus_one_and_resolution(self):
        r=compare(data(),permutations=7)
        self.assertEqual(r['p_value'],(r['exceedances']+1)/8)
        self.assertGreaterEqual(r['p_value'],r['p_resolution'])

    def test_unequal_group_sizes_and_swap(self):
        d={'reference':[[0],[1]],'current':[[2],[3],[4]]}
        a=compare(d,exact=True);b=compare({'reference':d['current'],'current':d['reference']},exact=True)
        self.assertAlmostEqual(a['mmd2_unbiased'],b['mmd2_unbiased'])
        self.assertEqual(a['p_value'],b['p_value'])
        self.assertEqual(a['draws'],10)

    def test_row_order_exact_invariance(self):
        a=compare(data(),exact=True);d=data();d['reference'].reverse()
        self.assertEqual(a['p_value'],compare(d,exact=True)['p_value'])

    def test_one_cached_kernel_per_test(self):
        with patch('driftlens.core.kernel_matrix',wraps=kernel_matrix) as cached:
            r=compare(data(),permutations=30)
            cached.assert_called_once()
        self.assertEqual(r['kernel_evaluations'],6)

    def test_bad_rows_and_options(self):
        for value in [[],[[1]],[[1],[True]],[[1],[float('nan')]],[[1],[1,2]]]:
            d=data();d['current']=value
            with self.subTest(value=value),self.assertRaises(ValueError):compare(d)
        for kwargs in [{'bandwidth':0},{'bandwidth':float('inf')},{'permutations':True},
                       {'permutations':0},{'alpha':1},{'seed':True},{'exact':1}]:
            with self.subTest(kwargs=kwargs),self.assertRaises(ValueError):compare(data(),**kwargs)

    def test_work_limits(self):
        for d,kwargs in [({'reference':[[1]]*251,'current':[[2]]*250},{}),
                         ({'reference':[[1]]*10,'current':[[2]]*10},{'exact':True}),
                         ({'reference':[[1]]*100,'current':[[2]]*100},{'permutations':1000})]:
            with self.assertRaises(ValueError):compare(d,**kwargs)

    def test_extreme_finite_coordinates(self):
        k=kernel_matrix([[1e308],[-1e308]],[[0],[0]],1e308)
        self.assertAlmostEqual(k[0][1],math.exp(-2))
        k=kernel_matrix([[1e154]*10,[-1e154]*10],[[0]*10,[0]*10],1)
        self.assertEqual(k[0][1],0)

    def test_demo_pinned(self):
        r=run_demo()
        self.assertEqual(r['same_generator']['p_value'],.55)
        self.assertEqual(r['mean_shift']['p_value'],.002)
        self.assertFalse(r['same_generator']['reject_at_alpha'])
        self.assertTrue(r['mean_shift']['reject_at_alpha'])
