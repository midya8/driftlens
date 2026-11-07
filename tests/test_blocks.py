import copy
import itertools
import unittest
from driftlens import compare, compare_family
from driftlens.allocations import partition, allocation_count, allocations
import random


def data():
    return {'reference':[[0],[5]],'current':[[.1],[5.1]],
            'blocks':{'reference':['a','b'],'current':['a','b']}}


class TestBlocks(unittest.TestCase):
    def test_exact_product_count_and_counts_preserved(self):
        groups, blocked = partition(data());self.assertTrue(blocked)
        self.assertEqual(allocation_count(groups),4)
        for exact in [True,False]:
            for chosen in allocations(groups,exact,20,random.Random(7)):
                for indices,n in groups:self.assertEqual(len(set(indices)&set(chosen)),n)
        self.assertEqual(compare(data(),exact=True)['draws'],4)

    def test_single_block_equals_pooled(self):
        d=data();d['blocks']={'reference':['one']*2,'current':['one']*2}
        a=compare(d,exact=True);del d['blocks'];b=compare(d,exact=True)
        self.assertEqual(a['p_value'],b['p_value'])
        self.assertEqual(a['mmd2_unbiased'],b['mmd2_unbiased'])

    def test_frozen_groups_are_reported_not_silently_dropped(self):
        d=data();d['blocks']={'reference':['a','a'],'current':['b','b']}
        r=compare(d,exact=True)
        self.assertEqual(r['total_label_allocations'],1)
        self.assertEqual(r['informative_blocks'],0)
        self.assertEqual(r['p_value'],1)

    def test_invalid_block_schema(self):
        for b in [{}, {'reference':['a'], 'current':['a','b']},
                  {'reference':[True,'b'],'current':['a','b']},
                  {'reference':['','b'],'current':['a','b']}]:
            d=data();d['blocks']=b
            with self.subTest(blocks=b),self.assertRaises(ValueError):compare(d)

    def test_seed_and_family_consistency_no_mutation(self):
        d=data();saved=copy.deepcopy(d)
        a=compare(d,permutations=31,seed=7)
        self.assertEqual(a,compare(d,permutations=31,seed=7))
        f=compare_family(d,[1.,2.],permutations=31,seed=7)
        self.assertEqual(f['tests'][0]['p_value'],a['p_value'])
        self.assertEqual(d,saved)

    def test_exact_within_block_null_rank(self):
        rows=[[0],[5],[.1],[5.1]];ps=[]
        for left in itertools.product([0,2],[1,3]):
            right=[i for i in range(4) if i not in left]
            d={'reference':[rows[i] for i in left], 'current':[rows[i] for i in right],
               'blocks':{'reference':['a','b'],'current':['a','b']}}
            ps.append(compare(d,exact=True)['p_value'])
        for alpha in [.1,.5,.75]:self.assertLessEqual(sum(p<=alpha for p in ps)/4,alpha)

    def test_mixture_fixture_conditional_vs_pooled(self):
        from driftlens.demo import run_demo
        r=run_demo();a=r['mixture_pooled'];b=r['mixture_blocked']
        self.assertEqual(a['mmd2_unbiased'],b['mmd2_unbiased'])
        self.assertLess(a['p_value'],.05)
        self.assertEqual(b['p_value'],1)
