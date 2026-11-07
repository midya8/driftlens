import itertools
import io
import json
import unittest
from contextlib import redirect_stdout
from unittest.mock import patch
from driftlens import compare, compare_family
from driftlens.family import holm_adjust
from driftlens.cli import main


def data():return {'reference':[[0],[.1]],'current':[[3],[3.1]]}


class TestFamily(unittest.TestCase):
    def test_holm_known_vector_and_original_order(self):
        raw=[.04,.001,.5,.01,.2]
        expected=[.12,.005,.5,.04,.4]
        for a,b in zip(holm_adjust(raw),expected):self.assertAlmostEqual(a,b)

    def test_ties_and_caps(self):
        self.assertEqual(holm_adjust([.03,.03,.8]),[.09,.09,.8])
        self.assertEqual(holm_adjust([.9,1]),[1.,1.])
        self.assertEqual(holm_adjust([0,0]),[0.,0.])

    def test_single_family_matches_single_test(self):
        a=compare(data(),exact=True)
        b=compare_family(data(),[1.],exact=True)['tests'][0]
        self.assertEqual(a['p_value'],b['holm_p_value'])
        self.assertEqual(a['reject_at_alpha'],b['reject_within_family'])
        self.assertNotIn('reject_at_alpha',b)

    def test_bandwidth_order_preserved_and_common_seed_results(self):
        result=compare_family(data(),[2.,1.],permutations=31,seed=7)
        self.assertEqual([t['bandwidth'] for t in result['tests']],[2.,1.])
        for t in result['tests']:
            self.assertEqual(t['p_value'],compare(data(),t['bandwidth'],31,7)['p_value'])
            self.assertGreaterEqual(t['holm_p_value'],t['p_value'])

    def test_invalid_values_and_family_work_limit(self):
        for bs in [[],[1,1],[0],[-1],[True],[float('nan')],list(range(1,22))]:
            with self.subTest(bs=bs),self.assertRaises(ValueError):compare_family(data(),bs)
        for ps in [[],[True],[-.1],[float('inf')]]:
            with self.assertRaises(ValueError):holm_adjust(ps)
        d={'reference':[[0]]*30,'current':[[1]]*30}
        with self.assertRaises(ValueError):compare_family(d,list(range(1,14)))

    def test_exact_family_null_rank_bound_all_allocations(self):
        rows=[[0],[.1],[3],[3.1]]; rejections=0
        for a in itertools.combinations(range(4),2):
            b=[i for i in range(4) if i not in a]
            r=compare_family({'reference':[rows[i] for i in a],'current':[rows[i] for i in b]},
                             [.2,1.,4.],alpha=.5,exact=True)
            rejections+=r['reject_any']
        self.assertLessEqual(rejections/6,.5)

    def test_cli_family(self):
        out=io.StringIO()
        with patch('sys.stdin',io.StringIO(json.dumps(data()))),redirect_stdout(out):
            self.assertEqual(main(['compare','--input','-','--bandwidths','1','2','--exact']),0)
        self.assertEqual(json.loads(out.getvalue())['adjustment'],'Holm step-down')

    def test_decision_uses_adjusted_not_raw_p(self):
        r=compare_family(data(),[.2,1.],alpha=.4,exact=True)
        self.assertTrue(all(t['p_value']<=.4 for t in r['tests']))
        self.assertFalse(r['reject_any'])
