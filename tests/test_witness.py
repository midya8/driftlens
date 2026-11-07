import io
import json
import unittest
from contextlib import redirect_stdout
from unittest.mock import patch
from driftlens import compare, compare_family
from driftlens.core import kernel_matrix
from driftlens.cli import main


def data():return {'reference':[[0],[1]],'current':[[2],[3],[4]]}


class TestWitness(unittest.TestCase):
    def test_manual_kernel_means_and_biased_mmd(self):
        d=data();k=kernel_matrix(d['reference'],d['current'],1.)
        r=compare(d,witness=5)['witness'];n=2;m=3
        expected=[sum(row[:n])/n-sum(row[n:])/m for row in k]
        for name,start,end in [('reference',0,n),('current',n,n+m)]:
            for entry in r['groups'][name]['top_absolute']:
                self.assertAlmostEqual(entry['witness'],expected[start+entry['row_index']])
        biased=sum(k[i][j] for i in range(n) for j in range(n))/n**2+sum(k[i][j] for i in range(n,n+m) for j in range(n,n+m))/m**2-2*sum(k[i][j] for i in range(n) for j in range(n,n+m))/(n*m)
        self.assertAlmostEqual(r['mmd2_biased'],biased)
        self.assertGreaterEqual(r['mmd2_biased'],-1e-12)

    def test_sign_flips_on_group_swap(self):
        d=data();a=compare(d,witness=5)['witness']
        b=compare({'reference':d['current'],'current':d['reference']},witness=5)['witness']
        av={e['row_index']:e['witness'] for e in a['groups']['current']['top_absolute']}
        bv={e['row_index']:e['witness'] for e in b['groups']['reference']['top_absolute']}
        for i in av:self.assertAlmostEqual(av[i],-bv[i])
        self.assertAlmostEqual(a['mmd2_biased'],b['mmd2_biased'])

    def test_top_absolute_ties_stable_and_no_vectors(self):
        d={'reference':[[0],[0]],'current':[[0],[0]]}
        w=compare(d,witness=1)['witness']
        for g in w['groups'].values():self.assertEqual(g['top_absolute'],[{'row_index':0,'witness':0.}])
        self.assertNotIn('coordinates',json.dumps(w))

    def test_one_kernel_no_p_or_stat_change(self):
        a=compare(data())
        with patch('driftlens.core.kernel_matrix',wraps=kernel_matrix) as cached:
            b=compare(data(),witness=3)
            cached.assert_called_once()
        del b['witness'];self.assertEqual(a,b)
        self.assertNotIn('witness',a)

    def test_invalid_count_and_family_propagation(self):
        for value in [True,-1,101,.5]:
            with self.subTest(value=value),self.assertRaises(ValueError):compare(data(),witness=value)
        f=compare_family(data(),[1.,2.],witness=1)
        self.assertTrue(all('witness' in t for t in f['tests']))

    def test_blocked_scope_is_not_conditional_attribution(self):
        d=data();d['blocks']={'reference':['a']*2,'current':['a']*3}
        w=compare(d,witness=2)['witness']
        self.assertIn('not a block-adjusted',w['scope'])

    def test_cli(self):
        out=io.StringIO()
        with patch('sys.stdin',io.StringIO(json.dumps(data()))),redirect_stdout(out):
            self.assertEqual(main(['compare','--input','-','--bandwidth','1','--witness','2']),0)
        self.assertEqual(json.loads(out.getvalue())['witness']['top_per_group'],2)
