"""Finite Gram transformations, separate from floating RBF comparisons."""
import copy
from fractions import Fraction as F
import unittest

from driftlens.orbit import audit_orbit
from driftlens.orbit_verify import verify_orbit


def requests():
    vectors = [[0, 2, 0, 1, -1, -3], [-1, 1, -3, 1, 0, 1],
               [1, 3, -3, 1, 2, 3]]
    for reference, blocks in [([0, 2, 3], None),
                              ([0, 2], ['a', 'a', 'b', 'b', 'c', 'c']),
                              ([0, 1, 4], ['a', 'a', 'b', 'b', 'c', 'c'])]:
        yield dict(row_names=['row' + str(i) for i in range(6)],
                   reference=reference, blocks=blocks, alpha='1/5',
                   kernels=[dict(name='view' + str(k),
                                 gram=[[str(x * y) for y in v] for x in v])
                            for k, v in enumerate(vectors)])


class OrbitTransformTests(unittest.TestCase):
    def assert_transformed_report(self, request, original, scales, offsets, diagonals):
        changed = copy.deepcopy(request)
        for k, kernel in enumerate(changed['kernels']):
            kernel['gram'] = [[str(scales[k] * F(x) + offsets[k] +
                                  (diagonals[k][i] if i == j else 0))
                               for j, x in enumerate(row)]
                              for i, row in enumerate(kernel['gram'])]
        saved = copy.deepcopy(changed)
        report = audit_orbit(changed)
        self.assertTrue(verify_orbit(report))
        self.assertEqual(changed, saved)
        expected = copy.deepcopy(original)
        expected['request'] = changed
        # The independent verifier checks the changed LDL witnesses. Their ranks
        # need not be invariant, even when all test statistics are unchanged.
        expected['psd_certificates'] = report['psd_certificates']
        for row in expected['allocations']:
            row['mmd2_unbiased'] = [str(scales[k] * F(x))
                                     for k, x in enumerate(row['mmd2_unbiased'])]
        self.assertEqual(report, expected)
        return report

    def test_positive_kernel_scales_and_constant_offsets(self):
        checks = 0
        for request in requests():
            original = audit_orbit(request)
            self.assertTrue(verify_orbit(original))
            for scale in [F(1, 3), F(2), F(7, 2)]:
                for offset in [F(0), F(5, 3)]:
                    self.assert_transformed_report(request, original,
                        [scale, scale * 2, scale / 2],
                        [offset, offset * 2, offset / 2], [[F(0)] * 6] * 3)
                    checks += 1
        self.assertEqual(checks, 18)

    def test_nonnegative_diagonal_changes_leave_complete_test_unchanged(self):
        checks = 0
        for request in requests():
            original = audit_orbit(request)
            for pattern in [[F(1)] * 6, [F(i, 3) for i in range(6)],
                            [F(0), F(0), F(0), F(0), F(0), F(7)]]:
                report = self.assert_transformed_report(request, original,
                    [F(1)] * 3, [F(0)] * 3,
                    [[x * (k + 1) for x in pattern] for k in range(3)])
                self.assertNotEqual(report['psd_certificates'], original['psd_certificates'])
                checks += 1
        self.assertEqual(checks, 9)
