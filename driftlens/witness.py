"""Descriptive unnormalized empirical RKHS witness, including self kernels."""
import math


def witness_report(matrix, reference_n, top):
    n = reference_n; size = len(matrix); m = size - n
    scores = [math.fsum(row[:n]) / n - math.fsum(row[n:]) / m for row in matrix]
    groups = {}
    means = {}
    for name, start, end in [('reference', 0, n), ('current', n, size)]:
        mean = math.fsum(scores[start:end]) / (end - start)
        means[name] = mean
        order = sorted(range(start, end), key=lambda i: (-abs(scores[i]), i))[:top]
        groups[name] = {'mean_witness': mean,
                        'top_absolute': [{'row_index': i-start, 'witness': scores[i]} for i in order]}
    return {'definition': 'mean_reference k(x,z) - mean_current k(y,z), unnormalized, includes self kernels',
            'top_per_group': top, 'groups': groups,
            'mmd2_biased': means['reference'] - means['current'],
            'scope': 'Descriptive in-sample row-index localization with pooled group means. '
                     'Positive means more reference-like kernel similarity; negative more current-like. '
                     'Not a density, row-level p-value, anomaly verdict, feature attribution or cause. '
                     'Self kernels and sample imbalance affect scores. With block permutations this '
                     'is still a pooled witness, not a block-adjusted or conditional-effect explanation. '
                     'Biased MMD squared here is not the unbiased statistic used for permutation testing.'}
