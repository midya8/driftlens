"""Original dependency-free implementation of standard RBF MMD permutation tests."""
import math
import random
from .allocations import partition, allocation_count, allocations


def _validate(x, y, bandwidth):
    if type(bandwidth) not in (int, float) or not math.isfinite(bandwidth) or bandwidth <= 0:
        raise ValueError('bandwidth must be finite and positive, chosen independently of group labels')
    if not isinstance(x, list) or not isinstance(y, list) or min(len(x), len(y)) < 2:
        raise ValueError('reference/current must each contain at least two rows')
    if not isinstance(x[0], list) or not x[0]:
        raise ValueError('rows must be nonempty numeric lists')
    width = len(x[0])
    if len(x) + len(y) > 500 or width > 4096:
        raise ValueError('limit is 500 pooled rows and 4096 dimensions')
    for row in x + y:
        if not isinstance(row, list) or len(row) != width or any(
            type(v) not in (int, float) or not math.isfinite(v) for v in row
        ):
            raise ValueError('rows must have equal dimension and finite numeric values')
    return width


def kernel_matrix(x, y, bandwidth):
    _validate(x, y, bandwidth)
    rows = x + y; n = len(rows)
    matrix = [[0.] * n for _ in rows]
    for i, a in enumerate(rows):
        matrix[i][i] = 1.
        for j in range(i):
            terms = []
            for u, v in zip(a, rows[j]):
                d = (u - v) / bandwidth
                if not math.isfinite(d) and bandwidth >= 1:
                    d = u / bandwidth - v / bandwidth
                # Large scaled differences have kernels indistinguishable from 0.
                terms.append(d * d)
            try:
                distance = math.fsum(terms) if all(math.isfinite(t) for t in terms) else float('inf')
            except OverflowError:
                distance = float('inf')
            value = math.exp(-.5 * distance)
            matrix[i][j] = matrix[j][i] = value
    return matrix


def statistic(matrix, left, right):
    """Unbiased squared MMD: off-diagonal within sums, all cross pairs.

    Finite-sample estimates can be negative. Do not clip before permutation.
    """
    n, m = len(left), len(right)
    a = math.fsum(matrix[i][j] for index, i in enumerate(left) for j in left[index + 1:])
    b = math.fsum(matrix[i][j] for index, i in enumerate(right) for j in right[index + 1:])
    c = math.fsum(matrix[i][j] for i in left for j in right)
    return 2 * a / (n * (n - 1)) + 2 * b / (m * (m - 1)) - 2 * c / (n * m)


def compare(payload, bandwidth=1., permutations=499, seed=13, alpha=.05, exact=False):
    groups, blocked = partition(payload)
    x, y = payload['reference'], payload['current']
    width = _validate(x, y, bandwidth)
    if type(permutations) is not int or not 1 <= permutations <= 10000:
        raise ValueError('permutations must be an integer from 1 to 10000')
    if type(seed) is not int or type(exact) is not bool:
        raise ValueError('seed must be integer and exact must be boolean')
    if type(alpha) not in (int, float) or not math.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError('alpha must be finite and strictly between zero and one')
    n, m = len(x), len(y); size = n + m; total = allocation_count(groups)
    draws = total if exact else permutations
    if exact and total > 10000:
        raise ValueError('exact mode exceeds 10000 label allocations')
    if draws * size * size > 20000000:
        raise ValueError('permutation work exceeds 20000000 pooled-row squares')
    matrix = kernel_matrix(x, y, bandwidth)
    observed = statistic(matrix, list(range(n)), list(range(n, size)))
    indices = list(range(size)); rng = random.Random(seed)
    exceedances = 0
    for chosen in allocations(groups, exact, permutations, rng):
        left = sorted(chosen); members = set(left)
        right = [i for i in indices if i not in members]
        value = statistic(matrix, left, right)
        # Conservative inclusion of floating near-ties, fixed before seeing data.
        exceedances += value >= observed - 1e-12
    p = exceedances / total if exact else (exceedances + 1) / (permutations + 1)
    return {'mmd2_unbiased': observed, 'p_value': p, 'reject_at_alpha': p <= alpha,
            'alpha': alpha, 'bandwidth': bandwidth, 'reference_n': n, 'current_n': m,
            'dimensions': width, 'mode': 'exact' if exact else 'monte_carlo',
            'draws': draws, 'exceedances': exceedances, 'seed': seed,
            'p_resolution': 1 / total if exact else 1 / (permutations + 1),
            'kernel_evaluations': size * (size - 1) // 2,
            'total_label_allocations': total,
            'permutation_scheme': 'within_blocks' if blocked else 'pooled',
            'blocks': [{'size': len(group), 'reference_n': count, 'current_n': len(group)-count,
                        'allocations': math.comb(len(group), count)} for group,count in groups] if blocked else [],
            'informative_blocks': sum(0 < count < len(group) for group,count in groups) if blocked else None,
            'scope': 'One upper-tail two-sample permutation test. Requires exchangeability '
                     'of group labels under the null (within each prespecified block if supplied), '
                     'joint invariance under allowed permutations, and prespecified or label-invariant '
                     'bandwidth/preprocessing. No cause, model-quality, practical-effect or '
                     'continuous-monitoring guarantee. Do not interpret non-rejection as '
                     'equivalence. Repeated/multiple tests need an explicit error-control plan. '
                     'Near ties within 1e-12 count conservatively as exceedances.'}
