"""Prespecified bandwidth families with standard Holm p-value adjustment."""
import math
from .core import compare, _validate
from .allocations import partition, allocation_count


def holm_adjust(p_values):
    if not isinstance(p_values, list) or not p_values or any(
        type(p) not in (int, float) or not math.isfinite(p) or not 0 <= p <= 1 for p in p_values
    ):
        raise ValueError('p_values must be a nonempty list of finite numbers in [0,1]')
    order = sorted(range(len(p_values)), key=p_values.__getitem__)
    adjusted = [0.] * len(order); running = 0.
    for rank, index in enumerate(order):
        running = max(running, min(1., (len(order) - rank) * p_values[index]))
        adjusted[index] = running
    return adjusted


def compare_family(payload, bandwidths, permutations=499, seed=13, alpha=.05, exact=False, witness=0):
    """Run a fixed kernel family, reusing the seed for common label allocations.

    Correction covers this supplied family only. It cannot repair adaptive
    selection, invalid input p-values, sequential looks or omitted tests.
    """
    if not isinstance(bandwidths, list) or not 1 <= len(bandwidths) <= 20 or any(
        type(b) not in (int, float) or not math.isfinite(b) or b <= 0 for b in bandwidths
    ) or len(set(bandwidths)) != len(bandwidths):
        raise ValueError('bandwidths must be 1..20 distinct finite positive numbers')
    groups, _ = partition(payload)
    x, y = payload['reference'], payload['current']
    _validate(x, y, bandwidths[0]); size = len(x) + len(y)
    if type(exact) is not bool or type(permutations) is not int or not 1 <= permutations <= 10000:
        raise ValueError('exact must be boolean; permutations must be integer 1..10000')
    draws = allocation_count(groups) if exact else permutations
    if len(bandwidths) * draws * size * size > 20000000:
        raise ValueError('family exceeds 20000000 draw-by-pooled-row-squared work units')
    reports = [compare(payload, b, permutations, seed, alpha, exact, witness) for b in bandwidths]
    adjusted = holm_adjust([r['p_value'] for r in reports])
    for report, p in zip(reports, adjusted):
        report['holm_p_value'] = p
        report['reject_within_family'] = p <= alpha
        # Raw p-values remain visible, but do not expose a competing decision.
        del report['reject_at_alpha']
    return {'adjustment': 'Holm step-down', 'alpha': alpha, 'family_size': len(reports),
            'reject_any': any(r['reject_within_family'] for r in reports),
            'tests': reports,
            'scope': 'Family-wise correction for this fixed, prespecified bandwidth list only. '
                     'Requires valid individual permutation p-values; arbitrary dependence '
                     'among them is allowed. Same seed uses common label allocations. '
                     'Does not cover adaptive bandwidth/preprocessing choices, omitted tests, '
                     'repeated windows, optional stopping or invalid exchangeability. '
                     'Not a posterior probability, causal diagnosis or performance guarantee.'}
