"""Prespecified within-stratum label allocations, not a cluster bootstrap."""
import itertools
import math


def partition(payload):
    if not isinstance(payload, dict) or set(payload) not in (
        {'reference', 'current'}, {'reference', 'current', 'blocks'}
    ):
        raise ValueError('input must have reference/current and optional blocks')
    x, y = payload['reference'], payload['current']
    if not isinstance(x, list) or not isinstance(y, list):
        raise ValueError('reference/current must be lists')
    n, m = len(x), len(y)
    if 'blocks' not in payload:
        return [(list(range(n + m)), n)], False
    blocks = payload['blocks']
    if not isinstance(blocks, dict) or set(blocks) != {'reference', 'current'}:
        raise ValueError('blocks must have exactly reference/current label lists')
    for key, length in [('reference', n), ('current', m)]:
        labels = blocks[key]
        if not isinstance(labels, list) or len(labels) != length or any(
            not isinstance(v, str) or not v.strip() for v in labels
        ):
            raise ValueError('block labels must be aligned nonempty strings')
    labels = blocks['reference'] + blocks['current']
    names = list(dict.fromkeys(labels))
    return [(indices, sum(i < n for i in indices)) for name in names
            for indices in [[i for i, value in enumerate(labels) if value == name]]], True


def allocation_count(groups):
    return math.prod(math.comb(len(indices), count) for indices, count in groups)


def allocations(groups, exact, draws, rng):
    if exact:
        choices = [list(itertools.combinations(indices, count)) for indices, count in groups]
        for combination in itertools.product(*choices):
            yield [i for group in combination for i in group]
    else:
        for _ in range(draws):
            yield [i for indices, count in groups for i in rng.sample(indices, count)]
