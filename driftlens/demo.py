"""Selected development fixtures, not a detection-power benchmark."""
import random
from .core import compare
from .family import compare_family


def run_demo():
    rng = random.Random(13)
    x = [[rng.gauss(0, 1) for _ in range(4)] for _ in range(30)]
    y = [[rng.gauss(0, 1) for _ in range(4)] for _ in range(30)]
    shifted = [[v + 2 for v in row] for row in y]
    confounded = {'reference': [[0]]*8 + [[5]]*2, 'current': [[0]]*2 + [[5]]*8}
    blocked = dict(confounded, blocks={'reference':['low']*8+['high']*2,
                                       'current':['low']*2+['high']*8})
    return {'mixture_pooled': compare(confounded),
            'mixture_blocked': compare(blocked),
            'shift_bandwidth_family': compare_family({'reference': x, 'current': shifted}, [1., 2., 4.]),
            'same_generator': compare({'reference': x, 'current': y}, bandwidth=2),
            'mean_shift': compare({'reference': x, 'current': shifted}, bandwidth=2),
            'exact_tiny': compare({'reference': [[0], [.1]], 'current': [[3], [3.1]]}, exact=True),
            'fixture_scope': 'Seed13 selected development fixture: two 30x4 standard-normal '
                             'samples and the same current rows shifted by +2 per coordinate. '
                             'Bandwidth2 and alpha.05. This is not a power/type-I error study.'}
