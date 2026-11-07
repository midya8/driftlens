"""Selected development fixtures, not a detection-power benchmark."""
import random
from .core import compare


def run_demo():
    rng = random.Random(13)
    x = [[rng.gauss(0, 1) for _ in range(4)] for _ in range(30)]
    y = [[rng.gauss(0, 1) for _ in range(4)] for _ in range(30)]
    shifted = [[v + 2 for v in row] for row in y]
    return {'same_generator': compare({'reference': x, 'current': y}, bandwidth=2),
            'mean_shift': compare({'reference': x, 'current': shifted}, bandwidth=2),
            'exact_tiny': compare({'reference': [[0], [.1]], 'current': [[3], [3.1]]}, exact=True),
            'fixture_scope': 'Seed13 selected development fixture: two 30x4 standard-normal '
                             'samples and the same current rows shifted by +2 per coordinate. '
                             'Bandwidth2 and alpha.05. This is not a power/type-I error study.'}
