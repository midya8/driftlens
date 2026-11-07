"""Complement-sum MMD scorer for imbalanced groups.

Keeps the direct statistic as reference. For near-balanced groups it uses
that reference to avoid extra complement cancellation without useful savings.
"""
import math


class CachedScorer:
    def __init__(self, matrix, reference_n):
        self.matrix = matrix
        self.n = reference_n
        self.size = len(matrix)
        self.m = self.size - self.n
        self.small = min(self.n, self.m)
        # Require at least a 2x per-allocation pair-read saving.
        direct = self.n * (self.n-1)//2 + self.m * (self.m-1)//2 + self.n*self.m
        fast = self.small * (self.small-1)//2
        self.use_complement = self.small * 4 <= self.size
        self.pair_reads = fast if self.use_complement else direct
        if self.use_complement:
            self.row_sums = [math.fsum(v for j,v in enumerate(row) if i != j)
                             for i,row in enumerate(matrix)]
            self.total = math.fsum(self.row_sums) / 2
        else:
            self.row_sums = None
            self.total = None

    def score(self, left, right):
        if not self.use_complement:
            from .core import statistic
            return statistic(self.matrix, left, right)
        group = left if self.n <= self.m else right
        a = math.fsum(self.matrix[i][j] for index,i in enumerate(group) for j in group[index+1:])
        cross = math.fsum(self.row_sums[i] for i in group) - 2*a
        other = self.total - a - cross
        small, large = self.small, self.size-self.small
        return 2*a/(small*(small-1)) + 2*other/(large*(large-1)) - 2*cross/(small*large)
