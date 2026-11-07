# DriftLens

Embedding-distribution comparison with an RBF kernel two-sample test, built
from scratch in dependency-free Python. Unbiased squared MMD, one cached
pooled kernel, exact enumeration for tiny samples and seeded Monte Carlo
permutations for larger ones. JSON CLI plus importable API. No paid APIs,
model downloads or private datasets.

## Measured development fixture

```bash
python -m driftlens demo
bash examples/demo.sh
```

| Arm | Reference/current | Unbiased MMD squared | p-value | Reject at .05? |
| --- | --- | ---: | ---: | --- |
| Independent same-generator samples | 30x4 / 30x4 | -0.004871 | 0.550 | No |
| Same current rows shifted +2 per coordinate | 30x4 / 30x4 | 0.634539 | 0.002 | Yes |
| Tiny separated samples, exact test | 2x1 / 2x1 | 1.967361 | 1/3 | No |

Seed13, bandwidth2 for the 30x4 arms, 499 permutations. The shifted p-value
is at the Monte Carlo resolution floor 1/500, not evidence of a smaller
measured p-value. The tiny test enumerates six allocations, so even visibly
separated samples cannot reject at .05. These fixtures/settings were selected
during development, not an independent power or false-positive benchmark.
No inference from one non-rejection that the two distributions are equal.

## Compare supplied embeddings

```bash
python -m driftlens compare --input examples/embeddings.json --bandwidth 1 --permutations 499
python -m driftlens compare --input examples/embeddings.json --bandwidth 1 --exact
```

```python
from driftlens import compare
report = compare({"reference": reference_rows, "current": current_rows},
                 bandwidth=2, permutations=499, seed=13, alpha=.05)
```

Input fields are exactly `reference` and `current`, each a numeric row matrix
with at least two observations and the same nonzero dimension. `--input -`
reads stdin. Invalid input exits2 with stderr and no partial JSON. Output
includes statistic, p-value, decision, exceedances, resolution, sample sizes,
permutation mode and kernel-evaluation count. Seeds use a local RNG.

Bandwidth is required on the CLI. Choose it before inspecting group-specific
results, or use a prespecified label-invariant pooled rule outside this tool.
Feature preprocessing must follow the same rule. Do not search bandwidths
and report the best uncorrected p-value. Different scales can dominate RBF
distance; this tool does not silently standardize, normalize or tune embeddings.

## Statistic and null test

For kernel `k(a,b)=exp(-||a-b||^2/(2*bandwidth^2))`, the statistic is
`sum(i!=j)k(x_i,x_j)/(n(n-1)) + sum(i!=j)k(y_i,y_j)/(m(m-1)) - 2sum(i,j)k(x_i,y_j)/(nm)`.
Within-sample diagonals are excluded. Finite-sample estimates can be negative
and are never clipped. With fixed kernel and IID samples it is an unbiased
estimate of population squared MMD. A data-dependent kernel changes that claim.

Exact mode enumerates every `choose(n+m,n)` label allocation and divides
upper-tail exceedances by that count. Monte Carlo mode samples independent
uniform allocations with replacement and reports `(exceedances+1)/(B+1)`.
The observed split is the additional reference allocation. Near ties within
1e-12 count as exceedances, conservatively. Seed/permutation count must be
chosen without optional stopping on interim results.

The pooled kernel is built once: N(N-1)/2 off-diagonal RBF evaluations,
O(N^2*d) construction and O(B*N^2) permutation scoring, O(N^2) memory.
This is not a streaming or million-vector approximate implementation.
Limits: 500 pooled rows, 4096 dimensions, 10000 draws/allocations and
20 million draw-by-pooled-row-squared work units. Exact mode rejects larger
allocation spaces rather than pretending a subsample is exact.

## What this does not establish

This is one two-sample upper-tail test. Null validity needs exchangeability
of group labels, such as independent samples from the same distribution.
Correlated time windows, duplicate/repeated users or group-dependent
preprocessing can violate that assumption. A rejection says neither why
features changed nor whether a deployed model's accuracy worsened. A large
sample can detect a tiny practically irrelevant change. Non-rejection is not
equivalence. Repeated/multiple tests need an explicit error-control plan;
there is no always-valid monitoring guarantee here.

## Verification and attribution

17 tests cover closed-form kernels, direct MMD arithmetic, negative estimates,
exact enumeration and rank superuniformity, ties, unequal groups, symmetry,
seed/no-mutation behavior, +1 resolution, one cached kernel, extremes,
work limits, CLI failures and pinned demo outputs. Tests run on Python
3.10-3.13; the offline demo asserts both gate directions.

Standard method, original small implementation, not a new statistical test:
- [Gretton et al., kernel two-sample test](https://www.gatsby.ucl.ac.uk/~gretton/mmd/mmd.htm)
- [Phipson and Smyth, permutation p-values should never be zero](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf)

See [DESIGN.md](DESIGN.md). MIT licensed.
