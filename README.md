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

## Prespecified bandwidth families

```bash
python -m driftlens compare --input examples/embeddings.json --bandwidths 1 2 4
```

```python
from driftlens import compare_family
family = compare_family(payload, [1., 2., 4.], permutations=499, seed=13)
```

The fixed supplied family gets standard Holm step-down adjusted p-values,
per-kernel `reject_within_family` and a family-level `reject_any`. Raw p-values
stay visible; the competing uncorrected decision field is omitted. The same
seed generates common label allocations across kernels. Arbitrary dependence
among valid input p-values is allowed by Holm. It cannot repair invalid
exchangeability, adaptively chosen bandwidths/preprocessing, omitted tests,
repeated windows or optional stopping. Prespecify the entire tested family.
This is family-wise error control, not a posterior error probability.

For the selected shifted demo, bandwidths [1,2,4] each have raw p=0.002 and
Holm p=0.006. A tiny exact example has raw p=1/3; testing two kernels at
alpha=0.4 does not reject after Holm, despite both raw p-values passing.
Families allow 1..20 distinct positive bandwidths. Work limits apply to the
whole family, not each kernel separately; kernels are built once each and
processed in sequence rather than all retained in memory.

Source: [Holm formula and family-wise scope](https://pingouin-stats.org/generated/pingouin.multicomp.html).
Standard correction, original small implementation.

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

25 tests cover closed-form kernels, direct MMD arithmetic, negative estimates,
exact enumeration and rank superuniformity, ties, unequal groups, symmetry,
seed/no-mutation behavior, +1 resolution, one cached kernel, extremes,
work limits, CLI failures and pinned demo outputs. Tests run on Python
3.10-3.13; the offline demo asserts both gate directions.

Standard method, original small implementation, not a new statistical test:
- [Gretton et al., kernel two-sample test](https://www.gatsby.ucl.ac.uk/~gretton/mmd/mmd.htm)
- [Phipson and Smyth, permutation p-values should never be zero](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf)

See [DESIGN.md](DESIGN.md). MIT licensed.
