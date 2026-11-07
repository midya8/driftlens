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

Required input fields are `reference` and `current`, each a numeric row matrix
with at least two observations and the same nonzero dimension. Optional `blocks`
selects the within-block scheme described below. `--input -`
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

Without blocks, exact mode enumerates every `choose(n+m,n)` label allocation and divides
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

## Prespecified within-block permutations

Optional input `blocks` contains aligned reference/current string labels:

```json
{"reference": [[0], [5]], "current": [[0.1], [5.1]],
 "blocks": {"reference": ["site-a", "site-b"], "current": ["site-a", "site-b"]}}
```

The same CLI/API accepts this schema, including exact and bandwidth-family
modes. Every allocation preserves the observed reference/current count in
each block. The exact space is the product of within-block combinations;
Monte Carlo draws independently within each block. Reports show the scheme,
per-block sizes/counts/allocation counts and informative-block count.
A block found on only one side is frozen. If all blocks are frozen there
is one allocation and p=1, explicitly uninformative, not evidence of equality.

This tests a different conditional null from unrestricted label permutations:
the joint data distribution must be invariant under the allowed within-block
label reallocations. Blocks must be prespecified, not discovered by searching
for a preferred result. It does not fix arbitrary correlation, permute whole
clusters, bootstrap repeated users, or implement nested exchangeability.
The pooled MMD formula is retained as the statistic, not a new unbiased
conditional-effect estimator or a causal adjustment.

Selected mixture demonstration: reference has eight zero rows/two five rows;
current has two zero rows/eight five rows. Pooled MMD squared=0.648886 in both
schemes. Unrestricted p=0.014; prespecified value-stratum blocks give p=1.000
because features are identical within each stratum. This artificial example
shows different null questions, not a claim that the blocked result is
universally more correct, a drift cause or practical validity on real data.
Unrestricted exact space=184756, blocked space=2025; reported p-values here
use499 Monte Carlo draws, seed13, bandwidth1.

Reference for within-group exchangeability:
[Berkeley StatLab stratified permutation testing](https://statlab.github.io/permute/api/stratified.html).
This implementation uses the same restriction principle with an MMD statistic,
not their mean statistic or code.

## Descriptive witness diagnostics

```bash
python -m driftlens compare --input examples/embeddings.json --bandwidth 1 --witness 3
```

API: `compare(payload, bandwidth=1., witness=3)`; family comparisons accept
`witness` too. Zero (default) omits diagnostics; integers1..100 request that
many highest-absolute-score rows per group, capped by available rows. Output
contains only group-local zero-based row indices and scores, not vectors.
Ties are broken by original row index. The existing pooled kernel is reused;
no new kernel evaluations or changes to p-values/test statistics occur.

At each observed row z, the unnormalized empirical witness is
`mean_reference k(x,z) - mean_current k(y,z)`, including self kernels.
Positive scores mean greater kernel similarity to reference; negative to
current. Group means and their difference, `mmd2_biased`, are also reported.
That biased squared MMD is distinct from the diagonal-excluding unbiased
statistic used for permutation tests. No norm normalization is performed.

This is descriptive in-sample localization, not row-level p-values, an anomaly
verdict, a density, feature attribution or a cause. Self kernels, sample
imbalance and bandwidth affect scores. With blocks supplied, the witness
still uses pooled means; it is not block-adjusted or a conditional-effect
explanation. Choosing rows after seeing these scores needs separate validation.

Source: [Gretton et al., Section2.3, empirical witness](https://www.jmlr.org/papers/volume13/gretton12a/gretton12a.pdf).

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

39 tests cover closed-form kernels, direct MMD arithmetic, negative estimates,
exact enumeration and rank superuniformity, ties, unequal groups, symmetry,
seed/no-mutation behavior, +1 resolution, one cached kernel, extremes,
work limits, CLI failures and pinned demo outputs. Tests run on Python
3.10-3.13; the offline demo asserts both gate directions.

Standard method, original small implementation, not a new statistical test:
- [Gretton et al., kernel two-sample test](https://www.gatsby.ucl.ac.uk/~gretton/mmd/mmd.htm)
- [Phipson and Smyth, permutation p-values should never be zero](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf)

See [DESIGN.md](DESIGN.md). MIT licensed.
