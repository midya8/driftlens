# DriftLens design

## A fixed kernel and explicit assumptions

Inputs are matrices of embedding-like numeric vectors. The library does not
train a model, choose semantic features, infer data lineage or prove that
rows are independent. Bandwidth/preprocessing must be specified independently
of group labels, or by a prespecified pooled label-invariant rule. The CLI
requires bandwidth to avoid silently presenting an arbitrary heuristic as
an optimal choice. Tuning it against group separation needs fresh data or
valid multiple-test correction, not the same nominal p-value.

## Correct statistic and cached computation

The unbiased squared MMD uses off-diagonal within-group averages and the
full cross-group average. It can be negative on finite samples. Pooled
kernel entries are computed once and then re-indexed for permutations.
RBF calculation divides coordinate differences by bandwidth before squaring;
opposite extreme floats use a division-first fallback when appropriate.
Unrepresentably large squared distances safely become zero kernel values.

Statistic sums use math.fsum. Comparisons count values within 1e-12 below
observed as exceedances, a prespecified conservative near-tie treatment.
Exact tail counts are identical under complement labels for equal sizes.
No clipping or data-dependent tie tolerance is used.

## Permutation conventions

Exact mode enumerates every allocation of n reference indices among N pooled
indices, including observed. Monte Carlo samples each allocation uniformly
and independently with replacement using a local seeded RNG. Thus duplicate
allocations and the observed allocation can recur. Its p-value uses the +1
formula, not a zero-capable tail fraction. Exact mode uses exceedances/total.
The finite resolution is reported and limits how small a p-value can be.
No adaptive permutation stopping or repeated-stream testing is implemented.

## Development evidence

Seed13 builds independent 30x4 standard-normal samples; the shifted arm adds2
to each current coordinate. Bandwidth2, B499, alpha.05. Observed same-generator
MMD squared -0.00487090068 and p.55; shifted .63453882456 and p.002. A tiny
2-vs-2 separated fixture has exact p1/3, illustrating low power/resolution.
These selected development fixtures are not a power/type-I error study.
They do not establish real deployment drift-detection performance.

## Limits and security

At least two rows per sample; finite numeric rectangular matrices. Up to500
pooled rows and4096dimensions. Draw count1..10000; exact allocation count<=10000;
draws*N^2<=20000000. JSON parsing and schema checks precede expensive kernel
work. Still, this is an offline local analysis package, not an internet-facing
service hardened against arbitrarily large request bodies. No eval, subprocess,
network request or file mutation is performed by the library. The CLI reads
only its input path and prints JSON; no observations are uploaded.

## Statistical interpretation

Exchangeability under the null is necessary for label permutation validity.
Sequentially overlapping windows, repeated users, sampling bias or changes
in preprocessing can violate it. Kernel tests do not diagnose a cause or
predict downstream performance. Non-rejection is not equality/equivalence;
rejection need not imply practical harm. A separate plan is needed for repeated
or multiple tests, alpha spending and investigation thresholds.

## Sources

- https://www.gatsby.ucl.ac.uk/~gretton/mmd/mmd.htm
- https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf

Standard algorithms with an original dependency-free implementation, no new
test, theorem or benchmark-superiority claim.

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
