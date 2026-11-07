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
