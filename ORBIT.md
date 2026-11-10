# Choosing is part of the test

DriftLens 0.2 adds a complete finite supplied-Gram randomization orbit. This is an additive dependency-free mode, separate from the original floating RBF `compare` and Holm family APIs. It does not estimate bandwidths, fit embeddings, train a model or upload data. Existing source, demos, scoring choices, statistical limits and API behavior stay unchanged.

## The original selected control

Six pooled rows, three fixed rational rank-one feature views, twenty balanced label allocations, alpha `1/5`:

| Rule | Rejected allowed allocations | Exact fraction |
|---|---:|---:|
| Uncorrected smallest raw p | 8/20 | 2/5 |
| Recalibrated global minimum-p test | 4/20 | 1/5 |
| Holm, any supplied kernel | 0/20 | 0 |

For observed reference indices `[0,2,3]`, the three exact unbiased MMD² values are `[-11/9,4/3,64/9]`, their upper-tail counts are `[14,8,2]`, raw p-values `[7/10,2/5,1/10]`, Holm p-values `[4/5,4/5,3/10]`, and the global minimum-p value is `1/5`. Its global tail contains allocation indices `[4,5,14,15]` in lexicographic reference-index order.

The feature vectors used to construct the three Gram matrices are `[0,2,0,1,-1,-3]`, `[-1,1,-3,1,0,1]`, and `[1,3,-3,1,2,3]`, using `K[i,j]=v[i]*v[j]`. These are original supplied development controls, searched and selected to make distinctions visible. They are not a random trial, holdout, practical-power comparison or real-data false-positive benchmark. The fractions above describe this entire finite orbit under a uniform allowed-allocation null. They do not establish that any observed real experiment has that null law. No estimator is claimed universally better.

## Exact finite construction

For pooled finite symmetric positive-semidefinite Gram matrices `K[k]`, and allowed reference-index sets `a` of size `n` with complement size `m`, compute the ordinary diagonal-excluding statistic:

```
t[k,a] = sum(i != j in a) K[k,i,j] / (n*(n-1))
       + sum(i != j in complement(a)) K[k,i,j] / (m*(m-1))
       - 2 * sum(i in a,j in complement(a)) K[k,i,j] / (n*m)
```

All sums and subsequent probabilities use rational arithmetic, including negative finite MMD² values. Native code uses `fractions.Fraction`; independent browser code uses reduced `BigInt` numerator/denominator pairs. No floating tolerance decides a rank. Entry diagonals are excluded from the MMD² statistic but are included in the PSD certificate.

Let `A` be the complete allowed allocation set and `N=|A|`. Every kernel uses the same ordered allocation list. For every possible observed allocation `a`:

```
r[k,a] = #{b in A: t[k,b] >= t[k,a]}
q[a] = min_k r[k,a]
global_tail[a] = {b in A: q[b] <= q[a]}
p_global[a] = |global_tail[a]| / N
```

All upper and lower tail ties count, including the allocation itself. Counts are integers, probabilities canonical reduced fractions. Complete enumeration uses `count/N`, not a Monte Carlo `+1` correction. There is no Monte Carlo option in this mode. The original RBF mode retains its existing Monte Carlo semantics.

This is one global randomization test under a justified uniform conditional allowed-allocation null and a fixed label-independent family. It is a bounded implementation of established rank/minimum-p randomization ideas, not a new statistical method, a full step-down Westfall-Young implementation, MMDAgg or an automatic bandwidth optimizer. It does not issue adjusted component-wise global-test decisions or claim strong family-wise control for arbitrary different component nulls. The separately shown Holm values use the usual step-down rule on each allocation's exact raw ranks. They remain comparison values, not the global procedure's per-kernel conclusion.

### Why the finite size count is bounded

For any finite scalar list `q[a]`, the lower-tail p-value `#{b:q[b]<=q[a]}/N` is superuniform when `a` is uniformly selected from `A`: ordering the distinct `q` levels partitions `A` into tied blocks. Each block receives its cumulative count divided by `N`; rejection at alpha includes only an initial collection with cumulative size at most `alpha*N`. This finite counting fact holds regardless of dependence among kernels and even if the original `q` list is a discrete statistic rather than a valid standalone p-value. It does not prove the experiment supplies a uniformly selected allowed allocation. That assumption belongs to the sampling/randomization design.

One kernel reduces to its usual tied upper-tail p-value. Duplicating a kernel does not change the minimum or its recalibration. Constant or zero geometries give p=1. There is no randomized tie breaking or anti-conservative mid-p option.

## Allowed labels and blocks

Input reference indices refer to the supplied fixed pooled row order. Without blocks, enumerate every reference subset of the observed reference size. With aligned prespecified block labels, each subset preserves the observed reference count inside each block. A block containing only one side is frozen. If all blocks are frozen there is one allocation and p=1: uninformative, not evidence of equality.

A blocked orbit and unrestricted orbit can use the same observed data and statistic while answering different null questions. Blocks are not learned or chosen for a preferred p-value. This restriction does not fix arbitrary repeated-user dependence, time correlation, hidden confounding, preprocessing differences or selection of the family after labels were inspected. Repeated windows require their own error-control plan. There is no continuous-monitoring guarantee.

## Finite PSD evidence

Each supplied matrix must be exactly symmetric and PSD. The native generator builds a unit-lower-triangular `L` and nonnegative diagonal `D` such that `K=L diag(D) L^T`. At each exact rational pivot, a negative Schur pivot is rejected; a zero pivot requires every residual entry below it to be zero. Singular and rank-zero matrices are supported, without a jitter or a fabricated positive pivot. Zero-pivot columns below the diagonal are canonically zero.

The native checker does not import either generator or its LDL routine. It verifies the certificate shape, canonical rational entries, unit lower triangle, nonnegative diagonal, zero-column convention, stated rank and every exact reconstructed matrix entry. PSD tests independently compare acceptance against all principal-minor signs on 729 symmetric small matrices, including negative and zero-pivot cases. General indefinite LDL factorizations may require pivoting; they are not accepted as PSD evidence here.

A finite PSD matrix certifies a finite inner-product geometry only. It does not show a characteristic population kernel, unbiasedness after adaptive kernel selection, label-independent preprocessing or a valid sampling design. The tool cannot authenticate those assumptions from a matrix. User-supplied row names are labels, not actual observations of a person or data source.

## Schema, limits and commands

```
python -m driftlens.orbit_cli audit --input examples/orbit-selection.json --output report.json
python -m driftlens.orbit_cli verify --input report.json
# Both commands also accept --input - for stdin.
# Installed console script: driftlens-orbit audit / driftlens-orbit verify
```

Request keys are exactly `row_names`, `reference`, `blocks`, `kernels`, `alpha`. Each kernel has exactly `name` and `gram`.

- Four through ten pooled rows, with at least two on each side. `reference` is a sorted list of distinct zero-based integer indices. Boolean, float and exponent-form index tokens are rejected.
- One through six kernels, aligned square matrices. Each Gram entry and alpha is a canonical reduced rational string, such as `-2/3`, `0` or `7`. Decimal strings, exponent strings, `-0`, nonreduced fractions, leading zeros and zero denominators are refused. Absolute reduced numerator and denominator are at most 1,000,000. Alpha is strictly between zero and one.
- At most `choose(10,5)=252` complete allocations. Block restrictions may reduce the count. Scores and rank tables are always fully enumerated, never silently sampled.
- Distinct nonempty row/kernel names, at most forty Unicode codepoints, no lone surrogates. Blocks are null or aligned nonempty bounded strings. No normalization or inferred grouping occurs.
- Audit input limit 65,536 UTF-8 bytes; verify input and serialized output limit 8 MiB. Duplicate keys, nonfinite JSON tokens, invalid UTF-8, oversized/deep malformed JSON, unknown keys and missing keys are refused. Report-certificate rational strings are bounded to 2,500 characters and 4,096 bits per numerator/denominator in the checker.
- Errors exit2 with stderr and no partial stdout report. `--output` uses a temporary file in the destination directory followed by atomic replacement; a failed command preserves an existing destination and removes its temporary file. Input and output paths cannot resolve to the same file. Concurrent readers may observe either complete version. No directory transaction, lock, disk-failure durability or adversarial filesystem-race guarantee is claimed.

The library validates the same structural/numeric limits but has no serialized-byte limit because it accepts an already materialized object. The browser parser enforces its request-byte limit and rejects duplicate keys/noninteger reference tokens before its independent solver runs. Downloads are generated locally from the current complete report. Editing clears results and disables download immediately.

## Report and independent replay

The report retains a deep copy of the request; exact finite PSD certificates; complete lexicographic reference/complement sets; every exact kernel statistic, tied upper-tail count and raw fraction; every selected-kernel tie; separate Holm values; complete global-tail allocation indices; global fractions and decisions; and exact whole-orbit rejection counts/rates for raw-minimum, global and Holm-any comparisons.

The native verifier independently visits all binary masks matching side and per-block counts, computes every ordered-pair contribution directly and rebuilds every score/rank/global tail. The generator instead constructs a product of within-block combinations and uses grouped within/cross sums. It shares strict request validation and scope wording, not allocation-generation or score-generation functions. The browser uses binary masks, ordered-pair sums and a Schur-update PSD calculation rather than the native generator's LDL recurrence. Native replay accepts only exact primitive types and exact rational strings; floats are not substituted for integer metadata. Replaying is consistency checking, not authentication, a signature or provenance proof. A fully self-consistent changed request is a different experiment, not necessarily a forgery.

Complexity is bounded but not scalable: full matrices, scores, tied ranks and tails are materialized. Native and browser work is roughly `O(K*Nalloc*Nrows² + K*Nalloc² + K*Nrows³)`, plus exact rational bit arithmetic and full output. This is a tiny finite teaching and audit tool, not a million-vector detector or streaming algorithm.

## Reproduction and original controls

```
python -m unittest discover -s tests -v
node --test tests/*.cjs
PYTHONPATH=. python examples/make_orbit_controls.py
git diff --exit-code
bash examples/check_orbit_wheel.sh
```

As checked in the private candidate before publication:

- 58 Python tests, retaining all 45 original tests. Fifty independently replayed rational orbits cover 1,410 allocations and every alpha grid threshold. 729 symmetric matrices match the independent all-principal-minor PSD acceptance oracle. 753 individually changed non-request primitive leaves are rejected. Maximum ten-row/six-kernel/252-allocation request and report sizes, group swap, row/kernel reorder, tied/duplicate/constant/zero/singular/full-rank controls and strict CLI behavior are checked.
- Four Node tests. Ninety-five complete independent BigInt/native Fraction reports match exactly: 3,289 allocations, 12,600 kernel scores and 179,225 global-tail memberships. All actual JavaScript serializations are accepted by the native independent checker. Tests include rational singular/full-rank matrices, blocks, varied unequal sides, six kernels and ten rows.
- Actual Chrome checks all eleven original controls at 1280/390/320: 33 complete reports, every 372 allocation view, 717 per-kernel display rows and 51 matrix/certificate views. Eleven trusted completed evidence downloads are read and independently native-verified. A custom ten-row/six-kernel request with Unicode names checks the 252-allocation boundary. Trusted keyboard input clears old evidence; duplicate/noncanonical/noninteger requests fail without stale output.
- Fresh isolated wheel reproduction verifies all fourteen module origins, eleven controls through both installed entrypoints, twenty-two forged report CLI rejections, the original seeded demo and original exact/Monte Carlo/family compare modes outside the checkout. No runtime third-party dependencies.

These finite tests are not a proof that every implementation input is bug-free, nor are synthetic orbits a practical power/type-I-error benchmark. Hosted feature/PR/main checks, deployed source-byte/behavior/pixel checks and actual portfolio/profile render checks are separate release gates and must not be inferred merely from this local checklist.

Eleven original controls include the selected multiplicity contrast, duplicate family, all ties, rank-zero PSD, within-site identical rows, unrestricted/restricted same mixture, completely frozen labels, unequal side sizes, nonlinear feature-map views and a full-rank Gram matrix.

## Attribution and source boundaries

Standard methods, original small implementation and original supplied controls. No copied source implementation or new statistical theorem is claimed.

1. [Gretton et al., A Kernel Two-Sample Test, JMLR2012](https://www.jmlr.org/papers/volume13/gretton12a/gretton12a.pdf): diagonal-excluding MMD statistic and kernel-test assumptions. A supplied finite Gram is not by itself a characteristic population kernel.
2. [Ge, Dudoit, Speed, Berkeley report633](https://statistics.berkeley.edu/sites/default/files/tech-reports/633.pdf), Sections4.4/4.4.1: shared permutation matrices, upper-tail ranks and tied ranks for minimum-p methods. This mode implements only one finite global minimum-p test, not their full step-down multiple-testing procedure.
3. [Westfall and Troendle, Multiple Testing with Minimal Assumptions](https://pmc.ncbi.nlm.nih.gov/articles/PMC3117234/): joint exchangeability and limitations of multiple-testing claims. No assumption-free license is inferred from its title.
4. [multtest minP documentation](https://web.mit.edu/r/current/arch/i386_linux26/lib/R/library/multtest/html/mt.maxT.html): shared resampling and tie handling; lower tie variants can be invalid. No implementation-equivalence claim.
5. [Berkeley StatLab stratified permutation documentation](https://statlab.github.io/permute/api/stratified.html): count-preserving within-block allocations and their conditional exchangeability question, not arbitrary correlation repair.
6. [Lanckriet et al., JMLR2004](http://www.jmlr.org/papers/volume5/lanckriet04a/lanckriet04a.pdf), Proposition2: symmetric PSD matrices and finite Gram geometry. This tool does not learn a kernel or implement their SDP.
7. [Schrab et al., MMD Aggregated Two-Sample Test, JMLR2023](https://jmlr.org/papers/v24/21-1289.html): established kernel aggregation context. This project is not MMDAgg and claims none of its adaptive minimax or bandwidth-calibration theory.
8. [Phipson and Smyth, permutation p-values, 2010/corrected2011](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf): distinguish randomly drawn permutation p-values from complete finite enumeration. There is no Monte Carlo path in the new mode.
9. [UCLA symmetric factorization notes](https://www.math.ucla.edu/~njhu/notes/nla/lin-direct/symfact.pdf): LDL and Schur-complement background. Singular PSD checks are independently verified rather than treating any unpivoted indefinite factorization as positive evidence.

MIT licensed. No paid API, model download, private dataset, causal effect, learned grouping, automatic scientific conclusion, deployment approval or sequential monitoring.
