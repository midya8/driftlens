# DriftLens

Compare two supplied embedding samples with an RBF kernel permutation test.
Inspect the statistic, allowed label allocations and p-value rather than
turning a shift into a claim about model accuracy or its cause.

Python 3.10+, no runtime dependencies, model downloads, paid APIs or private
datasets. A separate tiny supplied-Gram inspector uses exact rational arithmetic
in Python and BigInt JavaScript. It does not make the floating RBF API rational.

## Pick the question first

| Input and question | Entry point | Result |
| --- | --- | --- |
| Embedding rows, one prespecified bandwidth | `compare` / `driftlens.compare` | Unbiased MMD squared and exact-enumeration or Monte Carlo p-value |
| Embedding rows, a prespecified bandwidth family | `compare --bandwidths` / `compare_family` | Per-kernel Holm adjustment |
| Supplied rational PSD Gram matrices, a fixed family | `driftlens.orbit_cli audit` / browser orbit inspector | Every allocation and tail rank, recalibrated global p-value |

Optional blocks restrict label allocations; they do not establish
exchangeability. Witness scores describe rows, not drift causes. Choose
bandwidths, preprocessing, blocks and the tested family before examining
label-specific results. Neither mode provides always-valid repeated monitoring.

## Run a comparison

```sh
python -m driftlens compare --input examples/embeddings.json --bandwidth 1 --permutations 499
python -m driftlens compare --input examples/embeddings.json --bandwidths 1 2 4 --exact
python -m driftlens demo
bash examples/demo.sh
```

```python
from driftlens import compare
report = compare({"reference": [[0], [1]], "current": [[2], [3]]},
                 bandwidth=1, exact=True, witness=2)
```

Each group needs at least two numeric rows with the same nonzero dimension.
The CLI requires a bandwidth or family; it does not standardize or tune inputs.
`--input -` reads stdin. Invalid input exits 2 without partial JSON. Seeds use a
local RNG. [Schema, statistic, bounds and full API reference](REFERENCE.md).

## Small examples can expose limits

The retained seed-13 development demo uses 30x4 samples, bandwidth 2 and 499
permutations. The same-generator arm has MMD squared -.004871 and p=.550.
Shifting the same current rows by +2 per coordinate gives .634539 and p=.002.
That p-value is the Monte Carlo floor 1/500, not a measurement below that floor.
A visibly separated 2x1 exact example still has p=1/3: there are only six label
allocations. Non-rejection does not mean the distributions are equal.

Other retained controls show different questions rather than a preferred
answer: a selected mixture gives unrestricted p=.014 but prespecified
within-stratum p=1; the statistic is unchanged. A six-row fixed-Gram family has
raw-minimum rejection at 8/20 allocations, recalibrated global rejection at
4/20 and Holm-any at 0/20, using alpha=1/5. These are selected development
fixtures and finite conditional-orbit fractions, not real-data power or
false-positive estimates. [All settings and qualifications](REFERENCE.md).

## Inspect a finite kernel-family orbit

```sh
python -m driftlens.orbit_cli audit --input examples/orbit-selection.json --output report.json
python -m driftlens.orbit_cli verify --input report.json
bash examples/check_orbit_wheel.sh
```

This additive mode retains every exact rational MMD score, per-kernel tied tail,
PSD certificate and recalibrated minimum-p global tail. The independent checker
replays the report; the browser computes its own BigInt result. Limits are
4..10 pooled rows, 1..6 supplied PSD matrices and at most 252 allocations.
Geometry must be fixed independently of labels and the allowed uniform null
must be justified separately. A PSD certificate proves neither assumption.
[Formulas, editable inspector, controls and sources](ORBIT.md).

Positive per-kernel scale/constant-offset and nonnegative diagonal tests retain
all orbit decisions and tails while checking changed scores and PSD witnesses.
Those are finite consistency tests, not evidence of real-data accuracy.

## Check the implementation

```sh
python -m unittest discover -s tests -v
node --test --test-concurrency=1 tests/test_orbit_engine.cjs tests/test_orbit_browser.cjs
python tools/check_installed.py
```

The Python suite has 60 tests. The four Node tests include independent complete
BigInt/native report comparisons and actual headless Chrome checks of controls,
all allocation tables, PSD witnesses, native-verified downloads, invalid-state
clearing and narrow layouts. Node 22 and Google Chrome are needed for those
browser checks, not for the Python package.

The wheel check runs every Python test outside the source checkout, verifies
all fourteen module origins before and after, preserves the original orbit
controls/two-entrypoint/forgery checks and compares six complete floating CLI
reports, including family, witness and blocked cases. CI has separate Python
3.10..3.13, demo, installed-package and browser jobs. Source tests alone would
not prove that the wheel or browser engine works.

## What the result does not say

A rejection does not identify a drift cause or show that a deployed model got
worse. Correlated windows, repeated users or group-dependent preprocessing can
invalidate label exchangeability. A large sample can detect a tiny irrelevant
change. Witness scores are in-sample descriptions, not anomaly p-values.
Repeated tests require a separate error-control plan.

The cached floating implementation is bounded, not a streaming or million-row
system. Exact finite enumeration still uses the supplied null assumptions;
"exact" is not a real-data validity guarantee. [Full work counts, numerical
limits and family/block/witness semantics](REFERENCE.md).

Standard methods, original small implementations, not a new statistical test.
[DESIGN.md](DESIGN.md) and [ORBIT.md](ORBIT.md) retain sources and derivations;
[the full reference](REFERENCE.md#verification-and-attribution) includes the
original Gretton, Phipson-Smyth, Holm and stratified-permutation links. MIT
licensed.

### Browser download checks

The browser test now waits for the expected export's complete byte count and no
partial-download file before reading it. It still compares the full parsed report,
checks the directory has exactly the expected file and sends that file through
the native verifier. These semantic checks run once, without retrying assertions.
The transport wait is bounded at 30 seconds. Separate helper tests reject absent,
empty, partial, temporary and oversized downloads, then check completed files.
Chrome runs in a test-owned process group that is stopped before profile cleanup.

This addresses the transport failure observed in the PR browser job, where no
JSON file appeared within the old five-second wait and Chrome cleanup also raced.
Local passing checks do not establish that a hosted timing failure is fixed;
the updated hosted browser job must pass before this batch can be merged.
