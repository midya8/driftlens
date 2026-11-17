# driftlens

Checks whether two batches of embeddings look like they came from the same distribution. It's a kernel two-sample test (RBF MMD) with a permutation p-value, written in plain Python so I could see every step instead of trusting a library call.

The use case I had in mind: you have embeddings from last month and from this week, and you want a quick "did something move?" check before digging in. It won't tell you why things moved or whether your model got worse.

```sh
python -m driftlens compare --input examples/embeddings.json --bandwidth 1 --permutations 499
python -m driftlens compare --input examples/embeddings.json --bandwidths 1 2 4 --exact
python -m driftlens demo
```

Or from Python:

```python
from driftlens import compare
report = compare({"reference": [[0], [1]], "current": [[2], [3]]}, bandwidth=1, exact=True)
print(report["p_value"])
```

Input is JSON with a `reference` and a `current` list of rows (at least two each, same dimension). You pick the bandwidth yourself; nothing gets tuned or standardized for you. With several bandwidths you get a Holm adjustment across them.

## What the demo shows

With 30 rows of 4-d data per side, bandwidth 2 and 499 permutations, two samples from the same generator give MMD² of about -0.005 and p = 0.55. Shift the current rows by +2 in every coordinate and it jumps to 0.63 with p = 0.002, which is just the floor you get from 499 permutations.

Small samples have hard limits. Two points against two points, even if they're far apart, can't get below p = 1/3 because there are only six ways to split them. And a high p-value doesn't mean the distributions are the same.

## Extras

`python -m driftlens.orbit_cli` is a separate mode for tiny inputs where you supply the kernel matrices yourself and it enumerates every possible split. There's a browser version in `docs/`. REFERENCE.md and ORBIT.md have the long-form details.

Run the tests with `python -m unittest discover -s tests`.
