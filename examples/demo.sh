#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python -m driftlens demo
python - <<'PY'
from driftlens.demo import run_demo
r=run_demo()
assert not r['same_generator']['reject_at_alpha']
assert r['mean_shift']['reject_at_alpha']
assert r['exact_tiny']['p_value']==1/3
assert r['shift_bandwidth_family']['reject_any']
assert all(test['holm_p_value'] >= test['p_value'] for test in r['shift_bandwidth_family']['tests'])
try:
    assert r['same_generator']['p_value'] <= .05
except AssertionError:
    print('same-generator rejection gate fails as expected')
else:
    raise SystemExit('ERROR: invalid rejection gate passed')
PY
