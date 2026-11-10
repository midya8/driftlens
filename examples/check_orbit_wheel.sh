#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT HUP INT TERM
python -m pip wheel --no-deps "$ROOT" -w "$TMP/wheels"
python -m venv "$TMP/venv"
"$TMP/venv/bin/python" -m pip install --no-deps "$TMP"/wheels/driftlens-*.whl
mkdir "$TMP/outside"
cp "$ROOT/examples/verify_orbit_installed.py" "$TMP/outside/verify.py"
cd "$TMP/outside"
PATH="$TMP/venv/bin:$PATH" "$TMP/venv/bin/python" -I verify.py
"$TMP/venv/bin/python" -m driftlens demo > old-demo.json
"$TMP/venv/bin/python" -c 'import json; r=json.load(open("old-demo.json")); assert r'
echo 'all old demo and new exact orbit commands run from the installed wheel outside checkout'
