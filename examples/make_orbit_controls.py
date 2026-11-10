import json
from pathlib import Path
from driftlens.orbit_controls import definitions
from driftlens.orbit import audit_orbit
from driftlens.orbit_verify import verify_orbit
controls=definitions()
for c in controls:assert verify_orbit(audit_orbit(c['request']))
Path(__file__).resolve().parents[1].joinpath('docs/orbit-data.js').write_text('globalThis.ORBIT_CONTROLS='+json.dumps(controls,ensure_ascii=False,separators=(',',':'))+';\n')
