"""Run with an installed wheel from a working directory outside the checkout."""
import sys,json,subprocess,importlib,pkgutil,pathlib
import driftlens
from driftlens.orbit import audit_orbit
from driftlens.orbit_controls import definitions
from driftlens.orbit_verify import verify_orbit
root=pathlib.Path(driftlens.__file__).resolve().parent
origins=[]
for module in [driftlens]+[importlib.import_module('driftlens.'+x.name)for x in pkgutil.iter_modules(driftlens.__path__)if x.name!='__main__']:
    path=pathlib.Path(module.__file__).resolve();assert root in path.parents or path==root/'__init__.py';assert 'site-packages'in str(path);origins.append(str(path))
spec=importlib.util.find_spec('driftlens.__main__');assert 'site-packages'in spec.origin;origins.append(spec.origin)
count=0
for c in definitions():
    data=json.dumps(c['request']).encode();expected=audit_orbit(c['request']);assert verify_orbit(expected)
    for command in [[sys.executable,'-m','driftlens.orbit_cli'],['driftlens-orbit']]:
        a=subprocess.run(command+['audit','--input','-'],input=data,capture_output=True,check=True);assert json.loads(a.stdout)==expected
        b=subprocess.run(command+['verify','--input','-'],input=a.stdout,capture_output=True,check=True);assert json.loads(b.stdout)=={'valid':True,'schema':'driftlens-rational-orbit-v1'}
        forged=json.loads(a.stdout);forged['allocation_count']+=1;b=subprocess.run(command+['verify','--input','-'],input=json.dumps(forged).encode(),capture_output=True);assert b.returncode==2 and not b.stdout;count+=1
# Existing floating APIs, family, witness, blocks and seeded demo must remain available.
from driftlens import compare,compare_family
from driftlens.demo import run_demo
payload={'reference':[[0],[1]],'current':[[2],[3]]}
a=compare(payload,bandwidth=1,exact=True,witness=2);assert a['total_label_allocations']==6 and a['witness']
b=compare_family(payload,[1,2],exact=True);assert b
old=run_demo();out=subprocess.check_output([sys.executable,'-m','driftlens','demo']);assert json.loads(out)==old
for flags in [['--exact'],['--permutations','7'],['--exact','--bandwidths','1','2']]:
 command=[sys.executable,'-m','driftlens','compare','--input','-']
 if '--bandwidths'not in flags:command+=['--bandwidth','1']
 out=subprocess.run(command+flags,input=json.dumps(payload).encode(),capture_output=True,check=True);assert json.loads(out.stdout)
print(json.dumps({'installed_modules':len(origins),'module_origins':origins,'controls':len(definitions()),'entrypoints':2,'forged_cli_rejections':count}))
