"""Run the complete Python suite and retained CLI reports from a built wheel."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import venv

ROOT = Path(__file__).resolve().parents[1]
ENV = os.environ.copy()
ENV.pop('PYTHONPATH', None)
MODULES = ['driftlens', 'driftlens.__main__', 'driftlens.allocations',
           'driftlens.cli', 'driftlens.core', 'driftlens.demo', 'driftlens.family',
           'driftlens.orbit', 'driftlens.orbit_cli', 'driftlens.orbit_controls',
           'driftlens.orbit_input', 'driftlens.orbit_verify', 'driftlens.scoring',
           'driftlens.witness']


def call(args, cwd):
    return subprocess.check_output(args, cwd=cwd, env=ENV, text=True)


with tempfile.TemporaryDirectory() as directory:
    work = Path(directory)
    wheels = work / 'wheels'
    wheels.mkdir()
    call([sys.executable, '-m', 'pip', 'wheel', '--no-deps', '-w', str(wheels), str(ROOT)], work)
    venv.create(work / 'env', with_pip=True)
    py = str(work / 'env/bin/python')
    call([py, '-m', 'pip', 'install', '--no-deps', str(next(wheels.glob('*.whl')))], work)
    outside = work / 'outside'
    outside.mkdir()
    for name in ['tests', 'examples']:
        shutil.copytree(ROOT / name, outside / name)
    suite = '''
import importlib,importlib.util,json,pathlib,sys,unittest
names=json.loads(sys.argv[1])
for name in names:
    if name!='driftlens.__main__':importlib.import_module(name)
def check():
    for name in names:
        origin=pathlib.Path(importlib.util.find_spec(name).origin).resolve()
        assert 'site-packages' in origin.parts,(name,str(origin))
check()
suite=unittest.defaultTestLoader.discover(str(pathlib.Path.cwd()/'tests'))
result=unittest.TextTestRunner(verbosity=2).run(suite)
assert result.wasSuccessful()
check()
print(json.dumps(dict(tests=result.testsRun,installed_origins=len(names),checked_before_and_after=True)))
'''
    print(call([py, '-I', '-c', suite, json.dumps(MODULES)], outside))
    # Preserve the original eleven-control, two-entrypoint and forgery checks.
    installed_env = ENV.copy()
    installed_env['PATH'] = str(work / 'env/bin') + os.pathsep + installed_env['PATH']
    print(subprocess.check_output([py, '-I', str(outside / 'examples/verify_orbit_installed.py')],
                                  cwd=outside, env=installed_env, text=True))
    jobs = [(['demo'], None)]
    payload = {'reference': [[0],[1]], 'current': [[2],[3]]}
    for flags in [['--exact'], ['--permutations','7'],
                  ['--exact','--bandwidths','1','2'], ['--exact','--witness','2']]:
        args = ['compare','--input','-']
        if '--bandwidths' not in flags:
            args += ['--bandwidth','1']
        jobs.append((args+flags, json.dumps(payload)))
    blocked = {**payload, 'blocks': {'reference':['a','b'], 'current':['a','b']}}
    jobs.append((['compare','--input','-','--bandwidth','1','--exact'], json.dumps(blocked)))
    for args, data in jobs:
        expected = subprocess.check_output([sys.executable,'-m','driftlens']+args,
                                          cwd=ROOT,env=ENV,input=data,text=True)
        actual = subprocess.check_output([py,'-I','-m','driftlens']+args,
                                        cwd=outside,env=ENV,input=data,text=True)
        assert actual == expected,args
    print(json.dumps(dict(complete_floating_cli_reports=len(jobs),outside_checkout=True,
                          original_orbit_checks_preserved=True)))
