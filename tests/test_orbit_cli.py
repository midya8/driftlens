import unittest,subprocess,sys,tempfile,json,os
from pathlib import Path
from driftlens.orbit import audit_orbit
from test_orbit import request


class OrbitCLITests(unittest.TestCase):
    def run_cli(self,*args,data=b''):
        return subprocess.run([sys.executable,'-m','driftlens.orbit_cli',*args],input=data,capture_output=True)
    def test_stdin_roundtrip_and_atomic_file(self):
        r=request([[0,1,2,3]]);a=self.run_cli('audit','--input','-',data=json.dumps(r).encode());self.assertEqual(a.returncode,0,a.stderr);b=self.run_cli('verify','--input','-',data=a.stdout);self.assertEqual(b.returncode,0,b.stderr);self.assertTrue(json.loads(b.stdout)['valid'])
        with tempfile.TemporaryDirectory()as tmp:
            input=Path(tmp)/'input.json';output=Path(tmp)/'output.json';input.write_text(json.dumps(r));output.write_text('old');a=self.run_cli('audit','--input',str(input),'--output',str(output));self.assertEqual(a.returncode,0,a.stderr);self.assertEqual(json.loads(output.read_text()),audit_orbit(r));self.assertFalse(list(Path(tmp).glob('.driftlens-orbit-*')))
            old=output.read_bytes();input.write_text('bad');a=self.run_cli('audit','--input',str(input),'--output',str(output));self.assertEqual(a.returncode,2);self.assertEqual(a.stdout,b'');self.assertEqual(output.read_bytes(),old);self.assertFalse(list(Path(tmp).glob('.driftlens-orbit-*')))
            a=self.run_cli('audit','--input',str(input),'--output',str(input));self.assertEqual(a.returncode,2)
    def test_strict_parser_invalid_utf8_depth_size_and_no_output(self):
        for data in [b'{"alpha":"1/2","alpha":"1/3"}',b'NaN',b'Infinity',b'\xff',b'['*2000+b']'*2000,b' '*65537,b'null',b'{}']:
            a=self.run_cli('audit','--input','-',data=data);self.assertEqual(a.returncode,2,a.stdout);self.assertEqual(a.stdout,b'');self.assertTrue(a.stderr)
        r=request([[0,1,2,3]]);text=json.dumps(r).replace('"reference": [0, 1]','"reference": [0.0, 1]');a=self.run_cli('audit','--input','-',data=text.encode());self.assertEqual(a.returncode,2);self.assertEqual(a.stdout,b'')
    def test_forged_checker_cli_and_os_errors(self):
        a=audit_orbit(request([[0,1,2,3]]));a['allocations'][0]['global_p']='1';out=self.run_cli('verify','--input','-',data=json.dumps(a).encode());self.assertEqual(out.returncode,2);self.assertEqual(out.stdout,b'')
        a=self.run_cli('audit','--input','/does/not/exist');self.assertEqual(a.returncode,2);self.assertEqual(a.stdout,b'')
        with tempfile.TemporaryDirectory()as tmp:
            a=self.run_cli('audit','--input','-','--output',tmp,data=json.dumps(request([[0,1,2,3]])).encode());self.assertEqual(a.returncode,2);self.assertEqual(a.stdout,b'');self.assertFalse(list(Path(tmp).glob('.driftlens-orbit-*')))
