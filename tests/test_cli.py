import io
import json
import unittest
from contextlib import redirect_stdout,redirect_stderr
from unittest.mock import patch
from driftlens.cli import main


class TestCLI(unittest.TestCase):
    def test_stdin_json_exact(self):
        out=io.StringIO()
        d={'reference':[[0],[1]],'current':[[2],[3]]}
        with patch('sys.stdin',io.StringIO(json.dumps(d))),redirect_stdout(out):
            self.assertEqual(main(['compare','--input','-','--bandwidth','1','--exact']),0)
        self.assertEqual(json.loads(out.getvalue())['mode'],'exact')

    def test_invalid_input_no_partial_json(self):
        for body in ['{','{}']:
            out,err=io.StringIO(),io.StringIO()
            with patch('sys.stdin',io.StringIO(body)),redirect_stdout(out),redirect_stderr(err):
                self.assertEqual(main(['compare','--input','-','--bandwidth','1']),2)
            self.assertEqual(out.getvalue(),'');self.assertIn('input error',err.getvalue())
