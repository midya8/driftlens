"""Strict atomic command line for bounded rational orbit reports."""
import argparse,json,sys,os,tempfile
from .orbit import audit_orbit
from .orbit_verify import verify_orbit


def _pairs(pairs):
    out={}
    for k,v in pairs:
        if k in out:raise ValueError('duplicate JSON key')
        out[k]=v
    return out


def _constant(value):raise ValueError('nonfinite JSON token')


def main(argv=None):
    parser=argparse.ArgumentParser(prog='driftlens-orbit');parser.add_argument('command',choices=['audit','verify']);parser.add_argument('--input',required=True);parser.add_argument('--output')
    args=parser.parse_args(argv);path=None
    try:
        limit=65536 if args.command=='audit'else 8*1024*1024
        if args.output and args.input!='-'and os.path.realpath(args.input)==os.path.realpath(args.output):raise ValueError('input and output must differ')
        if args.input=='-':raw=sys.stdin.buffer.read(limit+1)
        else:
            with open(args.input,'rb')as f:raw=f.read(limit+1)
        if len(raw)>limit:raise ValueError('input size exceeds command limit')
        request=json.loads(raw.decode('utf-8'),object_pairs_hook=_pairs,parse_constant=_constant)
        if args.command=='audit':result=audit_orbit(request)
        else:verify_orbit(request);result={'valid':True,'schema':'driftlens-rational-orbit-v1'}
        output=(json.dumps(result,ensure_ascii=False,allow_nan=False,separators=(',',':'))+'\n').encode('utf-8')
        if len(output)>8*1024*1024:raise ValueError('report exceeds output size limit')
        if args.output:
            folder=os.path.dirname(os.path.abspath(args.output));fd,path=tempfile.mkstemp(prefix='.driftlens-orbit-',dir=folder)
            with os.fdopen(fd,'wb')as f:f.write(output);f.flush();os.fsync(f.fileno())
            os.replace(path,args.output);path=None
        else:sys.stdout.buffer.write(output)
        return 0
    except(ValueError,TypeError,OSError,UnicodeError,RecursionError,OverflowError)as exc:
        if path:
            try:os.unlink(path)
            except OSError:pass
        print('driftlens-orbit: '+str(exc),file=sys.stderr);return 2


if __name__=='__main__':raise SystemExit(main())
