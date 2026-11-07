import argparse
import json
import sys
from .core import compare
from .demo import run_demo
from .family import compare_family


def main(argv=None):
    parser = argparse.ArgumentParser(prog='driftlens')
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('demo')
    c = sub.add_parser('compare')
    c.add_argument('--input', required=True)
    bandwidth = c.add_mutually_exclusive_group(required=True)
    bandwidth.add_argument('--bandwidth', type=float)
    bandwidth.add_argument('--bandwidths', type=float, nargs='+', help='prespecified family, Holm adjusted')
    c.add_argument('--permutations', type=int, default=499)
    c.add_argument('--seed', type=int, default=13)
    c.add_argument('--alpha', type=float, default=.05)
    c.add_argument('--exact', action='store_true')
    c.add_argument('--witness', type=int, default=0, help='top absolute witness rows per group, 0 disables')
    try:
        args = parser.parse_args(argv)
        if args.command == 'demo': result = run_demo()
        else:
            if args.input == '-':payload = json.load(sys.stdin)
            else:
                with open(args.input, encoding='utf-8') as source:payload = json.load(source)
            if args.bandwidths is not None:
                result = compare_family(payload, args.bandwidths, args.permutations, args.seed, args.alpha, args.exact, args.witness)
            else:
                result = compare(payload, args.bandwidth, args.permutations, args.seed, args.alpha, args.exact, args.witness)
        body = json.dumps(result, indent=2, sort_keys=True, allow_nan=False)
    except (ValueError, OSError, TypeError, OverflowError) as exc:
        print('driftlens input error: ' + str(exc), file=sys.stderr)
        return 2
    print(body)
    return 0
