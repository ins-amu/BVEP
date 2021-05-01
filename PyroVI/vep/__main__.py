import sys
from .cli import main, build_parser

parser = build_parser()
args = parser.parse_args(sys.argv[1:])
main(args)