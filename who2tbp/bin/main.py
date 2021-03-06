"""
Convert WHO to TBprofiler format
"""

import argparse
import logging
import os
import sys

from who2tbp.lib.parse_excel_data import parse_file
from who2tbp.lib.gen_tbdb_csv import who2tbd

tool_name = os.path.basename(sys.argv[0])

logging.basicConfig(level=logging.DEBUG,
                    format=f"[{tool_name} – %(asctime)s]: %(message)s")
logger = logging


def arg_parser() -> argparse.Namespace:
    args = argparse.ArgumentParser("Convert WHO Excel sheet with MTB mutations to TBProfiler database format",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    args.add_argument("INFILE", help="The WHO Excel sheet", type=argparse.FileType("rb"), default=sys.stdin)
    args.add_argument("-f", "--filter", help="Limit to single category", choices=['assoc_resistance',
                                                                                  'no_assoc',
                                                                                  'assoc_resistance_interim',
                                                                                  'no_assoc_interim',
                                                                                  'combo',
                                                                                  'uncert_signif',
                                                                                  'all'],
                      default="assoc_resistance")
    args.add_argument("-o", "--outfile", default=sys.stdout, type=argparse.FileType('w'))
    return args.parse_args()


def main() -> int:
    args = arg_parser()
    logger.info(f"Welcome to {tool_name}")
    data = parse_file(args.INFILE, this_filter=args.filter)
    _ = who2tbd(data, args.outfile)
    return 0


if __name__ == "__main__":
    main()
