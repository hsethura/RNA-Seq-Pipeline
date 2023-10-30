#!/usr/bin/env python
#  The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2017 by the
# Broad Institute/Massachusetts Institute of Technology. All rights are
# reserved.

# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for its
# use, misuse, or functionality.
"""
Read groups of floating point numbers from stdin, sort each group, then write to stdout
"""
import argparse
import sys

def sort_and_print(lst):
    lst.sort()
    for f in lst: print str(f)

def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--chunk", "-c", type=int, default=100000,
                        help="Generate this many numbers.  Default: %(default)s")

    options = parser.parse_args(args)

    lst = []
    for strLine in sys.stdin:
        if len(lst) >= options.chunk:
            sort_and_print(lst)
            lst = []
        lst.append(float(strLine.rstrip("\n")))

    sort_and_print(lst)
    print >> sys.stderr, "Done sort_float_chunks"

if __name__ == "__main__":
    sys.exit(main())

