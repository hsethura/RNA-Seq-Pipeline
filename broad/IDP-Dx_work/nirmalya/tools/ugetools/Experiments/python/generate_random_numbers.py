#!/usr/bin/env python
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2017 by the
# Broad Institute/Massachusetts Institute of Technology. All rights are
# reserved.

# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for its
# use, misuse, or functionality.
"""
Generate random numbers in the range [0.0, 1.0) and write them to stdout
"""
import argparse
import random
import sys
import time

def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--count", "-c", type=int, default=None,
                        help="Generate this many numbers.  Default: no limit")
    parser.add_argument("--duration", "-d", type=float, default=None,
                        help="Run for this many minutes.  Default: no limit")
    parser.add_argument("--seed", "-s", type=int, default=None, help="RNG seed")

    options = parser.parse_args(args)

    if options.seed is not None:
        random.seed(options.seed)

    count = 0
    start_time = time.time()

    while True:
        if options.count is not None and count >= options.count:
            break
        if options.duration is not None and start_time + options.duration * 60 <= time.time():
            break
        print "%f" % random.random()
        count += 1

    print >> sys.stderr, "Done generating random numbers"

if __name__ == "__main__":
    sys.exit(main())
    
