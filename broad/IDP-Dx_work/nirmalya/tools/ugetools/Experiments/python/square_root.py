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
Read floating-point numbers from stdin, write sqrt(abs(number)) to stdout
"""
import math
import sys

def main(args=None):
    for strLine in sys.stdin:
        strLine = strLine.rstrip("\n")
        f = float(strLine)
        print "%f" % math.sqrt(abs(f))

    print >> sys.stderr, "Done square_root"


if __name__ == "__main__":
    sys.exit(main())

