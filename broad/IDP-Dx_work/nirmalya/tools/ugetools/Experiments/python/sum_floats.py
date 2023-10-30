#!/usr/bin/env python
#  The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2017 by the
# Broad Institute/Massachusetts Institute of Technology. All rights are
# reserved.

# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for its
# use, misuse, or functionality.
import sys


def main(args=None):
    sum = 0.0
    for strLine in sys.stdin:
        sum += float(strLine.rstrip("\n"))

    print str(sum)
    print >> sys.stderr, "Done sum_floats"

if __name__ == "__main__":
    sys.exit(main())

