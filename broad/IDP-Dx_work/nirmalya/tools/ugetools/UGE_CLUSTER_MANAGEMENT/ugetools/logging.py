# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2016 by the
# Broad Institute/Massachusetts Institute of Technology. All rights are
# reserved.

# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for its
# use, misuse, or functionality.
import sys
import time

def log_message(message, level=0, verbosity=0, file=sys.stderr):
    if level <= verbosity:
        print >> file, "\t".join([time.ctime(), message])