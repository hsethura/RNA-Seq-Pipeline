#!/usr/bin/env python

""" Trim the header of the transcript fasta files

Usage: 
    trim_header.py -i <infile> -o <outfile>
    trim_header.py (-h|--help)

    -i --infile input file
    -o --outfile output file

"""

from docopt import docopt
import re

def trim_header(infile, outfile):

    with open(outfile, "w") as outf:
        with open(infile) as inf:
            for line in inf:
                if line.startswith(">"):
                    regex = "\\.|\\s+"
                    parts = re.split(regex, line)
                    part = parts[0]
                    outf.write(part + "\n")
                else:
                    outf.write(line)



if __name__ == '__main__':
    args = docopt(__doc__, version = 'Trim Header 0.1')
    print(args)
    infile = args["--infile"]
    outfile = args["--outfile"]

    print("infile: " + infile)
    print("outfile: " + outfile)
    trim_header(infile, outfile)
