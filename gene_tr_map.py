#!/usr/bin/env python

""" Transcript to gene mapper

Usage: 
    gene_tr_map.py -i <infile> -o <outfile>
    gene_tr_map.py (-h|--help)

    -i --infile input file
    -o --outfile output file

"""

from docopt import docopt
import re

def map_single_line(line, outf):
    if "gene_symbol" in line: 
        # Old style
        #regex = ">(\\S+?)\\.\\d+\\s+(\\S+).*?gene:(\\S+).*gene_biotype:(\\S+).*gene_symbol:(\\S+)"
        regex = ">(\\S+?)\\.\\d+\\s+(\\S+)\\s+(\\S+)\\s+gene:(\\S+).*gene_biotype:(\\S+).*gene_symbol:(\\S+)"
        p = re.compile(regex)
        m = p.match(line)
        transcript = m.group(1)
        gene_type = m.group(2)
        gene_pos = m.group(3)
        geneid = m.group(4)
        gene_bio = m.group(5)
        gene_sym = m.group(6)
        outf.write(transcript + "\t" + geneid + ":" + gene_sym + ":" + gene_type+ ":" + gene_bio + "\n")
    else:
        regex = ">(\\S+?)\\..*\\((\\S+)\\)"
        p = re.compile(regex)
        m = p.match(line)
        transcript = m.group(1)
        desc = m.group(2)
        outf.write(transcript + "\t" + transcript + ":" + desc + "\n")

def tr_gene_mapper(infile, outfile):

    with open(outfile, "w") as outf:
        with open(infile) as inf:
            for line in inf:
                if line.startswith(">"):
                    map_single_line(line, outf)
                    


if __name__ == '__main__':
    args = docopt(__doc__, version = 'Transcript to gene mapper 0.1')
    print(args)
    infile = args["--infile"]
    outfile = args["--outfile"]

    print("infile: " + infile)
    print("outfile: " + outfile)
    tr_gene_mapper(infile, outfile)
