#!/usr/bin/env python

""" Transcript to gene mapper

Usage: 
    gene_tr_map.py -i <infile> --out_map <out_map_file> --out_fasta <out_fasta_file> 
    gene_tr_map.py (-h|--help)

    -i --infile input file
    --out_map <out_map_file>
    --out_fasta <out_fasta_file>

"""

from docopt import docopt
import re

def map_single_line(line, outf):
    retval = ""
    if "gene_symbol" in line: 
        regex = ">(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+gene:(\\S+)\\s+gene_biotype:(\\S+)\\s+transcript_biotype:(\\S+)\\s+gene_symbol:(\\S+)"
        p = re.compile(regex)
        m = p.match(line)
        transcript = m.group(1)
        gene_type = m.group(2)
        trans_pos_base = m.group(3)
        geneid = m.group(4)
        gene_bio = m.group(5)
        trans_bio = m.group(6)
        gene_sym = m.group(7)
        # extract part of the gene_pos
        trans_pos = re.sub('^\\S+?:\\S+?:', '', trans_pos_base)        
        # We have to extract the last part of trans_pos to get the orientation of the gene
        trans_pos_parts = trans_pos.split(':')
        gene_orient = trans_pos_parts[-1]        
        outf.write(transcript + ":" +  gene_type + ":" + trans_bio + ":" + trans_pos +  
            "\t" + geneid + ":" + gene_type + ":" + gene_sym + ":" + gene_bio + ":" + gene_orient + "\n")
        retval = ">" + transcript + ":" +  gene_type + ":" + trans_bio + ":" + trans_pos + "\n"
    else:
        regex = ">(\\S+?)\\..*\\((\\S+)\\)"
        print(line)
        p = re.compile(regex)
        m = p.match(line)
        transcript = m.group(1)
        desc = m.group(2)
        outf.write(transcript + ":" + desc + "\t" + transcript + ":" + desc + "\n")
        retval = ">" + transcript + ":" + desc + "\n"
    return retval

def tr_gene_mapper(infile, outfile, outfile_fasta):
    ip_count = 0
    entry_count = 0
    with open(outfile, "w") as outf, open(outfile_fasta, "w") as out_fasta:
        with open(infile) as inf:
            retval = ""
            for line in inf:
                if line.startswith(">"):
                    retval = map_single_line(line, outf)
                else:
                    retval = line
                out_fasta.write(retval)


if __name__ == '__main__':
    args = docopt(__doc__, version = 'Transcript to gene mapper 0.1')
    #print(args)
    infile = args["--infile"]
    out_map = args["--out_map"]
    out_fasta = args["--out_fasta"]

    print("infile: " + infile)
    print("out_map: " + out_map)
    print("out_fasta: " + out_fasta)
    tr_gene_mapper(infile, out_map, out_fasta)
