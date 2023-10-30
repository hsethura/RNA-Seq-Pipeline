#!/usr/bin/env python

""" Filter unique genes

Usage: 
    remove_duplicates.py --in_fa <in_fasta> --out_fa <out_fasta> --rm_fa <remove_fasta_file> --in_map <in_map_file> --out_map <out_map_file> --rm_map <remove_map_file>
    remove_duplicates.py (-h|--help)

    --in_fa <in_fasta_file>
    --out_fa <out_fasta_file>
    --rm_fa <remove_fasta_file>
    --in_map <in_map_file>
    --out_map <out_map_file>
    --rm_map <remove_map_file>

"""

from docopt import docopt
import re




def dup_remover_fa(in_fa, out_fa, rm_fa, chr_lst):
    in_seq = False
    with open(out_fa, "w") as outf, open(rm_fa, "w") as rmf:
        with open(in_fa) as inf:
            for line0 in inf:
                line1 = line0.rstrip()
                if line1.startswith(">"):
                    parts = line1.split(":")
                    if len(parts) < 4:
                        in_seq = True
                        outf.write(line1 + "\n")
                    elif (parts[3]).upper() in chr_lst:
                        in_seq = True
                        outf.write(line1 + "\n")
                    else:
                        in_seq = False
                        rmf.write(line1 + "\n")
                else:
                    if in_seq:
                        outf.write(line1 + "\n")
                    else:
                       rmf.write(line1 + "\n")
                    

def dup_remover_map(in_map, out_map, rm_map, chr_lst):
    with open(out_map, "w") as outf, open(rm_map, "w") as rmf:
        with open(in_map) as inf:
            for line0 in inf:
                line1 = line0.rstrip()
                parts = line1.split(":")
                if len(parts) < 4:
                    outf.write(line1 + "\n")
                elif (parts[3]).upper() in chr_lst:
                    outf.write(line1 + "\n")
                else:
                    rmf.write(line1 + "\n")
  

def dup_remover_main(in_fa, out_fa, rm_fa, in_map, out_map, rm_map):
    larr = range(1, 23)
    larr2 = list(map(str,larr))
    larr3 = ['X', 'Y', 'MT', '']
    chr_lst = larr2 + larr3
    
    dup_remover_fa(in_fa, out_fa, rm_fa, chr_lst)
    dup_remover_map(in_map, out_map, rm_map, chr_lst)


if __name__ == '__main__':
    args = docopt(__doc__, version = 'Remove duplicate genes')

    in_fa = args["--in_fa"]
    out_fa = args["--out_fa"]
    rm_fa = args["--rm_fa"]
    in_map = args["--in_map"]
    out_map = args["--out_map"]
    rm_map = args["--rm_map"]

    print("in_fa: " + in_fa)
    print("out_fa: " + out_fa)
    print("rm_fa: " + rm_fa)
    print("in_map: " + in_map)
    print("out_map: " + out_map)
    print("rm_map: " + rm_map)

    dup_remover_main(in_fa, out_fa, rm_fa, in_map, out_map, rm_map)


