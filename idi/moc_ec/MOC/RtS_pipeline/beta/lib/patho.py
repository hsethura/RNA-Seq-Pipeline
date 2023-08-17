#!/usr/bin/env python

import shutil
import re
import os.path
from subprocess import call

def copyLargeFile(src, dest, buffer_size=16000):
    print("src: " + src + " dest: " + dest)
    with open(src, 'rb') as fsrc:
        with open(dest, 'wb') as fdest:
            shutil.copyfileobj(fsrc, fdest, buffer_size)

def replace_refname(patho_gff, patho_fna, patho_fna_up, select_count):
    select_line = ''
    with open(patho_gff) as fp:
        for i, line in enumerate(fp):
            if i == select_count:
                select_line = line
    replace_part = re.split('\s*', select_line)[0]
    print("replace_part: " + replace_part)
    
    fna_first_line = ''
    with open(patho_fna) as f:
	    fna_first_line = f.readline()
    
    replace_part_s = ">" + replace_part
    with open(patho_fna_up, "w") as outf:
        with open(patho_fna, "r") as inf:
            lline = inf.readline()
            # Build the replacement line
            replace_line = re.sub('^\S+', replace_part_s, lline)
            print("Replace line: " + replace_line)
            outf.write(replace_line)
            shutil.copyfileobj(inf, outf)


def copy_and_update_patho(patho_id, gff_parser, add5, add3, inpath, outpath, \
        do_replace_refname, ldelim = '/', select_count = 10):
    """ sh gff_parse3.sh <NCBI_PATH> <OUT_DIR - should be patho_results?> 
    <TEMP_DIR - for intermediate files> <ADD5 - default 0> <ADD3 - default 0> 
    <list of refseq accessions - can be as many as you want - it will run a loop to make them all >
    """
   # gff_parser_cmd = "sh " + gff_parser + " " + inpath + " " + outpath + " " + outpath + " " + str(add5) + " " + str(add3) + " " + patho_id   
   # print("gff_parser_cmd: " + gff_parser_cmd)
   # call(gff_parser_cmd.split())
       
    
def touch(path):
    with open(path, 'a'):
        os.utime(path, None)


def files_concat(source_lst, dest_file):
    # Concatenate files in source_lst to des_file in order
    copy_len = 1024*1024*10

    for lpath in source_lst:
        touch(lpath)

    with open(dest_file,'wb') as wfd:
        for f in source_lst:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd, copy_len)

def patho_fna_concat(patho_lst, compo_name, ldir, ldelim = "/"):
    # Make path to input files
    lpath_lst = [ldir + ldelim + lpatho + ".fna" for lpatho in patho_lst]
    out_path = ldir + ldelim + compo_name + ".fna"
    files_concat(lpath_lst, out_path)


def make_index_bwa(lprefix, loutdir, bwa_path, ldelim = "/"):
    lfile_str = loutdir + ldelim + lprefix + ".fna"
    if os.path.isfile(lfile_str):
        bwa_index_cmd = bwa_path + " index " + lfile_str
        call(bwa_index_cmd.split())
    else:
        raise IOError("File not found: " + lfile_str)
    

if __name__ == "__main__":
    inpath = '/broad/IDP-Dx_storage/NCBI_files2/'
    outpath = '/broad/IDP-Dx_work/nirmalya/temp/temp2'
    patho_id = 'NC_007946'
    copy_and_update_patho(patho_id, inpath, outpath, )

    
    
