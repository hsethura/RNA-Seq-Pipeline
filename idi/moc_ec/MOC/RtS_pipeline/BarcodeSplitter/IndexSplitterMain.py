#!/usr/bin/env python

# This script would parse the input directory to idenify eligible flowcell ids, and then construct the 
# eligible input files and from there would submit index_splitter for each of the flowcell ids

import argparse
import shutil
import re
import os
import os.path
from subprocess import call

parser = argparse.ArgumentParser(description = "Process inputs for index splitter", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--indir', '-i', dest = 'indir', type = str, required = True, help = "Input directory")
parser.add_argument('--outdir', '-o', dest = 'outdir', type = str, required = True, help = "Output directory")
parser.add_argument('--dict_file', '-d', dest = 'dictfile', type = str, required = True, help = "p7 index dict file")
parser.add_argument('--no_qsub', dest = 'use_qsub', action = 'store_false', default = True, help = 'Does not submit qsub jobs.' )

args = parser.parse_args()

indir = args.indir
outdir = args.outdir
dictfile = args.dictfile
use_qsub = args.use_qsub

ldelim = '/'

print("Input dir: " + indir)
print("Output dir: " + outdir)
print("p7 index dict file" + dictfile)

# Get all the files that ends with fastq.gz

prefix_set = set()

for lfile in os.listdir(indir):
    if lfile.endswith("fastq.gz") or lfile.endswith("fastq"):
        #print("file: " + lfile)
        parts = lfile.split(".")
        prefix = parts[0]
        if prefix not in prefix_set:
            prefix_set.add(prefix)


# For each element in the prefix_set submit a job to the cluster

Script_dir = os.path.dirname(os.path.realpath(__file__))
index_split_cpp = Script_dir + ldelim + "index_splitter"

if not os.path.exists(outdir):
    os.makedirs(outdir)
joblist_path = outdir + ldelim + "index_bcsplit_joblist.txt"
jfile = open(joblist_path, "w")

def refine_gz_file(gz_file):
    if not os.path.isfile(gz_file):
        only_file = gz_file.replace(".gz", "")
        if not os.path.isfile(only_file):
            raise IOError('File not found: ' + gz_file)
        else:
            return only_file
    else: 
        return gz_file


for prefix in prefix_set:
    parts = prefix.split("_")
    lane = parts[0]

    index_file_str = indir + ldelim + prefix + "." + lane + ".barcode_1.fastq.gz"
    index_file = refine_gz_file(index_file_str)
    read1_file_str = indir + ldelim + prefix + "." + lane + ".1.fastq.gz"
    read1_file = refine_gz_file(read1_file_str)
    read2_file_str = indir + ldelim + prefix + "." + lane + ".2.fastq.gz"
    read2_file = refine_gz_file(read2_file_str)

    out_log = outdir + ldelim + prefix + "." + lane + "_out.txt"
    err_log = outdir + ldelim + prefix + "." + lane + "_err.txt"
    lprefix = prefix + "." + lane
    print("prefix: " + lprefix)
    #print("index: " + index_file)
    #print("read1_file: " + read1_file)
    #print("read2_file: " + read2_file)
    job_str = index_split_cpp + " -d " + dictfile + " -i " + index_file + " --file1 " + read1_file + \
        " --file2 " + read2_file + " -p " + lprefix + " -o " + outdir + " 1> " + out_log + " 2> " + err_log + "\n"
    print("job_str: " + job_str)
    jfile.write(job_str)
jfile.close()

print(joblist_path)

UGER_cbp = "/broad/IDP-Dx_work/nirmalya/tools/ugetools/UGE_SUBMISSIONS/UGER_cmd_batch_processor.py"
UGER_cbp_dir = outdir + ldelim + "UGER_cbp"

if use_qsub and os.path.exists(UGER_cbp_dir):
    shutil.rmtree(UGER_cbp_dir)

joblist_cmd = UGER_cbp + " --cmds_file " + joblist_path + \
                                " --batch_size 1" + \
                                " --queue long" + \
                                " --memory 16" + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header"

if use_qsub:
    call(joblist_cmd.split())    

