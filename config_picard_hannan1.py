#!/usr/bin/env python
import os
import re

bamdir = '/broad/hptmp/RNASeq_proj/nirmalya/Hannan1_picard'
jobfile_str = bamdir + "/joblist.txt"
jobfile = open(jobfile_str , "w")
cmd_str = "sh /broad/IDP-Dx_work/nirmalya/pipeline/beta/PICARD_metrics.sh"
delim = "/"
outdir = bamdir + delim + "outdir"
logdir = bamdir + delim + "logdir"

picard_path = '/broad/IDP-Dx_work/nirmalya/tools/picard/latest'
ref_file = '/broad/IDP-Dx_storage/Host_reference/data/Mouse_mm10.fna'
if not os.path.exists(outdir):
    os.makedirs(outdir)

if not os.path.exists(logdir):
    os.makedirs(logdir)

for lfile in os.listdir(bamdir):
    if lfile.endswith('.bam'):
        sample_id = re.sub("_pe.bam", "", lfile)
        outfile = logdir + delim + sample_id + "_out.txt"
        errfile = logdir + delim + sample_id + "_err.txt"
        bam_file = bamdir + delim + lfile
        cmd_path = cmd_str + " " + bam_file + " " + ref_file + " " + outdir + \
            " " + picard_path + " 1> " + outfile + " 2> " + errfile
        jobfile.write(cmd_path + "\n")
jobfile.close()

        
