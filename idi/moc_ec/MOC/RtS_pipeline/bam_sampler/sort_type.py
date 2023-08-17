import subprocess
import sys


infile = "/broad/hptmp/RNASeq_proj/MOC/nirmalya/gp_test/AbxTimeTrial/patho_result/AbxTimeTrial/T90_waterB_3_Pseudomonas_aeruginosa_UCBPP-PA14_pe.bam"

samtools_scr = "/broad/IDP-Dx_work/nirmalya/local/bin/samtools"

samtools_cmd = samtools_scr + " view -H " + infile

#result = subprocess.check_output(samtools_cmd, shell = True)


proc = subprocess.Popen(samtools_cmd.split(), stdout=subprocess.PIPE)
for line in proc.stdout.readlines():
    line1 = line.rstrip()
    if '@HD' in line1:
        parts = line1.split()
        for lpart in parts:
            if 'SO:queryname' in lpart:
                print("Type is queryname")
            elif 'SO:coordinate' in lpart:
                print("Type is coordinate.")






