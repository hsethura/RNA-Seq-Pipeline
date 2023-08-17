#!/usr/bin/env python

import os
import argparse

from subprocess import call

class SFCore:

    def __init__(self, args):
        self.infile = args.infile
        self.outfile = args.outfile
        self.top_seed = int(args.top_seed)
        self.sample_count = int(args.sample_count)

        self.sample_fastq_path = "/broad/IDP-Dx_work/nirmalya/research/fastq_sampler/sample_fastq"


    def exe_sample_fastq(self):
        infile = self.infile
        outfile = self.outfile
        top_seed = self.top_seed
        sample_count = self.sample_count

        
        sample_fastq_path = self.sample_fastq_path
        sample_fastq_cmd = sample_fastq_path + " -i " + infile + " -o " + outfile + " -t " + str(top_seed) + " -s " + str(sample_count)
        call(sample_fastq_cmd.split())
        print("sample_fastq_cmd: " + sample_fastq_cmd)



    def mainFunc(self):
        self.exe_sample_fastq() 



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Process command for a sample bam core instance", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--infile", dest = "infile", type = str, required = True, help = "Input fastq file")    
    parser.add_argument("--outfile", dest ="outfile", type = str, required = True, help = "Output fastq file")
    parser.add_argument("--top_seed", dest = "top_seed", type = str, required = True, help = "Top seed for the sampling")
    parser.add_argument("--sample_count", dest = "sample_count", type = str, required = True, help = "Number of counts to sample")

    args = parser.parse_args()
    sfco = SFCore(args)
    print("About to start sfco mainFunc")
    sfco.mainFunc()
    
    
