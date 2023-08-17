#!/usr/bin/env python

import os
import argparse
import subprocess
import ntpath
import shutil

from subprocess import call 
from shutil import copyfile

from sfcore import SFCore


class FastqSamplerMain:
    def __init__(self, args):
        self.count_file = args.count_file
        self.indir = args.indir
        self.outdir = args.outdir
        self.main_seed = args.main_seed
        self.use_qsub = args.use_qsub
        self.file_suffix = args.file_suffix

        basepath = os.path.dirname(os.path.realpath(__file__))
        self.Script_dir = basepath

        self.SFCore = self.Script_dir + "/" + "sfcore.py"


        self.seed_gen_path = "/broad/IDP-Dx_work/nirmalya/research/fastq_sampler/seed_fgen"
        self.sample_fastq_bath = "/broad/IDP-Dx_work/nirmalya/research/fastq_sampler/sample_fastq"
        self.UGER_cbp = "/broad/IDP-Dx_work/nirmalya/tools/ugetools/UGE_SUBMISSIONS/UGER_cmd_batch_processor.py"

        outdir = self.outdir
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        logdir = outdir + "/logdir"
        if not os.path.exists(logdir):
            os.makedirs(logdir)
        self.logdir = logdir

        mergedir = outdir + "/mergedir"
        if not os.path.exists(mergedir):
            os.makedirs(mergedir)
        self.mergedir = mergedir

        UGER_cbp_dir = outdir + "/UGER_cbp"

        if self.use_qsub and os.path.exists(UGER_cbp_dir):
            shutil.rmtree(UGER_cbp_dir)
        self.UGER_cbp_dir = UGER_cbp_dir


    def exe_seed_gen(self):
        seed_gen_path = self.seed_gen_path
        count_file = self.count_file
        main_seed = self.main_seed
        logdir = self.logdir

        count_leaf = os.path.basename(count_file)
        count_dir = os.path.dirname(count_file)
        seed_tab_file = count_leaf + ".seedtab"
        seed_table = logdir + "/" + seed_tab_file
        print("seed_table: " + seed_table)
        
        seed_gen_cmd = seed_gen_path + " -i " + count_file +\
             " -o " + seed_table + " -m " + str(main_seed)

        print("seed_gen_cmd: " + seed_gen_cmd)
        call(seed_gen_cmd.split())
        print("End of seed_gen")
        self.seed_table = seed_table
        return seed_table

    def get_sample_file_dict(self):

        indir = self.indir
        lsuffix = self.file_suffix
        sample_file_dict = dict()

        sample_set = self.sample_set

        for lfile in os.listdir(indir):
            fname = os.path.basename(lfile)
            if fname.endswith("gz"):
                lsuffix += ".gz"
            lprefix = fname.replace(lsuffix, '')

            if lprefix in sample_set:
                lpath = indir + "/" + fname
                sample_file_dict[lprefix] = lpath
        
        self.sample_file_dict = sample_file_dict
        print("sample_file_dict")
        print("........")
        print(sample_file_dict)

        return sample_file_dict
    

    def get_sample_set(self):
        seed_table = self.seed_table
        
        sample_set = set()
        with open(seed_table) as stab:
            for line1 in stab:
                line2 = line1.rstrip()
                parts = line2.split()
                lsample = parts[0]
                sample_set.add(lsample)

        self.sample_set = sample_set
        print("sample_set")
        print("........")
        print(sample_set)

        return sample_set

        
    def exe_sample_fastq_core(self):

        outdir = self.outdir
        mergedir = self.mergedir
        logdir = self.logdir
        use_qsub = self.use_qsub
        seed_table = self.seed_table
        sample_file_dict = self.sample_file_dict
        SFCore = self.SFCore

        sfcore_job_path = logdir + "/" + "sfcore_joblist_txt"
        num_cores = 1
        lmemory = 8
        UGER_cbp = self.UGER_cbp
        UGER_cbp_dir = self.UGER_cbp_dir
        lsuf_lst = list()

        jfile = open(sfcore_job_path, "w")
        with open(seed_table) as stab:
            stab.next()
            for line1 in stab:
                line2 = line1.rstrip()
                parts = line2.split()
                lsample = parts[0]
                lsample_count = int(parts[1])
                ltop_seed = int(parts[2])
                infile = sample_file_dict[lsample]
                inleaf = os.path.basename(infile)
                outfile = mergedir + "/" + inleaf
                

                cmd_str = SFCore + " --infile " + infile +\
                    " --outfile " + outfile + " --top_seed " + str(ltop_seed) +\
                     " --sample_count " + str(lsample_count)
                loutfile = logdir + "/" + lsample + "_out.txt"
                lerrfile = logdir + "/" + lsample + "_err.txt"
                cmd_str2 = cmd_str + " 1> " + loutfile  + " 2> " + lerrfile
                print(cmd_str2) 
                jfile.write(cmd_str2 + "\n")
        jfile.close()

        self.lsuf_lst = lsuf_lst

        joblist_cmd = UGER_cbp + " --cmds_file " + sfcore_job_path + \
                                " --batch_size 1" + \
                                " --num_cores " + str(num_cores) + \
                                " --memory " + str(lmemory) + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header"        
        if use_qsub:
            call(joblist_cmd.split())

 
    def mainFunc(self):

        self.exe_seed_gen()
        self.get_sample_set()
        self.get_sample_file_dict()
        self.exe_sample_fastq_core()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Process command for bam sampler.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--count_file", dest = "count_file", type = str, required = True, help = "Feagment counts in the pre-dir.")
    parser.add_argument("--indir", dest = "indir", type = str, required = True, help = "Location of the input fastq files.")
    parser.add_argument("--outdir", dest = "outdir", type = str, required = True, help = "Location of the output fastq files.")
    parser.add_argument("--main_seed", dest = "main_seed", type = int, default = 12345, help = "Value of main seed")
    parser.add_argument("--file_suffix", dest = "file_suffix", type = str, default = "_R2.fastq", help = "Suffix at the end of the input files")
    parser.add_argument('--no_qsub', dest = 'use_qsub', action = 'store_false', default = True, help = 'Does not submit qsub jobs.' )
    
    args = parser.parse_args()
    fsmo = FastqSamplerMain(args)
    print("About to start fsmo mainFunc")
    fsmo.mainFunc()
    


