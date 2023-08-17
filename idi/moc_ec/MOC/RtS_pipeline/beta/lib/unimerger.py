#!/usr/bin/env python

import argparse
import yaml
import os

from dictmap import DictMap
from merger import Merger

# This would execute the core part for RtS (RNATag-Seq) pipeline.

class UniMerger:
    def __init__(self, args):
        self.config_log_file = args.config_log_file
        
        self.sample_id = args.sample_id
        self.project_id = args.project_id
        self.prefix_set = args.prefix_set
        self.bc_set = args.bc_set

        cldict_d = yaml.load(open(self.config_log_file))
        cldict = DictMap(cldict_d)
        self.cldict = cldict

        sampd = dict()
        sampd['sample_id'] = self.sample_id
        sampd['project_id'] = self.project_id
        sampd['prefix_set'] = self.prefix_set
        sampd['bc_set'] = self.bc_set
        sampd_map = DictMap(sampd) 
        self.sampd = sampd_map
 
        mergo = Merger(cldict, sampd_map)
        self.mergo = mergo


    def unimerge_single(self):
        cldict = self.cldict
        sampd = self.sampd
        mergo = self.mergo
        project_id = self.project_id

        ldelim = cldict.ldelim
        Patho_dir = cldict.Patho_dir
        Merge_dir = cldict.Merge_dir
        outdir = Merge_dir
        sample_id = sampd.sample_id

        if not os.path.exists(outdir):
            os.makedirs(outdir)
        suffix = "R"

        trim_rs_5p =  cldict.trim_rs_5p
        trim_rs_3p = cldict.trim_rs_3p
        keep_rs_5p = cldict.keep_rs_5p
        keep_rs_3p = cldict.keep_rs_3p
        print("5P trim read-single count: " + str(trim_rs_5p))
        print("3P trim read-single count: " + str(trim_rs_3p))
        print("5P keep read-single count: " + str(keep_rs_5p))
        print("3P keep read-single count: " + str(keep_rs_3p))
        print("sample_id: " + sample_id)
        print("outdir: " + outdir)
        merged_file = mergo.merge_fastq_single(sample_id, suffix, outdir, False, False, trim_rs_5p, trim_rs_3p, keep_rs_5p, keep_rs_3p)
        return merged_file

    def unimerge_dropseq(self):
        # Unimerge would be different from the typical merge, since
        # 1. It would merge across possibly two sets of files in possibly
        # two directories, real and garbage
        # 2. It Would not use barcode. It would just merge across
        # multiple lanes and flocells
        cldict = self.cldict
        sampd = self.sampd
        sample_id = sampd.sample_id
        mergo = self.mergo
        print("from unimerge_dropseq sample_id: " + sample_id)
        Merge_dir = cldict.Merge_dir
        ldel = cldict.ldelim
        outdir_r = Merge_dir + ldel + "real"
        outdir_g = Merge_dir + ldel + "garbage"
        suf1_r = "R1_R" 
        suf2_r = "R2_R" 
        suf1_g = "R1_G"
        suf2_g = "R2_G"
        print("Merging R1 file for real.")
        read1_r_file = mergo.merge_fastq_single(sample_id, suf1_r, outdir_r, False, False)
        print("Merging R2 file for real.")
        read2_r_file = mergo.merge_fastq_single(sample_id, suf2_r, outdir_r, False, False)
        print("Merging R1 file for garbage.")
        read1_g_file = mergo.merge_fastq_single(sample_id, suf1_g, outdir_g, False, False)
        print("Merging R2 file for garbage.")
        read2_g_file = mergo.merge_fastq_single(sample_id, suf2_g, outdir_g, False, False)
        # We also have to merge the barcode files
        valid_bc_file = mergo.merge_valid_barcode_single(sample_id, outdir_r)
        

    def unimerge_paired(self):
        cldict = self.cldict
        sampd = self.sampd
        mergo = self.mergo

        ldelim = cldict.ldelim
        Merge_dir = cldict.Merge_dir
        project_id = sampd.project_id
        outdir = Merge_dir
        sample_id = sampd.sample_id
        do_R2_trim = cldict.do_R2_trim
        do_allseq_trim = cldict.do_allseq_trim


        if not os.path.exists(outdir):
            os.makedirs(outdir)
        suf1 = "R1"
        suf2 = "R2"

        print("sample_id: " + sample_id)
        print("outdir: " + outdir)
        LC_method_val = cldict.LC_method_val.lower()
        print("LC_method_val updated: " + LC_method_val)
        
        trim_r1_5p = cldict.trim_r1_5p
        trim_r1_3p = cldict.trim_r1_3p
        trim_r2_5p = cldict.trim_r2_5p
        trim_r2_3p = cldict.trim_r2_3p
        print("5P trim read1 count: " + str(trim_r1_5p))
        print("3P trim read1 count: " + str(trim_r1_3p))
        print("5P trim read2 count: " + str(trim_r2_5p))
        print("3P trim read2 count: " + str(trim_r2_3p))
        
        keep_r1_5p = cldict.keep_r1_5p
        keep_r1_3p = cldict.keep_r1_3p
        keep_r2_5p = cldict.keep_r2_5p
        keep_r2_3p = cldict.keep_r2_3p
        print("5P keep read1 count: " + str(keep_r1_5p))
        print("3P keep read1 count: " + str(keep_r1_3p))
        print("5P keep read2 count: " + str(keep_r2_5p))
        print("3P keep read2 count: " + str(keep_r2_3p))
        
        if LC_method_val == 'alignr1':
            print("Printing from alignr1")
            read1_file = mergo.merge_fastq_single(sample_id, suf1, outdir, False, False, trim_r1_5p, trim_r1_3p, keep_r1_5p, keep_r1_3p)
            print("alignr1- Merged read 1: " + read1_file)
        elif LC_method_val == "allseq" or LC_method_val == "alignr2":
            read2_file = mergo.merge_fastq_single(sample_id, suf2, outdir, do_R2_trim, do_allseq_trim, trim_r2_5p, trim_r2_3p, keep_r2_5p, keep_r2_3p)
            print("alignr1 or allseq - Merged read 2: " + read2_file)
        else:
            read1_file = mergo.merge_fastq_single(sample_id, suf1, outdir, False, False, trim_r1_5p, trim_r1_3p, keep_r1_5p, keep_r1_3p)
            read2_file = mergo.merge_fastq_single(sample_id, suf2, outdir, do_R2_trim, do_allseq_trim, trim_r2_5p, trim_r2_3p, keep_r2_5p, keep_r2_3p)


    def mainFunc(self):

        cldict = self.cldict
        Read_pairing_val = cldict.Read_pairing_val
        print("From unimerger: " +  Read_pairing_val)
        use_dropseq = cldict.use_dropseq
        
        if use_dropseq:
            self.unimerge_dropseq()
        elif Read_pairing_val == 'SINGLE':
            self.unimerge_single()
        elif Read_pairing_val == 'PAIRED':
            self.unimerge_paired()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Process command for a single sample which is essentially a row in the key file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--config_log_file", dest = "config_log_file", type = str, required = True, help = "Log file containig the config vars.")
    parser.add_argument("--sample_id", dest="sample_id", type=str, required=True, help="Id of the sample to be processed.")
    parser.add_argument("--project_id", dest = "project_id", required = True, type = str, help = "project id for this sample")
    parser.add_argument("--prefix_set", dest="prefix_set", type=str, required=True, help="Set of file prefixes for this sample")
    parser.add_argument("--bc_set", dest="bc_set", type=str, required=True, help = "Set of the barcodes used for this samples")

    args = parser.parse_args()
 
    unimo = UniMerger(args)
    print("About to start unimo mainFunc")
    unimo.mainFunc()
 
