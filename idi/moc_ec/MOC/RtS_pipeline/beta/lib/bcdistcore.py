#!/usr/bin/env python

import argparse
import yaml
from subprocess import call

from dictmap import DictMap

# For dropseq, this would create read distribution per barcode

class BCDistCore:
    def __init__(self, args):
        self.config_log_file = args.config_log_file

        # Loaded cldict/confd        
        self.sample_id = args.sample_id
        self.project_id = args.project_id
        #self.prefix_set = args.prefix_set
        #self.bc_set = args.bc_set
        self.ref_acc_str = args.ref_acc_str
        self.host_ref_str = args.host_ref_str
        cldict_d = yaml.load(open(self.config_log_file))
        cldict = DictMap(cldict_d)
        self.cldict = cldict

        # Loaded sampd
        sampd = dict()
        sampd['sample_id'] = self.sample_id
        sampd['project_id'] = self.project_id
        sampd['ref_acc_str'] = self.ref_acc_str
        sampd['host_ref_str'] = self.host_ref_str
        sampd_map = DictMap(sampd)
        self.sampd = sampd_map

   
    def exe_bc_counter(self):
        cldict = self.cldict
        ldel = cldict.ldelim
        Split_dir = cldict.Split_dir
        Merge_dir = cldict.Merge_dir
        Patho_dir = cldict.Patho_dir
        sampd = self.sampd
        sample_id = sampd.sample_id
        project_id = sampd.project_id
        ref_acc_str = sampd.ref_acc_str
        # bc files for valid reads
        merge_real = Merge_dir + ldel + "real"
        lbc_file = merge_real + ldel + sample_id + "_valid_bc.txt"
        real_sam_dir = Patho_dir + ldel + project_id + ldel + "real"
        bcdist_dir = Patho_dir + ldel + project_id + ldel + "bcdist"
        real_sam_file = real_sam_dir + ldel + sample_id + "_" + ref_acc_str + "_u.bam"
        aligned_outfile = bcdist_dir + ldel + sample_id + "_" + ref_acc_str + "_aligned.txt"
        total_outfile = bcdist_dir + ldel + sample_id + "_" + ref_acc_str + "_total.txt"
        ds_bc_counter_script = cldict.ds_bc_counter
        bc_counter_cmd = ds_bc_counter_script + " --bc_file " + lbc_file + " --sam_file " + real_sam_file + " --out_file_t " + total_outfile + " --out_file_a " + aligned_outfile
        print("bc_counter_cmd started: " + bc_counter_cmd)
        call(bc_counter_cmd.split()) 
        print("bc_counter_cmd done.")
        

    def mainFunc(self):

        cldict = self.cldict
        sampd = self.sampd
        self.exe_bc_counter() 

        


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Process command for a single sample which is essentially a row in the key file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--config_log_file", dest = "config_log_file", type = str, required = True, help = "Log file containig the config vars.")
    parser.add_argument("--sample_id", dest="sample_id", type=str, required=True, help="Id of the sample to be processed.")
    parser.add_argument("--project_id", dest = "project_id", required = True, type = str, help = "project id for this sample")
    parser.add_argument("--ref_acc_str", dest = "ref_acc_str", type = str, default = "none", help = "Contains the host accession ids")
    parser.add_argument("--host_ref_str", dest ="host_ref_str", default = "none", help = "Link to the host reference")

    args = parser.parse_args()
 
    bcdisto = BCDistCore(args)
    print("About to start uino mainFunc")
    bcdisto.mainFunc()
 
