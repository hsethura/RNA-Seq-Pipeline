#!/usr/bin/env python

import argparse
import yaml
import ntpath
import os
import re
from dictmap import DictMap
from subprocess import call
from alignerbwa import AlignerBwa
from alignerutils import bam_to_sam
from Parse_featureCounts import get_all_gene_counts

class UniDropCore:
    def __init__(self, args):
        self.config_log_file = args.config_log_file
        cldict_d = yaml.load(open(self.config_log_file))
        cldict = DictMap(cldict_d)
        self.confd = cldict
 
        self.sample_id = args.sample_id
        self.project_id = args.project_id
        self.bc_pre_str = args.bc_pre_str
        self.ref_acc_str = args.ref_acc_str

    def get_sample_prefix(self, lfile):
        sample_prefix = None
        if lfile.endswith("_u.sam"):
            sample_prefix = lfile.replace("_u.sam", "")
        elif lfile.endswith("_u.bam"):
            sample_prefix = lfile.replace("_u.bam", "")
        elif lfile.endswith("_se.bam"):
            sample_prefix = lfile.replace("_se.bam", "")
        elif lfile.endswith("_se.sam"):
            sample_prefix = lfile.replace("_se.sam", "")
        return (sample_prefix)


    def get_unsorted_bam_to_sam(self, unsorted_bam):
        confd = self.confd
        samtools_path = confd.samtools
        unsorted_sam = unsorted_bam.replace("_u.bam", "_u.sam")
        bam_to_sam(samtools_path, unsorted_bam, unsorted_sam)
        return unsorted_sam

    
    def path_leaf(self, path):
        # From https://stackoverflow.com/questions/8384737/extract-file-name-from-path-no-matter-what-the-os-path-format
        head, tail = ntpath.split(path)
        return tail or ntpath.basename(head)

    def umi_normalize(self, bam_sorted_path):
        confd = self.confd
        samtools_path = confd.samtools
        Patho_dir = confd.Patho_dir
        ldel = confd.ldelim
        umidir = self.umidir
        umi_log_dir = self.umi_logdir
        collapse_type = "coordinate"
        uminorm_path =  confd.uminorm_path
        bam_sorted_file = self.path_leaf(bam_sorted_path)
        sample_prefix = self.get_sample_prefix(bam_sorted_file)

        umi_cmd = uminorm_path + " -i " + bam_sorted_path + " -o " + umidir + \
            " -p " + sample_prefix + " -c " + collapse_type
        print("umi_cmd starting: " + umi_cmd)
        call(umi_cmd.split())
        print("umi_cmd ended.")

        out_bam_base_u = sample_prefix + '_u.bam'
        umi_bam_unsorted = umidir + ldel + out_bam_base_u
        print("returning umi_bam_unsorted: " + umi_bam_unsorted)
        return umi_bam_unsorted
      

    def create_dirs(self):
        confd = self.confd
        Patho_dir = confd.Patho_dir
        sample_id = self.sample_id
        project_id = self.project_id
        bc_pre_str = self.bc_pre_str
        ldel = confd.ldelim
        outdir = Patho_dir + ldel + project_id
        Patho_drops_dir = outdir + ldel + "drops"
        sample_drop_dir = Patho_drops_dir + ldel + sample_id
        bc_pre_dir = sample_drop_dir + ldel + bc_pre_str
        self.bc_pre_dir = bc_pre_dir

        bc_pre_logdir = bc_pre_dir + ldel + "logdir"
        self.bc_pre_logdir = bc_pre_logdir

        bc_pre_metrics_dir = bc_pre_dir + ldel + "metrics"
        self.bc_pre_metrics_dir = bc_pre_metrics_dir
        if not os.path.exists(bc_pre_metrics_dir):
            os.makedirs(bc_pre_metrics_dir)

        patho_temp_bamdir = bc_pre_dir + ldel + "temp_bamdir"
        self.patho_temp_bamdir = patho_temp_bamdir
        if not os.path.exists(patho_temp_bamdir):
            os.makedirs(patho_temp_bamdir)


        umidir = bc_pre_dir + ldel + "umidir"
        self.umidir = umidir
        if not os.path.exists(umidir):
            os.makedirs(umidir)

        umi_logdir = umidir + ldel + "logdir"
        self.umi_logdir = umi_logdir
        if not os.path.exists(umi_logdir):
            os.makedirs(umi_logdir)

        umi_metrics_dir = umidir + ldel + "metrics"
        self.umi_metrics_dir = umi_metrics_dir
        if not os.path.exists(umi_metrics_dir):
            os.makedirs(umi_metrics_dir)

        umi_temp_bamdir = umidir + ldel + "temp_bamdir"
        self.umi_temp_bamdir = umi_temp_bamdir
        if not os.path.exists(umi_temp_bamdir):
            os.makedirs(umi_temp_bamdir)
             

    def count_single(self, ref_acc, sampath, outdir, strand_rev):
        # sampath points to sorted or unsorted samfile
        cldict = self.confd
        JLCounter = cldict.JLCounter
        ldelim = cldict.ldelim
        Data_dir = cldict.Data_dir
        Script_dir = cldict.basepath
        Patho_dir = cldict.Patho_dir
        patho_temp_bamdir = outdir + ldelim + "temp_bamdir"
        patho_gff = Data_dir + ldelim + ref_acc + "_ALL.gff"
        shell_script_dir = Script_dir + ldelim + "shell_scripts"
        samfile = self.path_leaf(sampath)
        print("samfile from count_single: " + samfile)
        sample_prefix = self.get_sample_prefix(samfile)
        print("sample_prefix from count_single: " + sample_prefix)
        
        countfile_str = outdir + ldelim + sample_prefix + ".counts"
        l_strand_rev = None
        if strand_rev:
            l_strand_rev = "Y"
        else:
            l_strand_rev = "N"
        count_cmd = "sh " + JLCounter + " " + sampath + " " + patho_gff + \
            " " + shell_script_dir + " " + patho_temp_bamdir + " " + \
            countfile_str + " -STRAND_REV " + l_strand_rev
        print("count_cmd: " + count_cmd)
        call(count_cmd.split())

    def get_sam_to_bam_sorted(self, outsamfile):
        confd = self.confd
        samtools_path = confd.samtools
        outbamfile = outsamfile.replace('_u.sam', '_u.bam')
        print("outbamfile: " + outbamfile)
        print("samtools_path: " + samtools_path)
        AlignerBwa.sam_to_bam(samtools_path, outsamfile, outbamfile)
        suf = "se"
        suf_str = "_" + suf + ".bam"
        outsorted = outbamfile.replace('_u.bam', suf_str)
        AlignerBwa.sort_bam(samtools_path, outbamfile, outsorted)
        os.remove(outbamfile)
        print("Deleting outbamfile: " + outbamfile)
        return outsorted

    def exe_single_drop_path(self):
        confd = self.confd
        Patho_dir = confd.Patho_dir
        sample_id = self.sample_id
        project_id = self.project_id
        bc_pre_str = self.bc_pre_str
        ldel = confd.ldelim
        ref_acc = self.ref_acc_str
        outdir = Patho_dir + ldel + project_id
        Patho_drops_dir = outdir + ldel + "drops"
        sample_drop_dir = Patho_drops_dir + ldel + sample_id
        bc_pre_dir = sample_drop_dir + ldel + bc_pre_str
        sam_file_list_file = bc_pre_dir + ldel + "sam_file_list.txt"
        print("sam_file_list_file: " + sam_file_list_file)
        strand_rev = False
        with open(sam_file_list_file) as sf:
            for line in sf:
                file_name = line.rstrip()
                outsampath = bc_pre_dir + ldel + file_name 

                self.count_single(ref_acc, outsampath, bc_pre_dir, strand_rev) 
                bam_sorted = self.get_sam_to_bam_sorted(outsampath)
                umi_bam_unsorted = self.umi_normalize(bam_sorted)
                umi_sam_unsorted = self.get_unsorted_bam_to_sam(umi_bam_unsorted)
                umidir = self.umidir
                self.count_single(ref_acc, umi_sam_unsorted, umidir, strand_rev) 
                umi_bam_sorted = self.get_sam_to_bam_sorted(umi_sam_unsorted)
                os.remove(umi_sam_unsorted)
                os.remove(outsampath)

   
    def ref_RPG_drop_metrics(self, sample_prefix, tempdir, result_dir):
        # sh /broad/IDP-Dx_work/nirmalya/pipeline/alpha/shell_scripts/RPG_drop_metrics.sh 
        # -p 'F118.C2_postHC65_NC_009641_AAACAGGG_TGCCTCAC' -t . -r results
        confd = self.confd
        RPG_drop_metrics_script = confd.RPG_drop_metrics_script
        #RPG_drop_metrics_script = "/broad/IDP-Dx_work/nirmalya/pipeline/alpha/shell_scripts/RPG_drop_metrics.sh"
        RPG_drop_cmd = "sh " + RPG_drop_metrics_script + " -p " + sample_prefix + " -t " + tempdir + " -r " + result_dir
        print("RPG_drop_cmd: " + RPG_drop_cmd)
        call(RPG_drop_cmd.split())
        

    def combine_all_metrics(self):
        result_dir = self.bc_pre_dir
        outfile = result_dir
        # get all the prefix
        sample_prefix_lst = list()
        
        confd = self.confd
        Patho_dir = confd.Patho_dir
        sample_id = self.sample_id
        project_id = self.project_id
        bc_pre_str = self.bc_pre_str
        ldel = confd.ldelim
        ref_acc = self.ref_acc_str
        outdir = Patho_dir + ldel + project_id
        Patho_drops_dir = outdir + ldel + "drops"
        sample_drop_dir = Patho_drops_dir + ldel + sample_id
        bc_pre_dir = sample_drop_dir + ldel + bc_pre_str
        sam_file_list_file = bc_pre_dir + ldel + "sam_file_list.txt"
        print("sam_file_list_file: " + sam_file_list_file)
        strand_rev = False
        with open(sam_file_list_file) as sf:
            for line in sf:
                samfile = line.rstrip()
                sample_prefix = self.get_sample_prefix(samfile)
                sample_prefix_lst.append(sample_prefix) 

        # generate metrics for pre_pre_dir and be_pre_umi_dir
        
        print("sample_prefix_lst")
        print(sample_prefix_lst)
        sample_tag = sample_id + "_" + ref_acc + "_" + bc_pre_str
        print("Combining read counts in bc_pre_dir for sample: " + sample_tag)
        bc_pre_metrics_dir = self.bc_pre_metrics_dir
        loutfile = bc_pre_metrics_dir + ldel + sample_tag + ".counts"
        get_all_gene_counts(bc_pre_dir, loutfile, sample_prefix_lst, ref_acc = '', has_header = False)    
        print("Combining read counts for sample_tag: " + sample_tag)

        print("Combining read counts in umidir for sample_tag: " + sample_tag)
        umidir = self.umidir
        umi_metrics_dir = self.umi_metrics_dir
        lumi_outfile = umi_metrics_dir + ldel + sample_tag + ".counts"
        get_all_gene_counts(umidir, lumi_outfile, sample_prefix_lst, ref_acc = '', has_header = False)

        self.ref_RPG_drop_metrics(sample_tag, bc_pre_dir, bc_pre_metrics_dir)
        self.ref_RPG_drop_metrics(sample_tag, umidir, umi_metrics_dir)
                
        
    def mainFunc(self):
        # 1. Locate the folder for this bc_pre_str and create a log dir
        # 2. Create a UMI dir inside the folder and create a log dir inside it
        # 3. Locate all the bam files inside the folder
        # 4. Create a file with counting jobs for all the samples 
        # 5. Execute the jobs sequentially
        # 6. Create a file with jobs for UMI normalization for all the samples
        # 7. Execute the UMI jobs sequentially
        # 8. Crerate a file with counting jobs with UMI normalized drops

        # get the file containing the samfile
        self.create_dirs()
        self.exe_single_drop_path()
        self.combine_all_metrics()
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Process command for a single bc_pre dir", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--config_log_file", dest = "config_log_file", type = str, required = True, help = "Log file containig the config vars.")
    parser.add_argument("--sample_id", dest="sample_id", type=str, required=True, help="Id of the sample to be processed.")
    parser.add_argument("--bc_pre_str", dest = "bc_pre_str", required = True, type = str, help = "bc_pre_str for this run")
    parser.add_argument("--project_id", dest = "project_id", required = True, type = str, help = "project id for this sample")
    parser.add_argument("--ref_acc_str", dest = "ref_acc_str", type = str, default = "none", help = "Contains the host accession ids")
    args = parser.parse_args()
    usco = UniDropCore(args)
    
    usco.mainFunc()

