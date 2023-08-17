import os
import sys
import shutil
import gzip
import os.path
from subprocess import call
from alignerbwa import AlignerBwa

class UMI:
    def __init__(self, cldict, sampd):
        self.cldict = cldict
        self.sampd = sampd

    def get_sampleid_refacc(self):
        sampd = self.sampd
        cldict = self.cldict
        sample_id = sampd.sample_id
        ref_acc = sampd.ref_acc_str
        print("from get_sampleid_refacc, ref_acc: " + ref_acc)
        res = sample_id + "_" + ref_acc
        return res

    def get_sampleid_hostref(self):
        sampd = self.sampd
        cldict = self.cldict
        sample_id = sampd.sample_id
        host_ref = sampd.host_ref_str
        print("from get_sampleid_refacc, hostref: " + host_ref)
        res = sample_id + "_" + host_ref
        return res

    def exe_umi_host(self, sorted_bam):
        cldict = self.cldict
        sampd = self.sampd
        Host_dir = cldict.Host_dir
        samtools_path = cldict.samtools
        ldelim = cldict.ldelim
        project_id = sampd.project_id
        outdir = Host_dir + ldelim + project_id
        umidir = outdir + ldelim + "umidir"
        umi_log_dir = umidir + ldelim + "logdir"
        collapse_type = "feature"
        uminorm_path = cldict.uminorm_path
        sampleid_refacc = self.get_sampleid_hostref()
        umi_cmd = uminorm_path + " -i " + sorted_bam + " -o " + umidir + \
            " -p " + sampleid_refacc + " -c " + collapse_type
        print("umi_cmd starting: " + umi_cmd)
        call(umi_cmd.split())
        print("umi_cmd ended.")
        out_bam_base_u = sampleid_refacc + '_u.bam'
        umi_bam_u = umidir + ldelim + out_bam_base_u
        Bam_path = cldict.Bam_path
        bamdir = Bam_path + ldelim + project_id
        bamdir_umi = bamdir + ldelim + "umidir"
        suf_str = "_se.bam"
        umi_bam = umi_bam_u.replace('_u.bam', suf_str)
        AlignerBwa.sort_bam(samtools_path, umi_bam_u, umi_bam)
        AlignerBwa.copy_bam(samtools_path, bamdir_umi, umi_bam)
        umi_logfile = umi_log_dir + ldelim + sampleid_refacc + "_coll_len.txt"
        umi_coll_len_script = cldict.umi_coll_len_script
        umi_coll_len_cmd = "Rscript " + umi_coll_len_script + " -i " + \
        umi_logfile + " -o " + umi_log_dir + " -p " + sampleid_refacc
        print("umi_coll_len_cmd starting: " + umi_coll_len_cmd)
        call(umi_coll_len_cmd.split())
        print("umi_coll_len_cmd ended") 

        return umi_bam

    def exe_umi_patho(self, sorted_bam):
        cldict = self.cldict
        sampd = self.sampd
        Patho_dir = cldict.Patho_dir
        ldelim = cldict.ldelim
        project_id = sampd.project_id
        outdir = Patho_dir + ldelim + project_id
        umidir = outdir + ldelim + "umidir"
        umi_log_dir = umidir + ldelim + "logdir"
        collapse_type = "coordinate"
        uminorm_path = cldict.uminorm_path
        umi_coll_len_script = cldict.umi_coll_len_script
        sampleid_refacc = self.get_sampleid_refacc()
        print("sorted bam: " + sorted_bam)
        print("collapse_type: " + collapse_type)
        print("sampleid_refacc: " + sampleid_refacc)
        print("uminorm_path: " + uminorm_path)
        print("umidir: " + umidir)
        umi_cmd = uminorm_path + " -i " + sorted_bam + " -o " + umidir + \
            " -p " + sampleid_refacc + " -c " + collapse_type
        print("umi_cmd starting: " + umi_cmd)
        call(umi_cmd.split())
        print("umi_cmd ended.")
        umi_logfile = umi_log_dir + ldelim + sampleid_refacc + "_coll_len.txt"
        umi_coll_len_cmd = "Rscript " + umi_coll_len_script + " -i " + \
            umi_logfile + " -o " + umi_log_dir + " -p " + sampleid_refacc
        print("umi_coll_len_cmd starting: " + umi_coll_len_cmd)
        call(umi_coll_len_cmd.split())
        print("umi_coll_len_cmd ended") 
        samtools_path = cldict.samtools
        tdf_path = cldict.tdf_str
        out_bam_base_u = sampleid_refacc + '_u.bam'
        umi_bam_u = umidir + ldelim + out_bam_base_u
        Bam_path = cldict.Bam_path
        bamdir = Bam_path + ldelim + project_id
        bamdir_umi = bamdir + ldelim + "umidir"
        suf_str = "_se.bam"
        umi_bam = umi_bam_u.replace('_u.bam', suf_str)
        AlignerBwa.sort_bam(samtools_path, umi_bam_u, umi_bam)
        AlignerBwa.copy_bam_tdf(samtools_path, tdf_path, bamdir_umi, umi_bam)
        return umi_bam
