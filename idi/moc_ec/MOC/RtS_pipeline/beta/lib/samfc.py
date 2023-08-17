import os
import re
import os.path
from miscutils import exe_command
from miscutils import exe_command_stderr
from subprocess import call
from miscutils import copyLargeFile
import collections
import pandas as pd

class SamFC:
    def __init__(self, cldict, sampd):
        self.cldict = cldict
        self.sampd = sampd

    def exe_count_frags(self, sample_id, sorted_bam, outdir, proto = "allseq", end = "single", for_umi = False):
        # Get the location for the fragment count
        cldict = self.cldict
        ldelim = cldict.ldelim
        sam_fragcount = cldict.sam_fragcount
        infile = sorted_bam 
        sampd = self.sampd
        l_host_str = sampd.host_ref_str
        lc_lower = cldict.LC_method_val.lower()
        sample_id_s = sample_id + "_" + l_host_str
        outfile = outdir + ldelim + sample_id_s + "_read_count.txt"
        metrics = outdir + ldelim + sample_id_s + "_metrics.txt"
        count_frags_cmd = sam_fragcount + " -i " + infile + " -o " + \
            outfile + " -m " + metrics + " -p " + proto + " -e " + end
        if lc_lower == "allseq" and for_umi:
            count_frags_cmd += " --use_umi"
        print("Starting count_frags_cmd: " + count_frags_cmd)
        call(count_frags_cmd.split())
        return outfile, metrics

    def exe_frag_to_gene_count(self, sample_id, trans_count_file, outdir, for_umi = False):
        cldict = self.cldict
        ldelim = cldict.ldelim
        host_transcript_gene = cldict.host_transcript_gene 
        frag_to_gene_count = cldict.frag_to_gene_count
        sampd = self.sampd
        l_host_str = sampd.host_ref_str
        sample_id_s = sample_id + "_" + l_host_str
        lc_lower = cldict.LC_method_val.lower()
        out_file = outdir + ldelim + sample_id_s + ".counts"
        count_cmd = frag_to_gene_count + " -i " + trans_count_file + \
            " -m " + host_transcript_gene + " -o " + out_file
        if lc_lower == "allseq" and for_umi:
            count_cmd += " --use_umi"
            umilog_str = out_file + ".umilog"
            count_cmd += " -u " + umilog_str
        print("Starting frag_to_gene_count cmd: " + count_cmd)
        call(count_cmd.split())
        return out_file     

    def exe_detailed_metrics(self, sample_id, out_gene_count, basic_metrics, outdir, for_umi = False):
        cldict = self.cldict
        metrics_gen = cldict.metrics_gen
        ldelim = cldict.ldelim
        sampd = self.sampd
        l_host_str = sampd.host_ref_str
        sample_id_s = sample_id + "_" + l_host_str
        detailed_metrics = outdir + ldelim + sample_id_s + \
            ".metrics"
        met_cmd = metrics_gen + " -i " + out_gene_count + " -m " + \
            basic_metrics + " -o " + detailed_metrics
        lc_lower = cldict.LC_method_val.lower()
        if lc_lower == "allseq" and for_umi:
            met_cmd += " --use_umi"
            sample_id_s = sample_id + "_" + l_host_str
            out_file = outdir + ldelim + sample_id_s + ".counts"
            umilog_str = out_file + ".umilog"
            met_cmd += " -u " + umilog_str
        print("Starting detailed metrics file generation.")
        print("met_cmd: " + met_cmd)
        call(met_cmd.split())
        return detailed_metrics
        
        
