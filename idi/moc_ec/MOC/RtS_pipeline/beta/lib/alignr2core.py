#!/usr/bin/env python
import os
from merger import Merger
from unimerger import UniMerger
from alignerbwa import AlignerBwa
from counterfc import CounterFC
from samfc import SamFC
from counterjl import CounterJL
from metrics import Metrics

# This would execute the core part for AllSeq (RNATag-Seq) pipeline.
# With the present form allseq is single-ended on the second read.

class AlignR2Core:
    def __init__(self, cldict, sampd):
        self.cldict = cldict
        self.sampd = sampd
        self.mergo = Merger(cldict, sampd)
        self.meto = Metrics(cldict)
        
        lbwao = None
        lbbmapo = None
        lref_acc_str = sampd.ref_acc_str
        if lref_acc_str != "none":
            lbwao = AlignerBwa(cldict, sampd)
        self.bwao = lbwao
        self.samfco = SamFC(cldict, sampd)
         

    def merged_single_name(self, sample_id, suffix, merge_dir, ldelim):
        outfile_R = merge_dir + ldelim + sample_id + "_" + suffix + ".fastq"
        cldict = self.cldict
        gzip_merged = cldict.gzip_merged
        if gzip_merged:
            outfile_R += ".gz"
        return outfile_R

    def fastq_to_bam(self):
        cldict = self.cldict
        sampd = self.sampd
        mergo = self.mergo
        bwao = self.bwao

        ldelim = cldict.ldelim
        Patho_dir = cldict.Patho_dir
        project_id = sampd.project_id
        Bam_path = cldict.Bam_path
        bamdir = Bam_path + ldelim + project_id
        outdir = Patho_dir + ldelim + project_id
        sample_id = sampd.sample_id
        ref_acc = sampd.ref_acc_str
        merge_dir = cldict.Merge_dir

        suf2 = "R2"

        print("sample_id: " + sample_id)
        print("outdir: " + outdir)
        read2_file = self.merged_single_name(sample_id, suf2, merge_dir, ldelim)
        sorted_bam = bwao.exe_bwa_single(sample_id, ref_acc, read2_file,
            suf2, outdir, bamdir)
        return sorted_bam
    
    def exe_patho(self):
        cldict = self.cldict
        sampd = self.sampd
        bwao = self.bwao
        meto = self.meto
        Patho_dir = cldict.Patho_dir
        ldelim = cldict.ldelim
        project_id = sampd.project_id
        outdir = Patho_dir + ldelim + project_id
        sample_id = sampd.sample_id
        ref_acc = sampd.ref_acc_str
        do_align = cldict.do_align
        do_count = cldict.do_count
        sorted_bam = None
        Results_path = cldict.Results_path
        picard_outdir = Results_path + ldelim + project_id + ldelim + \
            "bam_metrics_patho"

        # For debug: start
        #do_align = False
        #do_count = False
        # For debug : end

        if do_align: 
            sorted_bam = self.fastq_to_bam()
        else:
            sorted_bam = self.sample_id_to_bam_patho()

        if do_count:
            read_counter = cldict.read_counter
            if read_counter == 'featureCounts':
                # Actually we do not use featureCounts for a long time. We have
                # switched to Jonathan's counter.
                fco = CounterFC(cldict, sampd)
                print ("Created featureCount object") 
                countfile_s_str = fco.exe_featureCounts(sample_id, ref_acc, sorted_bam, outdir, "s")  
                countfile_as_str = fco.exe_featureCounts(sample_id, ref_acc, sorted_bam, outdir, "as")  
                countfile_str = outdir + ldelim + sample_id + "_" + ref_acc + ".counts"
                CounterFC.combine_s_as(countfile_s_str, countfile_as_str, countfile_str)
                fco.get_fc_metrics(sample_id, ref_acc, outdir)
            elif read_counter == 'JL_counter':
                jlco = CounterJL(cldict, sampd)
                print("Created JLCounter object")
                count_strand_rev = cldict.count_strand_rev
                jlco.count_single(sample_id, ref_acc, sorted_bam, outdir, strand_rev = count_strand_rev)
                fnafile = bwao.get_patho_ref_path(ref_acc)
                meto.exe_picard_metrics(sorted_bam, fnafile, picard_outdir)
            else:
                raise StandardError('No read counter found for allseq pathogen side.')

        
    def sample_id_to_bam_patho(self):
        cldict = self.cldict
        sampd = self.sampd
        Patho_dir = cldict.Patho_dir
        ref_acc = sampd.ref_acc_str
        ldelim = cldict.ldelim
        project_id = sampd.project_id
        outdir = Patho_dir + ldelim + project_id
        sample_id = sampd.sample_id
        outbamfile = outdir + ldelim + sample_id + "_" + ref_acc + "_se.bam"
        return outbamfile
    
  
    def mainFunc(self):
        cldict = self.cldict
        do_patho = cldict.do_patho
        sampd = self.sampd
        ref_acc_str = sampd.ref_acc_str
        ref_acc_str_l = ref_acc_str.lower()

        if do_patho:
            if ref_acc_str_l and ref_acc_str_l != "unknown" and ref_acc_str_l != "none":
                self.exe_patho()
 
