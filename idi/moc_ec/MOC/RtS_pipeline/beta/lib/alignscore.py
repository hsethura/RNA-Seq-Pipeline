#!/usr/bin/env python
import os
from merger import Merger
from unimerger import UniMerger
from alignerbwa import AlignerBwa
from counterfc import CounterFC
from samfc import SamFC
from counterjl import CounterJL
from metrics import Metrics

# This would be one of the most flexible single read alignment core. It would 
# be able to align either R1 or R2 (specifically for pathogen) with any suffix
# and from and to any location. The design of this module was specifically 
# motivated by dropseq, as we need to align for real and garbage reads and
# possibly against either R1 or R2.

class AlignSCore:
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

    def fastq_to_bam(self, suf, subdir):
        
        cldict = self.cldict
        sampd = self.sampd
        mergo = self.mergo
        bwao = self.bwao

        ldelim = cldict.ldelim
        Patho_dir = cldict.Patho_dir
        project_id = sampd.project_id
        Bam_path = cldict.Bam_path
        bamdir = Bam_path + ldelim + project_id + ldelim + subdir
        outdir = Patho_dir + ldelim + project_id + ldelim + subdir
        sample_id = sampd.sample_id
        ref_acc = sampd.ref_acc_str
        merge_dir = cldict.Merge_dir + ldelim + subdir


        print("sample_id: " + sample_id)
        print("outdir: " + outdir)
        read_file = self.merged_single_name(sample_id, suf, merge_dir, ldelim)
        sorted_bam = bwao.exe_bwa_single(sample_id, ref_acc, read_file,
            suf, outdir, bamdir)
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
        use_dropseq = cldict.use_dropseq

        # For debug: start
        #do_align = False
        #do_count = False
        # For debug : end
        sorted_bam_real = None
        sorted_bam_garbage = None
        lc_lower = cldict.lc_lower

        if do_align: 
            if use_dropseq:
                align_readnum = cldict.align_readnum
                lsuf = 'R2'
                if align_readnum == 'alr1':
                    lsuf = 'R1'
                sorted_bam_garbage = self.fastq_to_bam(lsuf, 'garbage')
                sorted_bam_real = self.fastq_to_bam(lsuf, 'real')
            elif lc_lower == 'alignr2':
                sorted_bam = self.fastq_to_bam('R2', '')
            elif lc_lower == 'alignr1':
                sorted_bam = self.fastq_to_bam('R1', '')
            else:
                sorted_bam = self.fastq_to_bam('R2', '')
        else:
            if use_dropseq:
                sorted_bam_garbage = self.sample_id_to_bam_patho('garbage')
                sorted_bam_real = self.sample_id_to_bam_patho('real')
            else:
                sorted_bam = self.sample_id_to_bam_path('')

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
                if use_dropseq:
                    jlco = CounterJL(cldict, sampd)
                    print("Created JLCounter object")
                    count_strand_rev = cldict.count_strand_rev

                    outdir_r = outdir + ldelim + "real"
                    picard_outdir_r = Results_path + ldelim + project_id + ldelim + \
                        "real/bam_metrics_patho"
                    jlco.count_single(sample_id, ref_acc, sorted_bam_real, outdir_r, strand_rev = count_strand_rev)
                    fnafile = bwao.get_patho_ref_path(ref_acc)
                    meto.exe_picard_metrics(sorted_bam_real, fnafile, picard_outdir_r)

                    outdir_g = outdir + ldelim + "garbage"
                    picard_outdir_g = Results_path + ldelim + project_id + ldelim + \
                        "garbage/bam_metrics_patho"
                    jlco.count_single(sample_id, ref_acc, sorted_bam_garbage, outdir_g, strand_rev = count_strand_rev)
                    fnafile = bwao.get_patho_ref_path(ref_acc)
                    meto.exe_picard_metrics(sorted_bam_garbage, fnafile, picard_outdir_g)

                else:
                    jlco = CounterJL(cldict, sampd)
                    print("Created JLCounter object")
                    count_strand_rev = cldict.count_strand_rev
                    jlco.count_single(sample_id, ref_acc, sorted_bam, outdir, strand_rev = count_strand_rev)
                    fnafile = bwao.get_patho_ref_path(ref_acc)
                    meto.exe_picard_metrics(sorted_bam, fnafile, picard_outdir)
            else:
                raise StandardError('No read counter found for allseq pathogen side.')

        
    def sample_id_to_bam_patho(self, subdir):
        cldict = self.cldict
        sampd = self.sampd
        Patho_dir = cldict.Patho_dir
        ref_acc = sampd.ref_acc_str
        ldelim = cldict.ldelim
        project_id = sampd.project_id
        outdir = Patho_dir + ldelim + project_id + ldelim + subdir
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
 
