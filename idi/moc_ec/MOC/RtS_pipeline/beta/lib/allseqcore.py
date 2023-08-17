#!/usr/bin/env python
import os
from merger import Merger
from unimerger import UniMerger
from alignerbwa import AlignerBwa
from alignerbbmap import AlignerBBMap
from counterfc import CounterFC
from samfc import SamFC
from counterjl import CounterJL
from metrics import Metrics
from umi import UMI

# This would execute the core part for AllSeq (RNATag-Seq) pipeline.
# With the present form allseq is single-ended on the second read.

class AllSeqCore:
    def __init__(self, cldict, sampd):
        self.cldict = cldict
        self.sampd = sampd
        self.mergo = Merger(cldict, sampd)
        self.umio = UMI(cldict, sampd)
        self.meto = Metrics(cldict)
        
        lbwao = None
        lbbmapo = None
        lref_acc_str = sampd.ref_acc_str
        lhost_ref_str = sampd.host_ref_str 
        if lref_acc_str != "none":
            lbwao = AlignerBwa(cldict, sampd)
        if lhost_ref_str != "none":
            lbbmapo = AlignerBBMap(cldict, sampd)
        self.bwao = lbwao
        self.bbmapo = lbbmapo
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
        umio = self.umio
        Patho_dir = cldict.Patho_dir
        ldelim = cldict.ldelim
        project_id = sampd.project_id
        outdir = Patho_dir + ldelim + project_id
        sample_id = sampd.sample_id
        ref_acc = sampd.ref_acc_str
        do_align = cldict.do_align
        do_count = cldict.do_count
        do_umi_count = cldict.do_umi_count
        sorted_bam = None
        Results_path = cldict.Results_path
        picard_outdir = Results_path + ldelim + project_id + ldelim + \
            "bam_metrics_patho"

        # For debug: start
        #do_align = False
        #do_count = False
        print("do_umi_count: " + str(do_umi_count))
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
                jlco.count_single(sample_id, ref_acc, sorted_bam, outdir, strand_rev = "Y")
                fnafile = bwao.get_patho_ref_path(ref_acc)
                meto.exe_picard_metrics(sorted_bam, fnafile, picard_outdir)
            else:
                raise StandardError('No read counter found for allseq pathogen side.')
        if do_umi_count:
            umi_bam = umio.exe_umi_patho(sorted_bam)
            read_counter = cldict.read_counter
            if read_counter == 'JL_counter':
                jlco = CounterJL(cldict, sampd)
                print("Created JLCounter object for umi_bam")
                outdir_umi = outdir + ldelim + "umidir"
                picard_outdir_umi = Results_path + ldelim + project_id + ldelim + \
                    "umidir" + ldelim + "bam_metrics_patho"

                jlco.count_single(sample_id, ref_acc, umi_bam, outdir_umi, strand_rev = "Y")
                fnafile = bwao.get_patho_ref_path(ref_acc)
                meto.exe_picard_metrics(umi_bam, fnafile, picard_outdir_umi)


    def host_fastq_to_bam(self):
        """
        """
        sampd = self.sampd
        cldict = self.cldict
        Host_dir = cldict.Host_dir
        ldelim = cldict.ldelim
        project_id = sampd.project_id
        outdir = Host_dir + ldelim + project_id
        Bam_path = cldict.Bam_path
        bamdir = Bam_path + ldelim + project_id
        sample_id = sampd.sample_id
        bbmapo = self.bbmapo
        merge_dir = cldict.Merge_dir

        suf2 = "R2"

        print("sample_id: " + sample_id)
        print("outdir: " + outdir)
        read2_file = self.merged_single_name(sample_id, suf2, merge_dir, ldelim)
        # The assumption is that bbmap would behave similarly for single ended
        # allseq and single ended rts. We have to differentiate them
        # during subsequent read count step.
        outsorted, out_rpkm = bbmapo.exe_host_single(sample_id, read2_file, outdir, bamdir)
        #bbmapo.exe_fragcount(sample_id, outdir)
        return outsorted
      

    def host_bam_to_count(self, sorted_bam, for_umi = False): 
        """
        """
        sampd = self.sampd
        sample_id = sampd.sample_id
        cldict = self.cldict
        Host_dir = cldict.Host_dir
        ldelim = cldict.ldelim
        project_id = sampd.project_id
        outdir = None
        if for_umi:
            outdir = Host_dir + ldelim + project_id + ldelim + "umidir"
        else:
            outdir = Host_dir + ldelim + project_id 
        samfco = self.samfco
        out_count, basic_metrics = samfco.exe_count_frags(sample_id, sorted_bam, outdir, for_umi = for_umi)    
        out_gene_count = samfco.exe_frag_to_gene_count(sample_id, out_count, outdir, for_umi = for_umi)
        detailed_metrics_file = samfco.exe_detailed_metrics(sample_id, out_gene_count, basic_metrics, outdir, for_umi = for_umi)


    def sample_id_to_sorted_bam_path(self):
        sampd = self.sampd
        lhost_ref_str = sampd.host_ref_str
        cldict = self.cldict
        ldelim = cldict.ldelim
        Host_dir = cldict.Host_dir
        project_id = sampd.project_id
        outdir = Host_dir + ldelim + project_id
        sample_id = sampd.sample_id
        sorted_bam_path = outdir + ldelim + sample_id + "_" + lhost_ref_str + "_se.bam"
        return sorted_bam_path  
        
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
    

    def exe_host(self):
        cldict = self.cldict
        sampd = self.sampd
        do_align = cldict.do_align
        do_count = cldict.do_count
        Results_path = cldict.Results_path
        ldel = cldict.ldelim
        meto = self.meto
        bbmapo = self.bbmapo
        host_fna_path = bbmapo.get_host_fna_path()
        project_id = sampd.project_id
        do_umi_count = cldict.do_umi_count
        umio = self.umio
        Host_dir = cldict.Host_dir
        outdir = Host_dir + ldel + project_id

        if do_align:
            sorted_bam = self.host_fastq_to_bam()
        if do_count:
            sorted_bam = self.sample_id_to_sorted_bam_path()
            count_file = self.host_bam_to_count(sorted_bam)
            picard_outdir = Results_path + ldel + project_id + ldel + \
                "bam_metrics_host"
            print("sorted_bam: " + sorted_bam)
            print("host_fna_path: " + host_fna_path)
            print("picard_outdir: " + picard_outdir)
            meto.exe_picard_metrics(sorted_bam, host_fna_path, picard_outdir)
        if do_umi_count:
            sorted_bam = self.sample_id_to_sorted_bam_path()
            umi_bam = umio.exe_umi_host(sorted_bam)
            outdir_umi = outdir + ldel + "umidir"
            picard_outdir_umi = Results_path + ldel + project_id + ldel + \
                "umidir" + ldel + "bam_metrics_host"
            print("umi_bam: " + umi_bam)
            count_file_umi = self.host_bam_to_count(umi_bam, for_umi = True)
            print("umi_bam: " + umi_bam)
            print("host_fna_path: " + host_fna_path)
            print("umi_picard_outdir: " + picard_outdir_umi)
            meto.exe_picard_metrics(umi_bam, host_fna_path, picard_outdir_umi)

  
    def mainFunc(self):
        cldict = self.cldict
        do_host = cldict.do_host
        do_patho = cldict.do_patho
        sampd = self.sampd
        ref_acc_str = sampd.ref_acc_str
        ref_acc_str_l = ref_acc_str.lower()
        host_ref_str = sampd.host_ref_str
        host_ref_str_l = host_ref_str.lower()

        if do_patho:
            if ref_acc_str_l and ref_acc_str_l != "unknown" and ref_acc_str_l != "none":
                self.exe_patho()
        if do_host:
            if host_ref_str_l and host_ref_str_l != "unknown" and  host_ref_str_l != "none":
                self.exe_host()
 
