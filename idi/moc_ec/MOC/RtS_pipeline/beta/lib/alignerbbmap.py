import os
import re
import os.path
from miscutils import exe_command
from miscutils import exe_command_stderr
from subprocess import call
from miscutils import copyLargeFile
from alignerutils import sam_to_bam_sorted
import collections
import pandas as pd

class AlignerBBMap:
    def __init__(self, cldict, sampd):
        self.cldict = cldict
        self.sampd = sampd
        l_host_str = sampd.host_ref_str
        if re.search('human', l_host_str, re.IGNORECASE):
            self.host_transcript_gene = cldict.human_transcript_gene
        elif re.search('mouse', l_host_str, re.IGNORECASE):
            self.host_transcript_gene = cldict.mouse_transcript_gene
        elif re.search('rabbit', l_host_str, re.IGNORECASE):
            self.host_transcript_gene = cldict.rabbit_transcript_gene
        elif re.search('zebrafish', l_host_str, re.IGNORECASE):
            self.host_transcript_gene = cldict.zebrafish_transcript_gene
        elif re.search('fungal', l_host_str, re.IGNORECASE):
            self.host_transcript_gene = cldict.fungal_transcript_gene
        elif re.search('candida_albicans_SC5314', l_host_str, re.IGNORECASE):
            self.host_transcript_gene = cldict.candida_albicans_SC5314_transcript_gene
        else:
            raise ValueError('Wrong Host_reference: ' + l_host_str)
        print("host_transcript_gene: " + self.host_transcript_gene)

    def get_host_fna_path(self):
        cldict = self.cldict
        host_dbpath = cldict.host_dbpath
        host_aligner = cldict.host_aligner 
        sampd = self.sampd
        ldelim = cldict.ldelim
        
        l_host_str = sampd.host_ref_str
        if re.search('human', l_host_str, re.IGNORECASE):
            self.host_ref_str = cldict.human_ref_str
        elif re.search('mouse', l_host_str, re.IGNORECASE):
            self.host_ref_str = cldict.mouse_ref_str
        elif re.search('rabbit', l_host_str, re.IGNORECASE):
            self.host_ref_str = cldict.rabbit_ref_str
        elif re.search('zebrafish', l_host_str, re.IGNORECASE):
            self.host_ref_str = cldict.zebrafish_ref_str
        elif re.search('fungal', l_host_str, re.IGNORECASE):
            self.host_ref_str = cldict.fungal_ref_str
        elif re.search('candida_albicans_SC5314', l_host_str, re.IGNORECASE):
            self.host_ref_str = cldict.candida_albicans_ref_str
        else:
            raise ValueError('Wrong Host_reference: ' + l_host_str)
        host_fna_path = host_dbpath + ldelim + "data" + \
            ldelim + self.host_ref_str + ".fna"
        self.host_fna_path = host_fna_path
        return host_fna_path

    def get_host_ref_path(self):
        cldict = self.cldict
        host_dbpath = cldict.host_dbpath
        host_aligner = cldict.host_aligner 
        sampd = self.sampd
        ldelim = cldict.ldelim
        
        l_host_str = sampd.host_ref_str
        if re.search('human', l_host_str, re.IGNORECASE):
            self.host_ref_str = cldict.human_ref_str
        elif re.search('mouse', l_host_str, re.IGNORECASE):
            self.host_ref_str = cldict.mouse_ref_str
        elif re.search('rabbit', l_host_str, re.IGNORECASE):
            self.host_ref_str = cldict.rabbit_ref_str
        elif re.search('zebrafish', l_host_str, re.IGNORECASE):
            self.host_ref_str = cldict.zebrafish_ref_str
        elif re.search('fungal', l_host_str, re.IGNORECASE):
            self.host_ref_str = cldict.fungal_ref_str
        elif re.search('candida_albicans_SC5314', l_host_str, re.IGNORECASE):
            self.host_ref_str = cldict.candida_albicans_ref_str
        else:
            raise ValueError('Wrong Host_reference: ' + l_host_str)
        host_ref_path = host_dbpath + ldelim + host_aligner + \
            ldelim + self.host_ref_str
        self.host_ref_path = host_ref_path
        return host_ref_path


    @staticmethod
    def exe_bbmap_single(bbmap_path, read2_file, out_samfile, 
        host_ref_path, host_thread_count, host_total_memory_lim_str, 
        out_covstats, out_covhist, out_bincov, out_rpkm, out_summary):
        cmd_str = bbmap_path + " -Xmx" + host_total_memory_lim_str + "g" + \
            " ambiguous=random" + \
            " in=" + read2_file + " path=" + host_ref_path + \
            " out=" + out_samfile + " threads=" + str(host_thread_count) + \
            " covstats=" + out_covstats + " covhist=" + out_covhist + \
            " bincov=" + out_bincov + " rpkm=" + out_rpkm 
        print("bbmap_cmd: " + cmd_str)
        exe_command_stderr(cmd_str, out_summary)


    @staticmethod
    def exe_bbmap_paired(bbmap_path, read1_file, read2_file, out_samfile, 
        host_ref_path, host_thread_count, host_total_memory_lim_str, 
        out_covstats, out_covhist, out_bincov, out_rpkm, out_summary, paired_only_host):
        cmd_str = bbmap_path + " -Xmx" + host_total_memory_lim_str + "g" + \
            " ambiguous=random" + \
            " in=" + read1_file + " in2=" + read2_file + \
            " path=" + host_ref_path + " out=" + out_samfile + " threads=" + \
            str(host_thread_count)  + \
            " covstats=" + out_covstats + " covhist=" + out_covhist + \
            " bincov=" + out_bincov + " rpkm=" + out_rpkm 
        if paired_only_host:
            cmd_str += " pairedonly=t"
        print("bbmap_cmd: " + cmd_str)
        exe_command_stderr(cmd_str, out_summary)

    def get_total_memory_str(self, host_thread_count):
        cldict = self.cldict
        host_memory = cldict.host_memory
        host_thread_count_int = int(host_thread_count)
        host_memory_int = int(host_memory)
        host_total_memory = host_thread_count_int * host_memory_int
        # We allocate 2 GB for system usage. For example, for a total
        # allocated memory of 28GB, bbmap gets 26GB.
        host_total_memory_lim = host_total_memory -2
        host_total_memory_lim_str = str(host_total_memory_lim)
        return host_total_memory_lim_str

    def exe_host_paired(self, sample_id, read1_file, read2_file, outdir, bamdir):
        cldict = self.cldict    
        bbmap_path = cldict.bbmap_path
        sampd = self.sampd
        paired_only_host = cldict.paired_only_host
 
        host_thread_count = cldict.host_thread_count
        host_total_memory_lim_str = self.get_total_memory_str(host_thread_count)

        l_host_str = sampd.host_ref_str
        sample_id_s = sample_id + "_" + l_host_str
        ldelim = cldict.ldelim
        out_samfile = outdir + ldelim + sample_id_s + ".sam"
        out_covstats = outdir + ldelim + sample_id_s + "_covstats.txt"
        out_covhist = outdir + ldelim + sample_id_s + "_covhist.txt"
        out_bincov = outdir + ldelim + sample_id_s + "_bincov.txt"
        out_rpkm = outdir + ldelim + sample_id_s + "_rpkm.txt"
        out_summary = outdir + ldelim + sample_id_s + "_run_summary.txt"
        
        host_ref_path = self.get_host_ref_path()
        AlignerBBMap.exe_bbmap_paired(bbmap_path, read1_file, read2_file, 
            out_samfile, host_ref_path, host_thread_count, 
            host_total_memory_lim_str, out_covstats, out_covhist, out_bincov, 
            out_rpkm, out_summary, paired_only_host)
        # Now get the bbmap specific metric generation
        #exe_bbmap_metric()
        # Generate sorted bam 
        # move the bam and tdf files
        samtools_path = cldict.samtools
        tdf_path = cldict.tdf_str
        bam_suf = "pe"
        outsorted = sam_to_bam_sorted(out_samfile, samtools_path,\
            tdf_path, bamdir, bam_suf)
        os.remove(out_samfile)
        return outsorted, out_rpkm

    # Here the assumption is that, bbmap does not care about the 
    # orientation for single ended reads
    def exe_host_single(self, sample_id, read2_file, outdir, bamdir):
        cldict = self.cldict    
        sampd = self.sampd
        l_host_str = sampd.host_ref_str
        bbmap_path = cldict.bbmap_path

        host_thread_count = cldict.host_thread_count
        host_total_memory_lim_str = self.get_total_memory_str(host_thread_count)

        ldelim = cldict.ldelim
        sample_id_s = sample_id + "_" + l_host_str
        out_samfile = outdir + ldelim + sample_id_s + ".sam"
        out_covstats = outdir + ldelim + sample_id_s + "_covstats.txt"
        out_covhist = outdir + ldelim + sample_id_s + "_covhist.txt"
        out_bincov = outdir + ldelim + sample_id_s + "_bincov.txt"
        out_rpkm = outdir + ldelim + sample_id_s + "_rpkm.txt"
        out_summary = outdir + ldelim + sample_id_s + "_run_summary.txt"
        
        host_ref_path = self.get_host_ref_path()
        AlignerBBMap.exe_bbmap_single(bbmap_path, read2_file, 
            out_samfile, host_ref_path, host_thread_count, 
            host_total_memory_lim_str, out_covstats, out_covhist, 
            out_bincov, out_rpkm, out_summary)
        # Now get the bbmap specific metric generation
        #exe_bbmap_metric()
        # Generate sorted bam 
        # move the bam and tdf files
        samtools_path = cldict.samtools
        tdf_path = cldict.tdf_str
        bam_suf = "se"
        outsorted = sam_to_bam_sorted(out_samfile, samtools_path,\
            tdf_path, bamdir, bam_suf)
        os.remove(out_samfile)
        return outsorted, out_rpkm

    def exe_fragcount(self, sample_id, outdir):
        cldict = self.cldict
        ldelim = cldict.ldelim
        fragcount_file = outdir + ldelim + sample_id + "_fragcount.txt"
        rpkm_file = outdir + ldelim + sample_id + "_rpkm.txt"
        host_transcript_gene = self.host_transcript_gene
        trans_to_gene = self.get_trans_to_gene(host_transcript_gene)   
        self.get_gene_to_frag(rpkm_file, fragcount_file, trans_to_gene) 

    def get_trans_to_gene(self, trans_to_gene_file):
        trans_to_gene = dict()
        reg_str = "\t"
        with open(trans_to_gene_file) as tgf:
            for line0 in tgf:
                line = line0.rstrip()
                parts = re.split(reg_str, line)
                trans = parts[0]
                gene = parts[1]
                trans_to_gene[trans] = gene
                # print("transcript: " + trans + " gene: " + gene)    
        return trans_to_gene


    def get_gene_to_frag(self, rpkm_file, fragcount_file, trans_to_gene):
        print("rpkm_file: " + rpkm_file)
        print("fragcount_file: " + fragcount_file)
        skip_count = 5
        frag_count_pos = 6
        reg_str = "\s+"
        gene_to_frag = dict()
        with open(rpkm_file) as ff:
            for _ in xrange(skip_count):
                next(ff)
            for line0 in ff:
                line = line0.rstrip()
                parts = re.split(reg_str, line)
                trans = parts[0]
                frag_count_str = parts[frag_count_pos]
                frag_count = int(frag_count_str)
                gene = trans_to_gene[trans]
                if gene not in gene_to_frag:
                    gene_to_frag[gene] = frag_count
                elif gene in gene_to_frag:
                    gene_to_frag[gene] += frag_count
        #print(gene_to_frag)
        gene_to_frag_s = collections.OrderedDict(sorted(gene_to_frag.items()))
        gene_to_frag_df = pd.DataFrame(gene_to_frag_s.items())
        filename = os.path.basename(os.path.normpath(rpkm_file))
        file_pre = filename.replace("_rpkm.txt", "")
        gene_to_frag_df.columns = ['gene', file_pre]
        gene_series = gene_to_frag_df.iloc[:,0]
        fragc_series = gene_to_frag_df.iloc[:,1]
        gene_to_frag_df.to_csv(fragcount_file, sep = '\t', index = False)
        #return gene_series, fragc_series

 

    
    

    





