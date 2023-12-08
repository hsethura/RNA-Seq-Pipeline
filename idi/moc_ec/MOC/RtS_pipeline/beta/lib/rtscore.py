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
from subprocess import call
from alignerutils import sort_by_qname
from alignerutils import sort_bam

# This would execute the core part for RtS (RNATag-Seq) pipeline.

class RtSCore:
    def __init__(self, cldict, sampd):
        self.cldict = cldict
        self.sampd = sampd
        self.mergo = Merger(cldict, sampd)
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
        self.jlco = CounterJL(cldict, sampd)
        print("Created JLCounter object")


    def merged_single_name(self, sample_id, suffix, merge_dir, ldelim):
        outfile_R = merge_dir + ldelim + sample_id + "_" + suffix + ".fastq"
        cldict = self.cldict
        gzip_merged = cldict.gzip_merged
        if gzip_merged:
            outfile_R += ".gz"
        elif not os.path.isfile(outfile_R):
            outfile_R += ".gz"

        if not os.path.isfile(outfile_R):
            raise StandardError("Merged file does not exist: " + outfile_R)

        return outfile_R


    def fastq_to_bam_single(self):
        cldict = self.cldict
        sampd = self.sampd
        mergo = self.mergo
        bwao = self.bwao

        ldelim = cldict.ldelim
        Patho_dir = cldict.Patho_dir
        merge_dir = cldict.Merge_dir
        project_id = sampd.project_id
        Bam_path = cldict.Bam_path
        bamdir = Bam_path + ldelim + project_id
        outdir = Patho_dir + ldelim + project_id
        sample_id = sampd.sample_id
        ref_acc = sampd.ref_acc_str


        suffix = "R"

        print("sample_id: " + sample_id)
        print("outdir: " + outdir)
        read_file = self.merged_single_name(sample_id, suffix, merge_dir, ldelim)
        sorted_bam = bwao.exe_bwa_single(sample_id, ref_acc, read_file, suffix, outdir, bamdir)
        return sorted_bam

    def fastq_to_bam_paired(self):
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


        suf1 = "R1"
        suf2 = "R2"

        print("sample_id: " + sample_id)
        print("outdir: " + outdir)
        read1_file = self.merged_single_name(sample_id, suf1, merge_dir, ldelim)
        read2_file = self.merged_single_name(sample_id, suf2, merge_dir, ldelim)
        print ("just going to execute bwa paired")
        # sorted_bam = bwao.exe_bwa_paired(sample_id, ref_acc, read1_file, read2_file,
            # suf1, suf2, outdir, bamdir)
        print ("finished executing sorted bam")
        return sorted_bam
            
 
    def sample_id_to_bam_patho(self):
        cldict = self.cldict
        sampd = self.sampd
        Patho_dir = cldict.Patho_dir
        ref_acc = sampd.ref_acc_str
        ldelim = cldict.ldelim
        project_id = sampd.project_id
        outdir = Patho_dir + ldelim + project_id
        sample_id = sampd.sample_id
        pair_suffix = None
        Read_pairing_val = cldict.Read_pairing_val
        if Read_pairing_val == 'SINGLE':
            pair_suffix = 'se'
        elif Read_pairing_val == 'PAIRED':
            pair_suffix = 'pe'
        outbamfile = outdir + ldelim + sample_id + "_" + ref_acc + "_" + \
            pair_suffix + ".bam"
        return outbamfile

    def exe_patho(self):
        cldict = self.cldict
        meto = self.meto
        bwao = self.bwao
        sampd = self.sampd
        Patho_dir = cldict.Patho_dir
        ldelim = cldict.ldelim
        ldel = cldict.ldelim
        project_id = sampd.project_id
        samtools_path = cldict.samtools
        tdf_path = cldict.tdf_str

        Bam_path = cldict.Bam_path
        bamdir = Bam_path + ldelim + project_id

        outdir = Patho_dir + ldelim + project_id
        Results_path = cldict.Results_path
        picard_outdir = Results_path + ldel + project_id + ldel + "bam_metrics_patho" 
        nodup_picard_outdir = Results_path + ldel + project_id + ldel + "nodupdir" + ldel + "bam_metrics_patho"
        sample_id = sampd.sample_id
        ref_acc = sampd.ref_acc_str
        Read_pairing_val = cldict.Read_pairing_val
        rm_rts_dup = cldict.rm_rts_dup
        jlco = self.jlco

        sorted_bam = ''
        do_align = cldict.do_align
        if do_align:
            if Read_pairing_val == 'SINGLE':
                sorted_bam = self.fastq_to_bam_single()
            elif Read_pairing_val == 'PAIRED':
                print("Started fastq_to_bam_paired from exe_patho in RtSCore")
                sorted_bam = self.fastq_to_bam_paired()
            else:
                raise StandardError("Unknown Read_pairing: " + Read_pairing_val)
            
        else:
            sorted_bam = self.sample_id_to_bam_patho()

        do_count = cldict.do_count
        if do_count:
            read_counter = cldict.read_counter
            if read_counter == 'featureCounts':
                fco = CounterFC(cldict, sampd)
                print ("Created featureCount object") 
                countfile_s_str = fco.exe_featureCounts(sample_id, ref_acc, sorted_bam, outdir, "s")  
                countfile_as_str = fco.exe_featureCounts(sample_id, ref_acc, sorted_bam, outdir, "as")  
                countfile_str = outdir + ldelim + sample_id + "_" + ref_acc + ".counts"
                CounterFC.combine_s_as(countfile_s_str, countfile_as_str, countfile_str)
            elif read_counter == 'JL_counter':
                if Read_pairing_val == 'SINGLE':
                    jlco.count_single(sample_id, ref_acc, sorted_bam, outdir, strand_rev = "N")
                elif Read_pairing_val == 'PAIRED':
                    jlco.count_paired(sample_id, ref_acc, sorted_bam, outdir)
                fnafile = bwao.get_patho_ref_path(ref_acc)
                meto.exe_picard_metrics(sorted_bam, fnafile, picard_outdir)
            else:
                raise StandardError('No read counter found for rts pathogen side.')

            if rm_rts_dup:
                # Do the duplicate removal to a new directory called no_dup_dir
                no_dup_dir = outdir + ldelim + "nodupdir"
                # Create a new bam file in the no_dup_dir after removal of the 
                # pcr duplicates from pcr; ideally that folder should also 
                # contain the metrics related to dup removal
                dup_marked_sorted_bam = self.mark_dup_reads(sorted_bam) 
                nodup_sorted_bam = self.remove_dup_reads(dup_marked_sorted_bam)
                if Read_pairing_val == 'SINGLE':
                    jlco.count_single(sample_id, ref_acc, nodup_sorted_bam, no_dup_dir, strand_rev = "N")
                elif Read_pairing_val == 'PAIRED':
                    jlco.count_paired(sample_id, ref_acc, nodup_sorted_bam, no_dup_dir)
                fnafile = bwao.get_patho_ref_path(ref_acc)
                meto.exe_picard_metrics(nodup_sorted_bam, fnafile, nodup_picard_outdir)
                # We have to copy the bam file to appropriate place
                nodup_bamdir = bamdir + ldelim + "nodupdir"
                AlignerBwa.copy_bam_tdf(samtools_path, tdf_path, nodup_bamdir, nodup_sorted_bam)  


    def remove_dup_reads(self, dm_sorted_bam):
        cldict = self.cldict
        samtools_path = cldict.samtools
        remove_dup_script = cldict.remove_dup_script
        
        # 1. sort wrt queryname
        dm_unsorted_bam = dm_sorted_bam.replace("_dm.bam", "_u_dm.bam") 
        if not os.path.isfile(dm_unsorted_bam):
            sort_by_qname(samtools_path, dm_sorted_bam, dm_unsorted_bam)
            print("sort by qname: " + dm_sorted_bam + " to " + dm_unsorted_bam)
        
        # 2. remove the pair of marked reads 0x400
        nodup_unsorted_bam = dm_unsorted_bam.replace("_u_dm.bam", "_u.bam")
        dup_bam = dm_unsorted_bam.replace("_u_dm.bam", "_u_dup.bam")
        # Call the c++ program remove_dup
        remove_dup_cmd = remove_dup_script + " -i " + dm_unsorted_bam +\
            " -o " + nodup_unsorted_bam + " -d " + dup_bam
        print("remove_dup_cmd starts: " + remove_dup_cmd)
        call(remove_dup_cmd.split())
        print("remove_dup_cmd ends")
        
        # 3. sort the no_dup_unsorted_bam bam file
        nodup_sorted_bam = nodup_unsorted_bam.replace("_u.bam", "_pe.bam")
        if not os.path.isfile(nodup_sorted_bam):
            sort_bam(samtools_path, nodup_unsorted_bam, nodup_sorted_bam)
            print("sort by coordinate: " + nodup_unsorted_bam + " to " +\
                nodup_sorted_bam)
        os.remove(dm_unsorted_bam)
        print("Removed dm_unsorted_bam: " + dm_unsorted_bam)
        os.remove(dm_sorted_bam)
        print("Remove dm_sorted_bam: " + dm_sorted_bam)
        return nodup_sorted_bam

    def mark_dup_reads(self, input_bam):
        # Get the output directory which has been created before
        # Do not remove the duplicates, just let picard mark them using read is PCR or optical 
        # duplicate (0x400). We shall remove them later in pair
        cldict = self.cldict
        ldelim = cldict.ldelim 
        picard_bindir = cldict.picard_bindir
        picard_jar = picard_bindir + "/picard.jar"

        bamdir = os.path.dirname(input_bam)
        no_dup_dir = bamdir + ldelim + "nodupdir"
        bam_file = os.path.basename(input_bam)
        output_bam_ori = no_dup_dir + ldelim + bam_file
        # Output bam is the dup marked one
        output_bam = output_bam_ori.replace("_pe.bam", "_dm.bam")

        metrics_txt = output_bam.replace('.bam', '_rm_dup_met.txt')

        # Then rename the old bam file

        mem_str = "-Xmx6G"
        dup_str = "MarkDuplicates"
        picard_cmd = "java " + mem_str + " -jar " + picard_jar + " " + dup_str + " I=" + input_bam +\
            " O=" + output_bam + " M=" + metrics_txt + " REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT"
        print("start: picard_cmd: " + picard_cmd)
        call(picard_cmd.split())
        print("end: picard_cmd:")
        return (output_bam)

    
    def sample_id_to_unsorted_bam_path(self):
        sampd = self.sampd
        lhost_ref_str = sampd.host_ref_str
        cldict = self.cldict
        ldelim = cldict.ldelim
        Host_dir = cldict.Host_dir
        project_id = sampd.project_id
        outdir = Host_dir + ldelim + project_id
        sample_id = sampd.sample_id
        unsorted_bam_path = outdir + ldelim + sample_id + "_" + lhost_ref_str + "_u.bam"
        if not os.path.exists(unsorted_bam_path):
            raise StandardError("File not found: " + unsorted_bam_path)
        return unsorted_bam_path

    def sample_id_to_sorted_bam_path(self, pair_type):
        sampd = self.sampd
        lhost_ref_str = sampd.host_ref_str
        cldict = self.cldict
        ldelim = cldict.ldelim
        Host_dir = cldict.Host_dir
        project_id = sampd.project_id
        outdir = Host_dir + ldelim + project_id
        sample_id = sampd.sample_id
        unsorted_bam_path = outdir + ldelim + sample_id + "_" + lhost_ref_str + "_" + pair_type + ".bam"
        if not os.path.exists(unsorted_bam_path):
            raise StandardError("File not found: " + unsorted_bam_path)
        return unsorted_bam_path

    def host_fastq_to_bam_paired(self):
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
        do_align = cldict.do_align
        do_count = cldict.do_count

        suf1 = "R1"
        suf2 = "R2"

        print("sample_id: " + sample_id)
        print("outdir: " + outdir)
        read1_file = self.merged_single_name(sample_id, suf1, merge_dir, ldelim)
        read2_file = self.merged_single_name(sample_id, suf2, merge_dir, ldelim)
        if do_align:
            bbmapo.exe_host_paired(sample_id, read1_file, read2_file, outdir, bamdir)

        #bbmapo.exe_fragcount(sample_id, outdir)


    def host_bam_to_count(self, unsorted_bam):
        """
        """
        sampd = self.sampd
        sample_id = sampd.sample_id
        cldict = self.cldict
        Host_dir = cldict.Host_dir
        ldelim = cldict.ldelim
        project_id = sampd.project_id
        outdir = os.path.dirname(unsorted_bam)
        samfco = self.samfco
        proto = "rts"
        end = "paired"
        out_count, basic_metrics = samfco.exe_count_frags(sample_id, unsorted_bam, outdir, proto, end)
        out_gene_count = samfco.exe_frag_to_gene_count(sample_id, out_count, outdir)
        detailed_metrics_file = samfco.exe_detailed_metrics(sample_id, out_gene_count, basic_metrics, outdir)
       
    def get_unsorted_bam(self, sorted_bam):
        # Check if the unsorted bam exists and if not just generates it.
        cldict = self.cldict
        samtools_path = cldict.samtools
        unsorted_bam = None

        if not os.path.isfile(sorted_bam):
            lstr = "Bam file not found: " + sorted_bam
            raise StandardError(lstr)

        if sorted_bam.endswith("_pe.bam"):
            unsorted_bam = sorted_bam.replace("_pe.bam", "_u.bam")
        elif sorted_bam.endswith("_se.bam"):
            unsorted_bam = sorted_bam.replace("_se.bam", "_u.bam")

        if not os.path.isfile(unsorted_bam):
            sort_by_qname(samtools_path, sorted_bam, unsorted_bam)

        return (unsorted_bam)

    def exe_host(self):
        cldict = self.cldict
        sampd = self.sampd
        do_count = cldict.do_count
        Results_path = cldict.Results_path
        ldel = cldict.ldelim
        project_id = sampd.project_id
        ref_acc = sampd.ref_acc_str
        bbmapo = self.bbmapo
        rm_rts_dup = cldict.rm_rts_dup
        samtools_path = cldict.samtools
        tdf_path = cldict.tdf_str
        
        Bam_path = cldict.Bam_path
        bamdir = Bam_path + ldel + project_id
        
        Read_pairing_val = cldict.Read_pairing_val
        sorted_bam = None
        if Read_pairing_val == 'SINGLE':
            self.host_fastq_to_bam_single()
            sorted_bam = self.sample_id_to_sorted_bam_path("se")
        elif Read_pairing_val == 'PAIRED':
            self.host_fastq_to_bam_paired()
            sorted_bam = self.sample_id_to_sorted_bam_path("pe")
        else:
            raise ValueError('Incorrect value for Read_pairing_val: ' + Read_pairing_val) 
        
        host_fna_path = bbmapo.get_host_fna_path() 
        print("host_fna_path: " + host_fna_path)

        if do_count:
            unsorted_bam = self.sample_id_to_unsorted_bam_path()
            count_file = self.host_bam_to_count(unsorted_bam)
            picard_outdir = Results_path + ldel + project_id + ldel + \
                "bam_metrics_host" 
            meto = self.meto
            print("picard_outdir: " + picard_outdir)
            meto.exe_picard_metrics(sorted_bam, host_fna_path, picard_outdir)
        

            if rm_rts_dup:
                outdir = Results_path + ldel + project_id
                # Do the duplicate removal to a new directory called no_dup_dir
                no_dup_dir = outdir + ldel + "nodupdir"
                # Create a new bam file in the no_dup_dir after removal of the pcr duplicates from pcr; 
                # ideally that folder should also contain the metrics related to dup removal
                jlco = self.jlco
                dup_marked_sorted_bam = self.mark_dup_reads(sorted_bam) 
                nodup_sorted_bam = self.remove_dup_reads(dup_marked_sorted_bam)
                
                nodup_unsorted_bam = self.get_unsorted_bam(nodup_sorted_bam)
                count_file = self.host_bam_to_count(nodup_unsorted_bam)
                nodup_picard_outdir = no_dup_dir + ldel + "bam_metrics_host" 
                meto.exe_picard_metrics(nodup_sorted_bam, host_fna_path, nodup_picard_outdir)
                nodup_bamdir = bamdir + ldel + "nodupdir"
                AlignerBwa.copy_bam_tdf(samtools_path, tdf_path, nodup_bamdir, nodup_sorted_bam)  
   
  
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
                print("Started exe_patho from RtSCore")
                self.exe_patho()
        if do_host:
            if host_ref_str_l and host_ref_str_l != "unknown" and  host_ref_str_l != "none":
                self.exe_host()
 
