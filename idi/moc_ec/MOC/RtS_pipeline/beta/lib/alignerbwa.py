import os
import os.path
from miscutils import exe_command
from subprocess import call
from miscutils import copyLargeFile

class AlignerBwa:
    def __init__(self, cldict, sampd):
        self.cldict = cldict
        self.sampd = sampd

    def get_patho_ref_path(self, ref_acc_str):
        cldict = self.cldict
        ldelim = cldict.ldelim
        Data_dir = cldict.Data_dir

        patho_ref_path = Data_dir + ldelim + ref_acc_str + ".fna"
        return patho_ref_path

    def exe_bwa_paired(self, sample_id, ref_acc, read1_file, read2_file, 
            suf1, suf2, outdir, bamdir):

        cldict = self.cldict
        samtools_path = cldict.samtools
        tdf_path = cldict.tdf_str
        bwa_mem = cldict.bwa_mem
        picard_bindir = cldict.picard_bindir 
        picard_jar = picard_bindir + "/picard.jar"
        rm_rts_dup = cldict.rm_rts_dup 
        paired_only_patho = cldict.paired_only_patho
        paired_only_script = cldict.paired_only_script
        outsamfile = None
        
        if (bwa_mem):
            outsamfile = self.gen_sam_paired_bwamem(sample_id, ref_acc, read1_file,\
                 read2_file, outdir)
        else:
            outsamfile = self.gen_sam_paired(sample_id, ref_acc, read1_file,\
                 read2_file, suf1, suf2, outdir)
        bam_suf = "pe"
        outsorted = AlignerBwa.sam_to_bam_sorted(outsamfile, samtools_path,\
            picard_jar, tdf_path, bamdir, bam_suf, rm_rts_dup, \
            paired_only_patho = paired_only_patho, paired_only_script = paired_only_script)
        do_count = cldict.do_count
        if not do_count:
            if os.path.isfile(outsamfile):
                os.remove(outsamfile)
        return outsorted
 
    def exe_bwa_single(self, sample_id, ref_acc, read_file, suf, outdir, bamdir):

        cldict = self.cldict
        samtools_path = cldict.samtools
        tdf_path = cldict.tdf_str
        bwa_mem = cldict.bwa_mem
        outsamfile = None
        picard_bindir = cldict.picard_bindir
        picard_jar = picard_bindir + "/picard.jar"
        rm_rts_dup = cldict.rm_rts_dup 
        if bwa_mem:
            outsamfile = self.gen_sam_single_bwamem(sample_id, ref_acc, read_file, outdir)
        else:
            outsamfile = self.gen_sam_single(sample_id, ref_acc, read_file, suf, outdir)
        bam_suf = "se"
        outsorted = AlignerBwa.sam_to_bam_sorted(outsamfile, samtools_path,\
            picard_jar, tdf_path, bamdir, bam_suf, rm_rts_dup)
        do_count = cldict.do_count
        if not do_count:
            if os.path.isfile(outsamfile):
                os.remove(outsamfile)
        return outsorted

    @staticmethod
    def get_paired_only_bam(outbamfile, paired_only_script):
        outbamfile_ori = outbamfile.replace('_u.bam', '_u_ori.bam')
        os.rename(outbamfile, outbamfile_ori)
        outbamfile_broken = outbamfile.replace('_u.bam', '_u_br.bam') 
        paired_only_cmd = paired_only_script + " -i " + outbamfile_ori + " -o " + outbamfile + " -b " + outbamfile_broken
        print("paired_only_cmd started: " + paired_only_cmd)
        call(paired_only_cmd.split())
        print("paired_only_cmd end")
        os.remove(outbamfile_ori)
        return outbamfile
        

    @staticmethod
    def sam_to_bam_sorted(outsamfile, samtools_path, picard_jar, tdf_path, bamdir, suf, rm_rts_dup, paired_only_patho = False, paired_only_script = None):
        outbamfile_ori = outsamfile.replace('.sam', '_u.bam')
        AlignerBwa.sam_to_bam(samtools_path, outsamfile, outbamfile_ori)
        outbamfile = None
        # If paired_only_patho modify the bam
        if paired_only_patho:
            outbamfile = AlignerBwa.get_paired_only_bam(outbamfile_ori, paired_only_script)
            # Fix for bug_72: --paired_only_patho not working
            if os.path.isfile(outsamfile):
                print("Removing from sam_to_bam_sorted: " + outsamfile)
                os.remove(outsamfile)
        else:
            outbamfile = outbamfile_ori
        suf_str = "_" + suf + ".bam"
        outsorted = outbamfile.replace('_u.bam', suf_str)
        AlignerBwa.sort_bam(samtools_path, outbamfile, outsorted)
        AlignerBwa.copy_bam_tdf(samtools_path, tdf_path, bamdir, outsorted) 
        # Decided to keep both the sorted and unsorted bam files and delete 
        # sam files. This is a trade-off between time and space.
        #os.remove(outbamfile)
        return outsorted
    
    

    def sample_to_sai(self, outdir, sample_id, ref_acc, suf, read_file, patho_path):
     
        cldict = self.cldict
        ldelim = cldict.ldelim
        bwa_path = cldict.bwa_path
     
        R_sai = outdir + ldelim + sample_id + "_" + ref_acc + "_" + suf + ".sai"
        read_cmd =  bwa_path + " aln " + patho_path + " " + read_file
        print("bwa aln cmd: " + read_cmd)
        exe_command(read_cmd, R_sai)
        print("bwa aln cmd done.")
        return R_sai

     
    def gen_sam_paired(self, sample_id, ref_acc, read1_file, read2_file, suf1, suf2, outdir): 
        cldict = self.cldict
        ldelim = cldict.ldelim
        bwa_path = cldict.bwa_path
        samtools_path = cldict.samtools
        patho_path = self.get_patho_ref_path(ref_acc)

        R1_sai = self.sample_to_sai(outdir, sample_id, ref_acc, suf1, read1_file, patho_path)
        R2_sai = self.sample_to_sai(outdir, sample_id, ref_acc, suf2, read2_file, patho_path)

        sampecmd = bwa_path + " sampe " + patho_path + " " + R1_sai + " " + R2_sai +\
            " " + read1_file + " " + read2_file
        outsamfile = outdir + ldelim + sample_id + "_" + ref_acc + ".sam"
        print("bwa sampe cmd: " + sampecmd)
        exe_command(sampecmd, outsamfile)
        print("bwa sampe cmd done.")
        os.remove(R1_sai)
        os.remove(R2_sai)
        return outsamfile
   
    def gen_sam_paired_bwamem(self, sample_id, ref_acc, read1_file, read2_file, outdir): 
        cldict = self.cldict
        ldelim = cldict.ldelim
        bwa_path = cldict.bwa_path
        patho_path = self.get_patho_ref_path(ref_acc)
        # bwa mem ref.fa read1.fq read2.fq > aln-pe.sam

        bwamem_cmd = bwa_path + " mem " + patho_path + " " + read1_file + " " + read2_file   
        print("Aligning with bwa mem: " + bwamem_cmd)
        outsamfile = outdir + ldelim + sample_id + "_" + ref_acc + ".sam"
        exe_command(bwamem_cmd, outsamfile)
        return outsamfile
  
    
    def gen_sam_single(self, sample_id, ref_acc, read_file, suf, outdir): 
        cldict = self.cldict
        ldelim = cldict.ldelim
        bwa_path = cldict.bwa_path
        samtools_path = cldict.samtools
        patho_path = self.get_patho_ref_path(ref_acc)
        R_sai = self.sample_to_sai(outdir, sample_id, ref_acc, suf, read_file, patho_path)

        sampecmd = bwa_path + " samse " + patho_path + " " + R_sai + " " + read_file
        print("sampe cmd: " + sampecmd)
        outsamfile = outdir + ldelim + sample_id + "_" + ref_acc + ".sam"
        exe_command(sampecmd, outsamfile)
        print("sampe cmd done.")
        os.remove(R_sai)
        print("removed " + R_sai)
        return outsamfile

    def gen_sam_single_bwamem(self, sample_id, ref_acc, read_file, outdir): 
        cldict = self.cldict
        ldelim = cldict.ldelim
        bwa_path = cldict.bwa_path
        patho_path = self.get_patho_ref_path(ref_acc)
        # bwa mem ref.fa reads.fq > aln-se.sam

        bwamem_cmd = bwa_path + " mem " + patho_path + " " + read_file
        print("Aligning with bwa mem: " + bwamem_cmd)
        outsamfile = outdir + ldelim + sample_id + "_" + ref_acc + ".sam"
        exe_command(bwamem_cmd, outsamfile)
        return outsamfile

    @staticmethod
    def sam_to_bam(samtools_path, outsamfile, outbamfile):
        bamcmd = samtools_path + ' view -b -S ' + outsamfile
        print("sam to bam cmd: " + bamcmd)
        exe_command(bamcmd, outbamfile)
        print("sam to bam cmd done.")

    @staticmethod
    def sort_bam(samtools_path, outbamfile, outsorted):
        sortedcmd = samtools_path + ' sort -o ' + outsorted + " " + outbamfile
        print("sort bam cmd: " + sortedcmd)
        call(sortedcmd.split())
        print("sort bam cmd done.")
    
    @staticmethod
    def copy_bam_tdf(samtools_str, tdf_str, bam_path, outsorted):
        ldelim = "/"
        # Copy the bam file to bam_path

        out_bam_base = os.path.basename(os.path.normpath(outsorted))
        copy_bam_path = bam_path + ldelim + out_bam_base
        copyLargeFile(outsorted, copy_bam_path)
        bai_cmd = samtools_str + " index " + copy_bam_path
        print("bai_cmd: " + bai_cmd)
        call(bai_cmd.split())
        print("copy_bam_path: " + copy_bam_path)
        print("tdf_str: " + tdf_str)
        tdf_cmd = "java -jar " + tdf_str + " " + copy_bam_path
        call(tdf_cmd.split())

    @staticmethod
    def copy_bam(samtools_str, bam_path, outsorted):
        ldelim = "/"
        # Copy the bam file to bam_path

        out_bam_base = os.path.basename(os.path.normpath(outsorted))
        copy_bam_path = bam_path + ldelim + out_bam_base
        copyLargeFile(outsorted, copy_bam_path)
        bai_cmd = samtools_str + " index " + copy_bam_path
        print("bai_cmd: " + bai_cmd)
        call(bai_cmd.split())
        print("copy_bam_path: " + copy_bam_path)




