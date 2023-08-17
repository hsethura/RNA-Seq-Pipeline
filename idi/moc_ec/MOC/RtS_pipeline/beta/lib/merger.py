import re
import os
import gzip
import shutil
import yaml
from miscutils import exe_command
from miscutils import file_concat
from subprocess import call

class Merger:
    def __init__(self, cldict, sampd):
        self.cldict = cldict
        self.sampd = sampd

    def R_trim(self, outfile_R_ori, outfile_R):
        cldict = self.cldict
        trim_script = cldict.trim_script
        trim_cmd = trim_script + " " + outfile_R_ori
        print("trim_cmd: " + trim_cmd + ", outfile_R: " + outfile_R)
        exe_command(trim_cmd, outfile_R)
        return


    def trim_allseq(self, outfile_R_ori, outfile_R):
        cldict = self.cldict
        allseq_trim_script = cldict.allseq_trim_script
        allseq_trim_lim = cldict.allseq_trim_lim
        allseq_con_a = cldict.AllSeq_con_A
        allseq_trim_len = cldict.AllSeq_trim_len
        allseq_trim_cmd = allseq_trim_script + " -i " + outfile_R_ori + " -o " + outfile_R + " -c " + str(allseq_con_a) + " -r " + str(allseq_trim_len)
        print("allseq_trim_cmd: " + allseq_trim_cmd)
        call(allseq_trim_cmd.split())
        
    def get_out_file_suffix(self, suffix):
        cldict = self.cldict
        out_suffix = None
        if suffix == "R1":
            out_suffix = "R1"
        elif suffix == "R2":
            out_suffix = "R2"
        elif suffix == "R":
            out_suffix = "R"
        elif suffix == "R1_R":
            out_suffix = "R1"
        elif suffix == "R2_R":
            out_suffix = "R2"
        elif suffix == "R1_G":
            out_suffix = "R1"
        elif suffix == "R2_G":
            out_suffix = "R2"
        else:
            raise StandardError("Unknown suffix: " + suffix)
        print("from get_out_file_suffix, suffix: " + suffix + ", out_suffix: " + out_suffix)
        return out_suffix

 
    def get_seq_file_suffix(self, suffix):
        cldict = self.cldict
        seq_suffix = None
        if suffix == "R1":
            seq_suffix = cldict.suffix_s1
        elif suffix == "R2":
            seq_suffix = cldict.suffix_s2
        elif suffix == "R":
            seq_suffix = cldict.suffix_ne
        else:
            raise StandardError("Unknown suffix: " + suffix)
        print("from get_seq_file_suffix, suffix: " + suffix + ", seq_suffix: " + seq_suffix)
        return seq_suffix


    # Have to do tomorrow
    def merge_valid_barcode_single(self, lsample, Merge_dir):
        sampd = self.sampd
        cldict = self.cldict

        prefix_str = sampd.prefix_set
        bc_str = sampd.bc_set
        Split_dir = cldict.Split_dir
        ldelim = cldict.ldelim
        no_bc_split = cldict.no_bc_split

        ldelim1 = "\s*;\s*"
        prefix_lst =  re.split(ldelim1, prefix_str)

        ## Get bc_lst
        bc_lst = list()
        if bc_str and not no_bc_split:
            bc_lst = re.split(ldelim1, bc_str)

        templates = list()
        for lprefix in prefix_lst:
            if bc_lst:
                for lbc in bc_lst:
                    ltem = lprefix + "_" + lbc
                    templates.append(ltem)
            else:
                templates.append(lprefix)

        print("templates: " + ', '.join(templates))

        ## Get insuffix
        insuffix = "_valid_bc.txt"
        print("from merge_valid_barcode_single, templates: ")
        print(templates)
        print("insuffix: " + insuffix)

        ## Get bc_files
        bc_files = list()
        for lprefix in templates:
            lfile = Split_dir + ldelim + lprefix + insuffix
            print("lfile: " + lfile)
            if not os.path.isfile(lfile):
                lfile += ".gz"
                if not os.path.isfile(lfile):
                    raise OSError("File not found: " + lfile)
            bc_files.append(lfile)
        print("sample: " + lsample)
        print("Now print bc_files")
        print(bc_files)

        ## Get outfile_R_ori name
        print("bc_files: " + ', '.join(bc_files))
        out_suffix = "_valid_bc.txt"
        merged_bcfile = Merge_dir + ldelim + lsample + out_suffix
        print("merged_bcfile: " + merged_bcfile)

        ## Concatenamte R_files
        file_concat(bc_files, merged_bcfile)
        return (merged_bcfile)

    def exe_trim_5p_3p(self, outfile_R_ori, outfile_R, trim_5p, trim_3p):
        cldict = self.cldict
        read_trimmer_path = cldict.read_trimmer_path
        read_trimmer_cmd = read_trimmer_path + " -i " + outfile_R_ori + " -o " + outfile_R + " --trim_5p " + str(trim_5p) + " --trim_3p " + str(trim_3p)
        print("read_trimmer_cmd (trim_5p_3p) : " + read_trimmer_cmd)
        call(read_trimmer_cmd.split())
     
    def exe_keep_5p(self, outfile_R_ori, outfile_R, keep_5p):
        cldict = self.cldict
        read_trimmer_path = cldict.read_trimmer_path
        read_trimmer_cmd = read_trimmer_path + " -i " + outfile_R_ori + " -o " + outfile_R + " --keep_5p " + str(keep_5p)
        print("read_trimmer_cmd (keep_5p) : " + read_trimmer_cmd)
        call(read_trimmer_cmd.split())

    def exe_keep_3p(self, outfile_R_ori, outfile_R, keep_3p):  
        cldict = self.cldict
        read_trimmer_path = cldict.read_trimmer_path
        read_trimmer_cmd = read_trimmer_path + " -i " + outfile_R_ori + " -o " + outfile_R + " --keep_3p " + str(keep_3p)
        print("read_trimmer_cmd (keep_3p) : " + read_trimmer_cmd)
        call(read_trimmer_cmd.split())

    # If a trimmin is done then outfile_R would contained the trimmed version
    # and outfile_R_ori would be removed, otherwise outfile_R_ori would be
    # renamed to outfile_R

    def trim_fastq(self, outfile_R_ori, outfile_R,
            do_R_trim, 
            do_allseq_trim, 
            trim_5p, trim_3p, keep_5p, keep_3p):

        # The assumption is the trimming of one protocol will not mix with other
        # For example, R trimming is basically for RtS-TS
        # allseq trimming is obviously for allseq
        # trim_5p, trim_3p, keep_5p, keep_3p is mainly for experimental purpose
        # Will not interefere with one another
        
        command_count = 0
        if do_R_trim:
            command_count += 1
        if do_allseq_trim:
            command_count += 1
        if trim_5p > 0 or trim_3p > 0:
            command_count += 1
        if keep_5p >= 0:
            command_count += 1
        if keep_3p >= 0:
            command_count += 1

        if command_count > 1:
            raise StandardError("More than one command in trim_fastq: " + command_count)

        if trim_5p > 0 or trim_3p > 0:
            print("Value of trim_5p: " + str(trim_5p))
            print("Value of trim_3p: " + str(trim_3p))
            print("So going for 5p/3p trimming.")
            self.exe_trim_5p_3p(outfile_R_ori, outfile_R, trim_5p, trim_3p)

        if keep_5p >= 0:
            self.exe_keep_5p(outfile_R_ori, outfile_R, keep_5p)

        if keep_3p >= 0:
             self.exe_keep_3p(outfile_R_ori, outfile_R, keep_3p)
            

        # If R trimming is enabled, do the R trimming
        print("Value of do_R_trim: " + str(do_R_trim))
        if do_R_trim:
            self.R_trim(outfile_R_ori, outfile_R)

        # Trimming for allseq read2
        print("Value of do_allseq_trim: " + str(do_allseq_trim))
        if do_allseq_trim:
            self.trim_allseq(outfile_R_ori, outfile_R)

        # If none of the trimming is applied, then just rename
        if command_count == 0:
            os.rename(outfile_R_ori, outfile_R)
            print("Renamed: " + outfile_R_ori + " to " + outfile_R)
        elif command_count == 1:
            os.remove(outfile_R_ori)
            print("Removed: " + outfile_R_ori)
            


    def merge_all_files(self, lsample, suffix, Merge_dir):
        sampd = self.sampd
        cldict = self.cldict
        no_bc_split = cldict.no_bc_split
        soft_link_from_input = cldict.soft_link_from_input


        prefix_str = sampd.prefix_set
        bc_str = sampd.bc_set
        Split_dir = cldict.Split_dir
        ldelim = cldict.ldelim

        ldelim1 = "\s*;\s*"
        prefix_lst =  re.split(ldelim1, prefix_str)

        ## Get bc_lst
        bc_lst = list()
        if bc_str and not no_bc_split:
            bc_lst = re.split(ldelim1, bc_str)

        templates = list()
        for lprefix in prefix_lst:
            if bc_lst:
                for lbc in bc_lst:
                    ltem = lprefix + "_" + lbc
                    templates.append(ltem)
            else:
                templates.append(lprefix)

        print("templates: " + ', '.join(templates))

        ## Get insuffix
        insuffix = None
        if soft_link_from_input:
            insuffix = self.get_seq_file_suffix(suffix)
        else:
            insuffix = "_" + suffix + ".fastq"
        print("from merge_fastq_single, templates: ")
        print(templates)
        print("suffix: " + suffix + ", insuffix: " + insuffix)

        ## Get R_files
        R_files = list()
        for lprefix in templates:
            lfile = Split_dir + ldelim + lprefix + insuffix
            print("lfile: " + lfile)
            if not os.path.isfile(lfile):
                lfile += ".gz"
                if not os.path.isfile(lfile):
                    raise OSError("File not found: " + lfile)
            R_files.append(lfile)
        print("sample: " + lsample)
        print("Now print R_files")
        print(R_files)

        ## Get outfile_R_ori name
        print("R_files: " + ', '.join(R_files))
        out_suffix = self.get_out_file_suffix(suffix)
        outfile_R_ori = Merge_dir + ldelim + lsample + "_" + out_suffix + "_ori.fastq"
        print("outfile_R_ori: " + outfile_R_ori)

        ## Concatenamte R_files
        file_concat(R_files, outfile_R_ori)

        ## Remove the split files
        remove_splitted = cldict.remove_splitted
        if remove_splitted:
            for lfile in R_files:
                os.remove(lfile)

        return (outfile_R_ori)

    def merge_fastq_single(self, lsample, suffix, Merge_dir, 
            do_R_trim, 
            do_allseq_trim, 
            trim_5p = 0, trim_3p = 0, keep_5p = -1, keep_3p = -1):
        ## Get outfile_R_ori
        cldict = self.cldict
        ldelim = cldict.ldelim
        outfile_R_ori = self.merge_all_files(lsample, suffix, Merge_dir)

        ## Trim the outfile_R_ori and put in outfile_R
        out_suffix = self.get_out_file_suffix(suffix)
        outfile_R = Merge_dir + ldelim + lsample + "_" + out_suffix + ".fastq"
        print("outfile_R: " + outfile_R)
        self.trim_fastq(outfile_R_ori, outfile_R, do_R_trim, do_allseq_trim, trim_5p, trim_3p,
            keep_5p, keep_3p)


        ## Gzip the merged ones if requested
        gzip_merged = cldict.gzip_merged
        # Now compress the outfile for the next step
        outfile_R_final = None
        if not gzip_merged:
            outfile_R_final = outfile_R
        else:
            outfile_R_final = outfile_R + ".gz"
            print("unzipped file: " + outfile_R)
            print("zipped file: " + outfile_R_final)
            with open(outfile_R, 'rb') as f_in:
                with gzip.open(outfile_R_final, 'wb') as f_out:
                   shutil.copyfileobj(f_in, f_out)
            os.remove(outfile_R) 

        print("gzip_merged: " + str(gzip_merged)) 
        print("Returning final file: " + outfile_R_final)
        return outfile_R_final

