import os
import yaml
import os.path
import ntpath

class GarbageSorter:
    def __init__(self, confd):
        self.confd = confd

    def get_outfile(self, lprefix):
        confd = self.confd
        ldelim = confd.ldelim
        Log_dir = confd.Log_dir
        loutfile = Log_dir + ldelim + lprefix + "_out.txt"
        return loutfile

    def get_errfile(self, lprefix):
        confd = self.confd
        ldelim = confd.ldelim
        Log_dir = confd.Log_dir
        loutfile = Log_dir + ldelim + lprefix + "_err.txt"
        return loutfile


    def get_fastq_file(self, lprefix, suffix):
        confd = self.confd
        prefix_dict_out = confd.prefix_dict_out
        prefix_dict = yaml.load(open(prefix_dict_out))
        lfile = prefix_dict[lprefix] + suffix
        if not os.path.isfile(lfile):
            lfile += ".gz"
            if not os.path.isfile(lfile):
                raise OSError("File not found: " + lfile)
        return lfile


    def get_outread_name(self, lprefix, read_num, is_gz, out_type):
        confd = self.confd
        Split_dir = confd.Split_dir
        ldelim = confd.ldelim
        lfile = None
        if read_num == 0:
            lfile = Split_dir + ldelim + lprefix + "_R" + out_type + ".fastq"
        else:
            lfile = Split_dir + ldelim + lprefix + "_R" + str(read_num) + out_type + ".fastq"
        if is_gz:
            lfile += ".gz"
        return lfile

    def get_outlog_name(self, lprefix):
        confd = self.confd
        Split_dir = confd.Split_dir
        ldelim = confd.ldelim
        lfile = Split_dir + ldelim + lprefix + "_log.txt"
        return lfile

    def get_outbcfile_name(self, lprefix):
        confd = self.confd
        Split_dir = confd.Split_dir
        ldelim = confd.ldelim
        lfile = Split_dir + ldelim + lprefix + "_real_bc.txt"
        return lfile

    def get_bds_cmd(self, lprefix):
        confd = self.confd
        gs_script = confd.garbage_sorter_path
        suffix_s1 = confd.suffix_s1
        suffix_s2 = confd.suffix_s2
        lfile1 = self.get_fastq_file(lprefix, suffix_s1)
        lfile2 = self.get_fastq_file(lprefix, suffix_s2)
        is_gz = False
        lfile1_r = self.get_outread_name(lprefix, 1, is_gz, "_R")
        lfile2_r = self.get_outread_name(lprefix, 2, is_gz, "_R")
        lfile1_g = self.get_outread_name(lprefix, 1, is_gz, "_G")
        lfile2_g = self.get_outread_name(lprefix, 2, is_gz, "_G")
        log_file = self.get_outlog_name(lprefix)
        loutfile = self.get_outfile(lprefix)
        lerrfile = self.get_errfile(lprefix)
        loutbcfile = self.get_outbcfile_name(lprefix)
        gs_cmd = gs_script + " --in1 " + lfile1 + " --in2 " + lfile2 + \
            " --out1 " + lfile1_r + " --out2 " + lfile2_r + \
            " --out_g1 " + lfile1_g + " --out_g2 " + lfile2_g + \
            " --log " + log_file + " --b " + loutbcfile + " 1> " + loutfile + " 2> " + lerrfile + "\n"
        return gs_cmd


    def get_ind_cmd(self, lprefix):
        confd = self.confd
        bc1_file = confd.indrops_bc1_file
        bc2_file = confd.indrops_bc2_file
        Split_dir = confd.Split_dir
        gs_script = confd.garbage_sorter_ind_path
        suffix_s1 = confd.suffix_s1
        suffix_s2 = confd.suffix_s2
        lfile1 = self.get_fastq_file(lprefix, suffix_s1)
        lfile2 = self.get_fastq_file(lprefix, suffix_s2)
        log_file = self.get_outlog_name(lprefix)
        loutfile = self.get_outfile(lprefix)
        lerrfile = self.get_errfile(lprefix)
        loutbcfile = self.get_outbcfile_name(lprefix)
        gs_cmd = gs_script + " --bc1 " + bc1_file + " --bc2 " + bc2_file +\
            " --in1 " + lfile1 + " --in2 " + lfile2 + " --outdir " +\
            Split_dir + " --prefix " + lprefix + " 1> " + loutfile + " 2> " +\
            lerrfile + "\n"
        return gs_cmd

    def get_garbage_sorter_cmd(self, lprefix):
        confd = self.confd
        drop_type = confd.drop_type
        if drop_type == "bdropseq":
            return self.get_bds_cmd(lprefix)
        elif drop_type == "indrops":
            return self.get_ind_cmd(lprefix)
        else:
            raise StandardError("Illegal drop_type: " + drop_type)
        
        

        

