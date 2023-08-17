import os
import yaml
import os.path
import ntpath

class Splitter:
    def __init__(self, confd, ldict_file = None):
        self.confd = confd
        self.ldict_file = ldict_file

    def get_split_command(self, lprefix):
        confd = self.confd
        Read_pairing = confd.Read_pairing_val.upper()
        if Read_pairing == 'PAIRED':
            return self.get_split_command_paired(lprefix)
        elif Read_pairing == 'SINGLE':
            return self.get_split_command_single(lprefix)
        else:
            raise StandardError("Unknown Read_pairing: " + Read_pairing)

    def create_soft_links(self, lprefix):
        confd = self.confd
        Read_pairing = confd.Read_pairing_val.upper()
        if Read_pairing == 'PAIRED':
            return self.create_soft_links_paired(lprefix)
        elif Read_pairing == 'SINGLE':
            return self.create_soft_links_single(lprefix)

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

    def get_outread_name(self, lprefix, read_num, is_gz):
        confd = self.confd
        Split_dir = confd.Split_dir
        ldelim = confd.ldelim
        lfile = None
        if read_num == 0:
            lfile = Split_dir + ldelim + lprefix + "_R.fastq"
        else:
            lfile = Split_dir + ldelim + lprefix + "_R" + str(read_num) + ".fastq" 
        if is_gz:
            lfile += ".gz"
        return lfile

    def path_leaf(self, path):
        head, tail = ntpath.split(path)
        return tail or ntpath.basename(head)

    def get_softname_split(self, lfile):
        confd = self.confd
        Split_dir = confd.Split_dir
        ldelim = confd.ldelim
        lleaf = self.path_leaf(lfile)
        outfile = Split_dir + ldelim +  lleaf
        return outfile
        
    def create_soft_links_single(self, lprefix):
        confd = self.confd
        suffix_ne = confd.suffix_ne
        lfile_ne = self.get_fastq_file(lprefix, suffix_ne)

        is_gz_r_ne = False
        if lfile_ne.endswith(".gz"):
            is_gz_r_ne = True
        #soft_ne = self.get_outread_name(lprefix, 0, is_gz_r_ne)
        soft_ne = self.get_softname_split(lfile_ne)

        print("Create soft link, src: " + lfile_ne + " dest: " + soft_ne)
        if os.path.lexists(soft_ne):
            os.remove(soft_ne)
        os.symlink(lfile_ne, soft_ne) 

    def create_soft_links_paired(self, lprefix):
        confd = self.confd
        lbc_splitter = confd.lbc_splitter
        ldelim = confd.ldelim
        suffix_s1 = confd.suffix_s1
        suffix_s2 = confd.suffix_s2
        lfile1 = self.get_fastq_file(lprefix, suffix_s1)
        lfile2 = self.get_fastq_file(lprefix, suffix_s2)
        is_gz_r1 = False
        if lfile1.endswith(".gz"):
            is_gz_r1 = True
        #soft1 = self.get_outread_name(lprefix, 1, is_gz_r1)
        soft1 = self.get_softname_split(lfile1)

        is_gz_r2 = False
        if lfile2.endswith(".gz"):
            is_gz_r2 = True
        #soft2 = self.get_outread_name(lprefix, 2, is_gz_r2)
        soft2 = self.get_softname_split(lfile2)

        print("Create soft link, src: " + lfile1 + " dest: " + soft1)
        if os.path.lexists(soft1):
            os.remove(soft1)
        os.symlink(lfile1, soft1) 
        print("Create soft link, src: " + lfile2 + " dest: " + soft2)
        if os.path.lexists(soft2):
            os.remove(soft2)
        os.symlink(lfile2, soft2) 
        

    def get_split_command_paired(self, lprefix):
        confd = self.confd
        lbc_splitter = confd.lbc_splitter
        ldelim = confd.ldelim
        Split_dir = confd.Split_dir
        Log_dir = confd.Log_dir
        suffix_s1 = confd.suffix_s1
        suffix_s2 = confd.suffix_s2
        ldict_file = self.ldict_file
        lfile1 = self.get_fastq_file(lprefix, suffix_s1)
        lfile2 = self.get_fastq_file(lprefix, suffix_s2)
        loutfile = self.get_outfile(lprefix)
        lerrfile = self.get_errfile(lprefix)
        cmd_str = lbc_splitter + " -d " + ldict_file + " --file1 " +  lfile1 +\
            " --file2 " + lfile2 + " -p " + lprefix + " -o " + Split_dir
        cmd_str2 = cmd_str + " 1> " + loutfile + " 2> " + lerrfile + "\n"
        return cmd_str2
        
    def get_split_command_single(self, lprefix):
        confd = self.confd
        lbc_splitter = confd.lbc_splitter
        ldelim = confd.ldelim
        Split_dir = confd.Split_dir
        Log_dir = confd.Log_dir
        suffix_ne = confd.suffix_ne
        ldict_file = self.ldict_file

        lfile = self.get_fastq_file(lprefix, suffix_ne)
        loutfile = self.get_outfile(lprefix)
        lerrfile = self.get_errfile(lprefix)
        cmd_str = lbc_splitter + " -d " + ldict_file + " --file " +  lfile +\
            " -p " + lprefix + " -o " + Split_dir
        cmd_str2 = cmd_str + " 1> " + loutfile + " 2> " + lerrfile + "\n"
        return cmd_str2


