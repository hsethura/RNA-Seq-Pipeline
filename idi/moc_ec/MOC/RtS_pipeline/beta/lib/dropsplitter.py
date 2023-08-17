import os
import yaml
import os.path
import ntpath

class DropSplitter:
    def __init__(self, confd):
        self.confd = confd


    def get_split_command(self, sample_id, project_id, ref_acc_str, 
            sort_type = "aligned"):
        confd = self.confd
        sam_splitter_path = confd.sam_splitter_path
        Merge_dir = confd.Merge_dir
        Patho_dir = confd.Patho_dir
        ldel = confd.ldelim
        merge_real = Merge_dir + ldel + "real"
        bcdist_dir = Patho_dir + ldel + project_id + ldel + "bcdist"
        real_sam_dir = Patho_dir + ldel + project_id + ldel + "real"
        outdir = Patho_dir + ldel + project_id
        Patho_drops_dir = outdir + ldel + "drops"
        sample_drop_dir = Patho_drops_dir + ldel + sample_id
        if not os.path.exists(sample_drop_dir):
            os.makedirs(sample_drop_dir)
        prefix_str = sample_id + "_" + ref_acc_str
        lbc_sam_file = merge_real + ldel + sample_id + "_valid_bc.txt"
        lbc_dist_file = bcdist_dir + ldel + sample_id + "_" + ref_acc_str + "_" + sort_type + ".txt"
        real_sam_file = real_sam_dir + ldel + sample_id + "_" + ref_acc_str + "_u.bam"
        sam_splitter_cmd = sam_splitter_path + " --bc_sam " + lbc_sam_file  +\
            " --bc_dist " + lbc_dist_file + " --sam_file " + real_sam_file +\
            " -p " + prefix_str + " -o " + sample_drop_dir
        #print("sample_id: " + sample_id + ", sam_splitter_cmd: " + sam_splitter_cmd)
        return(sam_splitter_cmd)



