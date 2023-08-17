import os
import yaml
import os.path
import ntpath

class DropCounter:
    def __init__(self, confd):
        self.confd = confd

    def get_bc_pre_dir_list(self, sample_id, project_id):
        confd = self.confd
        Patho_dir = confd.Patho_dir
        ldel = confd.ldelim
        outdir = Patho_dir + ldel + project_id
        Patho_drops_dir = outdir + ldel + "drops"
        sample_drop_dir = Patho_drops_dir + ldel + sample_id
        bc_pre_dir_list_file = sample_drop_dir + ldel + "bc_pre_dir_list.txt"
        bc_pre_dir_list = list()
        with open(bc_pre_dir_list_file) as inf:
            for line in inf:
                bc_pre_dir_list.append(line.rstrip())
        return (bc_pre_dir_list)

    def get_bc_pre_dir(self, sample_id, project_id, bc_pre_str):
        confd = self.confd
        Patho_dir = confd.Patho_dir
        ldel = confd.ldelim
        outdir = Patho_dir + ldel + project_id
        Patho_drops_dir = outdir + ldel + "drops"
        sample_drop_dir = Patho_drops_dir + ldel + sample_id
        bc_pre_dir = sample_drop_dir + ldel + bc_pre_str
        return (bc_pre_dir)




