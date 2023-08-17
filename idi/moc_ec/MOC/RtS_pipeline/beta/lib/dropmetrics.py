#!/usr/bin/env python

import os
import argparse
import yaml
import os.path
import ntpath
import pandas as pd

from dictmap import DictMap


class DropMetrics:
    def __init__(self, args):
        self.config_log_file = args.config_log_file
        cldict_d = yaml.load(open(self.config_log_file))
        cldict = DictMap(cldict_d)
        self.confd = cldict

        self.sample_id = args.sample_id
        self.project_id = args.project_id
        self.ref_acc_str = args.ref_acc_str


    def get_sample_drop_dir(self, sample_id, project_id):
        confd = self.confd
        Patho_dir = confd.Patho_dir
        ldel = confd.ldelim
        outdir = Patho_dir + ldel + project_id
        Patho_drops_dir = outdir + ldel + "drops"
        sample_drop_dir = Patho_drops_dir + ldel + sample_id
        return (sample_drop_dir)

    def create_metrics_dir(self, ldir):
        confd = self.confd
        ldel = confd.ldelim
        metrics_dir = ldir + ldel + "metrics"
        if not os.path.exists(metrics_dir):
            os.makedirs(metrics_dir)
        return (metrics_dir)

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

    def get_count_path(self, bc_pre_dir, bc_pre_path, sample_id, ref_acc):
        confd = self.confd
        ldel = confd.ldelim
        metrics_dir = bc_pre_path + ldel + "metrics"
        count_tag = sample_id + "_" + ref_acc + "_" + bc_pre_dir 
        count_file = count_tag + ".counts"
        count_path = metrics_dir + ldel + count_file
        return (count_path)

    def get_metrics_path(self, bc_pre_dir, bc_pre_path, sample_id, ref_acc):
        confd = self.confd
        ldel = confd.ldelim
        metrics_dir = bc_pre_path + ldel + "metrics"
        metrics_tag = sample_id + "_" + ref_acc + "_" + bc_pre_dir 
        metrics_file = metrics_tag + "_metrics.txt"
        metrics_path = metrics_dir + ldel + metrics_file
        return (metrics_path)
              
    def combine_count_paths(self, count_path_lst, outfile): 
        # https://stackoverflow.com/questions/23668427/pandas-joining-multiple-dataframes-on-columns
        dfs = []
        name_tag = None
        started = False
        geneid_ori = None

        for count_path in count_path_lst:
            print("Processing count_path: " + count_path)
            ltab = pd.read_csv(count_path, sep = '\t', header = 0)
            if not started:
                started = True
                geneid_ori = ltab.iloc[:,0]
                name_tag = geneid_ori.name
            else: 
                geneid = ltab.iloc[:,0]
                if not geneid.equals(geneid_ori):
                    raise ValueError('Two geneids are not equal')    
            dfs.append(ltab)

        if dfs:
            df_final = reduce(lambda left,right: pd.merge(left,right,on= name_tag), dfs)
            df_final.to_csv(outfile, sep = ',', index = False)
        else:
            open(outfile, 'w')

    def combine_metrics_paths(self, metrics_path_lst, outfile): 
        # https://stackoverflow.com/questions/23668427/pandas-joining-multiple-dataframes-on-columns
        dfs = []

        for metrics_path in metrics_path_lst:
            print("Processing metrics_path: " + metrics_path)
            ltab = pd.read_csv(metrics_path, sep = '\t', header = 0)
            dfs.append(ltab)

        if dfs:
            df_final = pd.concat(dfs)
            df_final.to_csv(outfile, sep = ',', index = False)
        else:
            open(outfile, 'w')


    def summarize_counts(self, sample_id, project_id, ref_acc, do_umi = False):
        confd = self.confd
        ldel = confd.ldelim
        drop_dir = self.get_sample_drop_dir(sample_id, project_id)
        bc_pre_dir_lst = self.get_bc_pre_dir_list(sample_id, project_id)
        count_path_lst = list()
        umi_count_path_lst = list()
        for bc_pre_dir in bc_pre_dir_lst:
            bc_pre_path = None
            if do_umi:
                bc_pre_path = drop_dir + ldel + bc_pre_dir + ldel + "umidir"
            else:
                bc_pre_path = drop_dir + ldel + bc_pre_dir
            count_path = self.get_count_path(bc_pre_dir, bc_pre_path, sample_id, ref_acc)
            count_path_lst.append(count_path)
            
        print("for sample_id: " + sample_id + "count_path_lst: ")
        print(count_path_lst)
        final_summary_dir = drop_dir + ldel + "metrics"
        if do_umi:
            final_summary_dir += ldel + "umidir"
        if not os.path.exists(final_summary_dir):
            os.makedirs(final_summary_dir) 
        outfile = final_summary_dir + ldel + sample_id + "_" + ref_acc + "_counts.tsv"
            
        print("Combining all the count paths to outfile: " + outfile)
        self.combine_count_paths(count_path_lst, outfile)
 
   
    def summarize_metrics(self, sample_id, project_id, ref_acc, do_umi = False):
        confd = self.confd
        ldel = confd.ldelim
        drop_dir = self.get_sample_drop_dir(sample_id, project_id)
        bc_pre_dir_lst = self.get_bc_pre_dir_list(sample_id, project_id)
        metrics_path_lst = list()
        umi_count_path_lst = list()
        for bc_pre_dir in bc_pre_dir_lst:
            bc_pre_path = None
            if do_umi:
                bc_pre_path = drop_dir + ldel + bc_pre_dir + ldel + "umidir"
            else:
                bc_pre_path = drop_dir + ldel + bc_pre_dir
            metrics_path = self.get_metrics_path(bc_pre_dir, bc_pre_path, sample_id, ref_acc)
            metrics_path_lst.append(metrics_path)
            
        print("for sample_id: " + sample_id + "metrics_path_lst: ")
        print(metrics_path_lst)
        final_summary_dir = drop_dir + ldel + "metrics"
        if do_umi:
            final_summary_dir += ldel + "umidir"
        if not os.path.exists(final_summary_dir):
            os.makedirs(final_summary_dir) 
        outfile = final_summary_dir + ldel + sample_id + "_" + ref_acc + "_metrics.txt"
            
        print("Combining all the metrics paths to outfile: " + outfile)
        self.combine_metrics_paths(metrics_path_lst, outfile)
       

    def mainFunc(self):
        sample_id = self.sample_id
        project_id = self.project_id
        ref_acc_str = self.ref_acc_str
        self.summarize_counts(sample_id, project_id, ref_acc_str)
        self.summarize_counts(sample_id, project_id, ref_acc_str, do_umi = True)
        self.summarize_metrics(sample_id, project_id, ref_acc_str)
        self.summarize_metrics(sample_id, project_id, ref_acc_str, do_umi = True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Process command for a single sample_id", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--config_log_file", dest = "config_log_file", type = str, required = True, help = "Log file containig the config vars.")
    parser.add_argument("--sample_id", dest="sample_id", type=str, required=True, help="Id of the sample to be processed.")
    parser.add_argument("--project_id", dest = "project_id", required = True, type = str, help = "project id for this sample")
    parser.add_argument("--ref_acc_str", dest = "ref_acc_str", type = str, default = "none", help = "Contains the pathogen accession ids")
    args = parser.parse_args()
    dmeto = DropMetrics(args)

    dmeto.mainFunc()

