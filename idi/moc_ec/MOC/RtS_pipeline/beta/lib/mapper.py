#!/usr/bin/env python

import yaml
import argparse
import shutil
import pandas as pd
import csv
import re
import os
import ntpath
import codecs
import getpass
import json

from collections import defaultdict

class Mapper:

    def __init__(self, confd):
        self.confd = confd

    def get_prefix_set(self, ltab):
        confd = self.confd
        read_pairing = confd.Read_pairing_val.upper()
        print("read_pairing: " + read_pairing)
        if read_pairing == 'PAIRED':
            return self.get_prefix_set_paired(ltab) 
        elif read_pairing == 'SINGLE':
            return self.get_prefix_set_single(ltab)
        else:
            print("INFO: Returning prefix_set_paired just for best_acc old config!")
            return self.get_prefix_set_paired(ltab)
        

    # Inspired by 
    # https://stackoverflow.com/questions/9816816/relative-and-absolute-paths-of-all-files
    def absoluteFilePaths(self, directory):
        for dirpath, _, filenames in os.walk(directory):
            for f in filenames:
                yield os.path.abspath(os.path.join(dirpath, f))

    def get_raw_files(self, Input_dir):
        confd = self.confd
        raw_files = list()
        for Input_dir_l in Input_dir:
            raw_files_l = self.absoluteFilePaths(Input_dir_l)
            raw_files += raw_files_l
        return raw_files


    def get_prefix_set_single(self, ltab):
        confd = self.confd
        prefix_dict = dict()
        Input_dir = confd.Input_dir
        suffix_ne = confd.suffix_ne
        print("log: starting the get_prefix_set_single.")

        use_p7 = confd.use_p7
        use_p5 = confd.use_p5

        useful_p7_set = None
        useful_p5_set = None

        useful_p7 = ''
        if use_p7:
            useful_p7_set = self.get_useful_p7(ltab)
            print("useful_p7: " +  ', '.join(useful_p7_set))

        useful_p5 = ''
        if use_p5:
            useful_p5_set = self.get_useful_p5(ltab)
            print("useful_p5: " +  ', '.join(useful_p5_set))

        raw_files = self.get_raw_files(Input_dir)

        prefix_set = set()
        for lfile in raw_files:
            filename = os.path.basename(os.path.normpath(lfile))
            lsuffix_ne = suffix_ne
            if filename.endswith("gz"):
                lsuffix_ne += ".gz"
            lprefix = ''
            lprefix_ex = ''
            if filename.endswith(lsuffix_ne):
                lprefix = filename.replace(lsuffix_ne, '')
                lprefix_ex = lfile.replace(lsuffix_ne, '')

            if lprefix and lprefix not in prefix_set:
                valid_p7 = self.is_valid_p7(lprefix, useful_p7_set)
                valid_p5 = self.is_valid_p5(lprefix, useful_p5_set)
                if valid_p7 and valid_p5:
                    prefix_set.add(lprefix)
                    prefix_dict[lprefix] = lprefix_ex

        self.dump_prefix_dict(prefix_dict)
        return prefix_set


    def get_prefix_set_paired(self, ltab):
        confd = self.confd
        Input_dir = confd.Input_dir
        suffix_s1 = confd.suffix_s1
        suffix_s2 = confd.suffix_s2
        use_p7 = confd.use_p7
        use_p5 = confd.use_p5
        prefix_dict = dict()
      
        print("log: starting the get_prefix_set_paired.")

        useful_p7_set = None
        useful_p5_set = None
        if use_p7:
            useful_p7_set = self.get_useful_p7(ltab)
            print("useful_p7: " +  ', '.join(useful_p7_set))

        useful_p5 = ''
        if use_p5:
            useful_p5_set = self.get_useful_p5(ltab)
            print("useful_p5: " +  ', '.join(useful_p5_set))

        print("Input dir from get_prefix_set_paired: ")
        print(Input_dir)
        raw_files = self.get_raw_files(Input_dir)

        prefix_set = set()
        for lfile in raw_files:
            filename = os.path.basename(os.path.normpath(lfile))
            lsuffix_s1 = suffix_s1
            lsuffix_s2 = suffix_s2
            if filename.endswith("gz"):
                lsuffix_s1 += ".gz"
                lsuffix_s2 += ".gz"
            lprefix = ''
            lprefix_ex = ''
            if filename.endswith(lsuffix_s1):
                lprefix = filename.replace(lsuffix_s1, '')
                lprefix_ex = lfile.replace(lsuffix_s1, '')
            elif filename.endswith(lsuffix_s2):
                lprefix = filename.replace(lsuffix_s2, '')
                lprefix_ex = lfile.replace(lsuffix_s2, '')

            print("lprefix: " + lprefix)
            if lprefix and lprefix not in prefix_set:
                valid_p7 = self.is_valid_p7(lprefix, useful_p7_set)
                valid_p5 = self.is_valid_p5(lprefix, useful_p5_set)
                if valid_p7 and valid_p5:
                    prefix_set.add(lprefix)
                    prefix_dict[lprefix] = lprefix_ex

        self.dump_prefix_dict(prefix_dict)
        print("Done: dump_prefix_dict")
        print(prefix_dict)
        return prefix_set

    def is_valid_p7(self, lprefix, useful_p7_set):
        """ This function is written to include those prefix that are relevant
        with respect to useful_p7. useful_p7 is calculated based on the 
        project set relevant under current run.
        """
        confd = self.confd
        use_p7 = confd.use_p7
        use_p5 = confd.use_p5

        if not use_p7:
             return True
        else:

            lprefix_p7 = None
            if use_p5:
                lprefix_p7 = self.get_prefix_p7(lprefix)
            else:
                lprefix_p7 = lprefix

            for useful_p7 in useful_p7_set:
                if lprefix_p7.endswith(useful_p7):
                    return True
            return False

    def is_valid_p5(self, lprefix, useful_p5_set):
        """ This function is written to include those prefix that are relevant
        with respect to useful_p7. useful_p7 is calculated based on the 
        project set relevant under current run.
        """
        confd = self.confd
        use_p7 = confd.use_p7
        use_p5 = confd.use_p5

        if not use_p5:
             return True
        else:

            lprefix_p5 = None
            if use_p7:
                lprefix_p5 = self.get_prefix_p5(lprefix)
            else:
                lprefix_p5 = lprefix

            for useful_p5 in useful_p5_set:

                if lprefix_p5.endswith(useful_p5):
                    return True
            return False

    def create_sample_host_reference(self, ltab):
        confd = self.confd
        sample_id = confd.sample_id
        host_reference = confd.Host_reference
        sample_id_lst = ltab[sample_id].tolist()
        host_reference_lst = ltab[host_reference].tolist()
        sample_host_referenfce = dict(zip(sample_id_lst, host_reference_lst)) 
        return sample_host_referenfce

    def create_sample_refacc(self, ltab):
        confd = self.confd
        sample_id = confd.sample_id
        ref_accession = confd.Ref_accession
 
        sample_refacc = dict()
        for index, row in ltab.iterrows():
            lsample_id = row[sample_id]
            lref_str1 = row[ref_accession]
            regex = '\s*[;|,]\s*'
            lref_parts = re.split(regex, lref_str1)
            lref_str = lref_parts[-1]
            sample_refacc[lsample_id] = lref_str
        return sample_refacc
   
    def get_useful_p7(self, ltab):
        confd = self.confd
        project_id =   confd.project_id    
        P7_index =     confd.P7_index
        lproject_set = confd.project_set
        useful_p7 = set()
        for index, row in ltab.iterrows():
           lproject_id = row[project_id]
           lp7_index = row[P7_index] 
           print(project_id + " " + P7_index + " " + lproject_id + " " + \
               lp7_index + " " + ', '.join(lproject_set))
           if lproject_id in lproject_set:
               useful_p7.add(lp7_index)
        #print("useful_p7: " + ', '.join(useful_p7))
        return useful_p7

    def get_useful_p5(self, ltab):
        confd = self.confd
        project_id =   confd.project_id    
        P5_index =     confd.P5_index
        lproject_set = confd.project_set
        useful_p5 = set()
        for index, row in ltab.iterrows():
           lproject_id = row[project_id]
           lp5_index = row[P5_index] 
           print(project_id + " " + P5_index + " " + lproject_id + " " + \
               lp5_index + " " + ', '.join(lproject_set))
           if lproject_id in lproject_set:
               useful_p5.add(lp5_index)
        return useful_p5
  
    def get_prefix_p7(self, lprefix):
        parts = lprefix.rsplit("_", 1)
        lprefix_p7 = parts[0]
        return lprefix_p7

    def get_prefix_p5(self, lprefix):
        parts = lprefix.rsplit("_", 1)
        lprefix_p5 = parts[1]
        return lprefix_p5

    def process_p7(self, lprefix, lrow):
        if None == lprefix:
            return lprefix
        confd = self.confd
        use_p7 = confd.use_p7
        use_p5 = confd.use_p5
        if not use_p7:
            return lprefix
        else:
            P7_index = confd.P7_index
            lp7_index = lrow[P7_index]

            lprefix_p7 = None
            if use_p5:
                lprefix_p7 = self.get_prefix_p7(lprefix)
            else:
                lprefix_p7 = lprefix

            if lp7_index in lprefix_p7:
                return lprefix
            else:
                return None

    def process_p5(self, lprefix, lrow):
        if None == lprefix:
            return lprefix
        confd = self.confd
        use_p7 = confd.use_p7
        use_p5 = confd.use_p5
        if not use_p5:
            return lprefix
        else:
            P5_index = confd.P5_index
            lp5_index = lrow[P5_index]

            lprefix_p5 =  None
            if use_p7:
                lprefix_p5 = self.get_prefix_p5(lprefix)
            else:
                lprefix_p5 = lprefix

            if lp5_index in lprefix_p5:
                return lprefix
            else:
                return None

    def process_lane(self, lprefix, lrow):
        regex = '\s*[,|;]\s*'
        if None == lprefix:
            return lprefix
        confd = self.confd
        use_lane = confd.use_lane
        if not use_lane:
            return lprefix
        else:
            Lane_index = confd.Lane_index
            lane_val_all = lrow[Lane_index]
            l_lane_lst = re.split(regex, lane_val_all)
            for lane_val in l_lane_lst:
                lane_val_prefix = lane_val + "_"
                if lprefix.startswith(lane_val_prefix):
                    return lprefix
            return None

    def get_prefix_intersect(self, prefix_lst, l_seq_paths):
        l_int_lst = set()
        regex = '\s*[,|;]\s*'

        for lpre in prefix_lst:
            l_seq_lst = re.split(regex, l_seq_paths)
            for lseq in l_seq_lst:
                if lpre in lseq:
                    l_int_lst.add(lpre)
         
        return l_int_lst

    def get_prefix_equal(self, prefix_lst, lsample_id):
        l_int_lst = set()

        for lpre in prefix_lst:
            #if lpre == lsample_id:
            if lsample_id in lpre:
                l_int_lst.add(lpre)
        
        len_set = len(l_int_lst)
        if len_set == 1: 
            return l_int_lst
        elif len_set == 0:
            raise StandardError("From get_prefix_equal, l_int_lst empty")
        elif len_set > 1:
            raise StandardError("From get_prefix_equal, l_int_lst contains multiple prefix")
        else:
            raise StandardError("From get_prefix_equal, len_set has negative elements")


    def get_prefix_startswith(self, prefix_lst, lsample_id):
        l_int_lst = set()

        if not prefix_lst and not lsample_id:
            raise StandardError('From get_prefix_startswith, both args empty')
        if not prefix_lst:
            raise StandardError('From get_prefix_startswith, prefix_lst empty')
        if not lsample_id:
            raise StandardError('From get_prefix_startswith, lsample_id empty')
        for lpre in prefix_lst:
            if lpre.startswith(lsample_id):
                l_int_lst.add(lpre)
        
        len_set = len(l_int_lst)
        if len_set == 1: 
            return l_int_lst
        elif len_set == 0:
            raise StandardError("From get_prefix_startswith, l_int_lst empty")
        elif len_set > 1:
            raise StandardError("From get_prefix_startswith, l_int_lst contains multiple prefix")
        else:
            raise StandardError("From get_prefix_startswith, len_set has negative elements")
                  
       
    def build_prefix_map(self, prefix_lst, ltab):
        confd = self.confd
        use_seq_path = confd.use_seq_path
        use_sample_id = confd.use_sample_id
        sample_id = confd.sample_id
        Path_to_SeqFile_tag = confd.Path_to_SeqFile
        prefix_map = dict()

        for index, row in ltab.iterrows():
            lsample_id = row[sample_id]
            if use_seq_path:
                l_seq_paths = row[Path_to_SeqFile_tag]
                print("l_seq_paths: " + l_seq_paths)
                l_pre_lst = self.get_prefix_intersect(prefix_lst, l_seq_paths)
                if l_pre_lst:
                    prefix_map[lsample_id] = l_pre_lst
            elif use_sample_id:
                l_pre_lst = self.get_prefix_startswith(prefix_lst, lsample_id)
                if l_pre_lst:
                    prefix_map[lsample_id] = l_pre_lst
            else:
                prefix_map[lsample_id] = prefix_lst    
        print("prefix lst: ")
        print(prefix_lst)
        print("from build prefix map")
        print(prefix_map)
        return prefix_map

    
    def create_sample_map_prefix(self, ltab, prefix_map):
        # Ignore those samples that does not have relevant project ids.
        confd = self.confd
        lproject_set = confd.project_set
        use_p7 = confd.use_p7
        use_p5 = confd.use_p5
        sample_id = confd.sample_id
        P7_index = confd.P7_index
        project_id = confd.project_id

        print("from create_sample_map_prefix, prefix_map: ")
        print(prefix_map)
        sample_prefix = dict()
        for index, row in ltab.iterrows():
            lproject_id = row[project_id]
            if lproject_id in lproject_set: 
                lsample_id = row[sample_id]
                print("lsample_id: " + lsample_id)
                sample_prefix[lsample_id] = set()
                lp7_index = row[P7_index]
                for lprefix in prefix_map[lsample_id]:
                    lprefix_1_p7 = self.process_p7(lprefix, row)   
                    lprefix_1_p5 = self.process_p5(lprefix_1_p7, row)   
                    lprefix_2 = self.process_lane(lprefix_1_p5, row)
                    if None != lprefix_2:
                        sample_prefix[lsample_id].add(lprefix_2)
        return sample_prefix

    def create_sample_map_bc(self, ltab, bc_delim):
        confd = self.confd
        sample_id = confd.sample_id
        lbc = confd.lbc
        project_id = confd.project_id
        lproject_set = confd.project_set

        sample_bc = dict()
        for index, row in ltab.iterrows():
            lproject_id = row[project_id]
            if lproject_id in lproject_set: 
                lsample_id = row[sample_id]
                lbc_str = row[lbc]
                lbarcodes = re.split(bc_delim, lbc_str)
                lbc_set = set(lbarcodes) 
                sample_bc[lsample_id] = lbc_set
        return sample_bc
         
    def get_sample_project(self, ltab):
        confd = self.confd
        sample_id = confd.sample_id
        project_id = confd.project_id
        sample_id_lst = ltab[sample_id].tolist()
        project_id_lst = ltab[project_id].tolist()
        sample_project = dict(zip(sample_id_lst, project_id_lst))
        return sample_project
    

    def create_projid_sample_hostref(self, ltab):
        confd = self.confd
        sample_id = confd.sample_id
        project_id = confd.project_id
        host_reference = confd.Host_reference
        projid_sample_host_ref_d = defaultdict(list)
        for index, row in ltab.iterrows():
            lsample_id = row[sample_id]
            lproject_id = row[project_id]
            host_ref_str = row[host_reference]
            host_ref_str_l = host_ref_str.lower()
            if host_ref_str_l and host_ref_str_l != "unknown":
                proj_key = lproject_id + ";" + host_ref_str
                if lsample_id not in projid_sample_host_ref_d[proj_key]:
                    projid_sample_host_ref_d[proj_key].append(lsample_id)
        return projid_sample_host_ref_d 

    def create_projid_sample_refacc(self, ltab):
        """ This would create a set of dictionary between project_id to 
        sample_id suffix indexed by sample_id.
        """ 
        confd = self.confd
        sample_id = confd.sample_id
        project_id = confd.project_id
        ref_accession = confd.Ref_accession

        projid_sample_refacc_d = defaultdict(list)
        for index, row in ltab.iterrows():
            lsample_id = row[sample_id]
            lproject_id = row[project_id]
            lref_acc_str = row[ref_accession]
            lref_acc_str_l = lref_acc_str.lower()
            if lref_acc_str_l and lref_acc_str_l != "unknown":
                # Split ref_accession
                regex = '\s*[,|;]\s*'
                lref_acc_lst = re.split(regex, lref_acc_str)
                lref_item = lref_acc_lst[-1] 
                proj_key = lproject_id + ";" + lref_item
                if lsample_id not in projid_sample_refacc_d[proj_key]:
                    projid_sample_refacc_d[proj_key].append(lsample_id)
        return projid_sample_refacc_d

    def dump_prefix_dict(self, prefix_dict):
        confd = self.confd
        prefix_dict_out = confd.prefix_dict_out
        with open(prefix_dict_out, 'w') as outfile:
            yaml.dump(prefix_dict, outfile, default_flow_style=False)
        return prefix_dict_out
