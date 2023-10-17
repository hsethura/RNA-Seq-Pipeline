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
from itertools import chain
    
from subprocess import call
from collections import defaultdict
from lib.Parse_featureCounts import get_all_gene_counts
from lib.Parse_featureCounts import get_all_metrics_counts
#from lib.Parse_featureCounts import get_all_sample_summary
from lib.patho import copy_and_update_patho
from lib.patho import patho_fna_concat    
from lib.patho import make_index_bwa
from lib.patho import copyLargeFile
from lib.metrics import Metrics
from lib.mapper import Mapper
from lib.confdict import ConfDict
from lib.splitter import Splitter
from lib.expander import Expander
from lib.garbagesorter import GarbageSorter
from lib.dropsplitter import DropSplitter
from lib.dropcounter import DropCounter
from lib.dropmetrics import DropMetrics
    
class PipelineII:   
    def __init__(self, options):

        self.options = options
        basepath = os.path.dirname(os.path.realpath(__file__))
        confd = ConfDict(options, basepath)
        confd.loadConfig()
        self.config_log_out = confd.dumpConfigLog()
        self.confd = confd

        self.mapo  = Mapper(confd)
        self.meto = Metrics(confd)
        self.dcounto = DropCounter(confd)

    def exe_bestacc(self):
        confd = self.confd
        ldel = confd.ldelim
        Key_path = confd.Key_path
        config_file = confd.config_file
        Merge_dir = confd.Merge_dir
        bestacc_script = confd.bestacc_script
        basepath = confd.basepath
        Script_dir = basepath + ldel + "shell_scripts"
        bestacc_cmd = "sh " + bestacc_script + " " + Key_path + " " + Merge_dir + " " + Merge_dir + " " + config_file + " " + Script_dir
        print("Starting bestacc: " + bestacc_cmd)
        call(bestacc_cmd.split())

    
    def build_dict(self):
        confd = self.confd
        ldelim = confd.ldelim
        UGER_cbp = confd.UGER_cbp
        ldict_infile  = confd.Dict_infile
        UGER_cbp_dir = confd.UGER_cbp_dir
        Data_dir = confd.Data_dir
        Log_dir = confd.Log_dir
        use_qsub = confd.use_qsub
        dict_builder = confd.dict_builder

        ldict_infile2 = os.path.basename(os.path.normpath(ldict_infile))
        ldict_infile3 = ldict_infile2.replace(".txt", "")
        ldict_file1 = ldict_infile3 + ".dat"
        ldict_file = Data_dir + ldelim + ldict_file1 
        dict_builder_cmd = dict_builder + " -i " + ldict_infile + " -o " + \
            ldict_file 
        ldict_build_job_path = Log_dir + ldelim + "ldict_build_job.txt"
        ldict_build_jf = open(ldict_build_job_path, "w")
        ldict_build_jf.write(dict_builder_cmd + "\n")
        ldict_build_jf.close()
        ldict_build_cmd = UGER_cbp + " --cmds_file " + ldict_build_job_path + \
                                " --batch_size 1" + \
                                " --memory 1" + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header" 
        if use_qsub:
            call(ldict_build_cmd.split())
        print(ldict_build_cmd)
        return ldict_file
   
    def create_soft_links(self, prefix_lst):
        print("Creating soft links from input directory.")
        confd = self.confd
        splito = Splitter(confd)
        for lprefix in prefix_lst:
            splito.create_soft_links(lprefix) 
    
    def create_split_job_files(self, prefix_lst, ldict_file):
        confd = self.confd
        ldelim = confd.ldelim
        UGER_cbp = confd.UGER_cbp
        UGER_cbp_dir = confd.UGER_cbp_dir
        Log_dir = confd.Log_dir
        use_qsub = confd.use_qsub

        splito = Splitter(confd, ldict_file)
        joblist_path = Log_dir + ldelim + "bc_split_joblist.txt"
        jfile = open(joblist_path, "w")
        for lprefix in prefix_lst:
            cmd_str2 = splito.get_split_command(lprefix)
            print(cmd_str2)
            jfile.write(cmd_str2)

        jfile.close()
        joblist_cmd = UGER_cbp + " --cmds_file " + joblist_path + \
                                " --batch_size 1" + \
                                " --memory 8" + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header"
        if use_qsub:
            call(joblist_cmd.split())
 
        return joblist_path


    def create_drop_metrics_job_files(self, sample_prefix):

        confd = self.confd
        ldel = confd.ldelim
        UGER_cbp = confd.UGER_cbp
        UGER_cbp_dir = confd.UGER_cbp_dir
        Log_dir = confd.Log_dir
        use_qsub = confd.use_qsub
        config_log_out = self.config_log_out
        do_patho = confd.do_patho
        sample_project = self.sample_project
        sample_refacc = self.sample_refacc

        DropMetrics = confd.DropMetrics

        joblist_path = Log_dir + ldel + "dropmetrics_joblist.txt" 
        jfile = open(joblist_path, "w")
        for lsample in sample_prefix:
            lproject_id = sample_project[lsample]
            ref_acc_str = sample_refacc[lsample]
            print("sample_id: " + lsample)
            cmd_str = DropMetrics + " --config_log_file " + \
                config_log_out + " --sample_id " + lsample + \
                " --project_id " + lproject_id 
            if do_patho:
                ref_acc_str_l = None
                if ref_acc_str:
                    ref_acc_str_l = ref_acc_str
                else:
                    ref_acc_str_l = "none"
                cmd_str += " --ref_acc_str " + ref_acc_str_l
            loutfile = Log_dir + ldel + lsample + "_out.txt"
            lerrfile = Log_dir + ldel + lsample + "_err.txt"
            cmd_str2 = cmd_str + " 1> " + loutfile + " 2> " + lerrfile + "\n"
            print(cmd_str2)
            jfile.write(cmd_str2)

        jfile.close()
        joblist_cmd = UGER_cbp + " --cmds_file " + joblist_path + \
                                " --batch_size 1" + \
                                " --memory 8" + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header"
        if use_qsub:
            call(joblist_cmd.split())

        return joblist_path


            

    def create_drop_count_job_files(self, sample_prefix):
        confd = self.confd
        dcounto = self.dcounto
        sample_project = self.sample_project
        UniDropCore = confd.UniDropCore
        config_log_out = self.config_log_out
        ldel = confd.ldelim
        Log_dir = confd.Log_dir
        do_patho = confd.do_patho
        sample_refacc = self.sample_refacc
        UGER_cbp = confd.UGER_cbp
        UGER_cbp_dir = confd.UGER_cbp_dir
        use_qsub = confd.use_qsub

        unidropcore_job_path = Log_dir + ldel + "unidropcore_joblist.txt"
        jfile = open(unidropcore_job_path, "w")

        lcount = 0
        for lsample in sample_prefix:
            lproject_id = sample_project[lsample]
            ref_acc_str = sample_refacc[lsample]
            bc_pre_dir_list = dcounto.get_bc_pre_dir_list(lsample, lproject_id)
            for bc_pre_str in bc_pre_dir_list:
                bc_pre_dir = dcounto.get_bc_pre_dir(lsample, lproject_id, bc_pre_str)
                bc_pre_logdir = bc_pre_dir + ldel + "logdir"
                if not os.path.exists(bc_pre_logdir):
                    os.makedirs(bc_pre_logdir)
                lcount += 1
                cmd_str = UniDropCore + " --config_log_file " + \
                    config_log_out + " --sample_id " + lsample + \
                    " --bc_pre_str " + bc_pre_str + " --project_id " + \
                    lproject_id
                if do_patho:
                    ref_acc_str_l = None
                    if ref_acc_str:
                        ref_acc_str_l = ref_acc_str
                    else:
                        ref_acc_str_l = "none" 
                    cmd_str += " --ref_acc_str " + ref_acc_str_l 
                log_tag = lsample + "_" + bc_pre_str
                loutfile = Log_dir + ldel + log_tag + "_out.txt"
                lerrfile = Log_dir + ldel + log_tag + "_err.txt"
                cmd_str2 = cmd_str + " 1> " + loutfile + " 2> " + lerrfile + "\n"
                print(cmd_str2)
                jfile.write(cmd_str2)
        jfile.close()

        num_cores = confd.patho_thread_count
        lmemory = confd.patho_memory
        run_tim_h = 1
        joblist_cmd = UGER_cbp + " --cmds_file " + unidropcore_job_path + \
                                " --batch_size 1" + \
                                " --num_cores " + str(num_cores) + \
                                " --memory " + str(lmemory) + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --run_time " + str(run_tim_h) +  \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header"
        print("unidropcore_job: " + joblist_cmd)
        if use_qsub:
            call(joblist_cmd.split())
        return unidropcore_job_path

            
    def create_drop_split_job_files(self, sample_prefix):
        confd = self.confd
        ldelim = confd.ldelim
        UGER_cbp = confd.UGER_cbp
        UGER_cbp_dir = confd.UGER_cbp_dir
        Log_dir = confd.Log_dir
        use_qsub = confd.use_qsub

        sample_project = self.sample_project
        sample_refacc = self.sample_refacc 

        dsplito = DropSplitter(confd)
        joblist_path = Log_dir + ldelim + "drop_split_joblist.txt"
        jfile = open(joblist_path, "w")
        for lsample in sample_prefix:
            project_id = sample_project[lsample]
            ref_acc_str = sample_refacc[lsample]
            cmd_str = dsplito.get_split_command(lsample, project_id, ref_acc_str)
            loutfile = Log_dir + ldelim + lsample + "_out.txt"
            lerrfile = Log_dir + ldelim + lsample + "_err.txt"
            cmd_str2 = cmd_str + " 1> " + loutfile + " 2> " + lerrfile + "\n"
            print(cmd_str2)
            jfile.write(cmd_str2)

        jfile.close()
        joblist_cmd = UGER_cbp + " --cmds_file " + joblist_path + \
                                " --batch_size 1" + \
                                " --memory 8" + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header"
        if use_qsub:
            call(joblist_cmd.split())
 
        return joblist_path

    def exe_merge(self, ltab, sample_prefix, sample_bc):
        
        confd = self.confd
        mapo = self.mapo
        do_host = confd.do_host
 
        config_log_out = self.config_log_out
        UniMerger = confd.UniMerger
        Log_dir = confd.Log_dir
        ldelim = confd.ldelim
        project_set = confd.project_set
        use_qsub = confd.use_qsub
        UGER_cbp_dir = confd.UGER_cbp_dir
        UGER_cbp = confd.UGER_cbp

        sample_project = self.sample_project

        unicore_job_path = Log_dir + ldelim + "unimerger_joblist.txt"
        jfile = open(unicore_job_path, "w")

        for lsample in sample_prefix:
            lproject_id = sample_project[lsample]
            prefix_set = sample_prefix[lsample]
            prefix_str = ';'.join(prefix_set)
            bc_set = sample_bc[lsample]
            bc_str = ';'.join(bc_set)
            if not lproject_id in project_set:
                continue
            
            cmd_str = UniMerger + " --config_log_file " + config_log_out + \
              " --sample_id \"" + lsample + "\" --project_id " + lproject_id + \
              ' --prefix_set "' + prefix_str + '" --bc_set "' +  bc_str + '"'

            loutfile = Log_dir + ldelim + lsample + "_out.txt"
            lerrfile = Log_dir + ldelim + lsample + "_err.txt"
            cmd_str2 = cmd_str + " 1> " + loutfile + " 2> " + lerrfile + "\n"
            print(cmd_str2)
            jfile.write(cmd_str2)
        jfile.close()
        lmemory = "8"
        if do_host:
            lmemory = "16"
        joblist_cmd = UGER_cbp + " --cmds_file " + unicore_job_path + \
                                " --batch_size 1" + \
                                " --num_cores 1" + \
                                " --memory " + lmemory + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header"
        if use_qsub:
            call(joblist_cmd.split())
        return unicore_job_path

  
    def create_temp_bamdir(self, lproject_id, use_umi = False, subdir = ''):
        confd = self.confd
        ldelim = confd.ldelim
        Patho_dir = confd.Patho_dir
        patho_project_id = Patho_dir + ldelim + lproject_id
        if use_umi:
            patho_project_id += ldelim + "umidir"
        if subdir:
            patho_project_id += ldelim + subdir
        patho_temp_bamdir = patho_project_id + ldelim + "temp_bamdir"
        if not os.path.exists(patho_temp_bamdir):
            print("Creating temp_bam directory:" + patho_temp_bamdir)
            os.makedirs(patho_temp_bamdir)


    def create_picard_dir(self, lproject_id, species_type, use_umi = False):
        confd = self.confd
        ldel = confd.ldelim
        picard_outdir = None
        Results_path = confd.Results_path
        picard_project = Results_path + ldel + lproject_id
        if use_umi:
            picard_project += ldel + "umidir"
        picard_prefix = picard_project + ldel + "bam_metrics"
        if "patho" == species_type:
            picard_outdir = picard_prefix + "_patho"
        elif "host" == species_type:
            picard_outdir = picard_prefix + "_host"
        else:
            raise StandardError("Unknown species type: " + species_type)
            
        if not os.path.exists(picard_outdir):
            print("Creating picard metrics dir: " + picard_outdir)
            os.makedirs(picard_outdir)
        
    def submit_ds_bcdistcore(self, ltab, sample_prefix):
        confd = self.confd
        mapo = self.mapo
        ldelim = confd.ldelim
        do_patho = confd.do_patho
        do_host = confd.do_host
        config_log_out = self.config_log_out
        BCDistCore_p = confd.BCDistCore_path 
        use_qsub = confd.use_qsub
        UGER_cbp_dir = confd.UGER_cbp_dir
        UGER_cbp = confd.UGER_cbp
        Log_dir = confd.Log_dir
        project_set = confd.project_set

        sample_project = self.sample_project
        sample_refacc = self.sample_refacc
        sample_host_reference = self.sample_host_reference

        ds_bcdist_job_path = Log_dir + ldelim + "ds_bcdist_joblist.txt"
        jfile = open(ds_bcdist_job_path, "w")

        for lsample in sample_prefix:
            lproject_id = sample_project[lsample]
            ref_acc_str = sample_refacc[lsample]
            if not lproject_id in project_set:
                continue
            cmd_str = BCDistCore_p + " --config_log_file " + config_log_out + \
              " --sample_id " + lsample + " --project_id " + lproject_id
    
            if do_patho:
                ref_acc_str_l = None
                if ref_acc_str:
                    ref_acc_str_l = ref_acc_str
                else:
                    ref_acc_str_l = "none" 
                cmd_str += " --ref_acc_str " + ref_acc_str_l 

            if do_host:
                host_ref_str_l = None
                if host_ref_str:
                    host_ref_str_l = host_ref_str
                else:
                    host_ref_str_l = "none"
                cmd_str += " --host_ref_str " + host_ref_str_l

            loutfile = Log_dir + ldelim + lsample + "_out.txt"
            lerrfile = Log_dir + ldelim + lsample + "_err.txt"
            cmd_str2 = cmd_str + " 1> " + loutfile + " 2> " + lerrfile + "\n"
            print(cmd_str2)
            jfile.write(cmd_str2)
        jfile.close()

        num_cores = confd.patho_thread_count
        lmemory = confd.patho_memory
        joblist_cmd = UGER_cbp + " --cmds_file " + ds_bcdist_job_path + \
                                " --batch_size 1" + \
                                " --num_cores " + str(num_cores) + \
                                " --memory " + str(lmemory) + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header"
        if use_qsub:
            call(joblist_cmd.split())
        return ds_bcdist_job_path
              
    def submit_unicores(self, ltab, sample_prefix, sample_bc):
        
        confd = self.confd
        mapo = self.mapo
 
        config_log_out = self.config_log_out
        UniCore = confd.UniCore
        Log_dir = confd.Log_dir
        ldelim = confd.ldelim
        project_set = confd.project_set
        do_patho = confd.do_patho
        do_host = confd.do_host
        use_qsub = confd.use_qsub
        UGER_cbp_dir = confd.UGER_cbp_dir
        UGER_cbp = confd.UGER_cbp
        do_umi_count = confd.do_umi_count
        use_dropseq = confd.use_dropseq 

        sample_project = self.sample_project
        sample_refacc = self.sample_refacc
        sample_host_reference = self.sample_host_reference

        unicore_job_path = Log_dir + ldelim + "unicore_joblist.txt"
        jfile = open(unicore_job_path, "w")

        for lsample in sample_prefix:
            lproject_id = sample_project[lsample]
            #prefix_set = sample_prefix[lsample]
            #prefix_str = ';'.join(prefix_set)
            #bc_set = sample_bc[lsample]
            #bc_str = ';'.join(bc_set)
            ref_acc_str = sample_refacc[lsample]
            host_ref_str = sample_host_reference[lsample]
            if not lproject_id in project_set:
                continue
                # Create temp_bam file
            if do_patho:
                self.create_temp_bamdir(lproject_id)
                if use_dropseq:
                    self.create_temp_bamdir(lproject_id, subdir = 'real')
                    self.create_temp_bamdir(lproject_id, subdir = 'garbage')
                species_type = "patho"
                self.create_picard_dir(lproject_id, species_type)
                if do_umi_count:
                    self.create_temp_bamdir(lproject_id, use_umi = True)
                    self.create_picard_dir(lproject_id, species_type, use_umi = True)
            
            cmd_str = UniCore + " --config_log_file " + config_log_out + \
              " --sample_id " + lsample + " --project_id " + lproject_id
              #'\" --prefix_set "' + prefix_str + '" --bc_set "' +  bc_str + '"'

  
            if do_patho:
                ref_acc_str_l = None
                if ref_acc_str:
                    ref_acc_str_l = ref_acc_str
                else:
                    ref_acc_str_l = "none" 
                cmd_str += " --ref_acc_str " + ref_acc_str_l 
            if do_host:
                host_ref_str_l = None
                if host_ref_str:
                    host_ref_str_l = host_ref_str
                else:
                    host_ref_str_l = "none"
                cmd_str += " --host_ref_str " + host_ref_str_l

            loutfile = Log_dir + ldelim + lsample + "_out.txt"
            lerrfile = Log_dir + ldelim + lsample + "_err.txt"
            cmd_str2 = cmd_str + " 1> " + loutfile + " 2> " + lerrfile + "\n"
            print(cmd_str2)
            jfile.write(cmd_str2)
        jfile.close()
         
        num_cores = confd.patho_thread_count
        lmemory = confd.patho_memory
        min_resource = confd.min_resource
        ucore_time = int(confd.ucore_time)
        # When uses multiple cores total allocated memory is num_cores * lmemory
        if do_patho and do_umi_count:
            lmemory = 16
        if do_host:
            if not min_resource:
                num_cores = confd.host_thread_count
                lmemory = confd.host_memory
            else:
                num_cores = 1
                lmemory = 32
        joblist_cmd = UGER_cbp + " --cmds_file " + unicore_job_path + \
                                " --batch_size 1" + \
                                " --num_cores " + str(num_cores) + \
                                " --memory " + str(lmemory) + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header"
        if ucore_time > 0:
            joblist_cmd += " --run_time " + str(ucore_time)
        if do_host and not min_resource:
            joblist_cmd += " --reservations"
        
        if use_qsub:
            call(joblist_cmd.split())
        return unicore_job_path
 
          
    def prepare_patho_refs(self, ltab):
        
        """ This function would traverse over all the rows for the list of 
        pathogens in each row. If there are more than one NC ids in a row,
        the function would combine it in a single entity, and would create 
        an index out of that. That index would be used for subsequent steps.
        """
        confd = self.confd
        inpath = confd.patho_dbpath
        outpath = confd.Data_dir
        do_replace_refname = confd.do_replace_refname 
        sample_id = confd.sample_id
        bwa_path = confd.bwa_path
        lproject_set = confd.project_set 
        project_id = confd.project_id
        ref_accession = confd.Ref_accession
        gff_parser = confd.gff_parser
        add5 = confd.add5
        add3 = confd.add3

        ref_delim = "\s*[,|;]\s*";
        processed_acc_set = set()
        processed_compo_set = set()
        print(';'.join(lproject_set))

        for index, row in ltab.iterrows():
            lsample_id = row[sample_id]
            lproject_id = row[project_id]
            
            if lproject_id in lproject_set:
                ref_acc_str = row[ref_accession]
                # Split the ref accession string
                ref_acc_lst = re.split(ref_delim, ref_acc_str)
                ref_acc = ref_acc_lst[-1]
                if not ref_acc:
                   print('ref_acc is empty. sample_id: ' + lsample_id)
                elif "unknown" == ref_acc.lower():
                    print('ref_acc contains unknown: ' + lsample_id)
                else:
                    if ref_acc not in processed_acc_set:
                        print("Ref acc: " + ref_acc)
                        processed_acc_set.add(ref_acc)
                        copy_and_update_patho(ref_acc, gff_parser, add5, add3, \
                            inpath, outpath, do_replace_refname)
                        make_index_bwa(ref_acc, outpath, bwa_path)
                        print("Making bwa index: " + ref_acc + " " + outpath)
                    else:
                        print("Refacc {} with sample-id {} already \
                            processed".format(ref_acc, lsample_id))



    def exe_garbage_sorting(self, sample_prefix):
        # This sorts the garbage at least for the bDropSeq data
        confd = self.confd
        do_gs = confd.do_gs
        drop_type = confd.drop_type
        prefix_lst = set(chain.from_iterable(sample_prefix.values()))
        print("exe_garbage_sorting, prefix_lst: " + ', '.join(prefix_lst))

        joblist_path = self.create_gs_jobs(prefix_lst)

    def create_gs_jobs(self, prefix_lst):
        confd = self.confd
        gso = GarbageSorter(confd)
        use_qsub = confd.use_qsub
        Log_dir = confd.Log_dir
        ldelim = confd.ldelim
        UGER_cbp = confd.UGER_cbp
        UGER_cbp_dir = confd.UGER_cbp_dir

        joblist_path = Log_dir + ldelim + "garbage_sorter_joblist.txt"
        jfile = open(joblist_path, "w")
        for lprefix in prefix_lst:
            cmd_str = gso.get_garbage_sorter_cmd(lprefix)
            print(cmd_str)
            jfile.write(cmd_str)

        jfile.close()
        joblist_cmd = UGER_cbp + " --cmds_file " + joblist_path + \
                                " --batch_size 1" + \
                                " --memory 8" + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header"
        if use_qsub:
            call(joblist_cmd.split())

        return joblist_path

        

    def exe_splitting(self, KeyTblFinal, prefix_lst):
        confd = self.confd
        meto = self.meto
        do_split = confd.do_split
        use_qsub = confd.use_qsub
        soft_link_from_input = confd.soft_link_from_input

        if do_split:
            if soft_link_from_input:
                self.create_soft_links(prefix_lst)
            else:
                ldict_file = self.build_dict()
                joblist_path = self.create_split_job_files(prefix_lst, ldict_file)
                if use_qsub: 
                    meto.gen_all_bcLog()


    def exe_drop_count(self, KeyTblFinal, sample_prefix):
        confd = self.confd
        
        do_drop_count = confd.do_drop_count
        if do_drop_count:
            joblist_path = self.create_drop_count_job_files(sample_prefix)
    
    def exe_drop_splitting(self, KeyTblFinal, sample_prefix):
        confd = self.confd
        do_drop_split = confd.do_drop_split


        if do_drop_split:
            joblist_path = self.create_drop_split_job_files(sample_prefix)

    def exe_drop_metrics(self, KeyTblFinal, sample_prefix):
        confd = self.confd
        
        do_drop_metrics = confd.do_drop_metrics
        if do_drop_metrics:
            joblist_path = self.create_drop_metrics_job_files(sample_prefix)
                       

    def exe_alignment(self, KeyTblFinal, sample_prefix, sample_bc):
        confd = self.confd
        do_align = confd.do_align
        do_patho = confd.do_patho
        do_count = confd.do_count
        do_umi_count = confd.do_umi_count
        do_ref = confd.do_ref

        if do_align or do_count or do_umi_count:
            if do_align and do_patho and do_ref:
                self.prepare_patho_refs(KeyTblFinal)
            rtsjob_path = self.submit_unicores(KeyTblFinal, sample_prefix, sample_bc)

    def gen_host_metrics_umi(self, ltab):
        confd = self.confd
        mapo = self.mapo
        meto = self.meto
        Summary_dir = confd.Summary_dir
        Host_dir = confd.Host_dir
        project_set = confd.project_set
        ldelim = confd.ldelim
        print("----------------------------")
        print("Starting generation of host metrics")
        projid_sample_hostref_d = mapo.create_projid_sample_hostref(ltab)
        print(projid_sample_hostref_d)
        for project_hostref in projid_sample_hostref_d:
            lsample_lst = projid_sample_hostref_d[project_hostref]
            print(lsample_lst)
            print(project_set)
            regex = "\s*;\s*"
            project_hostref_parts = re.split(regex, project_hostref)
            lproject_id = project_hostref_parts[0]
            hostref = project_hostref_parts[-1]
            if lproject_id in project_set:
                lhost_summary_dir = Summary_dir + ldelim + lproject_id
                project_dir = Host_dir + ldelim + lproject_id
                print(hostref)
                meto.copy_fastq_etc_host(lproject_id, hostref, use_umi = True)
                meto.write_read_count_table(lproject_id, hostref, lsample_lst, \
                    species_type = "host", has_header = True, use_umi = True)
                meto.write_corr_file(lproject_id, hostref, species_type = "host", use_umi = True) 
                meto.write_host_metrics_table(lproject_id, hostref, lsample_lst, use_umi = True)
                view_mode = "gv"
                meto.write_bam_file_paths(lproject_id, hostref, lsample_lst, view_mode, append_mode = False, species_type = "host", use_umi = True)
                view_mode = "abs"
                meto.write_bam_file_paths(lproject_id, hostref, lsample_lst, view_mode, append_mode = False, species_type = "host", use_umi = True)
                meto.exe_picard_parser(lproject_id, "host", use_umi = True)
                meto.combine_umi_coll_pdfs(lproject_id, hostref, lsample_lst, species_type = "host")


    def gen_host_metrics(self, ltab, lsubdir = ''):
        confd = self.confd
        mapo = self.mapo
        meto = self.meto
        Summary_dir = confd.Summary_dir
        Host_dir = confd.Host_dir
        project_set = confd.project_set
        ldelim = confd.ldelim

        print("----------------------------")
        if not lsubdir:
            print("Starting generation of host metrics")
        elif lsubdir:
            print("Starting generation of host metrics with subdir: " + lsubdir)

        projid_sample_hostref_d = mapo.create_projid_sample_hostref(ltab)
        print(projid_sample_hostref_d)
        for project_hostref in projid_sample_hostref_d:
            lsample_lst = projid_sample_hostref_d[project_hostref]
            print(lsample_lst)
            print(project_set)
            regex = "\s*;\s*"
            project_hostref_parts = re.split(regex, project_hostref)
            lproject_id = project_hostref_parts[0]
            hostref = project_hostref_parts[-1]
            if lproject_id in project_set:
                lhost_summary_dir = Summary_dir + ldelim + lproject_id
                project_dir = Host_dir + ldelim + lproject_id
                print(hostref)
                meto.copy_fastq_etc_host(lproject_id, hostref, subdir = lsubdir)
                meto.write_read_count_table(lproject_id, hostref, lsample_lst, \
                    species_type = "host", has_header = True, subdir = lsubdir)
                meto.write_corr_file(lproject_id, hostref, species_type = "host", subdir = lsubdir) 
                meto.write_host_metrics_table(lproject_id, hostref, lsample_lst, subdir = lsubdir)
                view_mode = "gv"
                meto.write_bam_file_paths(lproject_id, hostref, lsample_lst, view_mode, append_mode = False, species_type = "host", subdir = lsubdir)
                view_mode = "abs"
                meto.write_bam_file_paths(lproject_id, hostref, lsample_lst, view_mode, append_mode = False, species_type = "host", subdir = lsubdir)
                meto.exe_picard_parser(lproject_id, "host", subdir = lsubdir)


    def exe_rpg_metrics(self, lproject_id, lsample_lst, ref_acc, use_umi = False, subdir = ''):
        confd = self.confd
        Summary_dir = confd.Summary_dir
        lc_lower = confd.LC_method_val.lower()
        read_counter = confd.read_counter
        if lc_lower == 'allseq':
            if read_counter == 'JL_counter':
                meto = self.meto
                meto.ref_RPG_metrics(lproject_id, use_umi)
                #meto.data_finish(lproject_id)
            elif read_counter == 'featureCounts':
                ldelim = confd.ldelim
                Patho_dir = confd.Patho_dir
                project_dir = Patho_dir + ldelim + lproject_id
                lpatho_summary_dir = Summary_dir + ldelim + lproject_id
                summary_tag = lproject_id + "_" + ref_acc + "_metrics.tsv"
                lpatho_summary_file = lpatho_summary_dir + ldelim + summary_tag
                get_all_metrics_counts(project_dir, lpatho_summary_file, 
                    lsample_lst, ref_acc = ref_acc, has_header = False)
                self.copy_summary_to_result(summary_tag, lproject_id)
        elif lc_lower == 'rts' or lc_lower == 'rts-ts' \
                or lc_lower == 'smarter' or lc_lower == 'alignr2' or lc_lower == 'alignr1' or lc_lower == 'scr':
            meto = self.meto
            meto.ref_RPG_metrics(lproject_id, subdir = subdir)
            #meto.data_finish(lproject_id)
        elif lc_lower == 'bdropseq' or lc_lower == 'indrops':
            meto = self.meto
            meto.ref_RPG_metrics(lproject_id, subdir = subdir)
        else:
            raise ValueError('Unknow LC_method_val: ' + lc_lower)

    def copy_summary_to_result(self, summary_tag, lproject_id):
        confd = self.confd
        ldelim = confd.ldelim
        Results_path = confd.Results_path
        Summary_dir = confd.Summary_dir
        lpatho_summary_dir = Summary_dir + ldelim + lproject_id
        results_dir = Results_path + ldelim + lproject_id
        src_file = lpatho_summary_dir + ldelim + summary_tag
        dest_file = results_dir + ldelim + summary_tag
        print("Copying metrics file")
        print("src_file: " + src_file)
        print("dest_file: " + dest_file)
        copyLargeFile(src_file, dest_file) 
               
    def gen_patho_metrics_umi(self, ltab):
        mapo = self.mapo
        meto = self.meto
        confd = self.confd
        project_set = confd.project_set

        write_path_set = set()
        projid_sample_refacc_d = mapo.create_projid_sample_refacc(ltab)
        for project_refacc in projid_sample_refacc_d:
            lsample_lst = projid_sample_refacc_d[project_refacc]
            regex = "\s*;\s*"
            project_refacc_parts = re.split(regex, project_refacc)
            lproject_id = project_refacc_parts[0]
            ref_acc = project_refacc_parts[-1] 
            if lproject_id in confd.project_set:
                meto.write_read_count_table(lproject_id, ref_acc, lsample_lst, \
                    species_type = "patho", has_header = False, use_umi = True)
                meto.write_fpkm_file(lproject_id, ref_acc, use_umi = True)
                meto.write_corr_file(lproject_id, ref_acc, species_type = "patho", use_umi = True) 
                meto.exe_picard_parser(lproject_id, "patho", use_umi = True)
                meto.copy_fastq_gff(lproject_id, ref_acc, use_umi = True)
                meto.combine_umi_coll_pdfs(lproject_id, ref_acc, lsample_lst, species_type = "patho")
                append_mode = '' 
                if lproject_id not in write_path_set:
                    append_mode = False
                    write_path_set.add(lproject_id)
                    self.exe_rpg_metrics(lproject_id, lsample_lst, ref_acc, use_umi = True)
                else:
                   append_mode = True
                view_mode = "gv"
                meto.write_bam_file_paths(lproject_id, ref_acc, lsample_lst, view_mode, append_mode, species_type = "patho", use_umi = True)
                view_mode = "abs"
                meto.write_bam_file_paths(lproject_id, ref_acc, lsample_lst, view_mode, append_mode, species_type = "patho", use_umi = True)


    def append_umi_paths(self, gv_path_file, combined_gv_path_file):
        with open(gv_path_file) as inf, open(combined_gv_path_file, "w") as outf:
            for inline in inf:
                inline_1 = inline.rstrip()
                parts = inline_1.split()
                filepath = parts[1]
                dirname = ntpath.dirname(filepath)
                filename = ntpath.basename(filepath)
                outf.write(inline)
                if filename.endswith("bam") or filename.endswith("tdf"):
                    umi_path = dirname + "/umidir/" + filename
                    outline = "EXTRA\t" + umi_path
                    outf.write(outline + "\n")


    def copy_moc_files(self):
        confd = self.confd
        Results_path = confd.Results_path
        ldelim = confd.ldelim
        # At this moment we shall copy a file to the results folder, this 
        # would be a hardcoded path.
        moc_dir = "/broad/IDP-Dx_storage/MOC/files"
        moc_file = "MOC_RNA-Seq_data_UserGuide_CURRENT.docx"
 
        moc_file_path = moc_dir + ldelim + moc_file
        result_dir = Results_path + ldelim + lproject_id
        moc_file_out_path = result_dir + ldelim + moc_file
        print("Copying: " + moc_file_path + " to " + moc_file_out_path)
        copyLargeFile(moc_file_path, moc_file_out_path)
         


    def write_patho_gv_session_info(self):
        confd = self.confd
        project_set = confd.project_set 
        ldelim = confd.ldelim
        Results_path = confd.Results_path
        Bam_path = confd.Bam_path
        do_umi_metrics = confd.do_umi_metrics

        for lproject_id in project_set:
            result_dir = Results_path + ldelim + lproject_id
            bam_path_file = result_dir + ldelim + lproject_id + "_patho_bam_paths_gv.txt"
            bam_dir = Bam_path + ldelim + lproject_id
            bam_path_file_out = bam_dir + ldelim + lproject_id + "_patho_bam_paths_gv.txt"
            print("Copying " + bam_path_file + " to " + bam_path_file_out)
            copyLargeFile(bam_path_file, bam_path_file_out)
            if do_umi_metrics:
                bam_path_file_combined = result_dir + ldelim + lproject_id + "_patho_bam_paths_gv_combined.txt"
                self.append_umi_paths(bam_path_file, bam_path_file_combined)
                bam_path_file_out_combined = bam_dir + ldelim + lproject_id + "_patho_bam_paths_gv_combined.txt"
                print("Copying " + bam_path_file_combined + " to " + bam_path_file_out_combined)
                copyLargeFile(bam_path_file_combined, bam_path_file_out_combined)
            
            # Create a file and write session info into that.
            session_file = result_dir + ldelim + lproject_id + "_gv_session_info.txt"
            session_f = open(session_file, "w")
            session_f.write("###gv_session_file_location " + "\n")
  
            mount_point = "/broad/IDP-Dx_storage/idpweb"
            idpweb_address = "http://idpweb.broadinstitute.org"
            bam_dir_http = bam_dir.replace(mount_point, idpweb_address)
            session_file_http = bam_dir_http + ldelim + lproject_id + "_patho_bam_paths_gv.txt"
            
            session_f.write("session_file_abs: " + bam_path_file_out + "\n")
            session_f.write("session_file_http: " + session_file_http + "\n")
            session_f.write("###gv_command\n")

            if do_umi_metrics:
                session_combined_http = bam_dir_http + ldelim + lproject_id + "_patho_bam_paths_gv_combined.txt"
                session_f.write("session_combined_http: " + session_combined_http + "\n\n")
                session_f.write("http://genomeview.org/start/launch.jnlp?--session "+ session_combined_http + "\n")
            else:
                session_f.write("http://genomeview.org/start/launch.jnlp?--session "+ session_file_http + "\n")
            
            session_f.close()
            session_file_out = bam_dir + ldelim + lproject_id + "_gv_session_info.txt"
            print("Copying " + session_file + " to " + session_file_out)
            copyLargeFile(session_file, session_file_out)


    def gen_patho_metrics(self, ltab, lsubdir = ''):
        mapo = self.mapo
        meto = self.meto
        confd = self.confd
        project_set = confd.project_set

        #self.copy_moc_files()

        write_path_set = set()
        projid_sample_refacc_d = mapo.create_projid_sample_refacc(ltab)
        for project_refacc in projid_sample_refacc_d:
            lsample_lst = projid_sample_refacc_d[project_refacc]
            regex = "\s*;\s*"
            project_refacc_parts = re.split(regex, project_refacc)
            lproject_id = project_refacc_parts[0]
            ref_acc = project_refacc_parts[-1] 
            if lproject_id in confd.project_set:
                meto.write_read_count_table(lproject_id, ref_acc, lsample_lst, \
                    species_type = "patho", has_header = False, subdir = lsubdir)
                meto.write_fpkm_file(lproject_id, ref_acc, subdir = lsubdir)
                meto.write_corr_file(lproject_id, ref_acc, species_type = "patho", subdir = lsubdir) 
                meto.exe_picard_parser(lproject_id, "patho", subdir = lsubdir)
                meto.copy_fastq_gff(lproject_id, ref_acc, subdir = lsubdir)
                append_mode = '' 
                if lproject_id not in write_path_set:
                    append_mode = False
                    write_path_set.add(lproject_id)
                    self.exe_rpg_metrics(lproject_id, lsample_lst, ref_acc, subdir = lsubdir)
                else:
                    append_mode = True
                view_mode = "gv"
                meto.write_bam_file_paths(lproject_id, ref_acc, lsample_lst, view_mode, append_mode, species_type = "patho", subdir = lsubdir)
                view_mode = "abs"
                meto.write_bam_file_paths(lproject_id, ref_acc, lsample_lst, view_mode, append_mode, species_type = "patho", subdir = lsubdir)

    def getKeyTblFinal(self):
        confd = self.confd
        KeyTblFinal_ori = confd.getKeyTblFinal()
        print("KeyTblFinal is loaded!")
        do_expand = confd.do_expand
        if not do_expand:
            return KeyTblFinal_ori
        else: 
            # Expand the table if required
            expo = Expander(confd)
            KeyTblFinal = expo.expandKeyTbl(KeyTblFinal_ori)
            return KeyTblFinal
       

    def mainFunc(self):
        confd = self.confd
        mapo = self.mapo

        do_split = confd.do_split
        do_merge = confd.do_merge
        do_patho = confd.do_patho
        do_host = confd.do_host
        do_align = confd.do_align
        do_count = confd.do_count
        do_metrics = confd.do_metrics
        do_umi_count = confd.do_umi_count
        do_umi_metrics = confd.do_umi_metrics
        do_bestacc = confd.do_bestacc
        do_gs = confd.do_gs
        use_dropseq = confd.use_dropseq
        do_bc_dist = confd.do_bc_dist
        do_drop_split = confd.do_drop_split
        do_drop_count = confd.do_drop_count
        do_drop_metrics = confd.do_drop_metrics
        rm_rts_dup = confd.rm_rts_dup

        KeyTblFinal = self.getKeyTblFinal()
        prefix_lst = mapo.get_prefix_set(KeyTblFinal)
        prefix_map = mapo.build_prefix_map(prefix_lst, KeyTblFinal)
        print(prefix_lst)
        #print("prefix_lst: " + ', '.join(prefix_lst))

        mapo = self.mapo
        sample_prefix = mapo.create_sample_map_prefix(KeyTblFinal, prefix_map)
        print("sample_prefix:")
        print(sample_prefix)
        bc_delim = '\s*;\s*'
        sample_bc = mapo.create_sample_map_bc(KeyTblFinal, bc_delim)

        sample_project = mapo.get_sample_project(KeyTblFinal)
        sample_refacc = mapo.create_sample_refacc(KeyTblFinal) 
        sample_host_reference = mapo.create_sample_host_reference(KeyTblFinal) 
        
        self.sample_project = sample_project
        self.sample_refacc = sample_refacc
        self.sample_host_reference = sample_host_reference
        
        print(sample_prefix)
        if do_gs and use_dropseq:
            self.exe_garbage_sorting(sample_prefix)
        if do_split and not use_dropseq:
            self.exe_splitting(KeyTblFinal, prefix_lst)
        if do_merge:
            self.exe_merge(KeyTblFinal, sample_prefix, sample_bc)
        if do_bestacc:
            self.exe_bestacc()
        if do_align or do_count or do_umi_count:
            self.exe_alignment(KeyTblFinal, sample_prefix, sample_bc)
        if do_patho and do_metrics:
            if use_dropseq:
                self.gen_patho_metrics(KeyTblFinal, lsubdir = 'real')
                self.gen_patho_metrics(KeyTblFinal, lsubdir = 'garbage')
            else:
                self.gen_patho_metrics(KeyTblFinal)
                self.write_patho_gv_session_info()
                if rm_rts_dup:
                    self.gen_patho_metrics(KeyTblFinal, lsubdir = 'nodupdir')
        if do_patho and do_umi_metrics:
            self.gen_patho_metrics_umi(KeyTblFinal)
        if use_dropseq and do_bc_dist:
            self.submit_ds_bcdistcore(KeyTblFinal, sample_prefix)
        if use_dropseq and do_drop_split:
            self.exe_drop_splitting(KeyTblFinal, sample_prefix)        
        if use_dropseq and do_drop_count:
            self.exe_drop_count(KeyTblFinal, sample_prefix)
        if use_dropseq and do_drop_metrics:
            self.exe_drop_metrics(KeyTblFinal, sample_prefix)
        if do_host and do_metrics:
            self.gen_host_metrics(KeyTblFinal)
            if rm_rts_dup:
                self.gen_host_metrics(KeyTblFinal, lsubdir = 'nodupdir')
        if do_host and do_umi_metrics:
            self.gen_host_metrics_umi(KeyTblFinal)
             

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process the options.')
    parser.add_argument('--config_file', '-c', type = str, required = True, help ='Path to main config file')
    parser.add_argument('--key_id', dest = 'key_id', required = False, default = 'none', help = 'Key file would be "key_id"_key.txt')
    parser.add_argument('--key_path', dest = 'key_path', required = False, default = 'none', help = 'Key file path (absolute)')
    parser.add_argument('--project_ids', required = False, default = 'none', help = 'Project ids')
    parser.add_argument('--seq_id', required = False, default = 'none', help = 'Sequencing id used to make raw_seq_path')
    parser.add_argument('--raw_seq_path', required = False, default = 'none', help = 'Directory for raw sequene files (absolute)')
    parser.add_argument('--temp_path', required = False, default = 'none', help = 'Will contain the temporary results')
    parser.add_argument('--bam_path', required = False, default = 'none', help = 'Will contain the aligned sorted bam files')
    parser.add_argument('--results_path', required = False, default = 'none', help = 'Will contain the path to the results')
    parser.add_argument("--do_patho", required=False, action="store_true", default=False, help="Do the patho")
    parser.add_argument("--do_host", required=False, action="store_true", default=False, help="Do the host")
    parser.add_argument("--remove_splitted", required=False, action="store_true", default=False, help="Remove the splitted files")
    parser.add_argument("--gzip_merged", required = False, action = "store_true", default = False, help = "Compress the merged files to gzip format")
    parser.add_argument('--host_aligner', dest = "host_aligner", type = str, \
            default = 'BBMap', help = 'Aligner for host (BBMap)')
    parser.add_argument('--read_counter', dest = "read_counter", type = str, \
        default = 'JL_counter', help = 'Tool for counting reads (default: Counter written by JLivny)')
    parser.add_argument('--no_p7', dest = 'use_p7', action = 'store_false', default = True, help = 'Use if P7 index is not used.' )
    parser.add_argument('--use_p5', dest = 'use_p5', action = 'store_true', default = False, help = 'Use if P5 index is used.' )
    parser.add_argument('--use_lane', dest = 'use_lane', action = 'store_true', default = False, help = 'Use if lane specific merging is required.' )
    parser.add_argument('--use_seq_path', dest = 'use_seq_path', action = 'store_true', default = False, help = 'Use if direct mapping of sample to raw seq path is required.' )
    parser.add_argument('--no_qsub', dest = 'use_qsub', action = 'store_false', default = True, help = 'Does not submit qsub jobs.' )
    parser.add_argument('--no_gs', dest = 'do_gs', action = 'store_false', default = True, help = 'Does not do the garbage sorting for dropseq' )
    parser.add_argument('--no_split', dest = 'do_split', action = 'store_false', default = True, help = 'Does not split the fastq files.' )
    parser.add_argument('--no_merge', dest = 'do_merge', action = 'store_false', default = True, help = 'Does not merge the split fastq files.' )
    parser.add_argument('--no_align', dest = 'do_align', action = 'store_false', default = True, help = 'Does not align.' )
    parser.add_argument('--no_count', dest = 'do_count', action = 'store_false', default = True, help = 'Does not count.' )
    parser.add_argument('--no_umi_count', dest = 'do_umi_count', action = 'store_false', default = True, help = 'Does not collapse/compute UMI.' )
    parser.add_argument('--no_metrics', dest = 'do_metrics', action = 'store_false', default = True, help = 'No metrics generation.' )
    parser.add_argument('--no_umi_metrics', dest = 'do_umi_metrics', action = 'store_false', default = True, help = 'Does not compute UMI metrics.' )
    parser.add_argument('--no_bc_dist', dest = 'do_bc_dist', action = 'store_false', default = True, help = 'Compute read distribution per barcode (dropseq only).' )
    parser.add_argument('--no_replace_refname', dest = 'do_replace_refname', action = 'store_false', default = True, help = 'Does not replace reference name in fna')
    parser.add_argument('--no_expand', dest = 'do_expand', action = 'store_false', default = True, help = 'Does not expand p7 and barcode entries in the key file if required.')
    parser.add_argument('--no_bc_split', dest = 'no_bc_split', action = 'store_true', default = False, help = 'Does not run the barcode splitter, instead create softlinks of the raw seq files to the split directory.')
    parser.add_argument('--no_drop_split', dest = 'do_drop_split', action = 'store_false', default = True, help = 'Does not split aligned sam files for dropseq.')
    parser.add_argument('--no_drop_count', dest = 'do_drop_count', action = 'store_false', default = True, help = 'Does not count reads for drops.')
    parser.add_argument('--no_drop_metrics', dest = 'do_drop_metrics', action = 'store_false', default = True, help = 'Does not create summarized dropseq count and metrics.')
    parser.add_argument('--no_ref', dest = 'do_ref', action = 'store_false', default = True, help = 'Does not generate patho ref.')
    parser.add_argument('--Suffix_s1', default = 'none', required = False, help = 'Update the value of Suffix_s1')
    parser.add_argument('--Suffix_s2', default = 'none', required = False, help = 'Update the value of Suffix_s2')
    parser.add_argument('--Suffix_ne', default = 'none', required = False, help = 'Update the value of Suffix_ne')
    parser.add_argument('--ADD5', dest = 'add5', type = int, default = 0, help = 'ADD5 for gff parser') 
    parser.add_argument('--ADD3', dest = 'add3', type = int, default = 0, help = 'ADD3 for gff parser') 
    parser.add_argument('--MOC_id', dest = 'MOC_id', type = str, default = 'none', help = 'Provide MOC string for adding MOC hierarchy')
    parser.add_argument('--trim_rs_5p', dest = 'trim_rs_5p', type = int, default = 0, help = '5p trim count for read-single') 
    parser.add_argument('--trim_rs_3p', dest = 'trim_rs_3p', type = int, default = 0, help = '3p trim count for read-single') 
    parser.add_argument('--keep_rs_5p', dest = 'keep_rs_5p', type = int, default = -1, help = '5p keep count for read-single') 
    parser.add_argument('--keep_rs_3p', dest = 'keep_rs_3p', type = int, default = -1, help = '3p keep count for read-single') 
    parser.add_argument('--trim_r1_5p', dest = 'trim_r1_5p', type = int, default = 0, help = '5p trim count for read1') 
    parser.add_argument('--trim_r1_3p', dest = 'trim_r1_3p', type = int, default = 0, help = '3p trim count for read1') 
    parser.add_argument('--trim_r2_5p', dest = 'trim_r2_5p', type = int, default = 0, help = '5p trim count for read2') 
    parser.add_argument('--trim_r2_3p', dest = 'trim_r2_3p', type = int, default = 0, help = '3p trim count for read2') 
    parser.add_argument('--keep_r1_5p', dest = 'keep_r1_5p', type = int, default = -1, help = '5p keep count for read1') 
    parser.add_argument('--keep_r1_3p', dest = 'keep_r1_3p', type = int, default = -1, help = '3p keep count for read1') 
    parser.add_argument('--keep_r2_5p', dest = 'keep_r2_5p', type = int, default = -1, help = '5p keep count for read2') 
    parser.add_argument('--keep_r2_3p', dest = 'keep_r2_3p', type = int, default = -1, help = '3p keep count for read2') 
    parser.add_argument('--MOC_id_ref', dest = 'MOC_id_ref', type = str, default = 'none', help = 'Adds MOC_id_ref to Bacterial_Ref_path')
    parser.add_argument('--AllSeq_con_A', dest = 'AllSeq_con_A', type = int, default = 7, help = 'Minimum number of consecutive A required to trim AllSeq read.') 
    parser.add_argument('--AllSeq_trim_len', dest = 'AllSeq_trim_len', type = int, default = 20, help = 'Minimum length of a read required to keep it after AllSeq trim') 
    parser.add_argument("--no_login_name", dest = "use_login_name", action = 'store_false', default = True, help = 'Generate results in a username specific directory')
    parser.add_argument("--min_resource", dest = "min_resource", action = 'store_true', default = False, help = 'Request for minimum resource during unicore allocation')
    parser.add_argument("--count_strand_rev", dest = "count_strand_rev", type = str, default = "N", help = 'Reverse the counting mechanism by making it "forward"')
    parser.add_argument("--do_bestacc", dest = "do_bestacc", action = "store_true", default = False, help = "Run bestacc after splitting and merging")
    parser.add_argument("--use_sample_id", dest = "use_sample_id", action = "store_true", default = False, help = "Use sample id for sample id mapping")
    parser.add_argument("--bwa_mem", dest = "bwa_mem", action = "store_true", default = False, help = "Use bwa mem instread of bwa backtrack")
    parser.add_argument("--rm_rts_dup", dest = "rm_rts_dup", action = "store_true", default = False, help = "Remove PCR duplicates from RNA-tagseq")
    parser.add_argument("--do_rerun", dest = "do_rerun", action = "store_true", default = False, help = "Rerun pipeline after a previous incomplete run")
    parser.add_argument("--paired_only_patho", dest = "paired_only_patho", action = "store_true", default = False, help = "Align/count only the reads when both ends are mapped on pathogen side")
    parser.add_argument("--paired_only_host", dest = "paired_only_host", action = "store_true", default = False, help = "Align/count only the reads when both ends are mapped on host side")
    parser.add_argument('--ucore_time', dest = 'ucore_time', type = int, default = 0, help = 'Timelimit of unicore for this run') 
    
    
    
    options = parser.parse_args()

    pipel = PipelineII(options) 
    pipel.mainFunc()

