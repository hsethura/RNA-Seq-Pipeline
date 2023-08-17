#!/usr/bin/env python

import yaml
import argparse
import shutil
import pandas as pd
import csv
import re
import os
import sys
import ntpath
import codecs
import getpass
from Parse_featureCounts import get_all_gene_counts
from Parse_featureCounts import get_all_metrics_counts
from patho import copyLargeFile

from subprocess import call

nirmalya_lib = "/broad/IDP-Dx_work/nirmalya/local/lib/python2.7/site-packages"
sys.path.append(nirmalya_lib)
import PyPDF2
from PyPDF2 import PdfFileMerger


class Metrics:
    def __init__(self, confd):
        self.confd = confd


    def combine_umi_coll_pdfs(self, lproject_id, ref_acc, lsample_lst, species_type):
        # Have to work on this; have to provide the proper directories
        confd = self.confd
        ldelim = confd.ldelim
        project_dir = None
        if "patho" == species_type:
            Patho_dir = confd.Patho_dir
            project_dir = Patho_dir + ldelim + lproject_id
        elif "host" == species_type:
            Host_dir = confd.Host_dir
            project_dir = Host_dir + ldelim + lproject_id
        else:
            raise StandardError("Unknown species type: " + species_type)

        umidir = project_dir + ldelim + "umidir"
        umi_logdir = umidir + ldelim + "logdir"

        merger = PdfFileMerger()

        for lsample in lsample_lst:
            # Check if the pdf file has been created 
            coll_len_file = umi_logdir + ldelim + lsample + "_" + ref_acc + "_coll_len.txt"
            coll_len_size = os.stat(coll_len_file).st_size
            if coll_len_size > 0:
                lpdf_file = umi_logdir + ldelim + lsample + "_" + ref_acc + "_coll_len_hist.pdf"
                if os.path.isfile(lpdf_file) == False:
                    raise StandardError("coll_len_pdf does not exists: " + lpdf_file)
                else:
                    merger.append(lpdf_file)
            elif coll_len_size == 0:
                print("coll_len.txt file is empty. no hist pdf: " + coll_len_file)
            else:
                raise StandardError("Illegal value of coll_len_size: " + coll_len_size + ", sample: " + lsample)

        Results_path = confd.Results_path
        result_dir = Results_path + ldelim + lproject_id
        result_file = result_dir + "/umidir/" + lproject_id + "_" + ref_acc + "_coll_len_hist_combined.pdf"
        merger.write(result_file)
        
    def exe_picard_metrics(self, bamfile, fnafile, outdir):
        """ This would run Jonathan's script for running picard metrics 
        generator.
        sh /broad/IDP-Dx_work/nirmalya/pipeline/beta/PICARD_metrics.sh 
        <path to bam file (no sorting needed)> <path to fna> <path to output 
        directory> <path to picard bin>
        """
        confd = self.confd
        picard_metrics = confd.picard_metrics
        picard_bindir = confd.picard_bindir
        picard_command = "sh " + picard_metrics + " " + bamfile + " " + fnafile + " " + outdir + " " + picard_bindir
        print("Picard metrics command: " + picard_command)
        call(picard_command.split())

    def exe_picard_parser(self, lproject_id, species_type, use_umi = False, subdir = ''):
        confd = self.confd
        ldel = confd.ldelim
        picard_parser_path = confd.picard_metrics_parse
        Results_path = confd.Results_path
        lresult_dir = Results_path + ldel + lproject_id
        if subdir:
            lresult_dir += ldel + subdir
        if use_umi:
            lresult_dir += ldel + "umidir"
        picard_outdir = None
        picard_prefix = lresult_dir + ldel + "bam_metrics"
        if "patho" == species_type:
            picard_outdir = picard_prefix + "_patho"
        elif "host" == species_type:
            picard_outdir =  picard_prefix + "_host"
        else:
            raise StandardError("Unknown species type: " + species_type)

        outfile_path =  lresult_dir + ldel + lproject_id + "_" + \
            species_type + "_AlignmentSummaryMetrics.txt"
        picard_parser_cmd = "sh " + picard_parser_path + " " + picard_outdir +\
            " " + outfile_path
        print("picard_parser_cmd: " + picard_parser_cmd)
        call(picard_parser_cmd.split())


    def write_corr_file(self, lproject_id, ref_acc, species_type, use_umi = False, subdir = ''):
        confd = self.confd
        Summary_dir = confd.Summary_dir
        ldelim = confd.ldelim
        Results_path = confd.Results_path
        corr_script = confd.corr_script

        loutdir = Results_path + ldelim + lproject_id
        if subdir:
            loutdir += ldelim + subdir
        if use_umi:
            loutdir += ldelim + "umidir"
        if species_type == "host":
            tsv_file = loutdir + ldelim +  lproject_id + "_" + ref_acc + "_counts.tsv"
        elif species_type == "patho":
            tsv_file = loutdir + ldelim +  lproject_id + "_" + ref_acc + "_fpkm.tsv"
        else:
            raise StandardError("Unknown species type: " + species_type)


        copy_dir = Results_path + ldelim + lproject_id
        if subdir:
            copy_dir += ldelim + subdir
        if use_umi:
            copy_dir += ldelim + "umidir"
        if not os.path.exists(copy_dir):
            os.makedirs(copy_dir)
        corr_file = copy_dir + ldelim +  lproject_id + "_" + ref_acc + "_corr.txt"
        corr_cmd = "Rscript " + corr_script + " " + tsv_file + " " + corr_file
        print("corr_cmd: " + corr_cmd)
        call(corr_cmd.split())
                       
    def write_read_count_table(self, lproject_id, ref_acc, lsample_lst, species_type, has_header, use_umi = False, subdir = ''):
        confd = self.confd
        Summary_dir = confd.Summary_dir
        ldelim = confd.ldelim 
        Results_path = confd.Results_path

        loutdir = Summary_dir + ldelim + lproject_id
        if subdir:
            loutdir += ldelim + subdir
        if use_umi:
            loutdir += ldelim + "umidir"

        if not os.path.exists(loutdir):
            os.makedirs(loutdir)
        # loutfile: location of the read count table in the tempdir
        loutfile = loutdir + ldelim +  lproject_id + "_" + ref_acc + "_counts.tsv"
        print(lproject_id + " :: " + ref_acc + " :: " + '; '.join(lsample_lst))

        project_dir = None
        if "patho" == species_type:
            Patho_dir = confd.Patho_dir
            project_dir = Patho_dir + ldelim + lproject_id
        elif "host" == species_type:
            Host_dir = confd.Host_dir
            project_dir = Host_dir + ldelim + lproject_id
        else:
            raise StandardError("Unknown species type: " + species_type)

        if subdir:
            project_dir += ldelim + subdir
        if use_umi:
            project_dir += ldelim + "umidir"
        get_all_gene_counts(project_dir, loutfile, lsample_lst, ref_acc, has_header = has_header)

        copy_dir = Results_path + ldelim + lproject_id
        if subdir:
            copy_dir += ldelim + subdir
        if use_umi:
            copy_dir += ldelim + "umidir"
        print("copy_dir: " + copy_dir)
        if not os.path.exists(copy_dir):
            os.makedirs(copy_dir)
        lout_copy = copy_dir + ldelim + lproject_id + "_" + ref_acc + "_counts.tsv"
        copyLargeFile(loutfile, lout_copy) 


    def write_fpkm_file(self, lproject_id, ref_acc, use_umi = False, subdir = ''):
        confd = self.confd
        ldel = confd.ldelim
        Summary_dir = confd.Summary_dir
        summary_path = Summary_dir + ldel + lproject_id
        if subdir:
            summary_path += ldel + subdir
        if use_umi:
            summary_path += ldel + "umidir"
        count_file = summary_path + ldel + lproject_id + "_" + ref_acc + \
            "_counts.tsv"
        fpkm_dir = summary_path + ldel + "fpkm_tempdir"
        if not os.path.exists(fpkm_dir):
            os.makedirs(fpkm_dir)
        Results_path = confd.Results_path
        copy_dir = Results_path + ldel + lproject_id
        if subdir:
            copy_dir += ldel + subdir
        if use_umi:
            copy_dir += ldel + "umidir"
        fpkm_file = fpkm_file = copy_dir + ldel + lproject_id + "_" + \
            ref_acc + "_fpkm.tsv"
        fpkm_script = confd.fpkm_script
        fpkm_cmd = "sh " + fpkm_script + " " + count_file + " " + \
            fpkm_file + " " + fpkm_dir
        print("fpkm_cmd: " + fpkm_cmd)
        call(fpkm_cmd.split())

    def copy_fastq_gff(self, lproject_id, ref_acc, use_umi = False, subdir = ''):
        confd = self.confd
        # Also copy the fastq files and gff files to the appropriate directory
        Data_dir = confd.Data_dir
        Bam_path = confd.Bam_path
        ldelim = confd.ldelim
        # Get the path to the gff and fna for the ref_acc
        fna_path = Data_dir + ldelim + ref_acc + ".fna"
        gff_path = Data_dir + ldelim + ref_acc + "_GENES.gff"
        bam_outpath_proj = Bam_path + ldelim + lproject_id

        if subdir:
            bam_outpath_proj += ldelim + subdir
        if use_umi:
            bam_outpath_proj += ldelim + "umidir"
        if not os.path.exists(bam_outpath_proj):
            os.makedirs(bam_outpath_proj)
        fna_outpath = bam_outpath_proj + ldelim + ref_acc + ".fna"
        gff_outpath = bam_outpath_proj + ldelim + ref_acc + "_GENES.gff" 
        copyLargeFile(fna_path, fna_outpath)
        copyLargeFile(gff_path, gff_outpath)

    def copy_fastq_etc_host(self, lproject_id, ref_acc, use_umi = False, subdir = ''):
        confd = self.confd
        Bam_path = confd.Bam_path
        ldel = confd.ldelim
        host_dbpath = confd.host_dbpath

        host_ref_str = confd.host_ref_str

       
        host_ref_path = host_dbpath + ldel + "data"
        fna_path = host_ref_path + ldel + host_ref_str + ".fna"
        trans_ge_path = host_ref_path + ldel + host_ref_str + "_transcript_gene.txt"

        bam_outpath_proj = Bam_path + ldel + lproject_id

        if subdir:
            bam_outpath_proj += ldel + subdir
        if use_umi:
            bam_outpath_proj += ldel + "umidir"
        if not os.path.exists(bam_outpath_proj):
            os.makedirs(bam_outpath_proj)

        fna_outpath = bam_outpath_proj + ldel + host_ref_str + ".fna"
        trans_ge_outpath = bam_outpath_proj + ldel + host_ref_str + "_transcript_gene.txt"
        print("fna_path: " + fna_path + ", fna_outpath: " + fna_outpath)
        copyLargeFile(fna_path, fna_outpath)
        copyLargeFile(trans_ge_path, trans_ge_outpath)


    def get_suffix(self):
        confd = self.confd
        suffix = None
        LC_method_val = confd.LC_method_val
        Read_pairing_val = confd.Read_pairing_val
        lc_lower = confd.LC_method_val.lower()
        if lc_lower == 'allseq' or lc_lower == 'alignr2' or lc_lower=='alignr1':
            suffix = "se"
        elif lc_lower == 'bdropseq' or lc_lower == 'indrops':
            suffix = "se"
        elif lc_lower == 'rts' or lc_lower == 'rts-ts':
            if Read_pairing_val == 'PAIRED':
                suffix = "pe"
            elif Read_pairing_val == 'SINGLE':
                suffix = "se"
        elif lc_lower == 'smarter':
            suffix = "pe"
        return suffix

   
    def write_bam_file_paths(self, lproject_id, ref_acc, lsample_lst, view_mode, append_mode, species_type, use_umi = False, subdir = ''):
        """ Create a file in the output dir containing the path to the bam 
        files.
        """
        confd = self.confd
        Results_path = confd.Results_path
        ldelim = confd.ldelim
        Bam_path = confd.Bam_path

        result_dir = Results_path + ldelim + lproject_id
        if subdir:
            result_dir += ldelim + subdir
        if use_umi:
            result_dir += ldelim + "umidir"

        bam_dir_ori = Bam_path + ldelim + lproject_id
        if view_mode == "gv":
            mount_point = "/broad/IDP-Dx_storage/idpweb"
            idpweb_address = "http://idpweb.broadinstitute.org"
            bam_dir = bam_dir_ori.replace(mount_point, idpweb_address)
        else:
            bam_dir = bam_dir_ori
        if subdir:
            bam_dir += ldelim + subdir
        if use_umi:
            bam_dir += ldelim + "umidir"

        bam_path_file = result_dir + ldelim + lproject_id + "_" + species_type + "_bam_paths_" + \
            view_mode + ".txt"
        write_bam_f = ''
        if not append_mode:
            write_bam_f = open(bam_path_file, 'w')
        else:
            write_bam_f = open(bam_path_file, 'a')

        if view_mode == "gv" and species_type == "patho":
            start_str = "##GenomeView session"
            write_bam_f.write(start_str + "\n")
        fna_outpath = bam_dir + ldelim + ref_acc + ".fna"
        if view_mode == "gv" and species_type == "patho":
            fna_outpath = "DATA\t" + fna_outpath
        write_bam_f.write(fna_outpath + '\n')

        if species_type == "patho":
            gff_outpath = bam_dir + ldelim + ref_acc + "_GENES.gff"
            if view_mode == "gv":
                gff_outpath = "DATA\t" + gff_outpath
            write_bam_f.write(gff_outpath + '\n')
        elif species_type == "host":
            trans_ge_file = bam_dir + ldelim + ref_acc + "_transcript_gene.txt"
            write_bam_f.write(trans_ge_file + '\n')

        suffix = self.get_suffix()
        for lsample in lsample_lst:
            print("bam_dir: " + bam_dir)
            print("lsample: " + lsample)
            print("ref_acc: " + ref_acc)
            print("suffix: " + suffix)
            lfile_path = bam_dir + ldelim + lsample + "_" + ref_acc + "_" + suffix + ".bam"
            if view_mode == "gv" and species_type == "patho":
                lfile_path = "EXTRA\t" + lfile_path
            write_bam_f.write(lfile_path + '\n')
            if species_type == "patho":
                tdf_lfile_path = bam_dir + ldelim + lsample + "_" + ref_acc + \
                    "_" + suffix + ".bam.tdf"
                if view_mode == "gv":
                    tdf_lfile_path = "EXTRA\t" + tdf_lfile_path
                write_bam_f.write(tdf_lfile_path + '\n')
        write_bam_f.close()

    def write_host_metrics_table(self, lproject_id, hostref, lsample_lst, use_umi = False, subdir = ''):
        confd = self.confd
        ldel = confd.ldelim
        Host_dir = confd.Host_dir
        project_dir = Host_dir + ldel + lproject_id
        Results_path = confd.Results_path
        results_dir = Results_path + ldel + lproject_id

        if subdir:
            project_dir += ldel + subdir
            results_dir += ldel + subdir
        if use_umi:
            project_dir += ldel + "umidir"
            results_dir += ldel + "umidir"
        lhost_metrics_file = results_dir + ldel + lproject_id + "_" + hostref + "_metrics.tsv"
        get_all_metrics_counts(project_dir, lhost_metrics_file, lsample_lst, ref_acc = hostref)

    def data_finish(self, lproject_id):
        confd = self.confd
        Out_dir = confd.Out_dir
        ldelim = confd.ldelim
        Results_path = confd.Results_path
        config_file = confd.config_file
        Script_dir = confd.basepath
        lscript_dir = Script_dir + ldelim + "shell_scripts"
        Data_finish_script = confd.Data_finish_script
        data_finish_cmd = "sh " + Data_finish_script + " " +  lproject_id + " " + Out_dir + " " + Results_path + " " + config_file + " " + lscript_dir
        print("Data_finish_cmd: " + data_finish_cmd)
        call(data_finish_cmd.split())
        print("Data_finish done")
 
    def ref_RPG_metrics(self, lproject_id, use_umi = False, subdir = ''):
        confd = self.confd
        Out_dir = confd.Out_dir
        ldelim = confd.ldelim
        
        Results_path = confd.Results_path
        #Run scripts for Jonathans metrics
        RPG_metrics_script = confd.RPG_metrics_script
      
        RPG_metrics_cmd = "sh " + RPG_metrics_script + " -p " +  lproject_id + " -t " + Out_dir + " -r " + Results_path
        if use_umi == True:
            RPG_metrics_cmd += " -u "
        if subdir:
            RPG_metrics_cmd += " -s " + subdir
        print("RPG_metrics_cmd: " + RPG_metrics_cmd)
        call(RPG_metrics_cmd.split())

    def gen_bcLog(self, lproject_id):
        confd = self.confd
        Split_dir = confd.Split_dir
        ldelim = confd.ldelim
        Results_path = confd.Results_path
        Pipeline_bcLog = confd.Pipeline_bcLog

        lresult_dir = Results_path + ldelim + lproject_id
        #if not os.path.exists(lresult_dir):
        #    os.makedirs(lresult_dir)
        bcLog_com = "sh " + Pipeline_bcLog + " " + Split_dir + " " + lresult_dir + " " + lproject_id
        call(bcLog_com.split())
        print("Call bcLog_command: " + bcLog_com)

    def gen_all_bcLog(self):
        confd = self.confd
        use_p7 = confd.use_p7
        project_set = confd.project_set

        if use_p7:
            for lproject_id in project_set:
                self.gen_bcLog(lproject_id)

