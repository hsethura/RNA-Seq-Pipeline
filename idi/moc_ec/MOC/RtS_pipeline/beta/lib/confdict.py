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

class ConfDict(object):

    def __init__(self, options, basepath):
        # Store all the properties to the dict "datastore"
        object.__setattr__(self, 'datastore', {})
        self.datastore['basepath'] = basepath

        self.options = options
        self.config_file = options.config_file
        self.mydict = yaml.load(open(self.config_file))
        self.ldelim = "/"
        self.inidict = dict()
        self.popu_inidict()

        self.file_path = os.path.abspath(__file__)
        self.project_root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(self.file_path)))))))

    def popu_inidict(self):
        inidict = self.inidict
        inidict['LC_method'] = 'LC_method'
        inidict['Read_pairing'] = 'Read_pairing'

    def __getattr__(self, key):
        if key in ['KeyTbl', 'options', 'mydict']:
            return super(ConfDict, self).__getattr__(key)
        else:
            return self.datastore[key]

    def __setattr__(self, key, value):
        if key in ['KeyTbl', 'options', 'mydict']:
            super(ConfDict, self).__setattr__(key, value)
        else:
            self.datastore[key] = value

    def get_project_set(self, ltab):
        project_id = self.project_id
        project_lst = ltab[project_id].tolist()
        project_set = list(set(project_lst))
        return project_set

    def loadKeyTable(self):
        Key_path = self.Key_path
        sample_id = self.sample_id
        lfile = codecs.open(Key_path, encoding='utf-8')
        tab0 = csv.reader(lfile, delimiter = '\t')
        tab1 = [row for row in tab0]
        tab2 = filter(lambda x: not x[0].startswith("###"), tab1)
        header = tab2.pop(0)	
        tab3 = pd.DataFrame(tab2, columns = header)
        for index, row in tab3.iterrows():
            lsample = row[sample_id]
            lsample2 = re.sub('\s+', '_', lsample)
            lsample3 = re.sub('/', '_', lsample2)
            tab3.loc[index, sample_id] = lsample3
        self.KeyTbl = tab3

    def storeConfigFromKeyTbl(self):
        KeyTbl = self.KeyTbl
        options = self.options
        mydict = self.mydict

        # Here we are modifying the code a little and still trying to keep it
        # backward compatible. If the LC_method_val obtained from the key file
        # ends with _alr1 ot _alr2, then set those flags on and keep that part
        # separate from the rest of LC_method_val

        LCmv_ori = KeyTbl[self.LC_method][0]
        self.align_readnum = ''
        LCmv_ori_u = LCmv_ori.lower()
        if LCmv_ori_u.endswith('_alr1') or LCmv_ori_u.endswith('_alr2'):
            lparts = LCmv_ori.rsplit('_', 1)
            self.LC_method_val = lparts[0]
            self.align_readnum = lparts[1].lower() 
        else:
            self.LC_method_val = LCmv_ori

        Read_pairing_val_ci = KeyTbl[self.Read_pairing][0]
        if not Read_pairing_val_ci:
            Read_pairing_val_ci = 'Paired'
            print("Rad_pairing_val empty. setting it: " + Read_pairing_val_ci)
        Read_pairing_val = Read_pairing_val_ci.upper()
        self.Read_pairing_val = Read_pairing_val.upper()


        if self.do_host:
            # Also get the name of the host
            l_host_str = KeyTbl[self.Host_reference][0]
            if re.search('human', l_host_str, re.IGNORECASE):
                self.host_ref_str = self.human_ref_str
                self.host_transcript_gene = self.human_transcript_gene
            elif re.search('mouse', l_host_str, re.IGNORECASE):
                self.host_ref_str = self.mouse_ref_str
                self.host_transcript_gene = self.mouse_transcript_gene
            elif re.search('rabbit', l_host_str, re.IGNORECASE):
                self.host_ref_str = self.rabbit_ref_str
                self.host_transcript_gene = self.rabbit_transcript_gene
            elif re.search('zebrafish', l_host_str, re.IGNORECASE):
                self.host_ref_str = self.zebrafish_ref_str
                self.host_transcript_gene = self.zebrafish_transcript_gene
            elif re.search('fungal', l_host_str, re.IGNORECASE):
                self.host_ref_str = self.fungal_ref_str
                self.host_transcript_gene = self.fungal_transcript_gene
            else:
                raise ValueError('Wrong Host_reference: ' + l_host_str)
         
        is_allseq = False
        self.use_dropseq = False
        self.drop_type = None

        self.soft_link_from_input = False
        lbc_splitter = ''
        lc_lower = self.LC_method_val.lower()
        if lc_lower == 'allseq':
            lbc_splitter = self.bc_splitter
            is_allseq = True
            self.Dict_infile = self.AllSeq_dict_file
            self.add5 = 30
            self.add3 = 20
        elif lc_lower == 'bdropseq':
            self.use_dropseq = True    
            self.drop_type = 'bdropseq'
            self.add5 = 20
            self.add3 = 30
        elif lc_lower == 'indrops':
            self.use_dropseq = True
            self.drop_type = 'indrops'
            self.add5 = 20
            self.add3 = 30
        elif lc_lower == 'rts' or lc_lower == 'rts-ts':
            if Read_pairing_val == 'PAIRED':
                lbc_splitter = self.bc_splitter_rts
            elif Read_pairing_val == 'SINGLE':
                lbc_splitter = self.bc_splitter_rts_se
            self.Dict_infile = self.RtS_dict_file
        elif lc_lower == 'scr':
            lbc_splitter = self.bc_splitter_scr
            self.Dict_infile = self.SCR_dict_file
        elif lc_lower == 'smarter':
            # We assume that for smarter there would be no split, but if there
            # is any split, we shall look into that in the future.
            self.soft_link_from_input = True
            lbc_splitter = 'None'
        elif lc_lower == 'alignr2':
            self.Dict_infile = self.RtS_dict_file
        elif lc_lower == 'alignr1':
            self.Dict_infile = self.RtS_dict_file

        umi_count_t = self.do_umi_count
        umi_metrics_t = self.do_umi_metrics
        self.do_umi_count = umi_count_t and is_allseq
        self.do_umi_metrics = umi_metrics_t and is_allseq

        no_bc_split = self.no_bc_split
        if no_bc_split:
            # This may not be required: Disturbing the dropseq
            #self.soft_link_from_input = True
            lbc_splitter = 'None'
                

        do_R2_trim = False
        if lc_lower == 'rts-ts':
            do_R2_trim = True

        do_allseq_trim = False
        if lc_lower == 'allseq':
            do_allseq_trim = True

        self.lbc_splitter = lbc_splitter
        self.do_R2_trim = do_R2_trim
        self.do_allseq_trim = do_allseq_trim
        self.lc_lower = lc_lower


    def get_from_mydict(self, key):
        mydict = self.mydict
        inidict = self.inidict

        if key in mydict:
            return mydict[key]
        elif key in inidict:
            return inidict[key]
        else:
            return '__NA__'

    def storeConfigFromConfig(self):
        mydict = self.mydict

        self.project_id      = self.get_from_mydict('Proj')
        self.P7_index        = self.get_from_mydict('P7')
        self.P5_index        = self.get_from_mydict('P5')
        self.Lane_index      = self.get_from_mydict('Lane')
        self.sample_id       = self.get_from_mydict('ID')
        self.LC_method       = self.get_from_mydict('LC_method')
        self.Path_to_SeqFile = self.get_from_mydict('Path_to_SeqFile')
        self.Ref_accession   = self.get_from_mydict('Ref_accession')
        self.Host_reference  = self.get_from_mydict('Host_reference')
        self.lbc             = self.get_from_mydict('bc')         
        self.ldelim          = self.get_from_mydict('delim')

        self.dict_builder    = self.get_from_mydict('dict_builder')
        if not os.path.isabs(self.dict_builder):
            self.dict_builder = os.path.join(self.project_root_dir, self.dict_builder)

        self.bash_header = self.get_from_mydict('bash_header')
        if not os.path.isabs(self.bash_header):
            self.bash_header = os.path.join(self.project_root_dir, self.bash_header)

        self.bc_splitter     = self.get_from_mydict('bc_splitter')
        if not os.path.isabs(self.bc_splitter):
            self.bc_splitter = os.path.join(self.project_root_dir, self.bc_splitter)

        self.bc_splitter_rts = self.get_from_mydict('bc_splitter_rts')
        if not os.path.isabs(self.bc_splitter_rts):
            self.bc_splitter_rts = os.path.join(self.project_root_dir, self.bc_splitter_rts) 

        self.bc_splitter_rts_se = self.get_from_mydict('bc_splitter_rts_se')
        if not os.path.isabs(self.bc_splitter_rts_se):
            self.bc_splitter_rts_se = os.path.join(self.project_root_dir, self.bc_splitter_rts_se)

        self.bc_splitter_scr = self.get_from_mydict('bc_splitter_scr')
        if not os.path.isabs(self.bc_splitter_scr):
            self.bc_splitter_scr = os.path.join(self.project_root_dir, self.bc_splitter_scr)

        self.sam_fragcount   = self.get_from_mydict('sam_fragcount')
        if not os.path.isabs(self.sam_fragcount):
            self.sam_fragcount = os.path.join(self.project_root_dir, self.sam_fragcount)

        self.paired_only_script = self.get_from_mydict('paired_only_script')
        if not os.path.isabs(self.paired_only_script):
            self.paired_only_script = os.path.join(self.project_root_dir, self.paired_only_script)

        self.frag_to_gene_count = self.get_from_mydict('frag_to_gene_count')
        if not os.path.isabs(self.frag_to_gene_count):
            self.frag_to_gene_count = os.path.join(self.project_root_dir, self.frag_to_gene_count)

        self.metrics_gen     =  self.get_from_mydict('metrics_gen')
        if not os.path.isabs(self.metrics_gen):
            self.metrics_gen = os.path.join(self.project_root_dir, self.metrics_gen)

        self.UGER_cbp        = self.get_from_mydict('UGER_cbp')
        if not os.path.isabs(self.UGER_cbp):
            self.UGER_cbp = os.path.join(self.project_root_dir, self.UGER_cbp)

        self.bwa_path = self.get_from_mydict('bwa')
        if not os.path.isabs(self.bwa_path):
            self.bwa_path = os.path.join(self.project_root_dir, self.bwa_path)

        self.split_plot_path = self.get_from_mydict('split_plot_path')
        if not os.path.isabs(self.split_plot_path):
            self.split_plot_path = os.path.join(self.project_root_dir, self.split_plot_path)

        self.patho_dbpath    = self.get_from_mydict('patho_dbpath')
        self.host_dbpath     = self.get_from_mydict('host_dbpath')
        self.Read_pairing    = self.get_from_mydict('Read_pairing')

        self.STAR            = self.get_from_mydict('STAR')
        if not os.path.isabs(self.STAR):
            self.STAR = os.path.join(self.project_root_dir, self.STAR)

        self.samtools       = self.get_from_mydict('samtools')
        if not os.path.isabs(self.samtools):
            self.samtools = os.path.join(self.project_root_dir, self.samtools)

        self.featureCounts  = self.get_from_mydict('featureCounts')
        if not os.path.isabs(self.featureCounts):
            self.featureCounts = os.path.join(self.project_root_dir, self.featureCounts)

        self.JLCounter      = self.get_from_mydict('JLCounter')
        if not os.path.isabs(self.JLCounter):
            self.JLCounter = os.path.join(self.project_root_dir, self.JLCounter)

        self.strand_dir     = self.get_from_mydict('strand_dir')
        self.tdf_str        = self.get_from_mydict('tdf_str')
        self.human_ref_str  = self.get_from_mydict('human_ref_str')
        self.mouse_ref_str  = self.get_from_mydict('mouse_ref_str')
        self.rabbit_ref_str  = self.get_from_mydict('rabbit_ref_str')
        self.zebrafish_ref_str = self.get_from_mydict('zebrafish_ref_str')
        self.fungal_ref_str = self.get_from_mydict('fungal_ref_str')

        self.bbmap_path     = self.get_from_mydict('bbmap_path')
        if not os.path.isabs(self.bbmap_path):
            self.bbmap_path = os.path.join(self.project_root_dir, self.bbmap_path)

        self.patho_thread_count = self.get_from_mydict('patho_thread_count')
        self.patho_memory   = self.get_from_mydict('patho_memory')
        self.host_thread_count = self.get_from_mydict('host_thread_count')
        self.host_memory    = self.get_from_mydict('host_memory')
        self.human_transcript_gene = self.get_from_mydict('human_transcript_gene')
        self.mouse_transcript_gene = self.get_from_mydict('mouse_transcript_gene')
        self.rabbit_transcript_gene = self.get_from_mydict('rabbit_transcript_gene')
        self.zebrafish_transcript_gene = self.get_from_mydict('zebrafish_transcript_gene')
        self.fungal_transcript_gene = self.get_from_mydict('fungal_transcript_gene')

        self.picard_bindir = self.get_from_mydict('picard_bindir')
        if not os.path.isabs(self.picard_bindir):
            self.picard_bindir = os.path.join(self.project_root_dir, self.picard_bindir)

        self.allseq_trim_script = self.get_from_mydict('AllSeq_read_trim')
        if not os.path.isabs(self.allseq_trim_script):
            self.allseq_trim_script = os.path.join(self.project_root_dir, self.allseq_trim_script)

        self.cutadapt       = self.get_from_mydict('cutadapt')   
        if not os.path.isabs(self.cutadapt):
            self.cutadapt = os.path.join(self.project_root_dir, self.cutadapt)

        self.AllSeq_dict_file = self.get_from_mydict('AllSeq_dict_file')
        if not os.path.isabs(self.AllSeq_dict_file):
            self.AllSeq_dict_file = os.path.join(self.project_root_dir, self.AllSeq_dict_file)

        self.RtS_dict_file = self.get_from_mydict('RtS_dict_file')
        if not os.path.isabs(self.RtS_dict_file):
            self.RtS_dict_file = self.project_root_dir + self.RtS_dict_file 

        self.SCR_dict_file = self.get_from_mydict('SCR_dict_file')
        if not os.path.isabs(self.SCR_dict_file):
            self.SCR_dict_file = self.project_root_dir + self.SCR_dict_file 

        self.Suffix_s1 = self.get_from_mydict('Suffix_s1')
        self.Suffix_s2 = self.get_from_mydict('Suffix_s2')
        self.Suffix_ne = self.get_from_mydict('Suffix_ne')
        self.uminorm_path = self.get_from_mydict('uminorm_path')
        self.garbage_sorter_path = self.get_from_mydict('garbage_sorter')
        self.garbage_sorter_ind_path = self.get_from_mydict('garbage_sorter_ind')
        self.indrops_bc1_file = self.get_from_mydict('indrops_bc1_file') 
        self.indrops_bc2_file = self.get_from_mydict('indrops_bc2_file') 
        self.ds_bc_counter = self.get_from_mydict('ds_bc_counter')
        self.read_trimmer_path = self.get_from_mydict('read_trimmer')
        self.sam_splitter_path = self.get_from_mydict('sam_splitter')
        self.Results_path = self.get_from_mydict('Results_path')
        
        # Added as a fix for issue 41 (https://github.com/broadinstitute/PipelineII/issues/41)
        self.Pipeline_bcLog = self.get_from_mydict('Pipeline_bcLog')
        if not os.path.isabs(self.Pipeline_bcLog):
            self.Pipeline_bcLog = os.path.join(self.project_root_dir, self.Pipeline_bcLog)

        self.corr_script = self.get_from_mydict('corr_script')
        if not os.path.isabs(self.corr_script):
            self.corr_script = os.path.join(self.project_root_dir, self.corr_script)

        self.gff_parser = self.get_from_mydict('gff_parser')
        if not os.path.isabs(self.gff_parser):
            self.gff_parser = os.path.join(self.project_root_dir, self.gff_parser)

        self.RPG_metrics_script = self.get_from_mydict('RPG_metrics_script')
        if not os.path.isabs(self.RPG_metrics_script):
            self.RPG_metrics_script = os.path.join(self.project_root_dir, self.RPG_metrics_script)

        self.RPG_drop_metrics_script = self.get_from_mydict('RPG_drop_metrics_script')
        if not os.path.isabs(self.RPG_drop_metrics_script):
            self.RPG_drop_metrics_script = os.path.join(self.project_root_dir, self.RPG_drop_metrics_script)

        self.Data_finish_script = self.get_from_mydict('Data_finish_script')
        if not os.path.isabs(self.Data_finish_script):
            self.Data_finish_script = os.path.join(self.project_root_dir, self.Data_finish_script)
            
        self.bcLog_metrics_script = self.get_from_mydict('bcLog_metrics_script')
        if not os.path.isabs(self.bcLog_metrics_script):
            self.bcLog_metrics_script = os.path.join(self.project_root_dir, self.bcLog_metrics_script)

        self.trim_script = self.get_from_mydict('trim_script')
        if not os.path.isabs(self.trim_script):
            self.trim_script = os.path.join(self.project_root_dir, self.trim_script)

        self.picard_metrics = self.get_from_mydict('picard_metrics')
        if not os.path.isabs(self.picard_metrics):
            self.picard_metrics = os.path.join(self.project_root_dir, self.picard_metrics)

        self.picard_metrics_parse = self.get_from_mydict('picard_metrics_parse')
        if not os.path.isabs(self.picard_metrics_parse):
            self.picard_metrics_parse = os.path.join(self.project_root_dir, self.picard_metrics_parse)

        self.fpkm_script = self.get_from_mydict('fpkm_script')
        if not os.path.isabs(self.fpkm_script):
            self.fpkm_script = os.path.join(self.project_root_dir, self.fpkm_script)

        self.bestacc_script = self.get_from_mydict('bestacc_script')
        if not os.path.isabs(self.bestacc_script):
            self.bestacc_script = os.path.join(self.project_root_dir, self.bestacc_script)

        self.remove_dup_script = self.get_from_mydict('remove_dup_script')
        if not os.path.isabs(self.remove_dup_script):
            self.remove_dup_script = os.path.join(self.project_root_dir, self.remove_dup_script)
        

    def storeDerivedPaths(self):
        ldelim = self.ldelim
        Script_dir = self.basepath 
        self.UniCore      = Script_dir + ldelim + "lib" + ldelim + "unicore.py" 
        self.UniDropCore = Script_dir + ldelim + "lib" + ldelim + "unidropcore.py"
        self.DropMetrics = Script_dir + ldelim + "lib" + ldelim + "dropmetrics.py"
        self.UniMerger      = Script_dir + ldelim + "lib" + ldelim + "unimerger.py" 
        self.BCDistCore_path = Script_dir + ldelim + "lib" + ldelim + "bcdistcore.py"
        self.umi_coll_len_script = Script_dir + ldelim + "rscripts" + ldelim + \
            "umi_coll_len.R"


    def storeConfigFromOptions(self):
        options = self.options
        self.Project_ids = options.project_ids
        self.do_replace_refname = options.do_replace_refname
        self.add5 = options.add5
        self.add3 = options.add3
        self.trim_rs_5p = options.trim_rs_5p
        self.trim_rs_3p = options.trim_rs_3p
        self.keep_rs_5p = options.keep_rs_5p
        self.keep_rs_3p = options.keep_rs_3p
        self.trim_r1_5p = options.trim_r1_5p
        self.trim_r1_3p = options.trim_r1_3p
        self.trim_r2_5p = options.trim_r2_5p
        self.trim_r2_3p = options.trim_r2_3p
        self.keep_r1_5p = options.keep_r1_5p
        self.keep_r1_3p = options.keep_r1_3p
        self.keep_r2_5p = options.keep_r2_5p
        self.keep_r2_3p = options.keep_r2_3p
        self.host_aligner = options.host_aligner
        self.read_counter = options.read_counter
        self.merge_dir = options.merge_dir
        self.use_p7 = options.use_p7
        self.use_p5 = options.use_p5
        self.use_lane = options.use_lane
        self.do_patho = options.do_patho
        self.do_host = options.do_host
        self.remove_splitted = options.remove_splitted
        self.gzip_merged = options.gzip_merged
        self.use_qsub = options.use_qsub
        self.do_gs    = options.do_gs
        self.do_split = options.do_split
        self.do_merge = options.do_merge
        self.do_align = options.do_align
        self.do_count = options.do_count
        self.do_metrics = options.do_metrics
        self.do_bc_dist = options.do_bc_dist
        self.min_resource = options.min_resource
        self.do_expand = options.do_expand
        self.AllSeq_con_A = options.AllSeq_con_A
        self.AllSeq_trim_len = options.AllSeq_trim_len
        self.no_bc_split = options.no_bc_split
        self.count_strand_rev = options.count_strand_rev
        self.use_seq_path = options.use_seq_path
        self.MOC_id = options.MOC_id
        self.MOC_id_ref = options.MOC_id_ref
        self.do_bestacc = options.do_bestacc
        self.do_umi_count = options.do_umi_count
        self.do_umi_metrics = options.do_umi_metrics
        self.use_sample_id = options.use_sample_id
        self.bwa_mem = options.bwa_mem
        self.do_ref = options.do_ref
        self.do_drop_split = options.do_drop_split
        self.do_drop_count = options.do_drop_count
        self.do_drop_metrics = options.do_drop_metrics
        self.rm_rts_dup = options.rm_rts_dup
        self.do_ref = options.do_ref
        self.do_rerun = options.do_rerun
        self.ucore_time = options.ucore_time
        self.paired_only_patho = options.paired_only_patho
        self.paired_only_host = options.paired_only_host


    def storeConfigMixed(self):
        options = self.options
        mydict = self.mydict
        ldelim = self.ldelim

        suffix_s1 = ''
        if options.Suffix_s1 != 'none':
            suffix_s1 = options.Suffix_s1
        else:
            suffix_s1 = self.Suffix_s1
        self.suffix_s1 = suffix_s1
        
        suffix_s2 = ''
        if options.Suffix_s2 != 'none':
            suffix_s2 = options.Suffix_s2
        else:
            suffix_s2 = self.Suffix_s2
        self.suffix_s2 = suffix_s2

        suffix_ne = ''
        if options.Suffix_ne != 'none':
            suffix_ne = options.Suffix_ne
        else:
            suffix_ne = self.Suffix_ne
        self.suffix_ne = suffix_ne

        use_dropseq = self.use_dropseq
        if self.use_sample_id:
            self.use_p7 = False
            self.no_bc_split = True
            if not use_dropseq:
                self.soft_link_from_input = True

        regex = '\s*[,|;]\s*'
        raw_seq_path = list()
        if options.seq_id != 'none':
            Seq_base = mydict['Seq_base']
            lseq_ids = re.split(regex, options.seq_id)
            for lseq_id in lseq_ids:
                raw_seq_path_l = Seq_base + self.ldelim + lseq_id
                raw_seq_path.append(raw_seq_path_l)
        elif options.raw_seq_path != 'none':
            raw_seq_path = re.split(regex, options.raw_seq_path)
        else:
            raise argparse.ArgumentTypeError('Either seq_id or raw_seq_path need to be provided.')
        print(raw_seq_path)
        self.Input_dir = raw_seq_path 

        if options.use_login_name:
            if options.login_name == 'Y':
                self.login_name = getpass.getuser()
            else:
                self.login_name = options.login_name

        temp_path = ''
        if options.temp_path != 'none':
            temp_path = options.temp_path
        else:
            temp_path = mydict['Temp_path']

        if options.use_login_name:
            self.Temp_dir = temp_path + self.ldelim + self.login_name
        else:
            self.Temp_dir = temp_path
        
        bam_path = ''
        if options.bam_path != 'none':
            bam_path = options.bam_path
        else:
            bam_path = mydict['Bam_path']

        if options.use_login_name:
            self.Bam_path = bam_path + self.ldelim + self.login_name
        else:
            self.Bam_path = bam_path

        results_path = ''
        if options.results_path != 'none':
            results_path = options.results_path 
        else:
            results_path = mydict['Results_path'] 
       
        if options.use_login_name:
            self.Results_path = results_path + self.ldelim + self.login_name
        else:
            self.Results_path = results_path

        # Add MOC _id hierarchy
        if self.MOC_id != 'none':
            self.Temp_dir = self.Temp_dir + self.ldelim + self.MOC_id
            self.Bam_path = self.Bam_path + self.ldelim + self.MOC_id
            self.Results_path = self.Results_path + self.ldelim + self.MOC_id

        if options.min_resource:
            self.host_thread_count = 1
            self.host_memory = 32


    def createSubPaths(self):
        Temp_dir = self.Temp_dir
        ldelim = self.ldelim
        use_qsub = self.use_qsub
        do_patho = self.do_patho
        do_host = self.do_host
        do_gs  = self.do_gs
        use_dropseq = self.use_dropseq
        do_rerun = self.do_rerun

        project_set = ''
        regex = '\s*[,|;]\s*'
        if self.Project_ids != 'none':
             project_set = set(re.split(regex, self.Project_ids))
        else:
            project_set = self.get_project_set(self.KeyTbl)
        self.project_set = sorted(project_set)
        out_dir_name = '_'.join(self.project_set)
        Out_dir = Temp_dir + ldelim + out_dir_name
        self.Out_dir = Out_dir
        
        UGER_cbp_dir = Out_dir + ldelim + "UGER_cbp"
        GS_dir = Out_dir + ldelim + "gsdir"
        Split_dir = Out_dir + ldelim + "splitdir"
        if self.merge_dir == 'none':
            Merge_dir = Out_dir + ldelim + "mergedir"
        else:
            Merge_dir = self.merge_dir
        Merge_real_dir = Out_dir + ldelim + "mergedir/real"
        Merge_garbage_dir = Out_dir + ldelim + "mergedir/garbage"
        Log_dir = Out_dir + ldelim + "logdir"
        Data_dir = Out_dir + ldelim + "datadir"
        Summary_dir = Out_dir + ldelim + "summary_dir"
        Host_dir = Out_dir + ldelim + "host_result"
        Patho_dir = Out_dir + ldelim + "patho_result"

        if use_qsub and os.path.exists(UGER_cbp_dir):
            if not do_rerun:
                shutil.rmtree(UGER_cbp_dir)
        if not os.path.exists(Out_dir):
            os.makedirs(Out_dir)
        #if do_gs and not os.path.exists(GS_dir):
        #    os.makedirs(GS_dir)
        if not os.path.exists(Split_dir):
            os.makedirs(Split_dir)
        if not os.path.exists(Merge_dir):
            os.makedirs(Merge_dir)
        if use_dropseq and not os.path.exists(Merge_real_dir):
            os.makedirs(Merge_real_dir)
        if use_dropseq and not os.path.exists(Merge_garbage_dir):
            os.makedirs(Merge_garbage_dir)
        if not os.path.exists(Log_dir):
            os.makedirs(Log_dir)
        if not os.path.exists(Data_dir):
            os.makedirs(Data_dir)
        if not os.path.exists(Summary_dir):
            os.makedirs(Summary_dir)
        if do_host and not os.path.exists(Host_dir):
            os.makedirs(Host_dir)
        if do_patho and not os.path.exists(Patho_dir):
            os.makedirs(Patho_dir)

        self.UGER_cbp_dir = UGER_cbp_dir
        self.Summary_dir = Summary_dir
        self.Patho_dir = Patho_dir
        self.Host_dir = Host_dir
        self.Data_dir = Data_dir
        self.Log_dir = Log_dir
        self.Split_dir = Split_dir
        self.Merge_dir = Merge_dir

        # Create the host_ref_path 
        if self.do_host:
            self.host_ref_path = self.host_dbpath + ldelim + \
                self.host_aligner + ldelim + self.host_ref_str
        self.prefix_dict_out = self.Log_dir + ldelim + "prefix_dict_out.txt"
       
    def initMain(self):
        options = self.options
        mydict = self.mydict
        key_path = ''
        if options.key_id != 'none':
            Key_base = mydict['Key_base']
            key_path = Key_base + self.ldelim + options.key_id + "_key.txt"
        elif key_path != 'none':
            key_path = options.key_path
        else:
            raise argparse.ArgumentTypeError('Either key_id or key_path need to be provided.')    
        self.Key_path = key_path

    def createUmiDirs(self, outdir_umi, logdir_umi):
        if not os.path.exists(outdir_umi):
            os.makedirs(outdir_umi)
        if not os.path.exists(logdir_umi):
            os.makedirs(logdir_umi)

    def createSubDirs(self):
        do_host = self.do_host
        do_patho = self.do_patho
        project_set = self.project_set
        Bam_path = self.Bam_path
        ldelim = self.ldelim
        use_dropseq = self.use_dropseq
        Summary_dir = self.Summary_dir
        ldel = ldelim

        do_umi_count = self.do_umi_count
        rm_rts_dup = self.rm_rts_dup
        Results_path = self.Results_path

        for project_id in project_set:
        
            result_dir = Results_path + ldel + project_id
            if not os.path.exists(result_dir):
                os.makedirs(result_dir)
    
            picard_outdir = result_dir + ldel + "bam_metrics_patho"
            if not os.path.exists(picard_outdir):
                os.makedirs(picard_outdir)
            
            picard_outdir_host = result_dir + ldel + "bam_metrics_host"
            if not os.path.exists(picard_outdir_host):
                os.makedirs(picard_outdir_host)

            if rm_rts_dup:
                no_dup_dir = result_dir + ldel + "nodupdir"
                no_dup_picard_outdir = no_dup_dir + ldel + "bam_metrics_patho"
                no_dup_picard_outdir_host = no_dup_dir + ldel + "bam_metrics_host"

                if not os.path.exists(no_dup_dir):
                    os.makedirs(no_dup_dir)

                if not os.path.exists(no_dup_picard_outdir):
                    os.makedirs(no_dup_picard_outdir)
            
                if not os.path.exists(no_dup_picard_outdir_host):
                    os.makedirs(no_dup_picard_outdir_host)

            if do_host:
                outdir = self.Host_dir + ldelim + project_id
                if not os.path.exists(outdir):
                    os.makedirs(outdir)

                if do_umi_count:
                    outdir_umi = self.Host_dir + ldelim + project_id + ldelim + "umidir"
                    logdir_umi = outdir_umi + ldelim + "logdir"
                    self.createUmiDirs(outdir_umi, logdir_umi)

                if rm_rts_dup:
                    no_dup_dir = outdir + ldelim + "nodupdir"
                    if not os.path.exists(no_dup_dir):
                        os.makedirs(no_dup_dir)

                    no_dup_picard_outdir = no_dup_dir + ldel + "bam_metrics_host"
                    if not os.path.exists(no_dup_picard_outdir):
                        os.makedirs(no_dup_picard_outdir)
            if do_patho:
                outdir = self.Patho_dir + ldelim + project_id
                summary_subdir = Summary_dir + ldelim + project_id
                temp_bamdir = outdir + ldelim + "temp_bamdir"

                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                if not os.path.exists(temp_bamdir):
                    os.makedirs(temp_bamdir)
                
                if do_umi_count:
                    outdir_umi = self.Patho_dir + ldelim + project_id + ldelim + "umidir"
                    logdir_umi = outdir_umi + ldelim + "logdir"
                    self.createUmiDirs(outdir_umi, logdir_umi)
                
                if use_dropseq:
                    Patho_real_dir = outdir + ldelim + "real"
                    Patho_garbage_dir = outdir + ldelim + "garbage"
                    Patho_bcdist_dir = outdir + ldelim + "bcdist"
                    Patho_drops_dir = outdir + ldelim + "drops"
                    Patho_drops_umidir = Patho_drops_dir + ldelim + "umidir"
                    summary_bcdist_dir = summary_subdir + ldelim + "bcdist"
                    if use_dropseq and not os.path.exists(Patho_real_dir):
                        os.makedirs(Patho_real_dir)
                    if use_dropseq and not os.path.exists(Patho_garbage_dir):
                        os.makedirs(Patho_garbage_dir)
                    if use_dropseq and not os.path.exists(Patho_bcdist_dir):
                        os.makedirs(Patho_bcdist_dir)
                    if use_dropseq and not os.path.exists(summary_bcdist_dir):
                        os.makedirs(summary_bcdist_dir)
                    if use_dropseq and not os.path.exists(Patho_drops_dir):
                        os.makedirs(Patho_drops_dir)
                    if use_dropseq and not os.path.exists(Patho_drops_umidir):
                        os.makedirs(Patho_drops_umidir)
               
                if rm_rts_dup:
                        no_dup_dir = outdir + ldelim + "nodupdir"
                        if not os.path.exists(no_dup_dir):
                            os.makedirs(no_dup_dir)
                        no_dup_temp_bamdir = no_dup_dir + ldelim + "temp_bamdir"
                        if not os.path.exists(no_dup_temp_bamdir):
                            os.makedirs(no_dup_temp_bamdir)

            bamdir = Bam_path + ldelim + project_id
            bamdir_real = Bam_path + ldelim + project_id + ldelim + "real"
            bamdir_garbage = Bam_path + ldelim + project_id + ldelim + "garbage"
            bamdir_nodupdir = Bam_path + ldelim + project_id + ldelim + "nodupdir"
            if not os.path.exists(bamdir):
                os.makedirs(bamdir)
            if use_dropseq and not os.path.exists(bamdir_real):
                os.makedirs(bamdir_real)
            if use_dropseq and not os.path.exists(bamdir_garbage):
                os.makedirs(bamdir_garbage)
            if do_umi_count:
                bamdir_umi = bamdir + ldelim + "umidir"
                if not os.path.exists(bamdir_umi):
                    os.makedirs(bamdir_umi)
            if rm_rts_dup and not os.path.exists(bamdir_nodupdir):
                os.makedirs(bamdir_nodupdir) 
                

    def doEndUpdate(self):
        if self.MOC_id_ref != 'none':
            ldel = self.ldelim
            self.patho_dbpath = self.patho_dbpath + ldel + self.MOC_id_ref
            
        
    def loadConfig(self):
        self.initMain()
        self.storeConfigFromConfig()
        self.storeConfigFromOptions()
        self.loadKeyTable()    
        self.storeConfigFromKeyTbl()
        self.storeDerivedPaths()
        self.storeConfigMixed()
        self.createSubPaths()
        self.createSubDirs()
        self.doEndUpdate()

    def dumpConfigLog(self):
        config_log_out = self.Log_dir + self.ldelim + "config_log_out.txt"
        with open(config_log_out, 'w') as outfile:
            yaml.dump(self.datastore, outfile, default_flow_style=False)
        return config_log_out

    def getKeyTblFinal(self):
        return self.KeyTbl

