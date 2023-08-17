#!/usr/bin/env python

import argparse
import csv
import codecs
import pandas as pd
import re
import os
import subprocess as sp
from subprocess import call
import ntpath

class revert_script:
    def __init__(self, options):
       self.key_path = options.key_path
       self.outdir = options.outdir
       self.use_qsub = options.use_qsub
       self.sample_id = 'Sample_ID'
       self.Path_to_SeqFile = 'Path_to_SeqFile'      
       self.UGER_cbp = "/broad/IDP-Dx_work/nirmalya/tools/ugetools/UGE_SUBMISSIONS/UGER_cmd_batch_processor.py" 
       self.basepath = os.path.dirname(os.path.realpath(__file__))

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

    def loadKeyTable(self):
        Key_path = self.key_path
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
        return tab3

    def getSampleidToSeqfile(self, KeyTbl):
        sample_id = self.sample_id
        Path_to_SeqFile = self.Path_to_SeqFile

        sample_id_lst = KeyTbl[sample_id].tolist()
        seqfile_lst = KeyTbl[Path_to_SeqFile].tolist()
        sample_seqfile = dict(zip(sample_id_lst, seqfile_lst))
        return sample_seqfile


    def is_sorted_qname(self, seqFile):
        return False
        """
        print("Checking if seqFile is sorted: " + seqFile)
        cmd1 = '/broad/IDP-Dx_work/nirmalya/local/bin/samtools view -H ' + seqFile 
        print("from is_sorted_qname: " + cmd1)
        cmd2 = 'grep @HD'
        p1 = sp.Popen(cmd1.split(), stdout = sp.PIPE, shell = True)
        p2 = sp.Popen(cmd2.split(), stdin = p1.stdout, stdout = sp.PIPE, shell = True)
        output,err = p2.communicate()
        qn_sort_str = 'SO:queryname'
        if qn_sort_str in output:
            return True
        else:
            return False
        """

    def get_sorted_seqFile(self, lsample, seqFile, tempdir):
        lfilename = ntpath.basename(seqFile)
        lfilename_u = lsample + "_u.bam"
        lpath_u = tempdir + "/" + lfilename_u
        sort_cmd = 'samtools sort -n -o ' + lpath_u + " " + seqFile
        return lpath_u, sort_cmd
  
    def exe_bam_revert(self, sampleid_to_seqFile):

        """
        sh ALN_BAM_to_FASTQ.sh <bam file> file1.fastq file2.fastq 1 1 <tempdir>
        """
        outdir = self.outdir
        fastq_dir = outdir + "/fastq"
        tempdir = outdir + "/tempdir"
        Log_dir = outdir + "/logdir"
        UGER_cbp_dir = outdir + "/UGER_cbp_dir"
        UGER_cbp = self.UGER_cbp
        use_qsub = self.use_qsub

        if not os.path.exists(outdir):
            os.makedirs(outdir)
        if not os.path.exists(fastq_dir):
            os.makedirs(fastq_dir)
        if not os.path.exists(tempdir):
            os.makedirs(tempdir)
        if not os.path.exists(Log_dir):
            os.makedirs(Log_dir)

        revert_script = "samtools fastq"

        sort_joblist_path = Log_dir + "/" + "bam_sort_joblist.txt"
        s_jfile = open(sort_joblist_path, "w")
        rev_joblist_path = Log_dir + "/" + "bam_revert_joblist.txt"
        r_jfile = open(rev_joblist_path, "w")

        for lsample in sampleid_to_seqFile:
            seqFile = sampleid_to_seqFile[lsample]
            if seqFile:
                loutfile = Log_dir + "/" + lsample + "_out.txt"
                lerrfile = Log_dir + "/" + lsample + "_err.txt"
                
                if self.is_sorted_qname(seqFile):
                    seqFile_s = seqFile
                else:
                    lseqFile, lsort_cmd = self.get_sorted_seqFile(lsample, seqFile, tempdir)
                    seqFile_s = lseqFile
                    lsort_cmd2 = lsort_cmd + " 1> " + loutfile + " 2> " + lerrfile + "\n"
                    print("sort cmd: " + lsort_cmd2)
                    s_jfile.write(lsort_cmd2)
     
                lfile1 = fastq_dir + "/" + lsample + "_R1.fastq.gz"
                lfile2 = fastq_dir + "/" + lsample + "_R2.fastq.gz"
                cmd_str = revert_script + " -1 " + lfile1 + " -2 " + lfile2 + " " + seqFile_s

                cmd_str2 = cmd_str + " 1> " + loutfile + " 2> " + lerrfile + "\n"
                r_jfile.write(cmd_str2)
                print("revert cmd: " + cmd_str)
        s_jfile.close()
        r_jfile.close()

        lmemory = "16"
        sort_joblist_cmd = UGER_cbp + " --cmds_file " + sort_joblist_path + \
                                " --batch_size 1" + \
                                " --num_cores 1" + \
                                " --memory " + lmemory + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header"
        if use_qsub:
            call(sort_joblist_cmd.split())

        rev_joblist_cmd = UGER_cbp + " --cmds_file " + rev_joblist_path + \
                                " --batch_size 1" + \
                                " --num_cores 1" + \
                                " --memory " + lmemory + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header"
        if use_qsub:
            call(rev_joblist_cmd.split())



    def mainFunc(self):
        KeyTbl = self.loadKeyTable()
        sampleid_to_seqfile = self.getSampleidToSeqfile(KeyTbl)
        #print(sampleid_to_seqfile)
        self.exe_bam_revert(sampleid_to_seqfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process the options.')
    parser.add_argument('--key_path', dest = 'key_path', required = True, help = 'Key file path (absolute)')
    parser.add_argument('--outdir', dest = 'outdir', required = True, help = 'Outdir path (absolute)')
    parser.add_argument('--no_qsub', dest = 'use_qsub', action = 'store_false', default = True, help = 'Does not submit qsub jobs.' )

    options = parser.parse_args()

    reverts = revert_script(options)
    reverts.mainFunc()  
