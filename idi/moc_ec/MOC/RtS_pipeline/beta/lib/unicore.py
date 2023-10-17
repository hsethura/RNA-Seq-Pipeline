#!/usr/bin/env python

import argparse
import yaml

from rtscore import RtSCore
from dictmap import DictMap
from allseqcore import AllSeqCore
from alignr2core import AlignR2Core
from alignscore import AlignSCore

# This would execute the core part for RtS (RNATag-Seq) pipeline.

class UniCore:
    def __init__(self, args):
        self.config_log_file = args.config_log_file
        
        self.sample_id = args.sample_id
        self.project_id = args.project_id
        #self.prefix_set = args.prefix_set
        #self.bc_set = args.bc_set
        self.ref_acc_str = args.ref_acc_str
        self.host_ref_str = args.host_ref_str
        cldict_d = yaml.load(open(self.config_log_file))
        cldict = DictMap(cldict_d)
        self.cldict = cldict

   

    def mainFunc(self):

        sampd = dict()
        sampd['sample_id'] = self.sample_id
        sampd['project_id'] = self.project_id
        #sampd['prefix_set'] = self.prefix_set
        #sampd['bc_set'] = self.bc_set
        sampd['ref_acc_str'] = self.ref_acc_str
        sampd['host_ref_str'] = self.host_ref_str
        sampd_map = DictMap(sampd) 
        cldict = self.cldict
        LC_method_val = cldict.LC_method_val
        lc_lower = LC_method_val.lower()
        print("lc_lower: " + lc_lower)
        if lc_lower == 'rts' or lc_lower == 'rts-ts':
            rtso = RtSCore(cldict, sampd_map)
            rtso.mainFunc()
        elif lc_lower == 'scr':
            rtso = RtSCore(cldict, sampd_map)
            rtso.mainFunc()
        elif lc_lower == 'allseq':
            allseqo = AllSeqCore(cldict, sampd_map)
            allseqo.mainFunc()
        elif lc_lower == 'alignr1':
            alignr1o = AlignSCore(cldict, sampd_map)
            alignr1o.mainFunc()
        elif lc_lower == 'alignr2':
            alignr2o = AlignSCore(cldict, sampd_map)
            alignr2o.mainFunc()
        elif lc_lower == 'bdropseq' or lc_lower == 'indrops':
            alignso = AlignSCore(cldict, sampd_map)
            alignso.mainFunc()
        elif lc_lower == 'smarter':
            print("Unicore for smarter protocol.")
            rtso = RtSCore(cldict, sampd_map)
            rtso.mainFunc()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Process command for a single sample which is essentially a row in the key file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--config_log_file", dest = "config_log_file", type = str, required = True, help = "Log file containig the config vars.")
    parser.add_argument("--sample_id", dest="sample_id", type=str, required=True, help="Id of the sample to be processed.")
    parser.add_argument("--project_id", dest = "project_id", required = True, type = str, help = "project id for this sample")
    #parser.add_argument("--prefix_set", dest="prefix_set", type=str, required=True, help="Set of file prefixes for this sample")
    #parser.add_argument("--bc_set", dest="bc_set", type=str, required=True, help = "Set of the barcodes used for this samples")
    parser.add_argument("--ref_acc_str", dest = "ref_acc_str", type = str, default = "none", help = "Contains the host accession ids")
    parser.add_argument("--host_ref_str", dest ="host_ref_str", default = "none", help = "Link to the host reference")

    args = parser.parse_args()
 
    unio = UniCore(args)
    print("About to start uino mainFunc")
    unio.mainFunc()
 
