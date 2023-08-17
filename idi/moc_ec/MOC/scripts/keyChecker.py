#!/usr/bin/env/python

'''
This takes in a MOC_ID and examines it's keyfile to see if there are any incomplete fields.

One Positional argument:
$1 = MOC_ID
'''

import pandas as pd
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('MOC_ID', help= 'the MOC_ID whose key will be assessed')
args= parser.parse_args()



KEY_DIR = "/broad/IDP-Dx_storage/MOC/Key_files/"
SAMPLE_KEY = KEY_DIR+"/"+args.MOC_ID+"_key.txt"


sampleKey = pd.read_csv(SAMPLE_KEY, sep="\t", skiprows=4, engine='python')
sampleKey = sampleKey[sampleKey["Unnamed: 0"] != "###"]
column_headers= sampleKey.columns


keyColNames = ["MOC_ID", "Sample_ID", "Read_pairing", "Requested_Seq_Platform", 
               "Total_cycles", "Desired_depth", "Bacterial_reference", "RiboZero_Request", "LC_method", "Inline_Name", 
              "Inline_Seq", "Pool_ID", "Index1_Name", "Index1_Seq", "Index2_Name", "Index2_Seq", "#_cycles_enrichment",
              "Seq_ID"]
              
              
for col in keyColNames:
    if sampleKey[col].hasnans is True:
        print("This Key is incomplete.  See:", col)
    
print("keychecker.py is done")        
        
        