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

from collections import defaultdict

class Expander:

    def __init__(self, confd):
        self.confd = confd

    def process_single_expand(self, lsample_id, ex_set, row, ex_type):
        confd = self.confd
        sample_id = confd.sample_id
        ex_id = None
        if ex_type == "p7":
            ex_id = confd.P7_index
        elif ex_type == "bc":
            ex_id = confd.lbc
        row_header = list(row.index)
        l_tbl = pd.DataFrame(columns = row_header)
        for lex in ex_set:
            ex_sample = lsample_id + "_" + ex_type + ":" + lex
            l_tbl.loc[ex_sample] = row
            l_tbl.loc[ex_sample, sample_id] = ex_sample
            l_tbl.loc[ex_sample, ex_id] = lex
        return l_tbl


    def process_double_expand(self, lsample_id, p7_set, bc_set, row):
        confd = self.confd
        sample_id = confd.sample_id
        P7_index = confd.P7_index
        lbc = confd.lbc
        row_header = list[row.index]
        l_tbl = pd.DataFrame(columns = row_header)
        for p7_single in p7_set:
            for bc_single in bc_set:
                ex_sample = lsample_id + "_p7:" + p7_single + "_bc:" + bc_single
                l_tbl.loc[ex_sample] = row
                l_tbl.loc[ex_sample, sample_id] = lsample_id
                l_tbl.loc[ex_sample, P7_index] = p7_single
                l_tbl.loc[ex_sample, lbc] = bc_single
        return l_tbl


    def getExpSubTbl(self, lsample_id, lp7_index, lbc_str, row):
        regex = '\s*[;|,]\s*'
        p7_lst  = re.split(regex, lp7_index)
        bc_lst = re.split(regex, lbc_str)

        if len(p7_lst) <= 1 and len(bc_lst) <= 1:
            return None
        elif len(p7_lst) > 1 and len(bc_lst) > 1:
            return self.process_double_expand(lsample_id, p7_lst, bc_lst, row)
        elif len(p7_lst) > 1:
            return self.process_single_expand(lsample_id, p7_lst, row, "p7")
        elif len(bc_lst) > 1:
            return self.process_single_expand(lsample_id, bc_lst, row, "bc")


    def expandKeyTbl(self, KeyTbl):
        confd = self.confd
        KeyTbl_out = None
        ltab = KeyTbl
        project_id =   confd.project_id
        lproject_set = confd.project_set
        sample_id = confd.sample_id
        P7_index =     confd.P7_index
        lbc = confd.lbc
        for index, row in ltab.iterrows():
            lproject_id = row[project_id]
            if lproject_id in lproject_set:
                lsample_id = row[sample_id]
                lp7_index = row[P7_index]
                lbc_str = row[lbc]
                expSubTbl = self.getExpSubTbl(lsample_id, lp7_index, lbc_str, \
                    row)
                if expSubTbl is not None:
                    if KeyTbl_out is None:
                        KeyTbl_out = KeyTbl.copy()
                        print("Copied the key table for expansion")
                    for ex_index, ex_row in expSubTbl.iterrows():
                        KeyTbl_out.loc[ex_index] = ex_row
        if KeyTbl_out is None:
            KeyTbl_out = KeyTbl
        return KeyTbl_out
