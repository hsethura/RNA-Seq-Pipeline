#!/usr/bin/env python
from subprocess import call
maindir = '/broad/hptmp/RNASeq_proj/viktoria_data/post'
joblist_path = maindir + '/joblist.txt'
UGER_cbp_dir = maindir + '/UGER_cbp_dir'
UGER_cbp = '/broad/IDP-Dx_work/nirmalya/tools/ugetools/UGE_SUBMISSIONS/UGER_cmd_batch_processor.py'

joblist_cmd = UGER_cbp + " --cmds_file " + joblist_path + \
                                " --batch_size 1" + \
                                " --queue long" + \
                                " --memory 8" + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header"
call(joblist_cmd.split())

