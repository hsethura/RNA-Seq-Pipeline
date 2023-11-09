#!/bin/bash

source /broad/software/scripts/useuse

reuse UGER 
reuse Anaconda3
reuse BLAST+

pip install -U kaleido



if [ -z "\$SGE_TASK_ID" ]; then
    echo "ERROR, env var SGE_TASK_ID is not set"
    exit 1
fi

echo "HOSTNAME ${HOSTNAME}, SGE_TASK_ID: ${SGE_TASK_ID}" > /home/unix/hsethura/RNA-Seq-Pipeline/Y/adapter_analysis/tmp.tracking/cmds/${SGE_TASK_ID}.cmd.err

cmd="/home/unix/hsethura/RNA-Seq-Pipeline/Y/adapter_analysis/tmp.tracking/cmds/${SGE_TASK_ID}.cmd 1>/home/unix/hsethura/RNA-Seq-Pipeline/Y/adapter_analysis/tmp.tracking/cmds/${SGE_TASK_ID}.cmd.out 2>>/home/unix/hsethura/RNA-Seq-Pipeline/Y/adapter_analysis/tmp.tracking/cmds/${SGE_TASK_ID}.cmd.err"

eval $cmd

exit 0
    