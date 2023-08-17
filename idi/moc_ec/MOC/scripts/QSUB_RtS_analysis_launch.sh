#!/bin/sh

MOC_ID=$1

source /idi/moc_ec/MOC/scripts/bash_header

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

### get options from command line
CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/PC_config.yaml" 1 $@`

# Get paths to dirs scripts from config file
read_config $CONFIG_FILE 

### run path_suff function to set RESPATH_SUFF.
##	If -moc_id N included in command line, do not moc_ID to RESPATH_SUFF.  
##	If -user_id included with Y or no string add USID to RESPATH_SUFF.  
##	If -user_id followed by userID, add userID to RESPATH_SUFF.

path_suff $@

echo $RtS_ANPIPE

RESULTS_DIR=$RESULTS_PATH"/"$RESPATH_SUFF"/UGER/"

echo "Results_dir: "$RESULTS_DIR

mkdir -p $RESULTS_DIR

TIME_STAMP=`date +%s` 

QSUB_ERR_FILE=$RESULTS_DIR"/"$MOC_ID"_"$TIME_STAMP"_qsub_err.txt"
QSUB_OUT_FILE=$RESULTS_DIR"/"$MOC_ID"_"$TIME_STAMP"_qsub_out.txt"

QSUB_FILE=~/$MOC_ID"_qsub.txt"

echo "source /idi/moc_ec/MOC/scripts/bash_header" > $QSUB_FILE
echo "sh $RtS_ANPIPE $@" >> $QSUB_FILE

echo "qsub -e $QSUB_ERR_FILE -o $QSUB_OUT_FILE -l h_rt=24:00:00 -l h_vmem=8g -l os=RedHat7 $QSUB_FILE" 
qsub -e $QSUB_ERR_FILE -o $QSUB_OUT_FILE -l h_rt=24:00:00 -l h_vmem=8g -l os=RedHat7 $QSUB_FILE

echo "Outfile and errfile directory is:"
echo $RESULTS_DIR