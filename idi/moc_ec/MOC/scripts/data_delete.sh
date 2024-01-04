#!/bin/sh

days=$1

# get path of the current file. if the file path is relative, convert it to absolute path
file_path="${BASH_SOURCE[0]}"
if [[ $file_path != /* ]]; then
  file_path="$PWD/${BASH_SOURCE[0]}"
fi

PROJECT_ROOT_DIR="$(dirname $(dirname $(dirname $(dirname $(dirname $file_path)))))"

# get parent directory
scripts_dir="$(dirname $file_path)"

# source /idi/moc_ec/MOC/scripts/bash_header
source "$scripts_dir/bash_header"

### source all functions 
# source "idi/moc_ec/MOC/scripts/MOC_functions.sh"
source "$scripts_dir/MOC_functions.sh"

### determining paths and headers 
### default config file is /broad/IDP-Dx_storage/MOC/config_files/PC_config.yaml
paths_and_headers $MOC_ID $@

USID=`USID`
MAIL_USID=`USID`

### set options
# -conf: sets path to config file (default /idi/moc_ec/MOC/config_files/PC_config.yaml)
# -MOVE_KEY: including or setting to Y skips moving google sheet to server (default N)

# CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/PC_config.yaml" 1 $@`
DEFAULT_CONFIG_PATH="$(dirname $(dirname $file_path))"/config_files/PC_config.yaml
CONFIG_FILE=`extract_option -conf $DEFAULT_CONFIG_PATH 1 $@`
if [[ $CONFIG_FILE != /* ]]; then
  CONFIG_FILE="$PROJECT_ROOT_DIR/$CONFIG_FILE"
fi
 
############## Import google sheet databases to server ###############=
if [ $IMPORT_GS == "Y" ];then
	echo "Running $GSIMPORT_SCRIPT to import google sheet databases to server"
	sh $GSIMPORT_SCRIPT
fi
########################################################

PC_DB=`ls -lrt $PCDB_DIR* | awk '{print $9}' | tail -1`

TEMP_OUTFILE=$TEMP_PATH"/data_delete_temp.txt"
OUTFILE=$TEMP_PATH"/data_delete.txt"

time_id $DATA_DB $days $TEMP_OUTFILE "Warning"

ALL_MOCPIDS=`cat $TEMP_OUTFILE | awk '{print $1}' | sort | uniq`

echo "Project_status	MOCP_ID	fastq_path	Days_since_last_warning" > $OUTFILE
echo $TEMP_OUTFILE


for MOCPID in $ALL_MOCPIDS
do
	PROJ_STATUS=`Index $PC_DB "MOCP_ID" "Project_Status" $MOCPID`
	cat $TEMP_OUTFILE | grep $MOCPID | awk -v PROJ_STATUS=$PROJ_STATUS '{print  PROJ_STATUS, $0}' >> $OUTFILE
	echo "" >> $OUTFILE
done 

ls -lrt $OUTFILE


exit



