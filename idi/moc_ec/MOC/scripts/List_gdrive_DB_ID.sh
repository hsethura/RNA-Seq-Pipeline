#!/bin/sh





# source /broad/software/scripts/useuse 2> /dev/null
# 
# use Python-2.7
# 
# source /idi/moc_ec/MOC/scripts/bash_header 2> /dev/null

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

### determining paths and headers 
paths_and_headers $ID $@

ALL_REFS_OPT=`extract_option -all_refs - 1 $@`
GREF=`extract_option -gref Y 1 $@`

### set paths to directories, files, and scripts from config file
### default config file is /idi/moc_ec/MOC/config_files/Universal_config.yaml

CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/Universal_config.yaml" 1 $@`

GPC_PATH=`config_read $CONFIG_FILE gdrivePC_path`
GPD_PATH=`config_read $CONFIG_FILE gdriveDEV_path`
GCID_PATH=`config_read $CONFIG_FILE gdriveGCIDTC_path`

GTRACK_PATH=`config_read $CONFIG_FILE gdrivetrack_path`
GPUNCH_PATH=`config_read $CONFIG_FILE gdrivepunch_path`
GLOCAL_PATH=`config_read $CONFIG_FILE gdrivelocal_path`
FILE_PATH=`config_read $CONFIG_FILE file_path`

GDRIVE_SCRIPT=`config_read $CONFIG_FILE gdrive_script`
TEMP_DIR=`config_read $CONFIG_FILE Temp_path`
GLOCAL_DIR=$GLOCAL_PATH"/"$GPD_PATH"/Combined_DBs/"
mkdir -p $GLOCAL_DIR

cd $GLOCAL_PATH
pwd 

echo $GPD_PATH
ALL_IDS=`$GDRIVE_SCRIPT file-id -depth 4 $GCID_PATH | grep "Expt_DB" | awk '{y=split($2, ar, "/");print $1";"ar[y]}' | sed 's/"//g'`

ALL_FILES=""
for ID in $ALL_IDS 
do
	GID=`echo $ID | cut -d";" -f1`
	FILE=`echo $ID | cut -d";" -f2`
	KEY_FILE=$TEMP_DIR$FILE"_key.txt"
	ALL_FILES=$ALL_FILES" "$KEY_FILE
	
	KEY_SHEET="Sheet1"
	echo "Moving key file to server..."
	echo "$KEY_SCRIPT -s $GID -t \"$KEY_SHEET\" -p $FILE --Key_dir $TEMP_DIR"

	$KEY_SCRIPT -s $GID -t "$KEY_SHEET" -p $FILE --Key_dir $TEMP_DIR
done

echo $ALL_FILES

COMB_DB=$GLOCAL_DIR"MOC_dev_combined_DB.tsv"

cat $ALL_FILES | grep -v "###"| head -1 > $COMB_DB
cat $ALL_FILES | grep -v "Status" | grep -v "###" >> $COMB_DB

ls -lrt $COMB_DB
cd $GLOCAL_DIR
pwd
ls -lrt 
echo $GLOCAL_DIR
$GDRIVE_SCRIPT push -convert -ignore-name-clashes -quiet

# cd $GLOCAL_PATH
# $GDRIVE_SCRIPT url MOC_-_Microbial_Omics_Core/IDMP_PC/Development//Combined_DBs/MOC_dev_combined_DB.tsv
