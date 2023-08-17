#!/bin/sh

MOC_ID="none"

source /broad/software/scripts/useuse

use Python-2.7

source /idi/moc_ec/MOC/scripts/bash_header

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"


CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/PC_config.yaml" 1 $@`
COMBINE=`extract_option -combine Y 1 $@`
LIMIT=`extract_option -limit MOCS 1 $@`
LIMIT_MOCS=`echo $LIMIT | sed 's/\-//g'`
OW_SYM=`extract_option -ow_sym N 1 $@`

### determining paths and headers 
### default config file is /idi/moc_ec/MOC/config_files/Universal_config.yaml
paths_and_headers $MOC_ID $@


### set paths to directories, files, and scripts from config file

TEMP_PATH=`config_read $CONFIG_FILE Temp_path`
GLOCAL_PATH=`config_read $CONFIG_FILE gdrivelocal_path`
GDRIVE_SCRIPT=`config_read $CONFIG_FILE gdrive_script`
RAWDATA_PATH=`config_read $CONFIG_FILE Seq_base`
TEMP_DIR=$TEMP_PATH"/fungal_move/"
mkdir -p $TEMP_DIR


DATA_DIR=`extract_option -data_dir $PARSE_DIR 1 $@`
PARSE=`extract_option -parse Y 1 $@`
IMPORT_GS=`extract_option -import_gs Y 1 $@`

GMOCF_PATH="MOC/GCID_TC/"
GWF_DIR=$GMOCF_PATH"Fungal_projects-Metrics/"
LOCAL_DIR=$GLOCAL_PATH"/"$GWF_DIR"/"

mkdir -p $LOCAL_DIR

echo "cd $GLOCAL_PATH"
cd $GLOCAL_PATH
pwd 


### pull desktop files for all MOCS submission wbs from $GWB_DIR
echo "Pulling all files from "$GWF_DIR
echo "$GDRIVE_SCRIPT pull -no-prompt -files $GWF_DIR"
$GDRIVE_SCRIPT pull -no-prompt -files $GWF_DIR

ls -lrt $LOCAL_DIR
echo $LOCAL_DIR

### pull MOCS submission wbs from $GWB_DIR
cd $LOCAL_DIR
pwd 


echo "" | sed 1d > ~/test.txt 

ALL_FILES=`grep -H $LOCAL_DIR*/* -e URL= -e URL= > ~/test.txt`


cat ~/test.txt

 

while read line 
do
	
	echo $line
	
	
	IN_FILE_NAME=`echo $line | awk -F":" '{y=split($1, ar, "/");print ar[y]}'`
	GID=`echo $line | awk -F":" '{print $3}' | cut -d"/" -f6`
	OUT_FILE_NAME=`echo $IN_FILE_NAME | sed -e 's/.desktop//g' -e 's/ /_/g'`
	echo "DFFD"
	echo $IN_FILE_NAME
	echo $OUT_FILE_NAME
	echo $GID
	
	
	OUT_FILE=$TEMP_DIR$OUT_FILE_NAME".txt"
	
	echo $NAME" "$GID
	wget --no-check-certificate --output-document=$OUT_FILE 'https://docs.google.com/spreadsheet/ccc?key='$GID'&output=txt'
echo wget --no-check-certificate --output-document=$OUT_FILE 'https://docs.google.com/spreadsheet/ccc?key='$GID'&output=txt'


# 	echo "$SCRIPTS_DIR"GS_import.py" -s $GID -t "Agg" -p $NAME --Key_dir $TEMP_DIR -S ".txt""
# 	$SCRIPTS_DIR"GS_import.py" -s $GID -t "Agg" -p $NAME --Key_dir $TEMP_DIR -S ".txt" 

done < ~/test.txt

cat $TEMP_DIR/*etrics* | sort | uniq | sort > $GLOCAL_PATH/BIG_DB.tsv

cd $LOCAL_DIR
pwd

$GDRIVE_SCRIPT push -no-prompt -destination $GWF_DIR $GLOCAL_PATH/BIG_DB.tsv


ls -lrt $TEMP_DIR

### change permissions for Results and temp dirs
change_perms $LOCAL_DIR 
