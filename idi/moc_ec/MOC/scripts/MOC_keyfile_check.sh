#!/bin/sh

MOC_ID=$1

### source all functions 
source "/broad/IDP-Dx_storage/MOC/scripts/MOC_functions.sh"


CONFIG_FILE=`extract_option -conf "/broad/IDP-Dx_storage/MOC/config_files/Universal_config.yaml" 1 $SCRIPT_OPTIONS`
Q_HEAD=`extract_option -q MOC_ID 1 $SCRIPT_OPTIONS`

### set paths to directories, files, and scripts from config file

KEY_DIR=`config_read $CONFIG_FILE Key_base`
RESULTS_PATH=`config_read $CONFIG_FILE Results_path`
TEMP_PATH=`config_read $CONFIG_FILE Temp_path`
SEQ_PATH=`config_read $CONFIG_FILE Seq_base`
I1_SPLIT_PATH=`config_read $CONFIG_FILE Index1_split_path`
DB_SCRIPT=`config_read $CONFIG_FILE DB_script`
KEY_SCRIPT=`config_read $CONFIG_FILE keyfile_importer`
BACT_REF_PATH=`config_read $CONFIG_FILE Bacterial_Ref_path`
BACT_REF_HEADER=`config_read $CONFIG_FILE Ref_accession`

KEY_FILE=$KEY_DIR$MOC_ID"_key.txt"
BACT_REF_FIELD=`FIELD_HEADER $KEY_FILE $BACT_REF_HEADER`

G_ID=`sh $DB_SCRIPT $Q_HEAD,$MOC_ID Google_ID | sed 1d | awk '{print $2}'`

if [ -z $G_ID ];then
	echo "ERROR: No GID found for" $MOC_ID
	exit
fi
echo $G_ID

echo "Moving key file to server..."

echo $KEY_SCRIPT -s $G_ID -t "Sample Information" -p $MOC_ID --Key_dir $KEY_DIR
$KEY_SCRIPT -s $G_ID -t "Sample Information" -p $MOC_ID --Key_dir $KEY_DIR

if [ ! -s $KEY_FILE ];then
	echo "ERROR: No key file moved for" $MOC_ID
	ls -lrt $KEY_FILE
	exit
fi

echo "Checking headers in key_file"

ALL_HEADERS=`cat $CONFIG_FILE | awk -v flag=N '{if(flag=="Y" && $0 !~ /KEY_HEADERS_END/) print $2;if($0 ~ /KEY_HEADERS_START/)flag="Y";if($0 ~ /KEY_HEADERS_END/) exit}'`

for HEADER in $ALL_HEADERS
do
	if [ `cat $KEY_FILE | grep -v "#" | head -1 | grep -w $HEADER | wc -l` == 0 ];then
	 	echo "ERROR: "$HEADER" header is missing!"
	fi 

done

ALL_BACT_REFS=`cat $KEY_FILE  | grep -v "#" | sed 1d | awk -F"\t" -v BACT_REF_FIELD=$BACT_REF_FIELD '{print $BACT_REF_FIELD}' | sort | uniq`

echo $ALL_BACT_REFS


for BACT_REFS in $ALL_BACT_REFS
do
	ls -lrt  $BACT_REF_PATH$BACT_REFS".gff"
	ls -lrt  $BACT_REF_PATH$BACT_REFS".fna"
	if [ ! -s $BACT_REF_PATH$BACT_REFS".gff" ];then
		echo "ERROR: "$BACT_REF_PATH$BACT_REFS".gff is missing!"
	fi
	if [ ! -s $BACT_REF_PATH$BACT_REFS".fna" ];then
		echo "ERROR: "$BACT_REF_PATH$BACT_REFS".fna is missing!"
	fi

done



exit

