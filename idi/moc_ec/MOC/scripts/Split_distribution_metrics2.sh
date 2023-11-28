#!/bin/sh

MOC_ID=$1

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

Q_HEAD="MOC_ID"

### source all functions 
# source "idi/moc_ec/MOC/scripts/MOC_functions.sh"
source "$scripts_dir/MOC_functions.sh"


### run path_suff function to set RESPATH_SUFF.
##	If -moc_id N included in command line, do not moc_ID to RESPATH_SUFF.  
##	If -user_id included with Y or no string add USID to RESPATH_SUFF.  
##	If -user_id followed by userID, add userID to RESPATH_SUFF.
path_suff $@

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

MOVE_KEY=`extract_option -move_key Y 1 $@`
USID=`extract_option -uid $USID 1 $@`
DUAL=`extract_option -dual Y 1 $@`

IN_DIR=$TEMP_PATH"/"$RESPATH_SUFF"/"
OUT_DIR=$TEMP_PATH"/"$RESPATH_SUFF"/inline_distribution/"

echo $IN_DIR

mkdir -p $OUT_DIR

# mkdir -p $BCSDB_DIR
# echo "$SCRIPTS_DIR"GS_import.py" -s $BCS_GID -t "NEW barcodes and indexes" -p $BCS_FILE_NAME --Key_dir $BCSDB_DIR"
# $SCRIPTS_DIR"GS_import.py" -s $BCS_GID -t "NEW barcodes and indexes" -p $BCS_FILE_NAME --Key_dir $BCSDB_DIR -S $BCS_FILE_SUFF

echo "$BC_FILE"
cat_bc_file=`cat $BC_FILE`
echo "$cat_bc_file"

ALL_INLINES=`cat $BC_FILE | awk '{printf "%s;%s ", $2,$1}'`

echo $ALL_INLINES


# find gid move keyfile to server

KEY_DIR=`config_read $CONFIG_FILE Key_base`
KEY_FILE=$KEY_DIR$MOC_ID"_key.txt"

echo "CONFIG_FILE: $CONFIG_FILE"
echo "KEY_DIR: $KEY_DIR"
### move key file to server
if [ $MOVE_KEY == "Y" ];then

	PROJ_TYPE=`project_type $MOC_ID`
	echo "Proj_type: "$PROJ_TYPE


	if [ $GID_OPT == "0" ];then
		echo "Getting GID from Gdrive using the path in $CONFIG_FILE"
		echo "gdrive_gid $MOC_ID -conf $CONFIG_FILE -proj_type $PROJ_TYPE | grep -v Dropping | grep -v Prepending | tail -1"
		G_ID=`gdrive_gid $MOC_ID -conf $CONFIG_FILE -proj_type $PROJ_TYPE | grep -v Dropping | grep -v Prepending | tail -1`
	else
		G_ID=$GID_OPT
	fi


	if [ -z $G_ID ];then

		echo "G_ID not found in gdrive or entered in command line with -gid option"
		exit
	fi

	echo "GID:" $G_ID 
	echo "Moving key file to server..."
	
	echo $KEY_SCRIPT -s $G_ID -t "Sample Information" -p $MOC_ID --Key_dir $KEY_DIR
	$KEY_SCRIPT -s $G_ID -t "Sample Information" -p $MOC_ID --Key_dir $KEY_DIR
fi




### get all project IDs

PROJID_F=`FIELD_HEADER $KEY_FILE Project_ID`
SAMPID_F=`FIELD_HEADER $KEY_FILE Sample_ID`
POOLID_F=`FIELD_HEADER $KEY_FILE Pool_ID`
INLINE_F=`FIELD_HEADER $KEY_FILE Inline_Seq`

ALL_PROJID=`FIND_VALUES_FOR_HEADER Project_ID $KEY_FILE | sort -u`
ALL_POOLID=`FIND_VALUES_FOR_HEADER Pool_ID $KEY_FILE | sort -u`


echo "ALL_POOLID:	"$ALL_POOLID
echo "ALL_PROJID:	"$ALL_PROJID


for PROJID in $ALL_PROJID
do
	TEMP_OUTFILE=$TEMP_PATH"/"$RESPATH_SUFF"/"$PROJID"_inline_distribution_temp.txt"
	OUTFILE1=$OUT_DIR"/"$PROJID"_inline_distribution.txt"

	if [ ! -s $SPLIT_DIR ];then
		echo $SPLIT_DIR" not found!"
	fi	
		 
	echo "" | sed 1d > $TEMP_OUTFILE 

	for INLINE_PAIR in $ALL_INLINES
	do
		INLINE_SEQ=`echo $INLINE_PAIR | cut -d";" -f1`
		INLINE_NAME=`echo $INLINE_PAIR | cut -d";" -f2`

		for POOL_ID in $ALL_POOLID
		do
			SAMPID=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v PROJID=$PROJID -v INLINE_SEQ=$INLINE_SEQ -v POOL_ID=$POOL_ID '{if($'$PROJID_F'==PROJID && $'$INLINE_F'==INLINE_SEQ && $'$POOLID_F'==POOL_ID) print $'$SAMPID_F'}' | sed 's/ /_/g'`
			#echo $POOL_ID, $SAMPID, $INLINE_SEQ
						
			if [ -z $SAMPID ];then
				SAMPID="-"
			fi
			if [ -s $SPLIT_DIR"/"*$POOL_ID"_"*combined.txt ];then
				data=`cat $SPLIT_DIR"/"*$POOL_ID"_"*combined.txt | grep $INLINE_SEQ | awk '{print $2, $3, $4}'` 
			else 
				if [ $SAMPID != "-" ];then
					data="DATA_MISSING!"
				else
					data="- - - -"
				fi
			fi
			echo $SAMPID $INLINE_SEQ $POOL_ID $data	$INLINE_NAME 
		done
	done >> $TEMP_OUTFILE
	
	
	echo $TEMP_OUTFILE
	
		
	MISSING=`cat $TEMP_OUTFILE | grep "DATA_MISSING" | wc -l`

	echo "Data missing for "$MISSING" samples" > $OUTFILE1
	echo "" >> $OUTFILE1
	echo "Sample_ID	Inline_Barcode_Seq	Index	%_of_reads_in_pool	total_reads %_no_mismatch	Inline_Barcode_Name"  >> $OUTFILE1
	cat $TEMP_OUTFILE | awk '{if($1 != "-") print $0}' | grep -v "DATA_MISSING" | sort -k4nr >> $OUTFILE1
	cat $TEMP_OUTFILE | awk '{if($1 != "-") print $0}' |  grep "DATA_MISSING"|  sort -k4n >> $OUTFILE1
	cat $TEMP_OUTFILE | awk '{if($1 == "-") print $0}' | sort -k4n >> $OUTFILE1
	echo "Metrics written to:"
	echo $OUTFILE1
	rm $TEMP_OUTFILE
	chmod 777 $OUTFILE1
	cp $OUTFILE1 $MERGE_DIR
	
	export TMPDIR=$TEMP_PATH

	cat $OUTFILE1 |  mailx -a $OUTFILE1 -s "Split metrics for $PROJID" $MAIL_USID"@broadinstitute.org"

done 



