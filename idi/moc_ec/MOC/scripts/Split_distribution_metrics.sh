#!/bin/sh

MOC_ID=$1

source /idi/moc_ec/MOC/scripts/bash_header

Q_HEAD="MOC_ID"

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"


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

CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/PC_config.yaml" 1 $@`
MOVE_KEY=`extract_option -move_key Y 1 $@`
USID=`extract_option -uid $USID 1 $@`
SPLIT_DIR=`extract_option -SPLIT_DIR N 1 $@`
DUAL=`extract_option -dual N 1 $@`


echo $USID

IN_DIR=$TEMP_PATH"/"$RESPATH_SUFF"/"
OUT_DIR=$TEMP_PATH"/"$RESPATH_SUFF"/inline_distribution/"
SPLIT_DIR=$IN_DIR"/*/splitdir/"
MERGE_DIR=$IN_DIR"/*/mergedir/"

mkdir -p $OUT_DIR


ALL_INLINES=`cat $BC_FILE | awk '{printf "%s;%s ", $2,$1}'`

echo $ALL_INLINES


# find gid move keyfile to server

G_ID=`sh $DB_SCRIPT $Q_HEAD,$MOC_ID Google_ID | sed 1d | awk '{print $2}'`
KEY_DIR=`config_read $CONFIG_FILE Key_base`
KEY_FILE=$KEY_DIR$MOC_ID"_key.txt"

### move key file to server
if [ $MOVE_KEY != "Y" ];then
	echo "Moving key file to server..."
	
	echo $KEY_SCRIPT -s $G_ID -t "Sample Information" -p $MOC_ID --Key_dir $KEY_DIR
	$KEY_SCRIPT -s $G_ID -t "Sample Information" -p $MOC_ID --Key_dir $KEY_DIR
fi

### determine if single or dual indexed

DUAL=`ls $SPLIT_DIR"/"*combined.txt | rev | cut -d"/" -f1 | rev | sed 's/_logfile_combined.txt//g' | awk '{if($1 ~ /_/) print "Y"; else print "N"}' | sort -u`

ls $SPLIT_DIR"/"*combined.txt
echo "Dual: "$DUAL


### get all project IDs

PROJID_F=`FIELD_HEADER $KEY_FILE Project_ID`
SAMPID_F=`FIELD_HEADER $KEY_FILE Sample_ID`
INDEX1_F=`FIELD_HEADER $KEY_FILE Index1_seq`
INDEX2_F=`FIELD_HEADER $KEY_FILE Index2_seq`
INLINE_F=`FIELD_HEADER $KEY_FILE Inline_Seq`

ALL_PROJID=`FIND_VALUES_FOR_HEADER Project_ID $KEY_FILE | sort -u`

if [ $DUAL == "N" ];then
	ALL_INDEXES=`FIND_VALUES_FOR_HEADER Index1_seq $KEY_FILE | sort -u`
else
	ALL_INDEXES=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" '{print $'$INDEX1_F'"_"$'$INDEX2_F'}' | sort -u`
fi

echo "ALL_INDEXES:	"$ALL_INDEXES
echo "ALL_PROJID:	"$ALL_PROJID


for PROJID in $ALL_PROJID
do
# 	ALL_INLINES=`FIND_VALUES_FOR_HEADER Inline_Seq $KEY_FILE | sort -u`
# 		
# 	echo "ALL_INLINES:"	$ALL_INLINES

	TEMP_OUTFILE=$TEMP_PATH"/"$RESPATH_SUFF"/"$PROJID"_inline_distribtion_temp.txt"
	OUTFILE1=$OUT_DIR"/"$PROJID"_inline_distribtion.txt"

	if [ ! -s $SPLIT_DIR ];then
		echo $SPLIT_DIR" not found!"
	fi	
		
	echo "" | sed 1d > $TEMP_OUTFILE 

	for INLINE_PAIR in $ALL_INLINES
	do
		INLINE_SEQ=`echo $INLINE_PAIR | cut -d";" -f1`
		INLINE_NAME=`echo $INLINE_PAIR | cut -d";" -f2`

		
		for INDEX in $ALL_INDEXES
		do
			if [ $DUAL == "N" ];then
				
				INDEX1=$INDEX
				SAMPID=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v PROJID=$PROJID -v INLINE_SEQ=$INLINE_SEQ -v INDEX1=$INDEX1 '{if($'$PROJID_F'==PROJID && $'$INLINE_F'==INLINE_SEQ && $'$INDEX1_F'==INDEX1) print $'$SAMPID_F'}' | sed 's/ /_/g'`
			else
				INDEX1=`echo $INDEX | cut -d"_" -f1`
				INDEX2=`echo $INDEX | cut -d"_" -f2`
				SAMPID=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v PROJID=$PROJID -v INLINE_SEQ=$INLINE_SEQ -v INDEX1=$INDEX1 -v INDEX2=$INDEX2 '{if($'$PROJID_F'==PROJID && $'$INLINE_F'==INLINE_SEQ && $'$INDEX1_F'==INDEX1 && $'$INDEX2_F'==INDEX2) print $'$SAMPID_F'}' | sed 's/ /_/g'`
				#echo $INDEX1, $INDEX2, $SAMPID 
			fi
			
			if [ -z $SAMPID ];then
				SAMPID="-"
			fi
			if [ -s $SPLIT_DIR"/"*$INDEX*combined.txt ];then
				data=`cat $SPLIT_DIR"/"*$INDEX*combined.txt | grep $INLINE_SEQ | awk '{print $2, $3, $4}'` 
			else 
				if [ $SAMPID != "-" ];then
					data="DATA_MISSING!"
				else
					data="- - - -"
				fi
			fi
			echo $SAMPID $INLINE_SEQ $INDEX $data	$INLINE_NAME 
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



