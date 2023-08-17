#!/bin/sh

MOC_ID=$1

#source /broad/IDP-Dx_storage/MOC/scripts/bash_header

### source all functions 
source "/broad/IDP-Dx_storage/MOC/scripts/MOC_functions.sh"

Q_HEAD="MOC_ID"
USID=`USID`

# identify options to pass to this script (as opposed to the pipeline)
SCRIPT_OPTIONS=`echo $@ | cut -d":" -f1,3`

### set options
# -conf: sets path to config file (default /broad/IDP-Dx_storage/MOC/config_files/Universal_config.yaml)
# -q: sets name of header for Q_VAL (default MOC_ID)
# -move_key: setting to N skips moving google sheet to server (default Y)
# -no_index1: indicates fastq file is not demultiplexed by P7 index (default N)
# -no_index2: indicates fastq file is not demultiplexed by P5 index (default Y)
# -no_pipe: including or setting to Y skips running pipeline (default N)
# -qsub: run pipeline in UGER (default N)
# -gid: enter google ID if missing from DB (default 0 => get it from DB)
# -smocid: enter smocid if missing from DB (default 0 => get it from DB)
# -import_gs: updates DBs from google drive on server (default Y)
# -user_id: add userID to temp and results paths (default N.  If Y, add login name, if value added after option that userID will be added)
# -moc_id: add MOC_ID to temp and results paths (default Y)

### determining paths and headers 
### default config file is /broad/IDP-Dx_storage/MOC/config_files/Universal_config.yaml
paths_and_headers $MOC_ID $@

### run path_suff function to set RESPATH_SUFF.
##	If -moc_id N included in command line, do not moc_ID to RESPATH_SUFF.  
##	If -user_id included with Y or no string add USID to RESPATH_SUFF.  
##	If -user_id followed by userID, add userID to RESPATH_SUFF.
 
RAW_SEQ_PATH=`extract_option -raw_seq_path $SEQ_PATH 1 $@`

path_suff $@ "-moc_id "$MOC_ID_OPT" -user_id "$USER_ID

PARSED_REFDIR=$TEMP_DIR"/parsed_ref/"


SAMP_ID_FIELD=`FIELD_HEADER $KEY_FILE Sample_ID`
PID_FIELD=`FIELD_HEADER $KEY_FILE Project_ID`
ACC_FIELD=`FIELD_HEADER $KEY_FILE Bacterial_reference`
I1_FIELD=`FIELD_HEADER $KEY_FILE Index1_Seq`
I2_FIELD=`FIELD_HEADER $KEY_FILE Index2_Seq`
IL_FIELD=`FIELD_HEADER $KEY_FILE Inline_Seq`
ILN_FIELD=`FIELD_HEADER $KEY_FILE Inline_Name`
POOL_FIELD=`FIELD_HEADER $KEY_FILE Pool_ID`

ALL_PROJID=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v PID_FIELD=$PID_FIELD '{print $PID_FIELD}' | sort | uniq`
echo $ALL_PROJID
echo $SAMP_ID_FIELD
echo $PID_FIELD
echo $ACC_FIELD

for PROJ in $ALL_PROJID
do

	RESULTS_DIR=$RESULTS_PATH"/"$RESPATH_SUFF"/"$PROJ"/"
	TEMP_DIR=$TEMP_PATH"/"$RESPATH_SUFF"/"$PROJ"/"
	BAM_DIR=$BAM_PATH"/"$RESPATH_SUFF"/"$PROJ"/"
	SPLIT_DIR=$TEMP_PATH"/"$RESPATH_SUFF"/"$PROJ"/splitdir"
	MERGE_DIR=$TEMP_PATH"/"$RESPATH_SUFF"/"$PROJ"/mergedir"
	PATHO_DIR=$TEMP_PATH"/"$RESPATH_SUFF"/"$PROJ"/patho_result/"$PROJ

	
# 	ls -lrt $TEMP_DIR
	echo $SPLIT_DIR
# 	ls -lrt $MERGE_DIR
 	echo  $PATHO_DIR
	ALL_SID=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v SID_FIELD=$SAMP_ID_FIELD -v PID_FIELD=$PID_FIELD -v PROJ=$PROJ '{if($PID_FIELD=PROJ) print $SID_FIELD}' | sort | uniq`
	#echo $ALL_SID
	
	for SAMP in $ALL_SID
	do
		I1=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v SID_FIELD=$SAMP_ID_FIELD -v I1_FIELD=$I1_FIELD -v SAMP=$SAMP '{if($SID_FIELD==SAMP) print $I1_FIELD}' `
		I2=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v SID_FIELD=$SAMP_ID_FIELD -v I2_FIELD=$I2_FIELD -v SAMP=$SAMP '{if($SID_FIELD==SAMP) print $I2_FIELD}' `
		IL=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v SID_FIELD=$SAMP_ID_FIELD -v IL_FIELD=$IL_FIELD -v SAMP=$SAMP '{if($SID_FIELD==SAMP) print $IL_FIELD}'`
		POOLID=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v SID_FIELD=$SAMP_ID_FIELD -v POOL_FIELD=$POOL_FIELD -v SAMP=$SAMP '{if($SID_FIELD==SAMP) print $POOL_FIELD}'`
		ILN=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v SID_FIELD=$SAMP_ID_FIELD -v ILN_FIELD=$ILN_FIELD -v SAMP=$SAMP '{if($SID_FIELD==SAMP) print $ILN_FIELD}'`
		
		NUM_FILES=`ls -lrt $SPLIT_DIR | grep $IL | grep $I1 | grep $I2 | wc -l`

		echo $SPLIT_DIR, $SAMP, $IL, $I1, $I2, $ILN, $POOLID, $NUM_FILES
		#ls -lrt $SPLIT_DIR | grep $IL | grep $I1 | grep $I2 
 
	done
	
	for SAMP in $ALL_SID
	do
		I1=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v SID_FIELD=$SAMP_ID_FIELD -v I1_FIELD=$I1_FIELD -v SAMP=$SAMP '{if($SID_FIELD==SAMP) print $I1_FIELD}' `
		I2=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v SID_FIELD=$SAMP_ID_FIELD -v I2_FIELD=$I2_FIELD -v SAMP=$SAMP '{if($SID_FIELD==SAMP) print $I2_FIELD}' `
		IL=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v SID_FIELD=$SAMP_ID_FIELD -v IL_FIELD=$IL_FIELD -v SAMP=$SAMP '{if($SID_FIELD==SAMP) print $IL_FIELD}'`
		POOLID=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v SID_FIELD=$SAMP_ID_FIELD -v POOL_FIELD=$POOL_FIELD -v SAMP=$SAMP '{if($SID_FIELD==SAMP) print $POOL_FIELD}'`
		ILN=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v SID_FIELD=$SAMP_ID_FIELD -v ILN_FIELD=$ILN_FIELD -v SAMP=$SAMP '{if($SID_FIELD==SAMP) print $ILN_FIELD}'`
		
		NUM_FILES=`ls -lrt $MERGE_DIR | grep $SAMP | wc -l`

		echo $MERGE_DIR, $SAMP, $POOLID, $NUM_FILES
		#ls -lrt $SPLIT_DIR | grep $IL | grep $I1 | grep $I2 
 
	done
	for SAMP in $ALL_SID
	do
		I1=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v SID_FIELD=$SAMP_ID_FIELD -v I1_FIELD=$I1_FIELD -v SAMP=$SAMP '{if($SID_FIELD==SAMP) print $I1_FIELD}' `
		I2=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v SID_FIELD=$SAMP_ID_FIELD -v I2_FIELD=$I2_FIELD -v SAMP=$SAMP '{if($SID_FIELD==SAMP) print $I2_FIELD}' `
		IL=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v SID_FIELD=$SAMP_ID_FIELD -v IL_FIELD=$IL_FIELD -v SAMP=$SAMP '{if($SID_FIELD==SAMP) print $IL_FIELD}'`
		POOLID=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v SID_FIELD=$SAMP_ID_FIELD -v POOL_FIELD=$POOL_FIELD -v SAMP=$SAMP '{if($SID_FIELD==SAMP) print $POOL_FIELD}'`
		ILN=`cat $KEY_FILE | grep -v "###" | sed 1d | awk -F"\t" -v SID_FIELD=$SAMP_ID_FIELD -v ILN_FIELD=$ILN_FIELD -v SAMP=$SAMP '{if($SID_FIELD==SAMP) print $ILN_FIELD}'`
		
		echo $SAMP, $POOLID, $PATHO_DIR
		ls -lrt $PATHO_DIR | grep $SAMP 

		#ls -lrt $SPLIT_DIR | grep $IL | grep $I1 | grep $I2 
 
	done

done

exit

ls -lrt  $BAM_DIR
ls -lrt  $RESULTS_DIR

