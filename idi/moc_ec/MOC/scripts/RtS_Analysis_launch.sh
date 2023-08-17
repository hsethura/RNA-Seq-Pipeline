#!/bin/sh

MOC_ID=$1

source /idi/moc_ec/MOC/scripts/bash_header

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

Q_HEAD="MOC_ID"
USID=`USID`

# identify options to pass to this script (as opposed to the pipeline)
SCRIPT_OPTIONS=`echo $@ | cut -d":" -f1,3`

### set options
# -conf: sets path to config file (default /idi/moc_ec/MOC/config_files/Universal_config.yaml)
# -q: sets name of header for Q_VAL (default MOC_ID)
# -move_key: setting to N skips moving google sheet to server (default Y)
# -use_p7: indicates fastq file is demultiplexed by P7 index (default Y)
# -use_p5: indicates fastq file is demultiplexed by P5 index (default N)
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

echo $RAW_SEQ_PATH
echo $RAW_SEQ_PATH

path_suff $@ "-moc_id "$MOC_ID_OPT" -user_id "$USER_ID
echo $RESPATH_SUFF

RESULTS_DIR=$RESULTS_PATH"/"$RESPATH_SUFF"/"
TEMP_DIR=$TEMP_PATH"/"$RESPATH_SUFF"/"
BAM_DIR=$BAM_PATH"/"$RESPATH_SUFF"/"
PARSED_REFDIR=$TEMP_DIR"/parsed_ref/"

echo "Results_dir: "$RESULTS_DIR
echo "Temp_dir: "$TEMP_DIR

mkdir -p $RESULTS_DIR
mkdir -p $TEMP_DIR

echo $KEY_SHEET

############## Set and assign pipeline options ###############=
### Set pipeline options
DEFAULT_PIPE_OPTIONS=" --ADD3 30 --ADD5 20 --gzip_merged --do_patho --MOC_id_ref $MOC_ID"
#DEFAULT_PIPE_OPTIONS=" --ADD3 30 --ADD5 20 --remove_splitted --gzip_merged --do_patho --MOC_id_ref $MOC_ID"
if [ $USER_ID == "N" ];then
	DEFAULT_PIPE_OPTIONS=$DEFAULT_PIPE_OPTIONS" --no_login_name"
fi
if [ $MOC_ID_OPT != "N" ];then
	DEFAULT_PIPE_OPTIONS=$DEFAULT_PIPE_OPTIONS" --MOC_id "$MOC_ID 
fi

# Assigns default pipeline options in case none are included in command line, replace default with 
# command line options, or add to default if --no_def_opt is not included
PIPE_COMMAND_OPTIONS=`echo $@ | cut -d':' -f2 | cut -d':' -f1`
PIPE_COMMAND_FLAG=`echo $@ | grep ":" | wc -w`
PIPE_COMMAND_DEF=`echo $@ | cut -d":" -f2 | cut -d":" -f1 | grep "\-no_def_opt" | wc -w`

if [ $PIPE_COMMAND_FLAG == "0" ];then
	PIPE_OPTIONS=$DEFAULT_PIPE_OPTIONS
else
	if [ $PIPE_COMMAND_DEF == "0" ];then
		PIPE_OPTIONS=$PIPE_COMMAND_OPTIONS" "$DEFAULT_PIPE_OPTIONS 
	else
		PIPE_OPTIONS=`echo $PIPE_COMMAND_OPTIONS" " | sed 's/\--no_def_opt//g'`
	fi
fi
if [ $USE_P5 == "Y" ];then
	PIPE_OPTIONS=$PIPE_OPTIONS" --use_p5"
fi
echo $PIPE_OPTIONS

########################################################

############## Import google sheet databases to server ###############=
if [ $IMPORT_GS == "Y" ];then
	echo "Running $GSIMPORT_SCRIPT to import google sheet databases to server"
	sh $GSIMPORT_SCRIPT
	echo "sh $GSIMPORT_SCRIPT"
fi
########################################################

### get GID from MOC_DB or use GID included as command line option
# if -GID_OPT not included in command line, get it from DB

if [ $GID_OPT == "0" ];then
	echo "sh $DB_SCRIPT $Q_HEAD,$MOC_ID Google_ID | grep "Google ID" | sed 's/ /_/g' | awk '{print $2}'"
	G_ID=`sh $DB_SCRIPT $Q_HEAD,$MOC_ID Google_ID | grep "Google ID" | sed 's/ /_/g' | awk '{print $2}'`
else
	G_ID=$GID_OPT
fi

if [ -z $G_ID ];then

	echo "GID not found in DB or entered in command line with -gid option"
	exit
fi

############## Running script to retrieve all entries from DB for MOC-ID ###############
echo " Running script to retrieve all entries from DB for $MOC_ID"
echo "sh $DB_SCRIPT $Q_HEAD,$MOC_ID"
sh $DB_SCRIPT $Q_HEAD,$MOC_ID
########################################################

############## moving key file to server ###############
if [ $MOVE_KEY == "Y" ];then 
	
	echo "Removing old key"
	rm -rf $KEY_FILE
	
	echo $KEY_SHEET
	
	echo "Moving key file to server..."
	echo "$KEY_SCRIPT -s $G_ID -t \"$KEY_SHEET\" -p $MOC_ID --Key_dir $KEY_DIR"
	
	$KEY_SCRIPT -s $G_ID -t "$KEY_SHEET" -p $MOC_ID --Key_dir $KEY_DIR
	### if key file not found or empty, stop pipeline
	if [ ! -s $KEY_FILE ];then
		ls -lrt $KEY_FILE
		exit
	fi
fi
########################################################


############## Move references from Gdrive or ref dir ###############
if [ $REFMOVE_OPT == "Y" ];then
	echo "Moving references to appropriate directory"
	if [ $PARSE_REF == "Y" ];then
	
		echo "...and parsing references, then moving them to $PARSED_REFDIR..."
	fi
	echo "sh $REF_MOVE_SCRIPT $MOC_ID -data_dir $PARSED_REFDIR -parse $PARSE_REF -move_key N"
	
	sh $REF_MOVE_SCRIPT $MOC_ID -data_dir $PARSED_REFDIR -parse $PARSE_REF -move_key N

	PROJID_HEAD=`config_read $CONFIG_FILE Proj`
	ALL_PROJIDS=`FIND_VALUES_FOR_HEADER $PROJID_HEAD $KEY_FILE | sort | uniq`
	echo "FIND_VALUES_FOR_HEADER $PROJID_HEAD $KEY_FILE | sort | uniq"
	
	if [ $PARSE_REF == "Y" ];then

		for PROJ_ID in $ALL_PROJIDS
		do
			REF_DIR=$TEMP_DIR"/"$PROJ_ID"/datadir/"
			mkdir -p $REF_DIR
			change_perms $REF_DIR 
		
			ALL_REFS=`FIND_HEADER_AND_FIELD $PROJ_ID $PROJID_HEAD Bacterial_reference $KEY_FILE | sort | uniq`
		
			for REF in $ALL_REFS
			do
				echo "cp $PARSED_REFDIR"/"*$REF* $REF_DIR"
				cp $PARSED_REFDIR"/"*$REF* $REF_DIR
			done
		done	
		
		echo $REF_DIR
		ls -lrt $REF_DIR
	fi
	
fi
#################################################

############## determining ALL_SMOCINDEX_SETS ###############
echo "Determining ALL_SMOCINDEX_SETS..."
all_smocindex_set $MOC_ID $@
echo "ALL_SMOCINDEX_SETS:" $ALL_SMOCINDEX_SETS
########################################################

############## make symlink in MOC_SPLIT_DIR fastq files in SEQ_DIR ###############

MOC_SPLIT_DIR=$INDXSPLIT_PATH"/"$RESPATH_SUFF"/"

rm -r $MOC_SPLIT_DIR
mkdir -p $MOC_SPLIT_DIR

###  

# if data has no index1, add --no_p7 as pipe option and set link to undemultiplexed file
if [ $USE_P7 == "N" ];then
	
	PIPE_OPTIONS=$PIPE_OPTIONS" --no_p7"
	for SMOCINDEX_SET in $ALL_SMOCINDEX_SETS
	do
		SMOC_ID=`echo $SMOCINDEX_SET | cut -d";" -f1 | sed 's/-//g'`
		SEQ_DIR=$RAW_SEQ_PATH$SMOC_ID"/"
		echo $SEQ_DIR
		
		ALL_FASTQS=`ls -lrt $SEQ_DIR*.*.fastq* | grep -v barcode | grep -v unmapped | grep -v unmatched | awk '{print $9}'`
		echo $ALL_FASTQS
		
		for FASTQ in $ALL_FASTQS
		do
			FILE=`echo $FASTQ | rev | cut -d"/" -f1 | rev`
			ln -s $FASTQ $MOC_SPLIT_DIR$FILE
		done
	done
	echo $MOC_SPLIT_DIR$FILE

# if data is indexed, sym link to demultiplexed files
else
	for SMOCINDEX_SET in $ALL_SMOCINDEX_SETS
	do
		echo "Making sym links to demultiplexed fastqs from $SEQ_DIR to $MOC_SPLIT_DIR..."
		
		SMOC_ID=`echo $SMOCINDEX_SET | cut -d";" -f1 | sed 's/-//g'`
		INDEX1=`echo $SMOCINDEX_SET | cut -d";" -f2`
		INDEX2=`echo $SMOCINDEX_SET | cut -d";" -f3`
	
		SEQ_DIR=$RAW_SEQ_PATH$SMOC_ID"/"

		### if no index2 is in SMOCINDEX_SET
		if [ $USE_P5 == "Y" ];then
			ALL_FASTQS=`ls -lrt $SEQ_DIR*unmapped.*.fastq* | grep -v barcode | grep $INDEX1 | grep $INDEX2 | awk '{print $9}'`
			NUM_FASTQS=`ls -lrt $SEQ_DIR*unmapped.*.fastq* | grep -v barcode | grep $INDEX1 | grep $INDEX2 | awk '{print $9}' | wc -l`
			echo $NUM_FASTQS, $INDEX1, $INDEX2
			if [ $NUM_FASTQS == 0 ];then
		
				echo "Missing fastq files for "$INDEX1" in "$SEQ_DIR"!"
			fi
		else 
			ALL_FASTQS=`ls -lrt $SEQ_DIR*unmapped.*.fastq* | grep -v barcode | grep $INDEX1 | awk '{print $9}'`
			NUM_FASTQS=`ls -lrt $SEQ_DIR*unmapped.*.fastq* | grep -v barcode | grep $INDEX1 | awk '{print $9}' | wc -l`
			echo $NUM_FASTQS, $INDEX2
			if [ $NUM_FASTQS == 0 ];then
		
				echo "Missing fastq files for "$INDEX1" and "$INDEX2" in "$SEQ_DIR"!"
			fi
		fi
			
		
		
		for FASTQ in $ALL_FASTQS
		do
			FILE=`echo $FASTQ | rev | cut -d"/" -f1 | rev`
			ln -s $FASTQ $MOC_SPLIT_DIR$FILE 2>/dev/null
		done
	done
fi


ls -lrt $MOC_SPLIT_DIR
echo $MOC_SPLIT_DIR
###########################################################################


############## Set command line and launch analysis pipeline ##############
if [ $NO_PIPE != "Y" ];then
	echo "Launching analysis pipeline..."
	if [ $QSUB == "Y" ];then
		echo "sh $QSUB_SCRIPT $PIPE_SCRIPT  -c $CONFIG_FILE --key_path $KEY_FILE --raw_seq_path $MOC_SPLIT_DIR $PIPE_OPTIONS"
		job_id=`sh $QSUB_SCRIPT $PIPE_SCRIPT -c $CONFIG_FILE --key_path $KEY_FILE --raw_seq_path $MOC_SPLIT_DIR $PIPE_OPTIONS | awk '{print $3}'`
		UGER_test $job_id
		
	else
		echo "$PIPE_SCRIPT -c $CONFIG_FILE --key_path $KEY_FILE --raw_seq_path $MOC_SPLIT_DIR $PIPE_OPTIONS"
		$PIPE_SCRIPT -c $CONFIG_FILE --key_path $KEY_FILE --raw_seq_path $MOC_SPLIT_DIR $PIPE_OPTIONS
	fi
fi
##########################################################################

######## Run metrics script to calculate distribution of reads per barcode

echo "sh $SPLITMET_SCRIPT $MOC_ID $@"
sh $SPLITMET_SCRIPT $MOC_ID $@


NUM_MET_FILES=`ls -lrt $RESULTS_DIR/* | grep metrics.txt | wc -l`
NUM_COUNT_FILES=`ls -lrt $RESULTS_DIR/* | grep _counts.tsv | wc -l`
NUM_CG_PAIRS=`FIND_VALUES_FOR_HEADER CG_ID $KEY_FILE | wc -l`
PROJID_HEAD=`config_read $CONFIG_FILE Proj`
ALL_PROJIDS=`FIND_VALUES_FOR_HEADER $PROJID_HEAD $KEY_FILE | sort | uniq`


echo "NUM_MET_FILES:"	$NUM_MET_FILES
echo "NUM_COUNT_FILES:"	$NUM_COUNT_FILES
echo "NUM_CG_PAIRS:"	$NUM_CG_PAIRS
echo "ALL_PROJIDS:"	$ALL_PROJIDS


if [ $NUM_MET_FILES -gt 0 ];then
	
	### run metric-key join script 
	echo "Running "$JOIN_SCRIPT" to combine metrics and key files"
	echo "sh $JOIN_SCRIPT $MOC_ID $SCRIPT_OPTIONS"
	sh $JOIN_SCRIPT $MOC_ID $SCRIPT_OPTIONS
	
	### test keymet file to ensure data has been generated for all samples
	for PROJ_ID in $ALL_PROJIDS
	do
		MET_FILE=`ls -lrt $RESULTS_DIR"/"$PROJ_ID"/"*"metrics.txt" | awk '{print $9}'`	
		JOIN_FILE=$JOIN_PATH"/"$MOC_ID"/"$PROJ_ID"_KeyMetrics.txt"
		ALL_FIELDS=`FIELD_HEADER $JOIN_FILE Sample_ID Pcnt_bc_in_pool Total_reads pcnt_aligned CDS_total_counts_for_replicon rRNA_pcnt_of_counted | tr '\n' ','`
		TEMP_FILE=$TEMP_DIR$MOC_ID"_"$PROJ_ID"_temp.txt"
		
		echo "Results directory: $RESULTS_DIR"/"$PROJ_ID"/"" > $TEMP_FILE
		echo "" >> $TEMP_FILE
		cat $JOIN_FILE | awk -F"\t" -v ALL_FIELDS=$ALL_FIELDS '{
															y=split(ALL_FIELDS,ar,",")
															for(i=1; i < y; i++)
															{
																if(i == 1)
																	printf"%-25s\t", $(ar[i])
																else
																	printf"%-10s\t", $(ar[i])
															}
															print ""
												
													}' | awk '{
																	if(NF != 6) 
																		print "*********DATA missing!***********\t" $0
																	else 
																		print $0
																	}' >> $TEMP_FILE
												
												
		MET_CHECK=`echo $TEMP_FILE	| awk '{print NF}' | sort | uniq | wc -l`
		echo $MET_CHECK	

	### send email with metrics attachment for each project alerting that analysis is complete
	
		cat $TEMP_FILE							
		cat $TEMP_FILE |  mailx -a $JOIN_FILE -s "Analysis of $MOC_ID $PROJ_ID is done" $USID"@broadinstitute.org"
	done


fi

if [ $NUM_CG_PAIRS -gt 0 ];then
	### run DE pipeline
	echo "Running "$DE_PIPE" to generate edgeR and DESeq2 comparisons"
	echo "sh $DE_PIPE $MOC_ID $SCRIPT_OPTIONS"
	sh $DE_PIPE $MOC_ID $SCRIPT_OPTIONS
fi

### change permissions for Results and temp dirs
change_perms $MOC_SPLIT_DIR $RESULTS_DIR $TEMP_DIR $KEY_FILE $BAM_DIR



