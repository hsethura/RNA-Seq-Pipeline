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
# -conf: sets path to config file (default /idi/moc_ec/MOC/config_files/PC_config.yaml)
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
# -PROJ_TYPE: project production (P) or development (D) (default P)

### determining paths and headers 
### default config file is /broad/IDP-Dx_storage/MOC/config_files/PC_config.yaml
paths_and_headers $MOC_ID $@

############################ MOVE KEY FILE TO SERVER ###################################
	
	############### moving key file to server #############################
	if [ $MOVE_KEY == "Y" ];then 
		echo "Running move_key..."
		move_key		
	fi
##########################################################################

### Repeating to reset path with updated proj name(s) in key
### determining paths and headers 
### default config file is /broad/IDP-Dx_storage/MOC/config_files/PC_config.yaml
paths_and_headers $MOC_ID $@


### run path_suff function to set RESPATH_SUFF.
##	If -moc_id N included in command line, do not moc_ID to RESPATH_SUFF.  
##	If -user_id included with Y or no string add USID to RESPATH_SUFF.  
##	If -user_id followed by userID, add userID to RESPATH_SUFF.


mkdir -p $RESULTS_DIR
mkdir -p $TEMP_DIR
mkdir -p $BAM_DIR

echo "RESULTS_DIR: "$RESULTS_DIR
echo "TEMP_DIR: "$TEMP_DIR

RAW_SYM_PATH=`extract_option -symlink $RAWSYM_PATH 1 $@`
RAW_SEQ_PATH=`extract_option -raw_seq_path $RAWSYM_PATH 1 $@`

echo "RAW_SYM_PATH: "$RAW_SYM_PATH


############################ SET ANALYSIS PIPELINE OPTIONS ###################################

DEFAULT_PIPE_OPTIONS=" --ADD3 30 --ADD5 20 --use_p5 --gzip_merged --do_patho --MOC_id_ref $MOC_ID"
if [ $USER_ID == "N" ];then
	DEFAULT_PIPE_OPTIONS=$DEFAULT_PIPE_OPTIONS" --no_login_name"
fi
if [ $MOC_ID_OPT != "N" ];then
	DEFAULT_PIPE_OPTIONS=$DEFAULT_PIPE_OPTIONS" --MOC_id "$MOC_ID 
fi
if [ $DO_HOST != "N" ];then
	DEFAULT_PIPE_OPTIONS=`echo $DEFAULT_PIPE_OPTIONS | sed 's/do_patho/do_host/g'`
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
echo $PIPE_OPTIONS


##########################################################################



############################ SET TEMP DIR PATHS ###################################

SAMP_KEY_ID="Sample_ID"
SID_F=`FIELD_HEADER $KEY_FILE $SAMP_KEY_ID`

echo "Results_dir: "$RESULTS_DIR
echo "Temp_dir: "$TEMP_DIR
echo "Data_dir: "$DATA_DIR
mkdir -p $RESULTS_DIR
mkdir -p $TEMP_DIR
rm -r $MOC_SYM_TEST_DIR
mkdir -p $MOC_SYM_TEST_DIR
echo $KEY_SHEET

############################ MOVE AND PARSE REF FILES ###################################

	############## Move references from Gdrive or ref dir and parse --> data dir ###############
	if [ $REFMOVE_OPT == "Y" ];then
		
		CHECK_DIR=$DATA_DIR
		ALL_SUFF="_ALL.gff,.fna" 

		echo "Moving references to appropriate directory"	
		echo "sh $REF_MOVE_SCRIPT $MOC_ID -data_dir $DATA_DIR -parse $PARSE_REF -move_key N -proj_type $PROJ_TYPE $@"
		sh $REF_MOVE_SCRIPT $MOC_ID -data_dir $DATA_DIR -parse $PARSE_REF -move_key N -proj_type $PROJ_TYPE $@
	


		############## Check if all gff files were parsed ##############
		echo "mod_check $MOC_ID $CHECK_DIR $ALL_SUFF -fail_exit Y -email Y -header "Bacterial_reference" -delimit none $@"
		mod_check $MOC_ID $CHECK_DIR $ALL_SUFF -fail_exit Y -email Y -header "Bacterial_reference" -delimit none $@
			
	fi
##########################################################################



######################## RUN MOD 1: SPLIT ###########################

	############## Launch pipeline to split by barcode ##############

	if [ $SPLIT == "Y" ] || [ $SPLIT == "O" ] || [ $SPLIT == "R" ];then

		MOD_KEY_FILE=$KEY_FILE

		######################  IF TEST_PIPE  #########################
		if [ $TEST_PIPE == "Y" ];then
			make_test_fastqs $NUM_TEST_READS
			MOC_SYM_DIR=$MOC_SYM_TEST_DIR
		fi
		
		############## Launch pipeline to split demultiplexed fastqs ##############
		MOD_PIPE_OPTIONS=$PIPE_OPTIONS" --no_merge --no_align --no_count"
		
		echo "Launching analysis pipeline to split fastqs by inline barcode..."
		echo "$PIPE_SCRIPT -c $CONFIG_FILE --key_path $MOD_KEY_FILE --raw_seq_path $MOC_SYM_DIR $MOD_PIPE_OPTIONS"
 		$PIPE_SCRIPT -c $CONFIG_FILE --key_path $MOD_KEY_FILE --raw_seq_path $MOC_SYM_DIR $MOD_PIPE_OPTIONS

		######## Run metrics script to calculate distribution of reads per barcode
		MOD_KEY_FILE=$KEY_FILE

		echo "sh $SPLITMET_SCRIPT $MOC_ID $@"
		sh $SPLITMET_SCRIPT $MOC_ID $@
			
	fi
##########################################################################


################################# RUN MOD 2: MERGE ###################################

	if [ $MERGE == "Y" ] || [ $MERGE == "O" ] || [ $MERGE == "R" ];then
		
		CHECK_DIR=$MERGE_DIR
		ALL_SUFF="R1.fastq.gz,R2.fastq.gz" 
		MOD_KEY_FILE=$KEY_FILE

		if [ $MERGE == "R" ];then
			redo_key $CHECK_DIR
			MOD_KEY_FILE=$REDO_KEY_FILE
		fi
		
		############## Launch pipeline to merge split fastqs ##############
		MOD_PIPE_OPTIONS=$PIPE_OPTIONS" --no_split --no_align --no_count"
		
		echo "Launching analysis pipeline to merge split fastqs by sample..."
		echo "$PIPE_SCRIPT -c $CONFIG_FILE --key_path $MOD_KEY_FILE --raw_seq_path $MOC_SYM_DIR $MOD_PIPE_OPTIONS"
		$PIPE_SCRIPT -c $CONFIG_FILE --key_path $MOD_KEY_FILE --raw_seq_path $MOC_SYM_DIR $MOD_PIPE_OPTIONS
		MOD_KEY_FILE=$KEY_FILE
		
		############## Check if all files were merged ##############
		mod_check $MOC_ID $CHECK_DIR $ALL_SUFF -fail_exit N -email Y $@
	fi
##########################################################################

	
################################# RUN MOD 3: ALIGN ###################################

	############## Launch pipeline to align merged fastqs ##############
	if [ $ALIGN == "Y" ] || [ $ALIGN == "O" ] || [ $ALIGN == "R" ];then

		CHECK_DIR=$ALIGN_DIR
		ALL_SUFF="_u.bam,_pe.bam" 
		MOD_KEY_FILE=$KEY_FILE

		if [ $ALIGN == "R" ];then
			redo_key 
			MOD_KEY_FILE=$REDO_KEY_FILE
		fi
		
		############## Launch pipeline to align merged fastqs ##############
		MOD_PIPE_OPTIONS=$PIPE_OPTIONS" --no_split --no_merge --no_count"

		echo "Launching analysis pipeline to align merged fastqs to reference genomes..."
		echo "$PIPE_SCRIPT -c $CONFIG_FILE --key_path $MOD_KEY_FILE --raw_seq_path $MOC_SYM_DIR $MOD_PIPE_OPTIONS"
		$PIPE_SCRIPT -c $CONFIG_FILE --key_path $MOD_KEY_FILE --raw_seq_path $MOC_SYM_DIR $MOD_PIPE_OPTIONS
		
		############## Check if all files were merged ##############
		mod_check $MOC_ID $CHECK_DIR $ALL_SUFF -fail_exit N -email Y -ref Y $@
	fi
##########################################################################

################################# RUN MOD 4: COUNT ###################################

	############## Launch pipeline to count aligned bams ##############
	if [ $COUNT == "Y" ] || [ $COUNT == "O" ] || [ $COUNT == "R" ];then
	
		CHECK_DIR=$ALIGN_DIR
		ALL_SUFF=".counts" 
		MOD_KEY_FILE=$KEY_FILE

		MOD_PIPE_OPTIONS=$PIPE_OPTIONS" --no_split --no_merge --no_align"

		echo "Launching analysis pipeline to count reads per genes in aligned bams..."
		echo "$PIPE_SCRIPT -c $CONFIG_FILE --key_path $MOD_KEY_FILE --raw_seq_path $MOC_SYM_DIR $MOD_PIPE_OPTIONS"
		$PIPE_SCRIPT -c $CONFIG_FILE --key_path $MOD_KEY_FILE --raw_seq_path $MOC_SYM_DIR $MOD_PIPE_OPTIONS
		
		############## Check if all files were merged ##############
		mod_check $MOC_ID $CHECK_DIR $ALL_SUFF -fail_exit N -email Y -ref Y $@
	fi
##########################################################################


NUM_CG_PAIRS=`FIND_VALUES_FOR_HEADER CG_ID $KEY_FILE | wc -l`
PROJID_HEAD=`config_read $CONFIG_FILE Proj`
ALL_PROJIDS=`FIND_VALUES_FOR_HEADER $PROJID_HEAD $KEY_FILE | sort | uniq`


echo "NUM_MET_FILES:"	$NUM_MET_FILES
echo "NUM_COUNT_FILES:"	$NUM_COUNT_FILES
echo "NUM_CG_PAIRS:"	$NUM_CG_PAIRS
echo "ALL_PROJIDS:"		$ALL_PROJIDS


#if [ $NUM_MET_FILES -gt 0 ];then
	
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
		cat $JOIN_FILE | sed '/^$/d' | awk -F"\t" -v ALL_FIELDS=$ALL_FIELDS '{
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
																	if(NF == 6) 
																	{
																		for(i=1; i < NF; i++)
																			printf"%-10s\t", $i
																		print $NF
																	
																	}	
																	else 
																	{
																		for(i=1; i < NF+1; i++)
																			printf"%-10s\t", $i
																		print "*********DATA missing!***********"

																	}
																}'  | sed 's/CDS_total_counts_for_replicon/total_CDS/g' >> $TEMP_FILE
												
												
		MET_CHECK=`echo $TEMP_FILE	| awk '{print NF}' | sort | uniq | wc -l`
		echo $MET_CHECK	

	### send email with metrics attachment for each project alerting that analysis is complete
		export TMPDIR=$TEMP_DIR
	    
		cat $TEMP_FILE							
		cat $TEMP_FILE |  mailx -a $JOIN_FILE -s "Analysis of $MOC_ID $PROJ_ID is done" $USID"@broadinstitute.org"
	done


#fi

if [ $NUM_CG_PAIRS -gt 0 ];then
	### run DE pipeline
	echo "Running "$DE_PIPE" to generate edgeR and DESeq2 comparisons"
	echo "sh $DE_PIPE $MOC_ID $SCRIPT_OPTIONS -move_key N"
	sh $DE_PIPE $MOC_ID $SCRIPT_OPTIONS -move_key N
fi

### change permissions for Results and temp dirs
change_perms $MOC_SYM_DIR $RESULTS_DIR $TEMP_DIR $KEY_FILE $BAM_DIR



