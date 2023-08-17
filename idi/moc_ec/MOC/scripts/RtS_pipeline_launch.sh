#!/bin/sh

### This scrip takes in a MOC_ID and RtS pipeline options and launches Nirmalya's pipeline

MOC_ID=$1

source /idi/moc_ec/MOC/scripts/bash_header

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

Q_HEAD="MOC_ID"
USID=`USID`

# identify options to pass to this script (as opposed to the pipeline)
LAUNCH_OPTIONS=`echo $@ | cut -d":" -f1,3 | awk '{for(i=2; i < NF+1; i++) print $i}'`

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



### identify options to pass to DE and metrics scripts etc (as opposed to the pipeline)
LAUNCH_OPTIONS=`echo $@ | cut -d":" -f1,3`

### identify options to pass to the pipeline
PIPE_OPTIONS=`echo $@ | cut -d':' -f2 | cut -d':' -f1`

RAW_SYM_PATH=`extract_option -symlink $RAWSYM_PATH 1 $@`
RAW_SEQ_PATH=`extract_option -raw_seq_path $RAWSYM_PATH 1 $@`


path_suff $@ "-moc_id "$MOC_ID_OPT" -user_id "$USER_ID
echo $RESPATH_SUFF

RESULTS_DIR=$RESULTS_PATH"/"$RESPATH_SUFF"/"
TEMP_DIR=$TEMP_PATH"/"$RESPATH_SUFF"/"
BAM_DIR=$BAM_PATH"/"$RESPATH_SUFF"/"
PARSED_REFDIR=$TEMP_DIR"/parsed_ref/"
SYM_DIR=$SYM_PATH"/"$MOC_ID"/"

echo "Results_dir: "$RESULTS_DIR
echo "Temp_dir: "$TEMP_DIR

mkdir -p $RESULTS_DIR
mkdir -p $TEMP_DIR


echo "$PIPE_SCRIPT $PIPE_OPTIONS"
$PIPE_SCRIPT $PIPE_OPTIONS
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
	echo "sh $JOIN_SCRIPT $MOC_ID $LAUNCH_OPTIONS"
	sh $JOIN_SCRIPT $MOC_ID $LAUNCH_OPTIONS
	
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
	echo "sh $DE_PIPE $MOC_ID $LAUNCH_OPTIONS -move_key N"
	sh $DE_PIPE $MOC_ID $LAUNCH_OPTIONS -move_key N
fi

### change permissions for Results and temp dirs
change_perms $MOC_SPLIT_DIR $RESULTS_DIR $TEMP_DIR $KEY_FILE $BAM_DIR



