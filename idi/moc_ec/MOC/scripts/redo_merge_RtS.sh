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
# -moc_id: add MOC_ID to temp and results paths (default Y)

### determining paths and headers 
### default config file is /broad/IDP-Dx_storage/MOC/config_files/PC_config.yaml
paths_and_headers $MOC_ID $@

### run path_suff function to set RESPATH_SUFF.
##	If -moc_id N included in command line, do not moc_ID to RESPATH_SUFF.  
##	If -user_id included with Y or no string add USID to RESPATH_SUFF.  
##	If -user_id followed by userID, add userID to RESPATH_SUFF.

path_suff $@ "-moc_id "$MOC_ID_OPT" -user_id "$USER_ID
echo $RESPATH_SUFF

SPLIT=`extract_option -split Y 1 $@`
COMBINE=`extract_option -combine Y 1 $@`
ALIGN=`extract_option -align Y 1 $@`
PROJ_TYPE=`project_type $MOC_ID`
INCL_ADD=`extract_option -incl_add N 1 $@`

echo "Proj_type: "$PROJ_TYPE


############## moving key file to server ###############
if [ $MOVE_KEY == "Y" ];then

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
	chmod 777 $KEY_FILE
fi


### set path to mergedir 
TEMP_DIR=$TEMP_PATH"/"$RESPATH_SUFF"/"
ALL_PROJ=`cat $KEY_FILE | grep -v "###" | sed 1d | awk '{print $2}' | sort | uniq`
ALL_PROJ_DIR=`echo $ALL_PROJ | sed 's/ /_/g'`
MERGE_DIR=$TEMP_DIR"/"$ALL_PROJ_DIR"/mergedir"
PREMERGE_DIR=$MERGE_DIR"/premerged"


### Running pipeline to split and merge with all samples
if [ $SPLIT == "Y" ];then
	
	rm -r $MERGE_DIR
	echo $PREMERGE_DIR
	echo "Running pipeline to split and merge with all samples"
	echo "/idi/moc_ec/MOC/scripts/MOC_RtS_pipe_v2.sh $MOC_ID $@ -move_key N -align N -count N"
	sh /idi/moc_ec/MOC/scripts/MOC_RtS_pipe_v2.sh $MOC_ID $@ -move_key N -align N -count N
	rm $PREMERGE_DIR
	
	echo "Making $PREMERGE_DIR and moving all files in $MERGE_DIR to it"
	mkdir -p $PREMERGE_DIR
	mv $MERGE_DIR/* $PREMERGE_DIR

fi



if [ $COMBINE == "Y" ];then

	if [ ! -s $PREMERGE_DIR ];then
		echo "Making $PREMERGE_DIR and moving all files in $MERGE_DIR to it"
		mkdir -p $PREMERGE_DIR
		mv $MERGE_DIR/* $PREMERGE_DIR
	fi

	echo $PREMERGE_DIR

	## merge all fastqs that have same ID +/- "_add", put all fastqs in $PREMERGE_DIR and rename merged files with ID w/o "add"
	echo "Merging fastqs in $PREMERGE_DIR"
	ALL_IDS=`ls -lrt $PREMERGE_DIR | awk '{print $9}'| grep fastq | sed 's/R[0-9].fastq.gz//g'| sed -e 's/_add//g' -e 's/_add1//g' -e 's/_add2//g' | sort | uniq`
	echo $ALL_IDS
	
	ALL_READS="R1 R2"
	for ID in $ALL_IDS
	do
		for READ in $ALL_READS
		do
			echo "Working on $READ of $ID"
			FIRST_FASTQ=$PREMERGE_DIR"/"$ID$READ".fastq.gz"
			ADD_FASTQ=$PREMERGE_DIR"/"$ID"add_"$READ".fastq.gz"
			OUT_ADD_FASTQ=$MERGE_DIR"/"$ID"add_"$READ".fastq.gz"
			OUT_FASTQ=$MERGE_DIR"/"$ID$READ".fastq.gz"
			
			if [ -f $ADD_FASTQ ];then
				
				OUT_FASTQ=$MERGE_DIR"/"$ID$READ".fastq.gz"
				cat $FIRST_FASTQ > $OUT_FASTQ
				cat $ADD_FASTQ >> $OUT_FASTQ
								
				ALL_FASTQ=`ls $ADD_FASTQ $FIRST_FASTQ`
	
				INSIZE=`ls -lrt $FIRST_FASTQ $ADD_FASTQ | awk '{print $5}' | awk -v total=0 '{for(i=1; i < NF+1; i++) total=total+$i; print total}' | tail -1`
				OUTSIZE=`ls -lrt $OUT_FASTQ | awk '{print $5}'`
				ls -lrt $ALL_FASTQ
				ls -lrt $OUT_FASTQ
				echo "Size diff between input and ouput is:"
				echo $INSIZE $OUTSIZE | awk '{print $2-$1}'
				echo ""
				
				OUT_FASTQ=$MERGE_DIR"/"$ID$READ".fastq.gz"
				echo "Making sym link between $ADD_FASTQ and $OUT_FASTQ "
				ln -s $ADD_FASTQ $OUT_ADD_FASTQ
			else
				echo "Making sym link between $FIRST_FASTQ and $OUT_FASTQ "
				ln -s $FIRST_FASTQ $OUT_FASTQ
			fi
		done
	
	done

	### change permissions for Results and temp dirs
	change_perms $MERGE_DIR
	ls -lrt $MERGE_DIR

	echo $MERGE_DIR
fi


if [ $ALIGN == "Y" ];then

	### Run pipeline to align and count merged fastq files. Key file is modified to remove samples with "redo"

	ls -lrt $KEY_FILE
	if [ $INCL_ADD == "N" ];then
		cat $KEY_FILE | awk -F"\t" '{if($4 !~ /_add/) print $0}' > $TEMP_DIR"/temp_keyfile.txt"
		ls -lrt  $TEMP_DIR"/temp_keyfile.txt"
		cat $TEMP_DIR"/temp_keyfile.txt" > $KEY_FILE
		MOVE_KEY="N"
	else 
	
		MOVE_KEY="Y"
	fi
	ls -lrt $KEY_FILE

	echo "sh /idi/moc_ec/MOC/scripts/QSUB_RtS_analysis_launch.sh $MOC_ID $@ -move_key $MOVE_KEY -split N -merge N"
	sh /idi/moc_ec/MOC/scripts/QSUB_RtS_analysis_launch.sh $MOC_ID $@ -move_key $MOVE_KEY -split N -merge N
	#sh /idi/moc_ec/MOC/scripts/MOC_RtS_pipe_v2.sh $MOC_ID $@ -move_key $MOVE_KEY -split N -merge N

fi

