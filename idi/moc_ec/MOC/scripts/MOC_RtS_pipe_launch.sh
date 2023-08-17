#!/bin/sh

MOC_ID=$1

source /idi/moc_ec/MOC/scripts/bash_header

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

### get options from command line
CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/Universal_config.yaml" 1 $@`

### Get paths to dirs scripts from config file
read_config $CONFIG_FILE 

### identify options to pass to this script (as opposed to the pipeline)
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

### run path_suff function to set RESPATH_SUFF.
#	If -moc_id N included in command line, do not moc_ID to RESPATH_SUFF.  
#	If -user_id included with Y or no string add USID to RESPATH_SUFF.  
#	If -user_id followed by userID, add userID to RESPATH_SUFF.
paths_and_headers $MOC_ID $@

path_suff $@ "-moc_id "$MOC_ID_OPT" -user_id "$USER_ID
echo $RESPATH_SUFF


### Set additional options and paths
Q_HEAD="MOC_ID"
USID=`USID`
RAW_SYM_PATH=`extract_option -symlink $RAWSYM_PATH 1 $@`
RAW_SEQ_PATH=`extract_option -raw_seq_path $RAWSYM_PATH 1 $@`

RESULTS_DIR=$RESULTS_PATH"/"$RESPATH_SUFF"/"
TEMP_DIR=$TEMP_PATH"/"$RESPATH_SUFF"/"
BAM_DIR=$BAM_PATH"/"$RESPATH_SUFF"/"
PARSED_REFDIR=$TEMP_DIR"/parsed_ref/"
SYM_DIR=$RAW_SYM_PATH"/"$MOC_ID"/"

mkdir -p $RESULTS_DIR
mkdir -p $TEMP_DIR

echo $KEY_SHEET


NO_SPLIT=`echo $PIPE_OPTIONS | grep no_split | wc -l`
MOC_SPLIT_DIR=$INDXSPLIT_PATH"/"$RESPATH_SUFF"/"

echo "RAW_SEQ_PATH: "$RAW_SEQ_PATH
echo "RAW_SYM_PATH: "$RAW_SYM_PATH
echo "Results_dir: "$RESULTS_DIR
echo "Temp_dir: "$TEMP_DIR
echo "MOC_SPLIT_DIR: "$MOC_SPLIT_DIR


############## Set and assign RtS pipeline options ###############=
### Set pipeline options
DEFAULT_PIPE_OPTIONS=" --ADD3 30 --ADD5 20 --remove_splitted --gzip_merged --do_patho --MOC_id_ref $MOC_ID"
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



############## Import google sheet databases to server ###############
if [ $IMPORT_GS == "Y" ];then
	echo "Running $GSIMPORT_SCRIPT to import google sheet databases to server"
	sh $GSIMPORT_SCRIPT
	echo "sh $GSIMPORT_SCRIPT"
fi

############# Get GID from MOC_DB or use GID included as command line option #######################
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

############## Move key file to server ###############
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
  
############## Symlink fastqs for this MOC from SEQ_DIR to MOC_SPLIT_DIR ###############

if [ $NO_SPLIT != 1 ];then

	############## determining ALL_SMOCINDEX_SETS ###############

	echo "Determining ALL_SMOCINDEX_SETS for $MOC_ID ..."
	all_smocindex_set $MOC_ID $@
	echo "ALL_SMOCINDEX_SETS:" $ALL_SMOCINDEX_SETS

	########################################################

	############## make symlink in MOC_SPLIT_DIR to MOC fastq files in SEQ_DIR ###############
	rm -r $MOC_SPLIT_DIR
	mkdir -p $MOC_SPLIT_DIR

	MISSING_FASTQ="N"

	for SMOCINDEX_SET in $ALL_SMOCINDEX_SETS
	do
		SMOC_ID=`echo $SMOCINDEX_SET | cut -d";" -f1 | sed 's/-//g'`
		INDEX1=`echo $SMOCINDEX_SET | cut -d";" -f2`
		INDEX2=`echo $SMOCINDEX_SET | cut -d";" -f3`

		SEQ_DIR=$RAW_SEQ_PATH$SMOC_ID"/"
	
		echo "Making sym links to demultiplexed fastqs from $SEQ_DIR to $MOC_SPLIT_DIR..."

		### if no index2 is in SMOCINDEX_SET
		if [ $USE_P5 == "Y" ];then
			ALL_FASTQS=`ls -lrt $SEQ_DIR*unmapped.*.fastq* | grep -v barcode | grep $INDEX1 | grep $INDEX2 | awk '{print $NF}'`
			NUM_FASTQS=`ls -lrt $SEQ_DIR*unmapped.*.fastq* | grep -v barcode | grep $INDEX1 | grep $INDEX2 | awk '{print $NF}' | wc -l`
			echo $NUM_FASTQS, $INDEX1, $INDEX2
			if [ $NUM_FASTQS == 0 ];then
	
				echo "Missing fastq files for "$INDEX1" in "$SEQ_DIR"!"
			fi
		else 
			ALL_FASTQS=`ls -lrt $SEQ_DIR*unmapped.*.fastq* | grep -v barcode | grep $INDEX1 | awk '{print $NF}'`
			NUM_FASTQS=`ls -lrt $SEQ_DIR*unmapped.*.fastq* | grep -v barcode | grep $INDEX1 | awk '{print $NF}' | wc -l`
			echo $NUM_FASTQS, $INDEX2
			if [ $NUM_FASTQS == 0 ];then
	
				echo "Missing fastq files for "$INDEX1" and "$INDEX2" in "$SEQ_DIR"!"
				MISSING_FASTQ="Y"
			fi
		fi
		
		for FASTQ in $ALL_FASTQS
		do
			FILE=`echo $FASTQ | rev | cut -d"/" -f1 | rev`
			ln -s $FASTQ $MOC_SPLIT_DIR$FILE 2>/dev/null
		done
	done

	if [ $MISSING_FASTQ == "Y" ];then
	
		echo "FASTQs missing!"
		exit
	
	fi

	echo $MOC_SPLIT_DIR

	###########################################################################

	typeset -i INDEX_NUM
	INDEX_NUM=`ls -lrt $MOC_SPLIT_DIR | awk '{print $NF}' | grep fastq | sed 's/.unmapped./ /g' | cut -d" " -f1 | rev | cut -d"/" -f1 | cut -d"." -f1 | rev | sed 's/_/ /g' | awk '{print NF}' | sort | uniq`

	if [ $INDEX_NUM == 2 ];then
		PIPE_OPTIONS=$PIPE_OPTIONS" --use_p5"
	fi
	if [ $INDEX_NUM != 1 ] && [ $INDEX_NUM != 2 ];then
		echo "Num indexes does not equal 1 or 2 (INDEX_NUM="$INDEX_NUM")"
		exit 
	fi
fi

PIPE_OPTIONS=":--config_file $CONFIG_FILE --key_path $KEY_FILE --raw_seq_path $MOC_SPLIT_DIR $PIPE_OPTIONS:"
echo $PIPE_OPTIONS
echo $MOC_ID
echo $RtS_ANPIPE
echo $LAUNCH_OPTIONS
 
############## Set command line and launch analysis pipeline ##############
if [ $NO_PIPE != "Y" ];then
	echo "Launching analysis pipeline..."
	if [ $QSUB == "Y" ];then
		UGER_DIR=$RESULTS_PATH"/UGER/"

		echo "UGER_dir: "$UGER_DIR

		mkdir -p $UGER_DIR

		TIME_STAMP=`date +%s`

		QSUB_ERR_FILE=$UGER_DIR"/"$MOC_ID"_"$TIME_STAMP"_qsub_err.txt"
		QSUB_OUT_FILE=$UGER_DIR"/"$MOC_ID"_"$TIME_STAMP"_qsub_out.txt"

		QSUB_FILE=~/$MOC_ID"_qsub.txt"

		echo "source /idi/moc_ec/MOC/scripts/bash_header" > $QSUB_FILE
		echo "sh $RtS_ANPIPE -c $CONFIG_FILE --key_path $KEY_FILE --raw_seq_path $MOC_SPLIT_DIR $PIPE_OPTIONS" >> $QSUB_FILE

		echo "qsub -e $QSUB_ERR_FILE -o $QSUB_OUT_FILE -l h_rt=24:00:00 -l h_vmem=8g -l os=RedHat7 $QSUB_FILE" 
		qsub -e $QSUB_ERR_FILE -o $QSUB_OUT_FILE -l h_rt=24:00:00 -l h_vmem=8g -l os=RedHat7 $QSUB_FILE

		echo "Outfile and errfile directory is:"
		echo $UGER_DIR
		
	else
		echo "$RtS_ANPIPE $MOC_ID $LAUNCH_OPTIONS $PIPE_OPTIONS" 
		sh $RtS_ANPIPE $MOC_ID $LAUNCH_OPTIONS $PIPE_OPTIONS
	fi 
fi
##########################################################################


