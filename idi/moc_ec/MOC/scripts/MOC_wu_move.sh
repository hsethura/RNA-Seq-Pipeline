#!/bin/sh

MOCS_ID=$1 
FC=$2

# get path of the current file. if the file path is relative, convert it to absolute path
file_path="${BASH_SOURCE[0]}"
if [[ $file_path != /* ]]; then
  file_path="$PWD/${BASH_SOURCE[0]}"
fi

# get parent directory
scripts_dir="$(dirname $file_path)"

# source /idi/moc_ec/MOC/scripts/bash_header
source "$scripts_dir/bash_header"

### source all functions 
# source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"
source "$scripts_dir/MOC_functions.sh"

### determining paths and headers 
### default config file is /idi/moc_ec//MOC/config_files/Universal_config.yaml
paths_and_headers $MOCS_ID $@

MOC_ID=`extract_option -moc_id N 1 $@`
DATA_TYPE=`extract_option -data_type RtS 1 $@`

# -symlink: set path to directory named $MOCS_ID containing symlinks to raw data (default $RAWSYM_PATH from config)
# -raw_seq_path: set path to directory $MOCS_ID containing raw data (default $SEQ_PATH from config)

echo $SEQ_PATH
echo $RAWSYM_PATH
echo $WU_PATH

echo $DATA_TYPE

if [ $DATA_TYPE != "RtS" ];then
	SEQ_PATH=`echo $SEQ_PATH | sed 's/RtS/'$DATA_TYPE'/g'`
fi


echo $RAWSYM_PATH
MOCS_PATH_ID=`echo $MOCS_ID | sed 's/-//g'`
SEQ_DIR=$SEQ_PATH"/"$MOCS_PATH_ID"/"
SYM_DIR=$RAWSYM_PATH"/"$MOCS_PATH_ID"/"
TEMP_DIR=$TEMP_PATH"/"$MOCS_PATH_ID"/"

echo "Making "$SEQ_DIR
echo "Making "$TEMP_DIR
mkdir -p $SEQ_DIR
mkdir -p $TEMP_DIR

echo $MOCS_PATH_ID

SEG_FILE=$TEMP_DIR$MOCS_PATH_ID"SGE_file_0.txt"

echo "" | sed 1d > $SEG_FILE

ls -lrt $SEG_FILE

echo $MOVE_DATA


if [ $MOVE_DATA == "Y" ];then

	cd $SEQ_DIR
	
	WU_SUBDIR=`ls -lrt $WU_PATH"/"$FC"/" | tail -1 | awk '{print $9}'`
	WU_DIR=$WU_PATH"/"$FC"/"$WU_SUBDIR"/"
	ALL_FILES=`ls -lrt $WU_DIR*/* | awk '{print $9}'`

	batch_size=1000
	i=1
	for FILE in $ALL_FILES
	do
		echo "qsub -e $SEQ_DIR"err.txt" -o $SEQ_DIR"out.txt" -l h_rt=24:00:00 -l h_vmem=8g -l os=RedHat7 -b Y cp $FILE $SEQ_DIR"
		qsub -e $SEQ_DIR"err.txt" -o $SEQ_DIR"out.txt" -l h_rt=24:00:00 -l h_vmem=8g -l os=RedHat7 -b Y cp $FILE $SEQ_DIR >> $SEG_FILE	

		if ((i % batch_size == 0)); then
			echo "Jobs submitted - waiting for them to complete...."
			echo $SEG_FILE

			# testing that all jobs launched are finished
			SGE_test $SEG_FILE	
			qstat

			seg_file_ix=$((i / batch_size))
			SEG_FILE=$TEMP_DIR$MOCS_PATH_ID"SGE_file_"$seg_file_ix".txt"

			echo "" | sed 1d > $SEG_FILE

			ls -lrt $SEG_FILE
		fi

		i=$((i+1))
	done

	echo "Jobs submitted - waiting for them to complete...."
	echo $SEG_FILE

	# testing that all jobs launched are finished
	SGE_test $SEG_FILE	
	qstat
	echo "All done copying demultiplexed data from $WU_DIR to $SEQ_DIR"


fi


# echo "Removing "$SYM_DIR
# rm -r $SYM_DIR
# echo "Making "$SYM_DIR
# mkdir -p $SYM_DIR

# echo "Making symlink to files in "$SEQ_DIR" in "$SYM_DIR
# echo "ln -s $SEQ_DIR* $SYM_DIR"
# ln -s $SEQ_DIR*  $SYM_DIR 2>/dev/null
# ls -lrt $SYM_DIR
# echo $SYM_DIR


# ### change permissions for Results and temp dirs
# change_perms $SYM_DIR $SEQ_DIR $TEMP_DIR

# ############## Run WB move script ###############=

# if [ $WB_MOVE == "Y" ] && [ $DATA_TYPE == "RtS" ];then
# 	echo "Running WB moving script"
# 	echo "sh $WBMOV_SCRIPT -limit $MOCS_ID"
# 	sh $WBMOV_SCRIPT -limit $MOCS_ID
# fi
# ########################################################

# ############## Run metrics script ###############=

# if [ $METRICS == "Y" ];then
# 	echo "Running metrics script"
# 	echo "sh $WU_SCRIPT $SYM_DIR $MOCS_ID"
# 	sh $WU_SCRIPT $SYM_DIR $MOCS_ID
# fi
# ########################################################
