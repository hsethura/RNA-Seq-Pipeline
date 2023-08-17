#!/bin/sh

SMOCID=$1 

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

### determining paths and headers 
### default config file is /idi/moc_ec//MOC/config_files/Universal_config.yaml
paths_and_headers $SMOCID $@

echo $USID

RAW_SEQ_PATH=`extract_option -raw_seq_path $SEQ_PATH 1 $@`
RAWSEQ_SYM_PATH=`extract_option -raw_sym_path $RAWSYM_PATH 1 $@`
MOC_ID=`extract_option -moc_id - 1 $@`

SEQ_PATH=$RAW_SEQ_PATH
SYM_PATH=$RAWSEQ_SYM_PATH

echo $SEQ_PATH

SMOCDIR_ID=`echo $SMOCID | sed 's/-//g'`

SEQ_DIR=$SEQ_PATH"/"$SMOCDIR_ID"/"
SYM_DIR=$SYM_PATH"/"$SMOCDIR_ID"/"

echo "Making "$SEQ_DIR
echo "Making "$SYM_DIR
mkdir -p $SEQ_DIR
mkdir -p $SYM_DIR

chmod 777 $SEQ_DIR

cd $SEQ_DIR

WGET_PATH=$SEQ_DIR"/get.broadinstitute.org/pkgs/"

if [ $USER != "N" ] && [ $PSSWD != "N" ];then
	MOVE_DATA="Y" 
fi	

echo $MOVE_DATA, $USER, $PSSWD

if [ $MOVE_DATA == "Y" ];then
	if [ $USER == "N" ] || [ $PSSWD == "N" ];then
		echo "Missing userID and/or psswd"
		exit 
	fi
	
	echo "wget --tries=10 --continue --mirror --user $USER --password $PSSWD --no-check-certificate https://get.broadinstitute.org/pkgs/$USER/"

	wget --tries=10 --continue --mirror --user $USER --password $PSSWD --no-check-certificate https://get.broadinstitute.org/pkgs/$USER/
	 
fi

	
echo "Making symlink to files in "$WGET_PATH" in "$SYM_DIR
echo "ln -s $WGET_PATH/*/* $SYM_DIR"
ln -s $WGET_PATH/*/* $SYM_DIR 2>/dev/null
ls -lrt $SYM_DIR
echo $SYM_DIR
	 


############## Run metrics script ###############=

if [ $METRICS == "Y" ];then
	echo "Running metrics script"
	echo "sh $WU_SCRIPT $SYM_DIR $SMOCID"
	sh $WU_SCRIPT $SYM_DIR $SMOCID
fi
########################################################


### change permissions for Results and temp dirs
change_perms $INDEX_SPLIT_OUTDIR $SEQ_DIR $WGET_DIR
