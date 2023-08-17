#!/bin/sh


### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

### determining paths and headers 
### default config file is /idi/moc_ec/MOC/config_files/Universal_config.yaml
paths_and_headers $SMOCID $@


MOVE_KEY=`extract_option -move_key Y 1 $SCRIPT_OPTIONS`
IMPORT_GS=`extract_option -import_gs N 1 $SCRIPT_OPTIONS`
RAW_SEQ_PATH=`extract_option -raw_seq_path $SEQ_PATH 1 $@`
CHECKFILES=`extract_option -checkfiles N 1 $@`
MOC_ID=`extract_option -moc_id - 1 $@`

SEQ_PATH=$RAW_SEQ_PATH

ALL_SMOCID=`echo $@ | awk '{
								for(i=1; i < NF+1; i++)
								{
									if($i ~ /^-/)
										i++
									else
										print $i
								}
							}'`

for SMOCID in $ALL_SMOCID
do

	SMOCDIR_ID=`echo $SMOCID | sed 's/-//g'`

	SEQ_DIR=$SEQ_PATH"/"$SMOCDIR_ID"/"
	RAWSYM_DIR=$RAWSYM_PATH"/"$SMOCDIR_ID"/"

	echo "Making "$SEQ_DIR
	echo "Making "$RAWSYM_DIR

	mkdir -p $RAWSYM_DIR
	chmod 777 $RAWSYM_DIR
	cd $SEQ_DIR

	ALL_WGET_DIRS=$SEQ_DIR"/get.broadinstitute.org/pkgs/*"

	echo $ALL_WGET_DIRS

	for WGET_DIR in $ALL_WGET_DIRS
	do
			echo "Making symlink to files in "$WGET_DIR" in "$RAWSYM_DIR
			echo "ln -s $WGET_DIR/* $RAWSYM_DIR"
			ln -s $WGET_DIR/* $RAWSYM_DIR 2>/dev/null
	done
 
done
