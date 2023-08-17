#!/bin/sh


### source all functions 
source "/broad/IDP-Dx_storage/MOC/scripts/MOC_functions.sh"

RAWSEQ_PATH=`extract_option -raw_path "/broad/IDP-Dx_storage/MOC/RawSeqData/" 1 $@`
TEMP_PATH="/broad/hptmp/MOC/"

mkdir -p $TEMP_PATH

for DIR in $@

do
	
	SEQ_DIR=$RAWSEQ_PATH$DIR"/"
	RM_FILES=`ls -lrt $SEQ_DIR* | grep -e unmapped -e counts -e unmatched | awk '{print $9}'`

	echo $RM_FILES

	for m in $RM_FILES
	do

		qsub -b y rm $m

	done 
done