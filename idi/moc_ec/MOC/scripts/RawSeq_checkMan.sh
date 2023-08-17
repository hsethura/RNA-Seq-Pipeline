#!/bin/sh


### source all functions 
source "/broad/IDP-Dx_storage/MOC/scripts/MOC_functions.sh"


DIR=$1

cd $DIR

MAN_FILE=$DIR"/MANIFEST"
VAL_FILE=$DIR"/Download_validation.txt"

ALL_SEQ_FILES=`ls -lrt $DIR | awk '{print $9}' | grep fastq.gz`

for SEQ_FILES in $ALL_SEQ_FILES
do
	NUM_MATCH=`cat $MAN_FILE | grep $SEQ_FILES | wc -l`
	
	
	if [ $NUM_MATCH != 0 ];then
		echo $SEQ_FILES" DOWNLOADED IN "$DIR
	else
		echo $SEQ_FILES" NOT DOWNLOADED IN "$DIR
	fi
done > $VAL_FILE

echo $VAL_FILE
