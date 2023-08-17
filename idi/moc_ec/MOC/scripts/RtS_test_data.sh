#!/bin/sh


SOURCE_DIR=$1
TARGET_PATH=$2
NUM_READS=$3

SMOC_ID=`basename $SOURCE_DIR | sed 's/-//g'`
echo $SMOC_ID

TARGET_DIR=$TARGET_PATH"/"$SMOC_ID"/"

mkdir -p $TARGET_DIR

ls -lrt $TARGET_DIR

ALL_FILES=`ls -lrt $SOURCE_DIR | grep fastq | awk '{print $NF}'`
echo $ALL_FILES

for FILE in $ALL_FILES
do

	NAME=`basename $FILE .gz`
	echo $NAME
	FASTQ_NUM=`echo $NUM_READS | awk '{print $1*4}'`
	OUTFILE=$TARGET_DIR"/"$NAME	
	
	echo $FASTQ_NUM
	
	zcat $FILE | head -$FASTQ_NUM > $OUTFILE



done

echo $TARGET_DIR