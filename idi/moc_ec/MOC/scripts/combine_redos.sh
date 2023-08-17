#!/bin/sh

OLD_DIR=$1
NEW_DIR=$2

mkdir -p $NEW_DIR

ALL_REDO_ID=`ls -lrt $OLD_DIR | grep fastq | grep redo | awk '{print $9}' | sed -e 's/_redo//g' -e 's/_R1.fastq.gz//g' -e 's/_R2.fastq.gz//g' | sort | uniq -c | sort -k1n | awk '{if($1==2) print $2}'`
ALL_SINGLE_ID=`ls -lrt $OLD_DIR | grep fastq | awk '{print $9}' | sed -e 's/_redo//g' -e 's/R1.fastq.gz//g' -e 's/R2.fastq.gz//g' | sort | uniq -c | sort -k1n | awk '{if($1==2) print $2}'`


for ID in $ALL_SINGLE_ID
do
	echo "single " $ID
	FILE=`ls -lrt $OLD_DIR* | grep $ID | awk '{print $9}'`
	echo $FILE
	mv $FILE $NEW_DIR

done


for ID in $ALL_REDO_ID
do
	
	
	echo "redo "$ID
	
	R1_OUTPUT_FILE=$NEW_DIR"/"$ID"_R1.fastq.gz"
	R2_OUTPUT_FILE=$NEW_DIR"/"$ID"_R2.fastq.gz"

	R1_FILES=`ls -lrt $OLD_DIR* | grep $ID | grep R1 | awk '{print $9}'`
	R2_FILES=`ls -lrt $OLD_DIR* | grep $ID | grep R2 | awk '{print $9}'`

	echo "cat $R1_FILES $R1_OUTPUT_FILE"
	echo "cat $R2_FILES $R2_OUTPUT_FILE"

	cat $R1_FILES > $R1_OUTPUT_FILE
	cat $R2_FILES > $R2_OUTPUT_FILE
done