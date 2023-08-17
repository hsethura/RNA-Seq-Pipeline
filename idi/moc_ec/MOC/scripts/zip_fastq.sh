#!/bin/sh

ALL_DIRS=$@
TEMP_DIR="/broad/hptmp/livny/"

mkdir -p $TEMP_DIR

for DIR in $ALL_DIRS
do
	echo $DIR
	
	FASTQ_FILES=`ls -lrt $DIR/*fastq | awk '{print $9}'`

	for FASTQ in $FASTQ_FILES
	do
	
		ZIP_NAME=$FASTQ".gz"
		QSUB_FILE=$TEMP_DIR"/fastq_zip_temp.txt"
		echo "Converting $FASTQ to $ZIP_NAME"
		echo "source /broad/IDP-Dx_work/nirmalya/bash_header" > $QSUB_FILE
		echo "gzip $ZIP_NAME $FASTQ" >> $QSUB_FILE
		echo "chmod 777 $ZIP_NAME" >> $QSUB_FILE
		echo "qsub -l h_vmem=8g -l os=RedHat7 $QSUB_FILE"
		qsub -l h_vmem=8g -l os=RedHat7 $QSUB_FILE
		
		ls -lrt $QSUB_FILE
		
	
	done
done
