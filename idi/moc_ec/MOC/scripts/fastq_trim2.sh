#!/bin/sh

source /broad/software/scripts/useuse
source /idi/moc_ec/MOC/scripts/bash_header

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

TRIM=`extract_option -trim Y 1 $@`

IN_DIR=$1
TRIM_DIR=$2
TRIM_1=$3
TRIM_2=$4

 
SCRIPTS_DIR="/idi/moc_ec/MOC/scripts/"
TEMP_DIR="/broad/hptmp/MOC/"

SEG_FILE=$TEMP_DIR$"trim_file.txt"
TEMP_FILE=$TEMP_DIR"temp.txt"

echo "" | sed 1d >  $SEG_FILE

ALL_R1=`ls -lrt $IN_DIR/*1.fastq* | grep -v "ndeterm" | awk '{print $9}'`
ALL_R2=`ls -lrt $IN_DIR/*2.fastq* | grep -v "ndeterm" | awk '{print $9}'`

echo $ALL_R1
echo $ALL_R2

mkdir -p $TRIM_DIR


TEMP_FILE=$TEMP_DIR"temp.txt"

if [ $TRIM == "Y" ];then
	echo "Sending all trim jobs to UGER..."


		for FILE in $ALL_R1
		do
			TRIM=$TRIM_1
			IN_FILE_NAME=`basename $FILE`
			TRIM_FILE_NAME=`basename $FILE .gz`
	
			INFILE=$IN_DIR"/"$IN_FILE_NAME
			TRIM_FILE=$TRIM_DIR"/"$TRIM_FILE_NAME
			
			echo sh $SCRIPTS_DIR"/trim_fastq_BD.sh" $INFILE $TRIM_FILE $TRIM > $TEMP_FILE
			cat $TEMP_FILE
			qsub $TEMP_FILE >> $SEG_FILE					
		done

		for FILE in $ALL_R2
		do
			TRIM=$TRIM_2
			IN_FILE_NAME=`basename $FILE`
			TRIM_FILE_NAME=`basename $FILE .gz`
	
			INFILE=$IN_DIR"/"$IN_FILE_NAME
			TRIM_FILE=$TRIM_DIR"/"$TRIM_FILE_NAME
			
			echo sh $SCRIPTS_DIR"/trim_fastq_BD.sh" $INFILE $TRIM_FILE $TRIM > $TEMP_FILE
			cat $TEMP_FILE
			qsub $TEMP_FILE >> $SEG_FILE					
		done

	
		echo "Jobs submitted - waiting for them to complete...."

		echo $SEG_FILE
		qstat

		###testing that all jobs launched are finished
		SGE_test $SEG_FILE	
		echo "All Jobs completed...."
fi

ls -lrt $TRIM_DIR
