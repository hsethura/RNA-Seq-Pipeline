#!/bin/sh


DIR1=$1
DIR2=$2
OUTDIR=$3



ALL_FASTQ=`ls -lrt $DIR1 $DIR2 | awk '{print $NF}' | grep fastq | sort | uniq`


for FASTQ in $ALL_FASTQ
do

	cat $DIR1"/"$FASTQ $DIR2"/"$FASTQ > $OUTDIR"/"$FASTQ
	echo ""
	ls -lrt $OUTDIR"/"$FASTQ
	echo ""
	SIZE1=`ls -lrt $DIR1/* | grep $FASTQ | awk '{print $5}'`
	SIZE2=`ls -lrt $DIR2/* | grep $FASTQ | awk '{print $5}'`
	OUTSIZE=`ls -lrt $OUTDIR/* | grep $FASTQ | awk '{print $5}'`
	echo $SIZE1 $SIZE2 | awk '{print $1+$2}'
	echo $OUTSIZE
done


