#!/bin/sh

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

### determining paths and headers 
paths_and_headers $ID $@


DIR=$1

PREMERGE_DIR=$DIR"/pre_redocombined/"

mkdir -p $PREMERGE_DIR

echo $PREMERGE_DIR

mv $DIR* $PREMERGE_DIR

ALL_IDS=`ls -lrt $PREMERGE_DIR | awk '{print $9}'| grep fastq | sed 's/R[0-9].fastq.gz//g'| sed 's/_redo//g' | sort | uniq`

echo $ALL_IDS

ALL_READS="R1 R2"
for ID in $ALL_IDS
do
	for READ in $ALL_READS
	do
		
		echo "Working on $READ of $ID"
		OUT_FASTQ=$DIR"/"$ID$READ".fastq.gz"
		ALL_FASTQ=`ls -lrt $PREMERGE_DIR* | awk '{print $9}' | grep $ID | grep "_"$READ`
		
		ls -lrt $ALL_FASTQ
		cat $ALL_FASTQ > $OUT_FASTQ

		
		ls -lrt $OUT_FASTQ

		INSIZE=`ls -lrt $ALL_FASTQ | awk '{print $5}' | awk -v total=0 '{for(i=1; i < NF+1; i++) total=total+$i; print total}' | tail -1`
		OUTSIZE=`ls -lrt $OUT_FASTQ | awk '{print $5}'`
		echo "Size diff between input and ouput is:"
		echo $INSIZE $OUTSIZE | awk '{print $2-$1}'
		echo ""

	done
	
done

change_perms $DIR

echo $DIR


