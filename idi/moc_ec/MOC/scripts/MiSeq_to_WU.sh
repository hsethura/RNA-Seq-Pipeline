#!/bin/sh

FASTQ_DIR=$1
SYM_DIR=$2
PROJ_ID=$3

mkdir -p $SYM_DIR"/"$PROJ_ID

ALL_FASTQ=`ls -lrt $FASTQ_DIR/*.fastq.gz | grep -v "ndeterm" | awk '{print $9}'`


echo $ALL_FASTQ

for FASTQ in $ALL_FASTQ
do
	
	echo $FASTQ
	POOL_ID=`echo $FASTQ | rev | cut -d"/" -f1 | rev | cut -d"_" -f1 | sed 's/pool-//g'`
	READ=`echo $FASTQ | rev | cut -d"/" -f1 | rev | awk '{
	
								y=split($1, ar, "_")
								
								for(i=1; i<y+1; i++)
								{
									if(ar[i] ~ /^R/)
										print substr(ar[i],2,1)
								}
					
								}'`


	
	echo ""
	echo $FASTQ
	echo $READ
	
	LN_FILE=$SYM_DIR"/"$PROJ_ID"/MiSeq-MOCS-NNNNNNNN-NNNNNNNN.1."$PROJ_ID"p"$POOL_ID"_"$PROJ_ID"p"$POOL_ID",.unmapped."$READ".fastq.gz"

	echo $POOL_ID, $READ
	echo "ln -s $FASTQ $LN_FILE"
	ln -s $FASTQ $LN_FILE

done
