#!/bin/sh

source /idi/moc_ec/MOC/scripts/bash_header

IN_DIR=$1
OUT_DIR=$2
REF_ACC=$3

#rm -r $OUT_DIR

mkdir -p $OUT_DIR

ALL_DIRS=`ls -lrtd $IN_DIR/* | awk '{print $9}' | grep -v param | grep -v Solexa`



for DIR in $ALL_DIRS
do
	OUT_NAME=`basename $DIR | sed 's/Jonathan_Livny_//g' `
	IN_BAM=`ls -lrt $DIR/*"marked.bam" | awk '{print $9}'`
	OUT_BAM=$OUT_DIR"/"$OUT_NAME"_"$REF_ACC"_u_alnsrtd.bam"
	
	QSUB_FILE="qsub.txt"

	echo "source /idi/moc_ec/MOC/scripts/bash_header" > $QSUB_FILE
	echo "cp $IN_BAM $OUT_BAM" >> $QSUB_FILE

	echo "qsub -e $QSUB_ERR_FILE -o $QSUB_OUT_FILE -l h_rt=24:00:00 -l h_vmem=8g -l os=RedHat7 $QSUB_FILE" 
	qsub -e $QSUB_ERR_FILE -o $QSUB_OUT_FILE -l h_rt=24:00:00 -l h_vmem=8g -l os=RedHat7 $QSUB_FILE
	
	cat $QSUB_FILE
	#ln -s $IN_BAM $OUT_BAM

# 	echo "samtools sort -n -T /broad/hptmp/RNASeq_proj/ -o $OUT_BAM $IN_BAM"
# 	samtools sort -n -T /broad/hptmp/RNASeq_proj/ -o $OUT_BAM $IN_BAM 
# 	

done