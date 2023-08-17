#!/bin/sh

rem_suff=$1
shift
OUT_DIR=$1
shift


TEMP_DIR="/broad/hptmp/RNASeq_proj/"

mkdir -p $TEMP_DIR

ALL_BAMS=$@


for BAM in $ALL_BAMS
do
	echo $BAM
	ROOT=`basename $BAM $rem_suff`
	echo $ROOT
	F1=$OUT_DIR"/"$ROOT"_R1.fastq"
	F2=$OUT_DIR"/"$ROOT"_R2.fastq"
	QSUB_FILE=$TEMP_DIR"File_"$ROOT"_temp.txt"
	
	echo "source /idi/moc_ec/MOC/scripts/bash_header" > $QSUB_FILE

	echo "sh /broad/IDP-Dx_work/nirmalya/pipeline/beta/shell_scripts/ALN_BAM_to_FASTQ.sh $BAM $ROOT $F1 $F2 1 1 $TEMP_DIR" >> $QSUB_FILE
	qsub -e err.txt -o out.txt -l h_vmem=8g -l h_rt=24:00:00 $QSUB_FILE
	echo $QSUB_FILE
		
done