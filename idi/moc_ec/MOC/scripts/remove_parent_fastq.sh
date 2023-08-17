#!/bin/sh


### source all functions 
source "/broad/IDP-Dx_storage/MOC/scripts/MOC_functions.sh"

RAWSEQ_PATH=`extract_option -raw_path "/broad/IDP-Dx_storage/MOC/RawSeqData/" 1 $@`
TEMP_PATH="/broad/hptmp/MOC/"

mkdir -p $TEMP_PATH

for DIR in $@

do
		
	RM_FILES=`ls -lrtd $RAWSEQ_PATH$DIR/get*/*/*/* | grep fastq | grep -v unmapped | grep -v counts | grep -v unmatched | awk '{print $NF}'`
	echo $RM_FILES
	
	for m in $RM_FILES
	do

		echo "qsub -b y rm $m"
		qsub -b y rm $m
 
	done 
done