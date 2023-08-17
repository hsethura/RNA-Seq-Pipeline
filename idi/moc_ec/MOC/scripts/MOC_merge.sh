#!/bin/sh

PROJID=$1

SCRIPT_PATH="/broad/IDP-Dx_storage/MOC/scripts/"
TEMP_DIR="/broad/hptmp/Dx/"$PROJID"/"
BC_PATH="/broad/IDP-Dx_work/nirmalya/research/BarcodeSplitter/"
ALLBC_FILE=$BC_PATH"RtS_bc_v1.txt"


ls -lrt $ALLBC_FILE

ALL_BC=`cat $ALLBC_FILE | awk '{print $2}' | sort | uniq`

for m in $ALL_BC
do
	echo $m
	reads="R1 R2"
	
	for p in $reads
	do
		INFILES=`ls -lrt $TEMP_DIR* | awk '{print $9}' | grep $m | grep $p  | tr '\n' ' '`
		OUTFILE=$TEMP_DIR"merged_"$m"_"$p".fastq"
		qsub -b y "cat $INFILES > $OUTFILE"
	done

done
