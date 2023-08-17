#!/bin/sh

Q=$1
shift

Q_VAL=`echo $Q | cut -d"," -f2`
Q_HEAD=`echo $Q | cut -d"," -f1`

SEQ_PATH="/broad/IDP-Dx_storage/MOC/RawSeqData/"


ALL_SMOC_ID=`sh /broad/IDP-Dx_storage/MOC/scripts/MOC_DB_find.sh $Q_HEAD,$Q_VAL SMOC_ID | sed 1d | awk '{print $2}' | sed 's/,/ /g'`
ALL_INDEX1=`sh /broad/IDP-Dx_storage/MOC/scripts/MOC_DB_find.sh $Q_HEAD,$Q_VAL Index1_Seq | sed 1d | awk '{print $2}' | sed 's/,/ /g'`

for SMOC_ID in $ALL_SMOC_ID
do
	
	NUM_FC=`sh /broad/IDP-Dx_storage/MOC/scripts/MOC_DB_find.sh SMOC_ID,$SMOC_ID Num_of_flowcells | sed 1d | awk '{print $2}' | sed 's/,/ /g'`
	NUM_LANES=`sh /broad/IDP-Dx_storage/MOC/scripts/MOC_DB_find.sh SMOC_ID,$SMOC_ID Num_of_lanes | sed 1d | awk '{print $2}' | sed 's/,/ /g'`
	
	SMOC_SEQ_ID=`echo $SMOC_ID | sed 's/-//g'`
	SEQ_DIR=$SEQ_PATH"/"$SMOC_SEQ_ID"/"
	
	EXPECTED_FILES=`echo $NUM_FC $NUM_LANES | awk '{print $1*$2}'`
	
	FLAG="N"
	
	echo "Looking for missing files in $SEQ_DIR"

	for INDEX1 in $ALL_INDEX1
	do
		NUM_READ1=`ls -lrt $SEQ_DIR*1.fastq.gz | grep $INDEX1 | grep -v barcode | wc -l`
		NUM_READ2=`ls -lrt $SEQ_DIR*2.fastq.gz | grep $INDEX1 | grep -v barcode | wc -l`
		if [ $NUM_READ1 != $EXPECTED_FILES ];then
			echo "Read1 files for index $INDEX1 missing!"
			FLAG="Y"
		fi 
		if [ $NUM_READ2 != $EXPECTED_FILES ];then
			echo "Read2 files for index $INDEX1 missing!"
			FLAG="Y"
		fi 
	done
	
	if [ $FLAG == "N" ];then
	
		echo "All files are there!"
		
	fi

done
