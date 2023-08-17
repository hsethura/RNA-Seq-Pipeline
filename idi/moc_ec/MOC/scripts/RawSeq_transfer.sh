#!/bin/sh


ORIGIN_PATH=$1
shift
TARGET_PATH=$1
shift


for SMOC_ID in $@ 
do
	ORIGIN_DIR=$ORIGIN_PATH"/"$SMOC_ID"/get.broadinstitute.org/pkgs/"
	TARGET_DIR=$TARGET_PATH"/"$SMOC_ID"/"
	echo "mkdir -p $TARGET_DIR"
	mkdir -p $TARGET_DIR

	#ls -lrt $ORIGIN_DIR*

	echo "qsub -l h_rt=24:00:00 -l os=RedHat7 -b Y mv $ORIGIN_DIR* $TARGET_DIR"
	qsub -l h_rt=24:00:00 -l os=RedHat7 -b Y mv $ORIGIN_DIR* $TARGET_DIR

done