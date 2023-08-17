#!/bin/sh


DIR=$1
TARGET_DIR=$2


ALL_FILES=`ls $DIR`
 
for FILE in $ALL_FILES
do


	echo "qsub -l h_rt=24:00:00 -l h_vmem=8g -l os=RedHat7 -b Y mv $DIR"/"$FILE $TARGET_DIR"
	qsub -l h_rt=24:00:00 -l h_vmem=8g -l os=RedHat7 -b Y mv $DIR"/"$FILE $TARGET_DIR

done

