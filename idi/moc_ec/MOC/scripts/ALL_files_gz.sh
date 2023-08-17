#!/bin/sh

DIR=$1

ALL_FILES=`ls -lrt $DIR | awk '{print $9}' | grep -v gz | grep -v txt`


for FILE in $ALL_FILES
do
	qsub -l h_rt=24:00:00 -l h_vmem=8g -b Y gzip $DIR$FILE

done



