#!/bin/sh

ALL_DIRS=$1
DAYS_OLD=$2

typeset -i FLAG_DIR
FLAG_DIR=`ls -lrt $ALL_DIRS | awk '{if($1 ~ /^d/) print $0}' | wc -l`

if [ $FLAG_DIR -gt 0 ];then

	ALL_BAMS=`find $ALL_DIRS/*/* -type f -ctime +$DAYS_OLD  -ls | awk '{print $11}'`
	rm $ALL_BAMS  
	#find $ALL_DIRS/*/* -mtime +$DAYS_OLD | awk '{if($1 ~ /bam$/) system ("ls -lrt "$1)}'
fi






