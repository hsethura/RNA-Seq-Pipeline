#!/bin/sh


PROJ1_DIR=$1"/"
PROJ2_DIR=$2"/"
TAG=$3
OUT_DIR=$4


mkdir -p $OUT_DIR

ALL_KEY="bcLogfile.txt:c+1 counts.tsv:p-1 fpkm.tsv:p-1 corr.txt:c+1 AlignmentSummaryMetrics.txt:c-h metrics.txt:c-h gv.txt:c abs.txt:c info.txt:c+1"

for KEY in $ALL_KEY
do

	SUFF=`echo $KEY | cut -d":" -f1`
	ACT=`echo $KEY | cut -d":" -f2`
	
	FILE1=`ls -lrt $PROJ1_DIR/*$SUFF | awk '{print $9}'`
	FILE2=`ls -lrt $PROJ2_DIR/*$SUFF | awk '{print $9}'`
	OUTFILE=$OUT_DIR"/"$TAG"_"$SUFF


	if [ $ACT == "c" ];then
		cat $FILE1 $FILE2 > $OUTFILE
	fi
	if [ $ACT == "c+1" ];then
		cat $FILE1 > $OUTFILE && echo "" >> $OUTFILE && cat $FILE2 >> $OUTFILE
	fi
	if [ $ACT == "c-h" ];then
		cat $FILE1 > $OUTFILE && cat $FILE2 | sed 1d >> $OUTFILE
	fi
	if [ $ACT == "p-1" ];then
		cat $FILE2 | awk '{for(i=2; i < NF; i++)printf ("%s\t"), $i; print $NF}'> temp.txt
		paste $FILE1 temp.txt > $OUTFILE
	fi
	if [ $ACT == "p" ];then
		paste $FILE1 $FILE2 > $OUTFILE
	fi

done 