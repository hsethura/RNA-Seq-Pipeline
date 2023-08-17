#!/bin/sh

BAMMET_OUTDIR=$1"/"
OUT_ALL_MET=$2


typeset -i i
i=1



f=`ls $BAMMET_OUTDIR | grep "_AlignmentSummaryMetrics" | grep -v ALL`
for m in $f
do
	ID=`echo $m | sed 's/_AlignmentSummaryMetrics.txt//g'`
	if [ $i -eq 1 ];then
	
		printf "\t" > $OUT_ALL_MET
		cat $BAMMET_OUTDIR"/"$m | grep "CATEGORY" | head -1 >> $OUT_ALL_MET
	fi
	
	cat $BAMMET_OUTDIR"/"$m | awk -v ID=$ID '{if(NR>7 && $0 != "" )print ID"\t"$0}' >> $OUT_ALL_MET
		
	i=`expr $i + 1`
done

echo $OUT_ALL_MET
   
