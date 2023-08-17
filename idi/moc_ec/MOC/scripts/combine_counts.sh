#!/bin/sh


source "/broad/IDP-Dx_storage/MOC/scripts/MOC_functions.sh"

FILE1=$1
FILE2=$2
OUTFILE=$3

ALL_HEADERS=`cat $FILE1 $FILE2 | sed 's/Geneid//g' | head -1 | sort | uniq `

 
cat $FILE1 | awk '{print $1}' > $OUTFILE

for HEADER in $ALL_HEADERS
do
	F1=`FIELD_HEADER $FILE1 $HEADER`
	F2=`FIELD_HEADER $FILE2 $HEADER`

	cat $FILE1 | awk -v i=$F1 '{print $i}' > temp1.txt
	cat $FILE2 | awk -v i=$F2 '{print $i}' > temp2.txt
	
	
	paste temp1.txt temp2.txt | awk '{if(NR==1) print $1; else print $1+$2}' > temp3.txt 

	cat $OUTFILE > temp4.txt
	paste temp4.txt temp3.txt > $OUTFILE

done

echo $OUTFILE