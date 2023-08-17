#!/bin/sh


DIR=$1

BC_FILE="/broad/IDP-Dx_storage/MOC/files/rts_bcs.txt"

ALL_INDEXES=`ls $DIR | grep fastq | grep unmapped | cut -d"." -f3 | sort | uniq`

for INDEX in $ALL_INDEXES
do

	FILE=`ls -lrt $DIR | grep $INDEX | grep unmapped.1.fastq | head -1 | awk '{print $9}'`

	zcat $DIR"/"$FILE | awk '{if(NR%4 ==2) print $1}' | head -100000 | cut -b1-9 | sort | uniq -c | awk '{print $1/100000*100, $2}' | sort | tail -200 > temp.txt

	while read line
	do

		bc=`echo $line | awk '{print $2}'`
		BC_FREQ=`echo $line | awk '{print $1}'`
		BC_NUM=`cat $BC_FILE | grep -w $bc | awk '{print $1}'`
		if [ ! -z $BC_NUM  ];then
			printf "%s\t%s\t%s\n" $BC_FREQ $bc $BC_NUM 
		fi
	done < temp.txt | sort -k3 > temp2.txt

	total=`cat temp2.txt | wc -l`
	total_pcnt=`awk -v total=0 '{total=total+$1; print total}' temp2.txt | tail -1`
	echo "$total_pcnt% reads assigned to $total bcs for $INDEX"
	cat temp2.txt
	echo ""

done