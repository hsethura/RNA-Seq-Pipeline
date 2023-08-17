#!/bin/sh


DB_FILE=$1
NCBI_DIR=$2



echo "" | sed 1d > $DB_FILE
f=`ls -lrt $NCBI_DIR"/"*".fna" | awk '{print $9}' `

echo "Pulling out data..."

for m in $f
do
	echo $m
	typeset -i FASTA_NUM
	FASTA_NUM=`cat $m | grep ">" | wc -l`
	
	if [ $FASTA_NUM == 1 ];then
		acc=`cat $m | grep ">" | sed 's/>//g' | cut -d"|" -f4`
		name=`cat $m | grep ">" | sed 's/>//g' | cut -d"|" -f5 | cut -d"," -f1 | sed 's/ /_/g'`
		size=`cat $m | grep -v ">" |  tr '\n' ' ' |  sed 's/ //g' | wc -c` 
		
		echo $name"	 "$acc"	"$size | sed 's/	_/ /g' | sed 's/^_//g' >> $DB_FILE
	fi
done


echo $DB_FILE