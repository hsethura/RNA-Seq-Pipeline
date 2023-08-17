#!/bin/sh

TRIM=$1
shift
DIR=$1
shift

TEMP_DIR="/broad/hptmp/MOC/trimmed_fastqs/"

mkdir -p $TEMP_DIR

INDEXES=$@

echo $INDEXES

for INDEX in $INDEXES
do

	FILES=`ls -lrt $DIR/* | awk '{print $9}' | grep "_"$INDEX"_R1"`
	
	echo $INDEX, $FILES
	
	for FILE in $FILES
	do
		ROOT=`basename $FILE .fastq`
		OUFTILE=$TEMP_DIR"/"$ROOT"_trimmed.fastq"
		echo $FILE, $ROOT, $OUFTILE
		cat $FILE | awk -v trim=$TRIM '{
							if(NR % 4 ==2 || NR % 4 ==0)
							{	
								y=substr($0, 1+trim, 10000)
								print y
							}
							else
								print $0
								
						}' > $OUFTILE
	
		echo "Moving $FILE $OUFTILE"
	
		mv $OUFTILE $FILE
	
	done

done


