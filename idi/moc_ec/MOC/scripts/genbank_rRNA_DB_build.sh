#!/bin/sh


GB_DIR=$1

ALL_DIRS=`ls -lrt $GB_DIR | sed 1d | grep -v "assembly_summary" | awk '{print $NF}' | sort`
ALL_COMBOS="23S;2000;3000 16S;1000;2000 5S;100;200"

pull_rRNAs ()
{
	GB_DIR=$1
	shift
	DIR=$1
	shift
	ALL_COMBOS=$@
	
	OUT_FILE=$TEMP_DIR$DIR"_GB_rRNAs.txt"
	echo "" | sed 1d > $OUT_FILE
	ALL_FILES=`ls -lrt $GB_DIR$DIR/*/*rna_from_genomic* | awk '{print $NF}' | head -20`


	for FILE in $ALL_FILES
	do
		#echo $FILE
		zcat $FILE | sed 's/ /_/g' | awk -v flag="N" '{ 
											if($1 ~ />/)
											{
												print ""
												printf "%s\t", $0
												flag="Y"
									
											}
											if($1 !~ />/ )
											{
												printf "%s", $1
											}
										}' | sed 's/^>/>'$DIR'/g' > temp2.txt
		
		for COMBO in $ALL_COMBOS
		do
			
			GENE=`echo $COMBO | cut -d";" -f1`
			MIN_SIZE=`echo $COMBO | cut -d";" -f2`
			MAX_SIZE=`echo $COMBO | cut -d";" -f3`

			cat temp2.txt | grep -v "NN" | grep $GENE |  awk -v MIN_SIZE=$MIN_SIZE -v MAX_SIZE=$MAX_SIZE '{if(length($2)>=MIN_SIZE && length($2)<=MAX_SIZE)print $1,$2}' >> $OUT_FILE
		done
	done
	ls -lrt $OUT_FILE


}

TEMP_DIR="/broad/hptmp/rRNA_DBs/"
mkdir -p $TEMP_DIR

echo $ALL_DIRS

for DIR in $ALL_DIRS
do
	pull_rRNAs $GB_DIR $DIR $ALL_COMBOS
done

echo $OUT_DB