#!/bin/sh


DIR=$1
FILE=$2
NUM_STRAINS=$3

ALL_GENES="23S 16S 5S"

ALL_ID=`cat $FILE`

for ID in $ALL_ID
do

	ALL_FILES=`ls -lrt $DIR* | grep $ID | awk '{print $NF}'`
	for FILE in $ALL_FILES
	do
		ALL_STRAINS=`cat $FILE | grep 23S | awk '{
										split($1, ar, "|")
										split(ar[2], pr, "_")
										print pr[1]
											
											
											
											
									}' | sort | uniq`
									
		for GENE in $ALL_GENES
		do 
			ALL_FASTA=`for STRAIN in $ALL_STRAINS
			do
				cat $FILE | grep $GENE | grep $STRAIN | awk '{
										split($1, ar, "|")
										print $1, length($2)
										}' | tail -1
			done | sort -k2n | tail -$NUM_STRAINS | cut -d"[" -f1`
			
			
			for FASTA in $ALL_FASTA
			do
				cat $FILE | grep $FASTA | awk '{print $1"\n"$2}'
				
			done
		done 
		
	done

done