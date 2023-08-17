#!/bin/sh

for m in $@
do
	
	DIR=$m
	
	ALL_GFFS=`ls -lrt $DIR/* | grep gff | awk '{print $9}'`
	
	for GFF_FILE in $ALL_GFFS
	do
		echo $GFF
		echo $DIR
	
		FNA_FILE=`echo $GFF_FILE | sed 's/_GENES.gff/.fna/g'`
	
		echo $GFF_FILE
		echo $FNA_FILE

		if [ -s $GFF_FILE ] && [ -s $FNA_FILE ];then
		
			GFF_ACC=`cat $GFF_FILE | grep -v "#" | head -10 | awk '{print $1}' | sort | uniq`
			FNA_ACC=`cat $FNA_FILE | head -1 | awk '{print $1}'| sed 's/>//g'`
	
	
			echo "Before GFF_ACC: "$GFF_ACC
			echo "Before FNA_ACC: "$FNA_ACC

			if [ $GFF_ACC != $FNA_ACC ];then
	
				echo "Switching ACCs..."
				cat $GFF_FILE | sed 's/'$GFF_ACC'/'$FNA_ACC'/g' > temp.txt
				cat temp.txt > $GFF_FILE
			fi
	
			GFF_ACC=`cat $GFF_FILE | grep -v "#" | head -10 | awk '{print $1}' | sort | uniq`
			FNA_ACC=`cat $FNA_FILE | head -1 | awk '{print $1}'| sed 's/>//g'`
	
		# 	echo $GFF_FILE
		# 	echo $FNA_FILE
	
			echo "After GFF_ACC: "$GFF_ACC
			echo "After FNA_ACC: "$FNA_ACC
		fi
	done
done

