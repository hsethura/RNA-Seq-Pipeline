#!/bin/sh


### source all functions 
source "/broad/IDP-Dx_storage/MOC/scripts/MOC_functions.sh"


DIR=$1

cd $DIR

MAN_FILE=$DIR"/MANIFEST"
RAW_FILE=$DIR"/Checksums.txt"
SEG_FILE=$DIR"/SGE_file.txt"
VAL_FILE=$DIR"/Download_validation.txt"

ALL_FILES=`cat $MAN_FILE | awk '{print $1}' | grep fastq.gz`

echo $ALL_FILES

echo "" | sed 1d > $RAW_FILE
echo "" | sed 1d > $SEG_FILE					

for FILE in $ALL_FILES
do

	echo "Launching md5sum for $FILE"
	echo "gzip -d -c $DIR"/"$FILE| md5sum | awk -v FILE=$FILE '{print FILE, \$1}' >> $RAW_FILE" > temp.txt
	qsub -l h_rt=12:00:00 temp.txt >> $SEG_FILE	
	echo $SEG_FILE				
done

qstat
echo "Waiting for UGER jobs to complete"
SGE_test $SEG_FILE	

echo "Comparing md5sums for files at server and from getsite" 
for FILE in $ALL_FILES
do
	DIR_SUM=`cat $RAW_FILE | grep $FILE | awk '{print $2}'`
	WU_SUM=`cat $MAN_FILE | grep $FILE | awk '{print $2}'`
	
	if [ $DIR_SUM == $WU_SUM ];then
		echo $FILE" DOWNLOADED" 
	else
		echo $FILE" ************** NOT FULLY DOWNLOADED *************" 
	fi
	
done | sort -k2 > $VAL_FILE

echo $VAL_FILE

ALL_FILES=`cat $VAL_FILE | grep "NOT FULLY" | awk '{print $1}'`

for FILE in $ALL_FILES
do
	
	echo $FILE
	ls -lrt $DIR | grep $FILE

done