#!/bin/sh


DIR=$1
shift

ALL_INDEXES=$@

TEMP_FILE="temp.txt"
OUTFILE_FILE="temp.txt"

ALL_OUT_FILES=`ls -lrt $DIR/* | grep "_out.txt" | awk '{print $9}'`

echo $ALL_OUT_FILES

ALL_P7=`cat $ALL_OUT_FILES | awk '{if(NF==2)print $1}' | sed 's/://g' | sort | uniq`

printf "Index\t" > $TEMP_FILE

for FILE in $ALL_OUT_FILES
do
	FILE_TAG=`echo $FILE | rev | cut -d"/" -f1 | rev | sed 's/_out.txt//g'`
	printf "%s\t" $FILE_TAG >> $TEMP_FILE
done

echo "Average" >> $TEMP_FILE

for P7 in $ALL_P7
do
	FLAG="N"
	for INDEX in $ALL_INDEXES
	do
		
		if [ $INDEX == $P7 ];then
			P7=$P7"*"
		fi
	
	done
	
	printf "%s:\t" $P7 >> $TEMP_FILE
	for FILE in $ALL_OUT_FILES
	do
		cat $FILE | grep $P7 | sed 's/%//g' | awk '{printf "%.3f\t", $2}' >> $TEMP_FILE
	done
		cat $ALL_OUT_FILES | grep $P7 | sed 's/%//g' | awk -v total=0 -v i=0 '{total=total+$2; i++;printf "%.3f\n", total/i}' | tail -1 >> $TEMP_FILE

done

cat $TEMP_FILE | head -1 | awk -F"\t" '{print $1, $6" | ", $2, $3, $4, $5}' 
cat $TEMP_FILE | sed 1d | awk -F"\t" '{print $1, $6" | ", $2, $3, $4, $5}' | sort -k6n
exit

echo $ALL_OUT_FILES

exit

for INDEX in $ALL_INDEXES
do



done

