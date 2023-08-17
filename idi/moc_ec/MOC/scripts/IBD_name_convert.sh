#!/bin/sh
 
ID_INDEX=$1
FILES=$2
OG_DIR=$3
NEW_DIR=$4

mkdir -p $NEW_DIR

cat $ID_INDEX | sed 's/(pooledLibrary=MOCW-00[0-9][1-9], lane=[1-4], mbs=//g' | sed 's/)//g' | sed 's/-/_/g' | awk '{
																													for(i=1; i<NF-1; i++)
																														printf "%s_", $i
																													printf "%s;", $(NF-1)
																													print $NF
																													
																													}' > temp_IDs.txt
																													
ALL_IDS=`cat temp_IDs.txt`

for ID in $ALL_IDS
do

	SID=`echo $ID | cut -d";" -f1`
	INDEX=`echo $ID | cut -d";" -f2`
	
	echo $SID, $INDEX
	
	ALL_FILES=`cat $FILES | grep -w $INDEX`
	
	for FILE in $ALL_FILES
	do
	
		 OUTPUT=`echo $FILE | sed -e s'*'$INDEX'*'$SID'_'$INDEX'*g' -e s'*'$OG_DIR'*'$NEW_DIR'*g'`
		 cp $FILE $OUTPUT
	done


done