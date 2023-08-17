#!/bin/sh


DIR=$1

ALL_FILES=`ls $DIR`

echo $ALL_FILES

for FILE in $ALL_FILES
do

	NEW_NAME=`echo $FILE | sed 's/+/_plus_/g'`
	
	mv $DIR$FILE $DIR$NEW_NAME

done