#!/bin/sh


FILE=$1

while read line
do
	echo "removing "$line
	rm -rf $line

done < $FILE