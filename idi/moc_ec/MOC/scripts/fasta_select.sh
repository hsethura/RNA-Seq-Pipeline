#!/bin/sh

FILE=$1
shift

for REGEXP in $@
do

	cat $FILE | grep -A1 $REGEXP 


done