#!/bin/sh

FILE1=$1
FILE2=$2

cat $FILE2 | awk -F"\t" '{for(i=2; i < NF; i++) printf "%s\t", $i; print $NF}' > temp.txt

paste $FILE1 temp.txt