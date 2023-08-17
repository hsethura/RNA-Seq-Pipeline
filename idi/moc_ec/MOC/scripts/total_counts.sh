#!/bin/sh

### Calculates the total reads per column from count table


FILE=$1

typeset -i NF
NF=`cat $FILE | head -1 | awk '{print NF}'`
i=2

while [ $i -lt $NF ]
do
	name=`cat $FILE | head -1 | awk '{print $'$i'}'`
	total=`cat $FILE | awk -v total=0 '{total=total+$'$i'; print total}' | tail -1`
	echo $name	$total
	i=`expr $i + 1`
done

