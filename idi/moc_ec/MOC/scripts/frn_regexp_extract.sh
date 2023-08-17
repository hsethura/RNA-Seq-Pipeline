#!/bin/sh

COMBINDED_REGEXP=$1
shift
FILES=$@


ALL_REGEXP=`echo $COMBINDED_REGEXP | sed 's/,/ /g' `


echo "" | sed 1d > temp.txt

for REGEXP in $ALL_REGEXP
do	
	cat $FILES | awk -F"\t" '{

						if($1 ~ />/)
						{
							print ""
						printf "%s ; ", $1
						
						}
						else 
							printf "%s", $1
					
					}' | grep "gene="$REGEXP >> temp.txt
done				
ALL_NCBI="~/all_NCBI_DB.txt"

while read line 
do

	ACC=`echo $line | cut -d"|" -f2`
	NAME=`cat ~/all_NCBI_DB.txt | grep $ACC | awk '{print $1}'`
	
	echo $line | sed 's/ref/ref|'$NAME'/g'

done < temp.txt | sed 's/;/\n/g'
