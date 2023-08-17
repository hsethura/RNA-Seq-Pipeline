#!/bin/sh

GENUS=$1
shift

source /idi/moc_ec/MOC/scripts/bash_header

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"
paths_and_headers $MOC_ID $@

for m in $@
do
	printf "%s\t" $m
	
	
	echo $ALN_DIR$m"_"$GENUS"_aligntest.txt"	
	cat $ALN_DIR$m"_"$GENUS"_aligntest.txt" | sort -k3nr | head -1 | awk '{print $1, $2, $3}'
done 


for m in $@
do
	cat $ALN_DIR$m"_"$GENUS"_aligntest.txt"
done > temp.txt


f=`cat temp.txt | awk '{print $2}' | sort | uniq`

for m in $f
do

	printf "%s\t" $m
	awk -v total=0  -v i=0 -v m=$m '{	
								if($2==m)
								{
									total=total+$3
									i++
									print total/i
								}
							}' temp.txt | tail -1

done | sort -k2n

ALL_BEST=`for m in $f
do

	printf "%s\t" $m
	awk -v total=0  -v i=0 -v m=$m '{	
								if($2==m)
								{
									total=total+$3
									i++
									print total/i
								}
							}' temp.txt | tail -1

done | sort -k2n  | tail -1 | awk '{print $1}'`

printf "strain\taccession\tpcnt_aligned\tdif_compared_to_best\trank\n"
for m in $@
do
	printf "%s\t" $m
	cat $ALN_DIR$m"_"$GENUS"_aligntest.txt" | sort -k3nr | awk -v m=$ALL_BEST '{	
														
														if(NR==1)
															best=$3*100
														if($2==m) 
															printf "%s\t%s\t%s\t%s\t%s\n", $1, $2, $3*100, best-($3*100), NR
													}'
done 
 