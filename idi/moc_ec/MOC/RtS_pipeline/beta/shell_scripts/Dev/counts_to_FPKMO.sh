#!/bin/sh


COUNTS_FILE=$1
OUTFILE=$2
TEMP_DIR=$3

TEMP_FILE1=$TEMP_DIR"/temp1.txt"
TEMP_FILE2=$TEMP_DIR"/temp2.txt"
TEMP_FILE3=$TEMP_DIR"/temp3.txt"

typeset -i i
typeset -i NF

NF=`cat $COUNTS_FILE | awk '{print NF}' | sort | uniq -c | sort | tail -1 | awk '{print $2}'`

i=2

typeset -i NL
NL=`cat $COUNTS_FILE | wc -l`


cat $COUNTS_FILE | awk '{print $1}' > $TEMP_FILE1

while [ $i -le $NF ]
do

	echo $i
	
	TOTALS=`awk -v i=$i -v total=0 -v total_CDS=0 '{
					
				if(NR > 1)
				{
					split($1, ar, ":")
					
					if(ar[1]=="CDS")
						total_CDS=total_CDS+$i
					total=total+$i
				}
				print total/1e6, total_CDS/1e6
	
	
		}' $COUNTS_FILE | tail -1`
	
	TOTAL=`echo $TOTALS | awk '{print $1}'`
	TOTAL_CDS=`echo $TOTALS | awk '{print $2}'`
	
	
	
	#echo $TOTAL, $TOTAL_CDS
	
	awk -v i=$i -v TOTAL=$TOTAL -v TOTAL_CDS=$TOTAL_CDS '{
					
				if(NR == 1)
					print $i
				if(NR > 1)
				{	
					y=split($1, ar, ":")
					len=ar[y-1]
				
					if(len==0)
					  exit
					if(NR == 1)
						print $i
					if(NR > 1)
					{
						print $i/(TOTAL_CDS+1)/(len+1/1000)
					}
				}				
				
		}' $COUNTS_FILE > $TEMP_FILE2
		
		cat $TEMP_FILE1 > $TEMP_FILE3
		paste $TEMP_FILE3 $TEMP_FILE2 > $TEMP_FILE1
	
	
	i=`expr $i + 1`

done


cat $TEMP_FILE1 > $OUTFILE

echo $OUTFILE




