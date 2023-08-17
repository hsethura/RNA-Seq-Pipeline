#!/bin/sh

MIN=$1
shift

ALL_FILES=$@

for FILE in $ALL_FILES
do
	typeset -i CDS_F
	echo $FILE
	CDS_F=`cat $FILE | head -3 | awk -F"\t" '{
									for(i=1; i < 200; i++)
										if($i ~ /CDS_total_counts/)
										{	
											print i
											exit
										}
									
									}'`
	CDS_F=`cat $FILE | head -3 | awk -F"\t" '{
								for(i=1; i < 200; i++)
									if($i ~ /CDS_total_counts/)
									{	
										print i
										exit
									}
									print "0"
								}' | tail -1`
		
	if [ $CDS_F != "0" ];then							
		cat $FILE | awk -F"\t" -v CDS_F=$CDS_F '{if($CDS_F < '$MIN') print $2, $3, $4, $CDS_F}'
	fi

done


