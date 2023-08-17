#!/bin/sh


DB_FILE=$1
NCBI_DIR=$2
TAG=$3




f=`cat $DB_FILE | grep $TAG | awk '{print $2}' | cut -d"." -f1`

ID=`cat $DB_FILE | grep $TAG | awk '{print $1}' | awk '{

															y=split($1, ar, "_")
															for(i=1; i < y+1; i++)
															{
																if(ar[i] == "chromosome" || ar[i] == "plasmid" )
																	break
																printf "%s_", ar[i]
															
															}
															print ""
														}' | sed 's/_$//g' | sort | uniq `
																

NUM=`cat $DB_FILE | grep $ID | awk '{print $1}' | awk '{

															y=split($1, ar, "_")
															for(i=1; i < y+1; i++)
															{
																if(ar[i] == "chromosome" || ar[i] == "plasmid" )
																	break
																printf "%s_", ar[i]
															
															}
															print ""
														}' | sed 's/_$//g' | sort | uniq `



if [ $NUM -gt 1 ];then

	echo "Too many matches..."
	
	exit

fi


OUTPUT_FNA=$NCBI_DIR"/"$ID"_All.fna"
OUTPUT_GFF=$NCBI_DIR"/"$ID"_All.gff"

echo "" | sed 1d > $OUTPUT_FNA
echo "" | sed 1d > $OUTPUT_GFF

echo $f

for m in $f
do
	echo "Working on $m..."
	ls -lrt $NCBI_DIR"/"$m".fna"
	ls -lrt $NCBI_DIR"/"$m".gff" 
	
	
	ACC=`cat $NCBI_DIR"/"$m".gff" | grep -v "#" | awk '{print $1}' | sort | uniq -c | sort | tail -1 | awk '{print $2}'`
	echo ">"$ACC  >> $OUTPUT_FNA
	cat $NCBI_DIR"/"$m".fna" | grep -v ">" >> $OUTPUT_FNA
	cat $NCBI_DIR"/"$m".gff" >> $OUTPUT_GFF
done


echo $OUTPUT_FNA
echo $OUTPUT_GFF