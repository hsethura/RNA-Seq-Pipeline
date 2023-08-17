#!/bin/sh

COMBINED_TAXA=$1
shift
NUM_MATCHES=$1
shift
FILES=$@

DIR="/idi/moc_ec/tmRNA"

ALL_TAXA=`echo $COMBINED_TAXA | sed 's/,/ /g' `
ALL_REFSEQ=~/"all_NCBI_DB.txt"
ALL_GENBANK="/idi/moc_ec/MOC/genbank/assembly_summary.txt"


# ALL_GENERA=`for TAXA in $ALL_TAXA 
# do
# 
# 	cat $ALL_REFSEQ | grep $TAXA | awk '{split($1, ar, "_"); print ar[1]}' | sort | uniq
# 
# done`
# 
# echo $ALL_GENERA

# ALL_SPEC=`for GENERA in $ALL_GENERA
# do	
# 
# 		cat $ALL_GENBANK | sed 's/ /_/g'| awk -F"\t" -v GENERA=$GENERA '{
# 														if($8 ~ GENERA)
# 													 	{
# 													 		split($8, ar, "_")
# 													 		print ar[1]"_"ar[2]
# 													 	}
# 												}' | sort -T /broad/hptmp | uniq
# 
# done`

echo "" | sed 1d > temp.txt

for TAXA in $ALL_TAXA
do
	ALL_SPEC=`cat $ALL_REFSEQ | grep $TAXA | awk '{split($1, ar, "_"); print ar[1]"_"ar[2]}' | sort | uniq`
	NUM_SPEC_TM=`cat $FILES | grep $TAXA | wc -l`
	echo $TAXA $NUM_SPEC_TM

	if [ $NUM_SPEC_TM -gt 0 ];then

		OUTFILE=$DIR"/"$TAXA"_tmRNA.fasta"

		echo "" | sed 1d > $OUTFILE

		ALL_SPEC=`cat $FILES | grep $TAXA | grep ">" | cut -d"|" -f2 | cut -d"_" -f1,2 | sort | uniq`

		for SPEC in $ALL_SPEC
		do	

			cat $FILES | awk -F"\t" -v SPEC=$SPEC -v flag="N" -v NUM_MATCHES=$NUM_MATCHES -v t=0 '{

				
								if($1 ~ /^>/ && $0 !~ SPEC )
									flag="N"
								if($1 ~ /^>/ && $0 ~ SPEC )
								{
									t++
									if(t <= NUM_MATCHES)
									{
										print ""
										flag="Y"
										print $0
									}
									else
										flag="N"
								}
								if($1 !~ /^>/ && flag=="Y")
									printf "%s", $1
						
			
							}' | sed '/^$/d'  >> $OUTFILE
							echo ""  >> $OUTFILE
		done	 
			echo $OUTFILE
	fi
done