#!/bin/sh


BAM=$1
GFF=$2

ALL_GENES="23S 16S 5S"
ROOT=`echo $BAM | rev | cut -d"/" -f1 | rev | sed 's/.bam//g'`
echo $ROOT

for GENE in $ALL_GENES
do
	
	OUTFILE=$ROOT"_"$GENE".txt"
	
	START=`cat $GFF | grep $GENE | awk '{if($3=="rRNA") print $0}' | head -1 | awk '{print $4}'`  
	END=`cat $GFF | grep $GENE | awk '{if($3=="rRNA") print $0}' | head -1 | awk '{print $5}'`  
	EXTEND=100
	#samtools view -h -f 0x0040 $BAM | awk -v start=$START -v end=$END '{if($4 > start-100  && $4 < end+100) print $0}'

	NUM_READS=`samtools view -h -f 0x0040 $BAM | awk -v start=$START -v end=$END '{if($4 > start-'$EXTEND'  && $4 < end+'$EXTEND') print $4}' | wc -l`

	echo $NUM_READS

	samtools view -h -f 0x0040 $BAM | awk -v start=$START -v end=$END '{if($4 > start-'$EXTEND' && $4 < end+'$EXTEND') print $4}' | sort | uniq -c | sort -k2n | awk -v NUM_READS=$NUM_READS '{print $1/NUM_READS*100, $2}' > temp.txt


	typeset -i NL
	typeset -i S
	typeset -i E

	inc=25

	S=`cat temp.txt | head -1 | awk '{print $2}'`
	E=`cat temp.txt | tail -1 | awk '{print $2}'`
	i=$S

	echo $S $E

	while [ $i -lt $E ]
	do
		printf "%s\t" $i
		awk -v i=$i -v inc=$inc -v total=0 '{
												if($2 >= i && $2 <= i+inc) 
													total=total+$1
												if($2 > i+inc)
												{
													print total
													exit
											
												}
										
										
											}' temp.txt
	
		i=`expr $i + 1`
	done  > temp2.txt

	cat temp2.txt | awk '{if(NF==2)print $0}'| sort -k2n > $OUTFILE

	cat $OUTFILE
	ls -lrt $OUTFILE

done