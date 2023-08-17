#!/bin/sh

#### 1st argument is path to directory containing BAMs
BAM_DIR=$1
shift
#### arguments 2-n are # of reads you want to pull
ALL_NUM_READS=$@

#number of bams to pull
NUM_BAMS=1000

echo $@



#set default for 2nd argument
if [ $# == 0 ];then
	ALL_NUM_READS=4	
fi



OUT_DIR=$BAM_DIR"/unaligned_fasta/"

mkdir -p $OUT_DIR


echo $BAM_DIR

### set array containing paths to all BAMs
ALL_BAMS=`ls $BAM_DIR*bam | head -$NUM_BAMS` 

### set num of BAMs
BAM_NUM=`echo $ALL_BAMS | wc -w`

echo $BAM_NUM

i=1

while [ $i -le $NUM_BAMS ] && [ $i -le $BAM_NUM ]
do

	### loop through BAMs
	for BAM in $ALL_BAMS
	do
	
		#loop through numbers in array
		for NUM_READS in $ALL_NUM_READS 
		do


			### set rootname for each bam, stripping suffix
			ROOT=`basename $BAM bam`
			OUTFILE=$OUT_DIR"/"$ROOT"_"$NUM_READS"_unaligned.fasta"
	
			samtools view -f4 $BAM | head -$NUM_READS | awk '{print ">"$1"\n"$10}' > $OUTFILE
	
			ls -lrt $OUTFILE
		
		done
		
	
	i=`expr $i + 1`
	done
done