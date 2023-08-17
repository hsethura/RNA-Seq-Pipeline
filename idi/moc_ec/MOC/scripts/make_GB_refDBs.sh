#!/bin/sh

FILE=$1
OUT_TAG=$2

ALL_SPEC=`cat $FILE`

echo $ALL_SPEC

OUTPUT_PATH="/idi/moc_ec/MOC/genbank/"

OUT_FNA=$OUTPUT_PATH"/"$OUT_TAG".fna"
OUT_GFF=$OUTPUT_PATH"/"$OUT_TAG".gff"

echo "" | sed 1d > $OUT_FNA
echo "" | sed 1d > $OUT_GFF

for SPEC in $ALL_SPEC
do
	echo $OUTPUT_PATH"/"$SPEC
	FNA_FILE=`ls -lrt $OUTPUT_PATH"/"$SPEC/*/*fna.gz | grep -v rna | sort -k5n | tail -1 | awk '{print $NF}'`
	ACC=`echo $FNA_FILE | sed 's/_genomic.fna.gz//g' | rev | cut -d"/" -f1 | rev`
	GFF_FILE=`ls -lrt $OUTPUT_PATH"/"$SPEC/*/*gff.gz | grep $ACC | awk '{print $NF}'`
	if [ ! -z $FNA_FILE ];then
		ls -lrt $FNA_FILE
		zcat $FNA_FILE >> $OUT_FNA
	fi
	if [ ! -z $GFF_FILE ];then
		ls -lrt $GFF_FILE
		zcat $GFF_FILE >> $OUT_GFF
	fi
done


ls -lrt $OUT_GFF
ls -lrt $OUT_FNA