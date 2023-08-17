#!/bin/sh


DIR=$1

if [ $# -eq 1 ];then
	ALL_BAMS=`ls -lrt $DIR/*bam | awk '{print $9}'`
fi

if [ $# -gt 2 ];then
	ALL_BAMS=`ls -lrt $DIR/*/*bam | awk '{print $9}'`
fi

echo $ALL_BAMS


for BAM in $ALL_BAMS
do
	echo $BAM
	ROOT=`basename $BAM .bam`
	OUT_ROOT="/broad/hptmp/MOC/"$ROOT
	samtools view -u -f 1 -F 12 $BAM > $OUT_ROOT"_map_map.bam"

	# R1 unmapped, R2 mapped
	samtools view -u -f 4 -F 264 $BAM > $OUT_ROOT"_unmap_map.bam"
	# R1 mapped, R2 unmapped
	samtools view -u -f 8 -F 260 $BAM > $OUT_ROOT"_map_unmap.bam"
	# R1 & R2 unmapped
	samtools view -u -f 12 -F 256 $BAM > $OUT_ROOT"_unmap_unmap.bam"

	samtools merge -u $OUT_ROOT"_unmapped.bam" $OUT_ROOT"_unmap_map.bam" $OUT_ROOT"_map_unmap.bam" $OUT_ROOT"_unmap_unmap.bam"

	samtools sort -n $OUT_ROOT"_map_map.bam" $OUT_ROOT"_mapped.sort"
	samtools sort -n $OUT_ROOT"_unmapped.bam" $OUT_ROOT"_unmapped.sort"
	
	samtools flagstat lib_002.sorted.md.bam

done