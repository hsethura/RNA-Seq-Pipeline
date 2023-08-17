#! /bin/bash
#################################################################################################
## this script is for aligning illumina reads to reference genome(s) given a FASTQ file and the reference genome(s) in a multi-fasta file
## $1: FASTQ file
## $2: genome sequence fasta file
## $3: the output prefix for all the output files


# Input params
FASTQ=$1
FNA=$2
TAG=$3
ALN_SAM=$4
TEMP_DIR=$5


# OUTPUT FILES
ALIGNED_SAI=$TEMP_DIR"/"$TAG".aln_sai"
ALIGNED_BAM="$OUTFILE_PREFIX.aligned.bam"

# Must preload Java, BWA and Samtools
# source ~/.bashrc # in case the paths are not set
# reuse -q Java-1.6
# reuse -q BWA
# reuse -q Samtools


echo "Outfile is $SAM"

echo "bwa index -a is $FNA"
bwa index -a is $FNA
echo "bwa aln $FNA $FASTQ > $ALIGNED_SAI"
bwa aln $FNA $FASTQ > $ALIGNED_SAI
echo "bwa samse $FNA $ALIGNED_SAI $FASTQ > $ALN_SAM"
bwa samse $FNA $ALIGNED_SAI $FASTQ > $ALN_SAM