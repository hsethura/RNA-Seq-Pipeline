#! /bin/bash
#################################################################################################
## this script is for aligning illumina reads to reference genome(s) given a FASTQ file and the reference genome(s) in a multi-fasta file
## $1: FASTQ file
## $2: genome sequence fasta file
## $3: the output prefix for all the output files


# Input params
INFILE=$1
REF_GENOME=$2
OUTFILE_PREFIX=$3
UNIQ_REF_GENOME=$OUTFILE_PREFIX"_ref.fna"



# OUTPUT FILES
UNALIGNED_BAM="$OUTFILE_PREFIX.unaligned.bam"
FASTQ="$INFILE"
ALIGNED_SAI="$OUTFILE_PREFIX.aln_sai"
SAM="$OUTFILE_PREFIX.aligned.sam"
ALIGNED_BAM="$OUTFILE_PREFIX.aligned.bam"

# Must preload Java, BWA and Samtools
# source ~/.bashrc # in case the paths are not set
# reuse -q Java-1.6
# reuse -q BWA
# reuse -q Samtools

# Make uniq ref file to avoid conflicts when running in LSF
echo "Making unique reference file for run"
cp $REF_GENOME $UNIQ_REF_GENOME

echo "Outfile is $SAM"

bwa index -a is $UNIQ_REF_GENOME
bwa aln $UNIQ_REF_GENOME $FASTQ > $ALIGNED_SAI
echo "bwa samse $UNIQ_REF_GENOME $ALIGNED_SAI $FASTQ > $SAM"
bwa samse $UNIQ_REF_GENOME $ALIGNED_SAI $FASTQ > $SAM