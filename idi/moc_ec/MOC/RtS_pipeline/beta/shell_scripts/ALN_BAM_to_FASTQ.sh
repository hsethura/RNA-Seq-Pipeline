#! /bin/bash
#################################################################################################
## this script is for aligning illumina reads to reference genome(s) given a BAM file and the reference genome(s) in a multi-fasta file
## $1: BAM file
## $2: the output prefix for all the output files
## $3: 1, if BAM is pre-aligned; 0, is un-aligned.
## $4: 1, if reads are pair-ended; 0, not pair-ended
## $5: path to Temp_dir


## um, need some code cleanup and also need more error checking and write functions!
#################################################################################################

usage="Usage: $0 <BAM file> <genome sequence fasta> <output prefix> <isPreAligned> <is pair-ended> <remove PCR dup>"

if [ -z "$1" ]; then
        echo "Must specify a BAM file"
        echo $usage
        exit 1
elif [ -z "$2" ]; then
        echo "Must specify the prefix for all the output files"
        echo $usage
        exit 3
elif [ -z "$3" ]; then
        echo "Please indicate whether the BAM file is pre-aligned or not"
        echo $usage
        exit 3
elif [ -z "$4" ]; then
        echo "Please indicate whether the reads are pair-ended or not" 
        echo $usage
        exit 3
fi

# Must preload Java, BWA and Samtools
source ~/.bashrc # in case the paths are not set
reuse -q Java-1.6
reuse -q BWA
reuse -q Samtools

# Input params
BAM_INFILE=$1
ROOT=$2
FASTQ1=$3
FASTQ2=$4

typeset -i ALIGN_CHOICE
ALIGN_CHOICE=$5
typeset -i PAIRED_CHOICE
PAIRED_CHOICE=$6


# Pre-defined variables
PICARD_BIN="/seq/software/picard/1.782/bin/"
TMP_DIR=$7
bwa="/home/radon00/khhuang/apps/bwa-0.5.9/bwa"


# OUTPUT FILES
UNALIGNED_BAM=$TMP_DIR"/"$ROOT".unaligned.bam"
FASTQ=$TMP_DIR"/"$ROOT".fastq"
ALIGNED_SAI=$TMP_DIR"/"$ROOT".aln_sai"
SAM=$TMP_DIR"/"$ROOT".aligned.sam"
ALIGNED_BAM=$TMP_DIR"/"$ROOT".aligned.bam"
SUMMARY_FILE=$TMP_DIR"/"$ROOT"_raw.summary.txt"

function revertBAM {
   
	if [ $ALIGN_CHOICE -eq 1 ]; then
			# Generate a summary file:
			echo "Generate summary report of the original BAM: $BAM_INFILE"
			samtools flagstat $BAM_INFILE > $SUMMARY_FILE
			echo "Revert BAM to unaligned BAM: $UNALIGNED_BAM"
			echo "java -Xmx2g -jar $PICARD_BIN/RevertSam.jar I=$BAM_INFILE O=$UNALIGNED_BAM TMP_DIR=$TMP_DIR"
			java -Xmx2g -jar $PICARD_BIN/RevertSam.jar I=$BAM_INFILE O=$UNALIGNED_BAM TMP_DIR=$TMP_DIR
	else
		cp $BAM_INFILE $UNALIGNED_BAM
	fi

}

function getFastQ {
    
   	echo "Running getFastQ..."
   	# NON_PF: include not passing filter reads, i.e. include all reads

    if [ $PAIRED_CHOICE -eq 0 ]; then
		echo "single end"
		echo "Convert unaligned BAM to FASTQ: $FASTQ"
		# it is important to unalign BAM first, otherwise the read seq might be reverse complemented
		java -Xmx2g -jar $PICARD_BIN/SamToFastq.jar I=$UNALIGNED_BAM F=$FASTQ TMP_DIR=$TMP_DIR NON_PF=true
    else
    # pair end reads
		echo "pair end"
		echo "Convert unaligned BAM to FASTQs $FASTQ1 and $FASTQ2"
		echo "java -Xmx2g -jar $PICARD_BIN/SamToFastq.jar I=$UNALIGNED_BAM F=$FASTQ1 F2=$FASTQ2 TMP_DIR=$TMP_DIR NON_PF=true"
		java -Xmx2g -jar $PICARD_BIN/SamToFastq.jar I=$UNALIGNED_BAM F=$FASTQ1 F2=$FASTQ2 TMP_DIR=$TMP_DIR NON_PF=true
    fi
	
	echo "###########################################################################################"

}



# Pipeline #
# additional note: this pipeline is essentially what I learned from Sean Sykes
# -Xmx2g: run in 2GB of JVM
# OQ: restore original qualities

echo $UNALIGNED_BAM
revertBAM
getFastQ

echo $FASTQ1
echo $FASTQ2

echo "finished"
