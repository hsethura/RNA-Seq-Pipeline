#!/bin/sh

#source /broad/IDP-Dx_work/nirmalya/bash_header

FASTQ_FILE=$1
NUM=$2

DB="/home/unix/livny/work/all_rRNA_wID.txt"

FASTA_FILE="fasta_temp.txt"

typeset -i GZ
GZ=`echo $FASTQ_FILE | grep ".gz" | wc -l`

if [ $GZ -gt 0 ];then

	zcat $FASTQ_FILE | awk -v i=0 -v NUM=$NUM '{
										if(NR % 4 == 1) 
											print ">"$1
										if(NR % 4 == 2) 
										{	
											y=substr($1, 5, 100)
											print y
											i++
										}
										if(i >= NUM)
											exit
										}' > $FASTA_FILE
else

	cat $FASTQ_FILE | awk -v i=0 -v NUM=$NUM '{
										if(NR % 4 == 1) 
											print ">"$1
										if(NR % 4 == 2) 
										{	
											print $1
											i++
										}
										if(i >= NUM)
											exit
										}' > $FASTA_FILE

fi

echo $FASTA_FILE
echo $DB

echo "blastn -query fasta_temp.txt -subject /broad/IDP-Dx_storage/NCBI_files2/NC_000913.fna -outfmt 6 | awk '{print $1}' | sort | uniq | wc -l"


exit

echo "Blasting $FASTA_FILE against $DB"
echo "blastn -query $FASTA_FILE -subject $DB -outfmt 6 -max_hsps 1"

NUM_HITS=`blastn -query $FASTA_FILE -subject $DB -outfmt 6 -max_hsps 1 | wc -l`

echo $NUM_HITS  " out of "$NUM" reads aligned to rRNA"