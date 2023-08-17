#!/bin/sh


FASTA_FILE=$1
FASTQ_FILE=$2
TRIM_5=$3
TRIM_3=$4


cat $FASTA_FILE | awk -v T5=$TRIM_5 -v T3=$TRIM_3 '{
														
														if($1 ~ />/)
															print $0
														if($1 !~ />/)
														{	
															len=length($1)
															print substr($1, (T5+1), (len-(T3+T5+1)))
															print "+"
															for(i=1; i < len-(T3+T5); i++)
																printf "A"
															print "" 

														}
												}' | sed 's/>/@/g' | sed 's/ /_/g' > $FASTQ_FILE
												

echo $FASTQ_FILE