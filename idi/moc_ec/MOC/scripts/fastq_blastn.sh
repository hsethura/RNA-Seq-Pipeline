#!/bin/sh


SPEC_GB_DIR=$1
shift
OUT_DIR=$1
shift
NUM_READS=$1
shift
ALL_FASTQ_FILES=$@


GFF_FILE_DIR=`ls -lrt $SPEC_GB_DIR*/* | grep genomic.gff.gz | awk '{print $NF}' | tail -1`
ACC=`echo $GFF_FILE_DIR | rev | cut -d"/" -f2 | rev`
REF_DIR=$SPEC_GB_DIR"/"$ACC"/"

FNA=`ls -lrt $REF_DIR* | grep genomic | grep -v rna | awk '{print $NF}'`
GFF=`ls -lrt $REF_DIR* | grep gff | grep -v RNA | awk '{print $NF}'`

echo $FNA


source /broad/software/scripts/useuse

use Python-2.7

source /idi/moc_ec/MOC/scripts/bash_header

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"
CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/Universal_config.yaml" 1 $@`
TEMP_PATH=`config_read $CONFIG_FILE Temp_path`


### determining paths and headers 
### default config file is /idi/moc_ec/MOC/config_files/Universal_config.yaml
paths_and_headers $MOC_ID $@


OUT_FNA=$TEMP_PATH"/"$ACC".fna"
OUT_GFF=$TEMP_PATH"/"$ACC".gff"

ls -lrt  $OUT_FNA
ls -lrt  $OUT_GFF

zcat $FNA > $OUT_FNA
zcat $GFF > $OUT_GFF

for FASTQ in $ALL_FASTQ_FILES
do
	ROOT=`basename $FASTQ .fastq.gz`

	S_FILE=$TEMP_PATH$ROOT".fasta"
	OUT_FILE=$TEMP_PATH$ROOT"_output.txt"
	COORD_FILE=$TEMP_PATH$ROOT".coords"

	zcat $FASTQ | awk '{if(NR%4==2)print ">"NR,$1}' | head -$NUM_READS | awk '{print $1"\n"$2}' > $S_FILE

	#cat $S_FILE
	echo "blastn -subject $Q_FILE -query $OUT_DIR"/"$ACC".fna"  -outfmt 7 -evalue 1e-10 > $OUT_FILE"
	blastn -subject $OUT_DIR"/"$ACC".fna" -query $S_FILE -outfmt 7 -evalue 1e-10 2>/dev/null > $OUT_FILE

ls -lrt  $OUT_DIR"/"$ACC".fna"
ls -lrt temp.fna
ls -lrt $S_FILE
# 	echo $OUT_FILE
# 	echo $Q_FILE

	NUM_HITS=`cat $OUT_FILE | awk '{print $1}' | sort -T $TEMP_PATH | uniq -c | wc -l`
	echo $ROOT $NUM_HITS
	echo $OUT_FILE
	
	sh /idi/moc_ec/MOC/scripts/blast_to_coords.sh $OUT_FILE $COORD_FILE
	
	echo $COORD_FILE
	ACC_FILE=$OUT_DIR"/"$ACC".txt"
	
	cat $OUT_DIR"/"$ACC".fna" | grep ">" | sed 's/>//g' > $ACC_FILE
	
	PARSED_GFF=$TEMP_PATH"/"$ACC"_parsed.gff"
	
	cat $OUT_GFF | grep -v "#" |  sed 's/ /_/g' |  awk -F"\t" '{
	
							if($3 != "gene")
							{
								$NF="type="$3";"$NF
								print $0
							}
						}' > $PARSED_GFF
	
	
	
	NUM_PARSED_GFF=`ls -lrt $PARSED_GFF | wc -l`	
	echo $NUM_PARSED_GFF
	ls -lrt $PARSED_GFF
	
# 	if [ $NUM_PARSED_GFF -eq 0 ];then
# 		echo "sh /broad/IDP-Dx_work/nirmalya/pipeline/beta/shell_scripts/gff_parse4.sh $TEMP_PATH $TEMP_PATH  $TEMP_PATH $TEMP_DIR 20 30 $ACC"
# 		sh /broad/IDP-Dx_work/nirmalya/pipeline/beta/shell_scripts/gff_parse4.sh $TEMP_PATH $TEMP_PATH  $TEMP_PATH $TEMP_DIR 20 30 $ACC
# 	fi

	/broad/IDP-Dx_work/nirmalya/pipeline/beta/shell_scripts/FEATURE_COUNT_COORDS $COORD_FILE $PARSED_GFF $OUT_DIR"test.mets" $ACC_FILE N > temp_counts.txt
	
	cat temp_counts.txt | grep "type=CDS" | awk -v total=0 '{total=total+$(NF-1); print total}' | tail -1
	cat temp_counts.txt | grep "type=rRNA" | awk -v total=0 '{total=total+$(NF-1); print total}' | tail -1

	
 done
