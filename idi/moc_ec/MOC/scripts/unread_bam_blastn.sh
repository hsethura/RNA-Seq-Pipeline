#!/bin/sh


S_FILE=$1
shift
NUM_READS=$1
shift
ALL_BAM_FILES=$@

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


for BAM_FILE in $ALL_BAM_FILES
do
	ROOT=`basename $BAM_FILE .bam`

	Q_FILE=$TEMP_PATH$ROOT".fasta"
	OUT_FILE=$TEMP_PATH$ROOT"_output.txt"

	samtools view -f4 $BAM_FILE | head -$NUM_READS | awk '{print ">"$1"_"$2"_"NR"\n"$10}' > $Q_FILE

	#cat $Q_FILE
	echo "blastn -subject $S_FILE -query $Q_FILE -outfmt 7 | grep -v "#" > $OUT_FILE"
	blastn -subject $S_FILE -query $Q_FILE -outfmt 7 -evalue 1e-10 > $OUT_FILE

# 	echo $OUT_FILE
# 	echo $Q_FILE

	NUM_HITS=`cat $OUT_FILE | awk '{print $1}' | sort -T $TEMP_PATH | uniq -c | wc -l`
	echo $ROOT $NUM_HITS
	echo $OUT_FILE
 done



