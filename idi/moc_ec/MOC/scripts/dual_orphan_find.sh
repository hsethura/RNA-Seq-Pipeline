#!/bin/sh


DIR=$1

source /broad/software/scripts/useuse

use Python-2.7

source /idi/moc_ec/MOC/scripts/bash_header

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/Universal_config.yaml" 1 $@`

### determining paths and headers 
### default config file is /idi/moc_ec/MOC/config_files/Universal_config.yaml
paths_and_headers $MOC_ID $@

BC1_FILE=`ls -lrt $DIR"/"*"1.unmatched.barcode_1.fastq.gz" | awk '{print $NF}'`
BC2_FILE=`echo $BC1_FILE | sed 's/1.fastq/2.fastq/g'`
echo $BC1_FILE
echo $BC2_FILE

NUM_READS=100000

TEMP_DIR="/broad/hptmp/livny/"
mkdir -p $TEMP_DIR 

TEMP_FILE1=$TEMP_DIR"ind1.txt"
TEMP_FILE2=$TEMP_DIR"ind2.txt"
TEMP_FILE3=$TEMP_DIR"temp.txt"
TEMP_BCS=$TEMP_DIR"BCS_temp.txt"

ls -lrt $BC1_FILE
ls -lrt $BC2_FILE
zcat $BC1_FILE | awk '{if(NR % 4 ==2) print $1}' | head -$NUM_READS > $TEMP_FILE1
zcat $BC2_FILE | awk '{if(NR % 4 ==2) print $1}' | head -$NUM_READS > $TEMP_FILE2

paste $TEMP_FILE1 $TEMP_FILE2 | sort -T $TEMP_DIR | uniq -c | sort | awk -v NUM=$NUM_READS '{print ($1/NUM)*100, $2, $3}' | tail -100 > $TEMP_FILE3

cat $TEMP_FILE3
ls -lrt $BCS_FILE

cat $BCS_FILE | sed 's/ /_/g' > $TEMP_BCS
while read line
do

		PCNT=`echo $line | awk '{print $1}'`
		IS1=`echo $line | awk '{print $2}'`
		IS2=`echo $line | awk '{print $3}'`
		ID1=`FIND_HEADER_AND_FIELD $IS1 "P7_Reverse_Complement_Index_Sequence" "P7_Index_Code" $TEMP_BCS | sort | uniq`	
		ID2=`FIND_HEADER_AND_FIELD $IS2 "P5_Reverse_Complement_Index_Sequence" "P5_Index_Code" $TEMP_BCS | sort | uniq`	
		
		if [ -z $ID1 ];then
			ID1="-"
		fi
		if [ -z $ID2 ];then
			ID2="-"
		fi

		echo $PCNT $IS1 $IS2 $ID1 $ID2
done < $TEMP_FILE3

ls -lrt $TEMP_BCS
ls -lrt $TEMP_FILE3


chmod 777 $TEMP_FILE1 $TEMP_FILE2