#!/bin/sh

SMOC=$1
shift

FIELD_HEADER ()
{

	KEY_FILE=$1
	HEADER_NAME=$2

	cat $KEY_FILE | grep -v "#" | head -1 | awk -F"\t" -v HEADER_NAME=$HEADER_NAME '{	
																for(i=1; i < NF+1; i++)															
																	if($i==HEADER_NAME)
																		print i
															}' 
}



KEY_DIR="/broad/IDP-Dx_storage/MOC/Key_files/"


f=`ls -lrt $KEY_DIR/* | awk '{print $9}'`

ALL_ID=""
ALL_KEY=""
ALL_PID=""

for m in $f 
do
	
	NUM_MATCH=`cat $m | grep -v "#" | grep $SMOC | wc -l`
	
	if [ $NUM_MATCH -gt 1 ];then	
		PID_FIELD=`FIELD_HEADER $m Project_ID`
		ALL_KEY=$ALL_KEY" "$m
		ID=`basename $m _key.txt`
		ALL_ID=$ALL_ID" "$ID
		PID=`cat $m | grep -v "#" | sed 1d |  awk -F"\t" -v PID_FIELD=$PID_FIELD '{print $PID_FIELD}'`
		ALL_PID=$ALL_PID" "$PID
	fi

done 


echo $ALL_PID | tr ' ' '\n' | sort | uniq 
echo $ALL_ID
echo $ALL_KEY

for ID in $ALL_ID
do

	QSUB_OUT=$ID"_qsubFile.txt"
	echo source "/broad/software/scripts/useuse" > $QSUB_OUT
	echo reuse -q UGER >> $QSUB_OUT
	echo reuse -q Java-1.8 >> $QSUB_OUT
	echo reuse -q Samtools >> $QSUB_OUT
	echo reuse -q GCC-5.2 >> $QSUB_OUT
	echo reuse -q R-3.2 >> $QSUB_OUT
	echo reuse -q Python-2.7 >> $QSUB_OUT

	echo "/broad/IDP-Dx_work/nirmalya/pipeline/beta/PipelineMain.py -c /broad/IDP-Dx_storage/MOC/config_files/Default_config.yaml --key_path /broad/IDP-Dx_storage/MOC/Key_files/"$ID"_key.txt --seq_id "$SMOC" --do_patho --ADD3 30 --ADD5 20 $@" >> $QSUB_OUT

	#qsub -q long $QSUB_OUT

done 

