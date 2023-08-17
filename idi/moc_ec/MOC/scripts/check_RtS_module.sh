#!/bin/sh

MOC_ID=$1
DIR=$2
ALL_SUFF=$3
CHECK_MOD_FILE=$4

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"
paths_and_headers $MOC_ID $@

DIR_NAME=`basename $DIR`

ALL_SUFF_ARRAY=`echo $ALL_SUFF | sed 's/,/ /g'`

FAIL_EXIT=`extract_option -fail_exit N 1 $@`
EMAIL=`extract_option -email N 1 $@`
HEADER=`extract_option -header Sample_ID 1 $@`
DELIMIT=`extract_option -delimit _ 1 $@`
REF_OPT=`extract_option -ref N 1 $@`

rm $CHECK_MOD_FILE
touch $CHECK_MOD_FILE

ALL_ID=`sh $DB_SCRIPT $Q_HEAD,$MOC_ID -key_only Y $HEADER | tail -1 | sed 's/,/ /g' | awk -F":" '{print $2}' | awk '{for(i=1; i < NF+1; i++)print $i}'`


ls -lrt $KEY_FILE


echo $ALL_ID

NUM_ALL=`echo $ALL_ID | wc -w`
NUM_SUFF=`echo $ALL_SUFF_ARRAY | wc -w`

typeset -i NUM_MISSING
ALL_MISSING=""
NUM_MISSING=0

for ID in $ALL_ID
do
	NUM_FOUND=0
	for SUFF in $ALL_SUFF_ARRAY
	do
		if [ $REF_OPT == "Y" ];then
			REF=`FIND_HEADER_AND_FIELD $ID "Sample_ID" "Bacterial_reference" $KEY_FILE | awk '{print $1}'`	
		else
			REF=""
		fi		
		echo $REF
		if [ $DELIMIT == "none" ];then
			FILE=$DIR/$ID$REF$SUFF
		else	
			FILE=$DIR/$ID$DELIMIT$REF$DELIMI$SUFF
		fi
		echo "Looking for "$FILE
		ls -lrt $FILE
		
		if [ -s $FILE ];then
		  	NUM_FOUND=`expr $NUM_FOUND + 1`
		else
			echo $FILE" not found"
		fi
	done
	echo "NUM_FOUND:" $NUM_FOUND
	echo "NUM_SUFF:" $NUM_SUFF
	if [ $NUM_FOUND -lt $NUM_SUFF ];then
		echo $ID "missing ("$NUM_FOUND" of "$NUM_SUFF")" >> $CHECK_MOD_FILE
		NUM_MISSING=`expr $NUM_MISSING + 1`
	fi

done

echo $NUM_MISSING" missing from "$NUM_ALL" files with suffixes $ALL_SUFF in "$DIR >> $CHECK_MOD_FILE

export TMPDIR=$TEMP_PATH

echo $USID

if [ $NUM_MISSING -ne 0 ];then
	if [ $EMAIL == "Y" ] || [ $EMAIL == "F" ];then	
		echo "$NUM_MISSING files missing from $DIR of $MOC_ID"
		echo "sending email..."
		cat $CHECK_MOD_FILE |  mailx -a $CHECK_MOD_FILE -s "Files missing from $DIR_NAME of $MOC_ID" $USID"@broadinstitute.org"
	fi
	
	if [ $FAIL_EXIT == "Y" ];then
		exit
	fi
fi
if [ $NUM_MISSING -eq 0 ];then
	if [ $EMAIL == "Y" ] || [ $EMAIL == "F" ];then	
		echo "All files found in $DIR of $MOC_ID"
		echo "sending email..."
		cat $CHECK_MOD_FILE |  mailx -s "All files found in $DIR_NAME of $MOC_ID" $USID"@broadinstitute.org"
	fi
fi
exit_script="P"
	
	
### change permissions for DIR
change_perms $DIR

	
