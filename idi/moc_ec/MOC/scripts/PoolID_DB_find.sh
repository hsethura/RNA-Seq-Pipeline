#!/bin/sh

MOC_ID_STRING=$1

#####################################################################################

ALL_HEADERS="Sample_type Low_Input_Day_1 RiboZero_Request Total_cycles Desired_depth"
POOLID_HEADER="Pool_ID"
SAMPID_HEADER="Sample_ID"

#####################################################################################

ALL_MOC_IDS=`echo $MOC_ID_STRING | sed 's/,/ /g'`

### source all functions 
source "/broad/IDP-Dx_storage/MOC/scripts/MOC_functions.sh"

CONFIG_FILE=`extract_option -conf "/broad/IDP-Dx_storage/MOC/config_files/Universal_config.yaml" 1 $SCRIPT_OPTIONS`
Q_HEAD=`extract_option -q MOC_ID 1 $SCRIPT_OPTIONS`

DB_SCRIPT=`config_read $CONFIG_FILE DB_script`
KEY_SCRIPT=`config_read $CONFIG_FILE keyfile_importer`
KEY_DIR=`config_read $CONFIG_FILE Key_base`

for MOC_ID in $ALL_MOC_IDS
do

	G_ID=`sh $DB_SCRIPT $Q_HEAD,$MOC_ID Google_ID | sed 1d | awk '{print $2}'`
	echo $G_ID
	if [ -z $G_ID ];then

		echo "GID not found in DB for "$MOC_ID
		exit
	fi

	echo "Moving key file to server..."
	
	echo $KEY_SCRIPT -s $G_ID -t "Sample Information" -p $MOC_ID --Key_dir $KEY_DIR
	#$KEY_SCRIPT -s $G_ID -t "Sample Information" -p $MOC_ID --Key_dir $KEY_DIR
	KEY_FILE=$KEY_DIR$MOC_ID"_key.txt"

	POOLID_F=`FIELD_HEADER $KEY_FILE $POOLID_HEADER`
	SAMPID_F=`FIELD_HEADER $KEY_FILE $SAMPID_HEADER`
	
	if [ -z $POOLID_F ];then
		echo $KEY_FILE" does not contain" \"$POOLID_HEADER\" "in the header row"
		continue
	fi
	
	
	ALL_POOLIDS=`cat $KEY_FILE | grep -v "#" | sed 1d | awk -F"\t" '{print $'$POOLID_F'}' | sort | uniq`
	
	printf "MOC_ID\tPool_ID\t"

	for HEADER in $ALL_HEADERS
	do
		printf "%s\t" $HEADER
	done
	echo "Num_Samples"
	
	for POOL_ID in $ALL_POOLIDS
	do
		printf "%s\t%s\t" $MOC_ID $POOL_ID
		for HEADER in $ALL_HEADERS
		do
			H_F=`FIELD_HEADER $KEY_FILE $HEADER`
			printf "%s\t" `FIND_FIELD $KEY_FILE $POOLID_F $POOL_ID $H_F | tr '\n' ' \t'`

		
		done
		printf "%s\n" `FIND_FIELD $KEY_FILE $POOLID_F $POOL_ID $SAMPID_F | sed 's/,/ /g' | wc -w`
			
	
	done
done

