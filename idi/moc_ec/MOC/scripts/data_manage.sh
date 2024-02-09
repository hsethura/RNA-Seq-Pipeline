#!/bin/sh

days=$1

# get path of the current file
	file_path="${BASH_SOURCE[0]}"
# if the file path is relative, convert it to absolute path
	if [[ $file_path != /* ]]; then
	  file_path="$PWD/${BASH_SOURCE[0]}"
	fi

scripts_dir="$(dirname $file_path)"

# source idi/moc_ec/MOC/scripts/bash_header
	source "$scripts_dir/bash_header"


### source all functions 
# source "idi/moc_ec/MOC/scripts/MOC_functions.sh"
	source "$scripts_dir/MOC_functions.sh"


### determining paths and headers 
	paths_and_headers $ID $@


### set options
	DEFAULT_CONFIG_PATH="$(dirname $(dirname $file_path))"/config_files/PC_config.yaml
	CONFIG_FILE=`extract_option -conf $DEFAULT_CONFIG_PATH 1 $@`
	if [[ $CONFIG_FILE != /* ]]; then
	  CONFIG_FILE="$PROJECT_ROOT_DIR/$CONFIG_FILE"
	fi

### Identify projects that have been transferred in > $ $days
	echo "Identify projects that have been transferred in > $ $days"
	TEMP_FILE=$TEMP_PATH"/data_email_temp.txt"
	TEMP_FILE1=$TEMP_PATH"/data_email_temp1.txt"
	TEMP_FILE2=$TEMP_PATH"/data_email_temp2.txt"

	ALL_IDs=`cat $DATA_DB | sed 1d | awk '{print $2}' | sort | uniq`

	echo "" | sed 1d > $TEMP_FILE
	for ID in $ALL_IDs
	do
		cat $DATA_DB | grep $ID | grep "Transfer"| tail -1 >> $TEMP_FILE
	done

	stamp=`date +%s`
	echo $stamp

	cat $TEMP_FILE | sed 1d | awk '{ printf "%.0f\t", (('$stamp'-$4)/60000)+1; print $0 }'
	cat $TEMP_FILE | sed 1d | awk '{ if((('$stamp'-$4)/60000)+1 > '$days') 
										{									
											printf "%.0f\t", (('$stamp'-$4)/60000)+1
											print $0
										}
									}' > $TEMP_FILE1							


	echo "" | sed 1d > $TEMP_FILE2
	ALL_IDs=`cat $TEMP_FILE1 | awk '{print $3}' | sort | uniq`

### Update DATA_DB and send warning emails
	echo "Updating DATA_DB and sending warning emails"

	for ID in $ALL_IDs
	do
		cat $TEMP_FILE1 | grep $ID | tail -1 >> $TEMP_FILE2
	done

	while read line
	do 
		typeset -i TIME_GAP
		TIME_GAP=`echo $line | awk '{print $1}'`
		MOC_ID=`echo $line | awk '{print $3}'`
		DATE=`echo $line | awk '{print $4}'`
		COLLAB_EMAIL=`echo $line | awk '{print $6}'`
		DATA_EMAIL=`echo $line | awk '{print $7}'`
		INTERNAL=`echo $line | awk '{print $8}'`

		TODAY=`date "+%D"`
		edate=`date +%s`

		IN_MESSAGE=$FILE_PATH"/Data_deletion_warning.txt"
		OUT_MESSAGE=$TEMP_PATH"/"$MOC_ID"_"$USID"_OUT_MESSAGE.txt"
	
		touch $OUT_MESSAGE
		if [ $TIME_GAP -gt 100 ];then
			GAP=">100"
  			cat $IN_MESSAGE | sed 's* on DATE**g' | sed 's/TIME_GAP/'$GAP'/g' | sed 's/DATA_EMAIL/'$DATA_EMAIL'/g' | sed 's/COLLAB_EMAIL/'$COLLAB_EMAIL'/g' | sed 's/MOC_ID/'$MOC_ID'/g' | sed 's/PROJ_ID/'$PROJ_ID'/g' | sed 's*RES_DIR*'$RES_DIR'*g' | sed 's*G_DIR_ADDRESS*'$G_DIR_ADDRESS'*g' |  sed 's*TAR_DIR*'$TAR_DIR'*g' > $OUT_MESSAGE
		
		else
			GAP=$TIME_GAP
  			echo $DATE
  			cat $IN_MESSAGE | sed 's*DATE*'$DATE'*g' | sed 's/TIME_GAP/'$GAP'/g' | sed 's/DATA_EMAIL/'$DATA_EMAIL'/g' | sed 's/COLLAB_EMAIL/'$COLLAB_EMAIL'/g' | sed 's/MOC_ID/'$MOC_ID'/g' | sed 's/PROJ_ID/'$PROJ_ID'/g' | sed 's*RES_DIR*'$RES_DIR'*g' | sed 's*G_DIR_ADDRESS*'$G_DIR_ADDRESS'*g' |  sed 's*TAR_DIR*'$TAR_DIR'*g' > $OUT_MESSAGE
		fi
		
		export TMPDIR=$TEMP_PATH
		echo ""
		echo "Sending email for $MOC_ID..."
		cat $OUT_MESSAGE |  mailx -s "*** IMPORTANT: Data deletion warning for $MOC_ID ***" $USID"@broadinstitute.org"

		echo "Warning"	$MOC_ID	$TODAY	$edate	$COLLAB_EMAIL	$DATA_EMAIL	$INTERNAL	$USID >> $DATA_DB
	done < $TEMP_FILE2

### Identify all files to delete 
	
	echo "Identifying files to delete..."
	############## Import google sheet databases to server ###############=
	if [ $IMPORT_GS == "Y" ];then
		echo "Running $GSIMPORT_SCRIPT to import google sheet databases to server"
		sh $GSIMPORT_SCRIPT
	fi
	########################################################

	PC_DB=`ls -lrt $PCDB_DIR* | awk '{print $9}' | tail -1`

	TEMP_OUTFILE=$TEMP_PATH"/data_delete_temp.txt"
	OUTFILE=$TEMP_PATH"/data_delete.txt"

	time_id $DATA_DB $days $TEMP_OUTFILE "Warning"

	ALL_MOCPIDS=`cat $TEMP_OUTFILE | awk '{print $1}' | sort | uniq`

	echo "Project_status	MOCP_ID	fastq_path	Days_since_last_warning" > $OUTFILE
	echo $TEMP_OUTFILE


	for MOCPID in $ALL_MOCPIDS
	do
		PROJ_STATUS=`Index $PC_DB "MOCP_ID" "Project_Status" $MOCPID`
		cat $TEMP_OUTFILE | grep $MOCPID | awk -v PROJ_STATUS=$PROJ_STATUS '{print  PROJ_STATUS, $0}' >> $OUTFILE
		echo "" >> $OUTFILE
	done 

	ls -lrt $OUTFILE
	ls -lrt $DATA_DB