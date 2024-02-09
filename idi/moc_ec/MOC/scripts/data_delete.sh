
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
	
	MOCP_IDS=`extract_option -m - 1 $@`
	DAYS=`extract_option -d - 1 $@`
	
	if [ $MOCP_IDS == "-" ] && [ $DAYS == "-" ];then
		echo "Please enter MOCP_ID using -m or # of days using -d"
		exit
	fi
	
	
### Identify projects that have been transferred in > $ $days
	
	if [ $MOCP_IDS == "-" ];then
		echo "Identify projects that have been transferred in > $DAYS"
		TEMP_FILE=$TEMP_PATH"/data_email_temp.txt"
		TEMP_FILE1=$TEMP_PATH"/data_email_temp1.txt"
		TEMP_FILE2=$TEMP_PATH"/data_email_temp2.txt"
	
		ALL_IDs=`cat $DATA_DB | sed 1d | awk '{print $2}' | sort | uniq`
	
		echo "" | sed 1d > $TEMP_FILE
		for ID in $ALL_IDs
		do
			cat $DATA_DB | grep $ID | grep "Warning" | grep -v "Deleted" | tail -1 >> $TEMP_FILE
		done
	
		stamp=`date +%s`
		echo $stamp
	
		cat $TEMP_FILE | sed 1d | awk '{ printf "%.0f\t", (('$stamp'-$4)/60000)+1; print $0 }'
		cat $TEMP_FILE | sed 1d | awk '{ if((('$stamp'-$4)/60000)+1 > '$DAYS') 
											{									
												printf "%.0f\t", (('$stamp'-$4)/60000)+1
												print $0
											}
										}' > $TEMP_FILE1							
	
	
		echo "" | sed 1d > $TEMP_FILE2
		ALL_IDs=`cat $TEMP_FILE1 | awk '{print $3}' | sort | uniq`
	else 
	
		ALL_IDs=`echo $MOCP_IDS | sed 's/,/ /g' `
	
	fi
	
	
	
	echo $ALL_IDs
	
	for ID in $ALL_IDs
	do
		echo $ID
		REMOVE_FILE=$TEMP_PATH"/"$ID"_remove.txt"
		

		du -sh $RAWSYM_PATH"/"$ID
		echo "df -h /idi/moc_ec" > $REMOVE_FILE
		ls -lrt $RAWSYM_PATH"/"$ID"/"* | grep -e $RAWSYM_PATH | sed 's/://g'| awk '{print "rm -r " $NF}' >> $REMOVE_FILE
 		echo "df -h /idi/moc_ec" >> $REMOVE_FILE
		
		ls -lrt $REMOVE_FILE
		sh $REMOVE_FILE
		
		
		TODAY=`date "+%D"`
		stamp=`date +%s`
		echo $stamp $TODAY
		INFO=`cat $DATA_DB | grep $ID | grep "Warning" | awk '{ print $5,$6,$7}'`
		
		echo "Deleted	"$ID		$TODAY	$stamp	$INFO	$USID >> $DATA_DB
		cat $DATA_DB | grep $ID

	done
	

	
	
	
	
	
	
