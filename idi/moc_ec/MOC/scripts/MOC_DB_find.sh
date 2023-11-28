#!/bin/sh

Q=$1

### source all functions 
# source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

# get path of the current file. if the file path is relative, convert it to absolute path
file_path="${BASH_SOURCE[0]}"
if [[ $file_path != /* ]]; then
  file_path="$PWD/${BASH_SOURCE[0]}"
fi

PROJECT_ROOT_DIR="$(dirname $(dirname $(dirname $(dirname $(dirname $file_path)))))"

# get parent directory
scripts_dir="$(dirname $file_path)"

source "$scripts_dir/MOC_functions.sh"

### determining paths and headers 
### default config file is /idi/moc_ec/MOC/config_files/Universal_config.yaml
paths_and_headers $MOC_ID $@

# CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/Universal_config.yaml" 1 $@`
DEFAULT_CONFIG_PATH="$(dirname $(dirname $file_path))"/config_files/PC_config.yaml
CONFIG_FILE=`extract_option -conf $DEFAULT_CONFIG_PATH 1 $@`
if [[ $CONFIG_FILE != /* ]]; then
  CONFIG_FILE="$PROJECT_ROOT_DIR/$CONFIG_FILE"
fi

KEY_DIR=`config_read $CONFIG_FILE Key_base`
KEY_ONLY=`extract_option -key_only N 1 $@`
IMPORT_GS=`extract_option -import_gs N 1 $@`
DB_CHOICE=`extract_option -db - 1 $@`

NUM_Q=`echo $Q | sed 's/,/ /g' | wc -w`

if [ $NUM_Q == 1 ];then
	Q_HEAD="MOC_ID"
	ALL_QVAL=$Q
	MOC_ID=$Q
else
	Q_HEAD=`echo $Q | awk '{split($1, ar, ","); print ar[1]}'`
	MOC_ID=`echo $Q | awk '{split($1, ar, ","); if(ar[1]=="MOC_ID") print ar[2]}'`
	ALL_QVAL=`echo $Q | awk '{y=split($1, ar, ","); for(i=2; i <y; i++) printf "%s ", ar[i];print ar[y]}'`
fi

SUB_FILE=`ls -lrt $SUB_DIR/*$SUB_SUFF | tail -1 | awk '{print $9}'`
POOL_FILE=`ls -lrt $POOL_DIR/*$POOL_SUFF | tail -1 | awk '{print $9}'`
DB_FILE=`ls -lrt $MOCDB_DIR/*$MOCDB_SUFF | tail -1 | awk '{print $9}'`
RTS_FILE=`ls -lrt $RTSDB_DIR/*$RTSDB_SUFF | tail -1 | awk '{print $9}'`

ALL_RETURN_HEADERS=""

############## Import google sheet databases to server ###############=
if [ $IMPORT_GS == "Y" ];then
	echo "Running $GSIMPORT_SCRIPT to import google sheet databases to server"
	sh $GSIMPORT_SCRIPT
fi
########################################################


LAST="N"
for m in $@
do	
	HEADER=`echo $m | awk -v LAST=$LAST '{if($1 !~ /^-/ && LAST == "N" ) print $1}'`
	LAST=`echo $m | awk '{if($1 ~ /^-/) print "Y"; else print "N"}'`
	ALL_RETURN_HEADERS=$ALL_RETURN_HEADERS" "$HEADER
done

NUM_HEADERS=`echo $ALL_RETURN_HEADERS | awk '{print NF}'`

if [ $DB_CHOICE == "-" ];then
	ALL_DBS=$DB_FILE" "$RTS_FILE
else

	ALL_DBS=$DB_CHOICE
fi

ALL_DB_HEADERS=" "
if [ $NUM_HEADERS == 0  ];then
	for DB in $ALL_DBS
	do
		DB_HEADERS=`cat $DB | sed '/^$/d' | grep -v "###" | head -1 | sed 's/ /_/g'`
		ALL_DB_HEADERS=$ALL_DB_HEADERS" "$DB_HEADERS
	done
		
	ALL_RETURN_HEADERS="$ALL_DB_HEADERS Index1_Seq Index2_Seq SMOC_ID"
fi


ALL_QVALS=`echo $Q_VAL | sed 's/;/ /g'`

for Q_VAL in $ALL_QVAL
do
	printf "%s:\t" $Q_HEAD
	echo $Q_VAL
			 
	for R_HEAD in $ALL_RETURN_HEADERS
	do
		printf "%s:\t" $R_HEAD | sed 's/_/ /g'		
		if [ $KEY_ONLY == "N" ];then
			KEY_FILE=$KEY_DIR"/"`FIND_HEADER_AND_FIELD $Q_VAL $Q_HEAD Key_file $ALL_DBS | tr ',' '\n' | sort | uniq | tr '\n' ',' | sed 's/^,//g' | sed 's/,$//g' | sed -e 's/Key/key/g' -e 's/KEY/key/g'`".txt"


			if [ ! -s $KEY_FILE ] && [ $KEY_ONLY == "Y" ];then

				echo "$KEY_FILE not found!"
				exit
			fi
			FIND_HEADER_AND_FIELD $Q_VAL $Q_HEAD $R_HEAD $ALL_DBS | tr ';' ',' | tr ',' '\n' | sort | uniq | tr '\n' ',' | sed 's/^,//g' | sed 's/,$//g' | sed 's/_/ /g'
		else
			KEY_FILE=$KEY_DIR"/"$MOC_ID"_key.txt"
			FIND_VALUES_FOR_HEADER $R_HEAD $KEY_FILE | sort | uniq | tr '\n' ',' | sed 's/,$//g' 

		fi
		echo ""
	done
done


