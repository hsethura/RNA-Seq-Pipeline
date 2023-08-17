#!/bin/sh

MOC_ID=$1
shift

### source all functions 
source "/broad/IDP-Dx_storage/MOC/scripts/MOC_functions.sh"


##### FUNCTIONS ######

### check DBs to make sure there are values for required headers
check_DB()
{
	Q_HEAD=$1
	Q_VAL=$2
	R_HEAD=$3
	
	R_VAL=`sh $DB_SCRIPT $Q_HEAD,$MOC_ID $R_HEAD | sed 1d | awk '{print $2}' | tr ',' '\n' | sort | uniq | tr '\n' ',' | sed 's/,$//g'`
	
if [ -z $R_VAL ];then
	echo ""
	echo "	ERROR: No "$R_HEAD" values found for "$Q_VAL  
	echo "########"
else
	echo ""
	echo "	"$R_HEAD "values found: "$R_VAL
fi
}

#### check key file to make sure the required headers are there


key_head_check ()
{
	CONFIG_FILE=$1
	KEY_FILE=$2

	ALL_HEADERS=`cat $CONFIG_FILE | awk -v flag=N '{if(flag=="Y" && $0 !~ /KEY_HEADERS_END/) print $2;if($0 ~ /KEY_HEADERS_START/)flag="Y";if($0 ~ /KEY_HEADERS_END/) exit}'`

	for HEADER in $ALL_HEADERS
	do
		if [ `cat $KEY_FILE | grep -v "#" | head -1 | grep -w $HEADER | wc -l` == 0 ];then
			echo "ERROR: "$HEADER" header is missing!"
		fi 
	done
}

###########################


CONFIG_FILE=`extract_option -conf "/broad/IDP-Dx_storage/MOC/config_files/Universal_config.yaml" 1 $SCRIPT_OPTIONS`
Q_HEAD=`extract_option -q MOC_ID 1 $SCRIPT_OPTIONS`

DB_SCRIPT=`config_read $CONFIG_FILE DB_script`
KEY_SCRIPT=`config_read $CONFIG_FILE keyfile_importer`
KEY_DIR=`config_read $CONFIG_FILE Key_base`

KEY_FILE=$KEY_DIR$MOC_ID"_key.txt"

BACT_REF_PATH=`config_read $CONFIG_FILE Bacterial_Ref_path`
BACT_REF_HEADER=`config_read $CONFIG_FILE Ref_accession`
BACT_REF_FIELD=`FIELD_HEADER $KEY_FILE $BACT_REF_HEADER`


# values required for 
ALL_DB_VALUES="Google_ID SMOC_ID Index1_Seq Index2_Seq"

for DB_VALUE in $ALL_DB_VALUES
do
	check_DB $Q_HEAD $MOC_ID $DB_VALUE
done 


#### Getting G_ID for MOC_ID

G_ID=`sh $DB_SCRIPT $Q_HEAD,$MOC_ID Google_ID | sed 1d | awk '{print $2}'`

if [ -z $G_ID ];then
	echo "ERROR: No GID found for" $MOC_ID
	exit
fi
echo $G_ID

#### moving key file to server

echo "Moving key file to server..."

$KEY_SCRIPT -s $G_ID -t "Sample Information" -p $MOC_ID --Key_dir $KEY_DIR

if [ ! -s $KEY_FILE ];then
	echo "ERROR: No key file moved for" $MOC_ID
	ls -lrt $KEY_FILE
	exit
fi


#### check for presence of all headers listed in $CONFIG_FILE 

key_head_check $CONFIG_FILE $KEY_FILE 

#### check for existence of all Bact_reference sequences and annotations

ALL_BACT_REFS=`cat $KEY_FILE  | grep -v "#" | sed 1d | awk -F"\t" -v BACT_REF_FIELD=$BACT_REF_FIELD '{print $BACT_REF_FIELD}' | sort | uniq`


for BACT_REFS in $ALL_BACT_REFS
do
	if [ ! -s $BACT_REF_PATH$BACT_REFS".gff" ];then
		echo "ERROR: "$BACT_REF_PATH$BACT_REFS".gff is missing!"
	else
		ls -lrt $BACT_REF_PATH$BACT_REFS".gff"

	fi
	if [ ! -s $BACT_REF_PATH$BACT_REFS".fna" ];then
		echo "ERROR: "$BACT_REF_PATH$BACT_REFS".fna is missing!"
	else
		ls -lrt $BACT_REF_PATH$BACT_REFS".fna"
	fi

done

exit


