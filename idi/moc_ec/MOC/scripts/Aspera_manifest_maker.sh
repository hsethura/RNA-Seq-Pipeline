#!/bin/sh

MOC_ID=$1
shift
BAM_DIR=$1
shift
TAR_DIR=$1
shift
MERGE_DIR=$1
shift
SUFF_STR=$1
shift

ALL_SUFF=`echo $SUFF_STR | sed 's/,/ /g'`

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

### source all functions 
# source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"
source "$scripts_dir/MOC_functions.sh"


USID=`USID`
PID=`echo $BAM_DIR | rev | cut -d"/" -f2 | rev`


### get options from command line
# -fastq: include fastq files along with bams in Aspera (default N)

# CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/Universal_config.yaml" 1 $@`
DEFAULT_CONFIG_PATH="$(dirname $(dirname $file_path))"/config_files/PC_config.yaml
CONFIG_FILE=`extract_option -conf $DEFAULT_CONFIG_PATH 1 $@`
if [[ $CONFIG_FILE != /* ]]; then
  CONFIG_FILE="$PROJECT_ROOT_DIR/$CONFIG_FILE"
fi
TRANS_FASTQ=`extract_option -fastq Y 1 $@`
TRANS_BAMS=`extract_option -bams Y 1 $@`
USER_ID=`extract_option -user_id N 1 $@`

### get paths etc. from config file
TEMP_PATH=`config_read $CONFIG_FILE Temp_path`
MAN_DIR=`config_read $CONFIG_FILE manifest_dir`

echo $PID
FLAG="Y"

mkdir -p $MAN_DIR

OUT_FILE=$MAN_DIR"/"$PID"_manifest.txt"
TEMP_FILE=$MAN_DIR"/"$PID"_manifest_temp.txt"
MESSAGE_FILE=~/"message.txt"

echo "" | sed 1d > $OUT_FILE
echo "" | sed 1d > $TEMP_FILE

echo $BAM_DIR

if [ $TRANS_BAMS == "Y" ];then
	for SUFF in $ALL_SUFF
	do
		ls -lrt $BAM_DIR*$SUFF | grep -v manifest | awk '{print $9}' | sed '/^$/d' >> $TEMP_FILE
	done
fi

if [ $TRANS_FASTQ == "Y" ];then
	ls -lrt $MERGE_DIR"/"*fastq* | awk '{print $9}' >> $TEMP_FILE
fi
  
f=`cat $TEMP_FILE | awk '{print $1}'`

for m in $f 
do
	if [ ! -s $m ];then
	
		echo "File $m does not exist!" 
		#FLAG="N" 
	fi 
done

if [ $FLAG == "Y" ];then
	echo ""  >  $OUT_FILE
	cat $TEMP_FILE | sed -e 's*//*/*g' -e 's*///*/*g' >> $OUT_FILE
	echo $TAR_DIR >> $OUT_FILE
fi
echo "Sending email to $USID"

echo "Login at the link below, select NEW, enter your username for Requester and Owner, and MOCP-ID for Request Name " > $MESSAGE_FILE
echo "Then copy/paste paths (below and in attachment) into dialogue box and submit" >> $MESSAGE_FILE
echo "You'll get a bunch of emails, just wait for the one that says it is ready to transfer and follow the instructions'" >> $MESSAGE_FILE
echo "You'll then get a bunch more emails, forward the one that has a hashtag box containing the text 'Please download and read the instructions at'" >> $MESSAGE_FILE
echo "Forward that to the collab and enter the date in the data tab of the MOC DB" >> $MESSAGE_FILE
echo "Link to aspera site: https://transfer.broadinstitute.org/auth/login" >> $MESSAGE_FILE

export TMPDIR=$TEMP_PATH

cat $MESSAGE_FILE $OUT_FILE |  mailx -a $OUT_FILE -s "Aspera request for $MOC_ID $PID " $USID"@broadinstitute.org"

ls -lrt $OUT_FILE
ls -lrt $TEMP_FILE
