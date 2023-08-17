#!/bin/sh

MOC_ID=$1 
KEY=$2


source /idi/moc_ec/MOC/scripts/bash_header

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

### determining paths and headers 
### default config file is /idi/moc_ec//MOC/config_files/PC_config.yaml
paths_and_headers $MOC_ID $@

GID=`extract_option -gid - 1 $@`
echo $GID $KEY_DIR

KEY_FILE=$KEY_DIR$MOC_ID"_key.txt"

if [ $GID != "-" ];then
	echo 
	"wget --no-check-certificate --output-document=$KEY_FILE 'https://docs.google.com/spreadsheet/ccc?key='$KEY'&output=txt&gid='$GID''"
else 
	wget --no-check-certificate --output-document=$KEY_FILE 'https://docs.google.com/spreadsheet/ccc?key='$KEY'&output=txt'
fi


ls -lrt $KEY_FILE