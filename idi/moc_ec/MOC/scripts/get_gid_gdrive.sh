#!/bin/sh


MOC_ID=$1


### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh" 

CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/Universal_config.yaml" 1 $@`

### determining paths and headers 
### default config file is /idi/moc_ec/MOC/config_files/Universal_config.yaml
paths_and_headers $MOC_ID $@

### set paths to directories, files, and scripts from config file

GMOC_PATH=`config_read $CONFIG_FILE gdrivemoc_path`
GLOCAL_PATH=`config_read $CONFIG_FILE gdrivelocal_path`
GDRIVE_SCRIPT=`config_read $CONFIG_FILE gdrive_script`

GKEY_FILE=$GMOC_PATH"/"$MOC_ID"/"$MOC_ID"_Key"

cd $GLOCAL_PATH

$GDRIVE_SCRIPT id $GKEY_FILE | cut -d'"' -f2 | awk '{print $NF}' | tail -1


