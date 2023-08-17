!#/bin/sh


PROJ_ID=$1
TEMP_PATH=$2
RESULT_PATH=$3
CONFIG_FILE=$4
SCRIPT_DIR=$5

FILE_DIR="/broad/IDP-Dx_storage/MOC/files/"

## Setting up data for transfer
sh $SCRIPT_DIR"/MOC_data_transfer.sh" $PROJ_ID $CONFIG_FILE Y

echo "Changing permissions in $TEMP_PATH"
chmod -R 777 $TEMP_PATH

echo "Changing permissions in $RESULT_PATH"
chmod -R 777 $RESULT_PATH


