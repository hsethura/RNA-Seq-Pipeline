#!/bin/sh


ID=$1


### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

### determining paths and headers 
paths_and_headers $ID $@

ALL_REFS_OPT=`extract_option -all_refs - 1 $@`
GREF=`extract_option -gref Y 1 $@`

### set paths to directories, files, and scripts from config file
### default config file is /idi/moc_ec/MOC/config_files/Universal_config.yaml

CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/Universal_config.yaml" 1 $@`

GPC_PATH=`config_read $CONFIG_FILE gdrivePC_path`
GPD_PATH=`config_read $CONFIG_FILE gdriveDEV_path`
GTRACK_PATH=`config_read $CONFIG_FILE gdrivetrack_path`
GLOCAL_PATH=`config_read $CONFIG_FILE gdrivelocal_path`
GDRIVE_SCRIPT=`config_read $CONFIG_FILE gdrive_script`
TEMP_DIR=`config_read $CONFIG_FILE temp_dir`
Q_PATH=`config_read $CONFIG_FILE qiime_path`
Q_DIR=$Q_PATH"/"$ID

if [ ! -d $Q_DIR ];then
	echo "********"
	echo $Q_DIR" does not exist!"
	echo "********"
	exit
fi

TYPE=`echo $ID | awk '{if($1 ~ /MOC/) print "P"; else print "16S" }'`

### use -gd_path option to set path to gdrive dir
GD_PATH=`extract_option -gd_path - 1 $@`
MOC_DIR=`echo $GPD_PATH | cut -d"/" -f1,2`

### set local gdrive path
if [ $GD_PATH == "-" ];then
	if [ $TYPE == "P" ];then
		DIR_PATH=$ID"/qiime_results"
		GLOCAL_DIR=$GLOCAL_PATH"/"$GPC_PATH"/"$DIR_PATH
		GTRANS_DIR=$GPC_PATH"/"$DIR_PATH
	else
		
		DIR_PATH="16S/16S_data/"$ID"/qiime_results"
		GLOCAL_DIR=$GLOCAL_PATH"/"$GPD_PATH"/"$DIR_PATH
		GTRANS_DIR=$GPD_PATH"/"$DIR_PATH
	fi
else
	TYPE="D"
	DIR_PATH=$ID"/16S_data/qiime_results"
	GLOCAL_DIR=$GLOCAL_PATH"/"$GD_PATH"/"$DIR_PATH
	GTRANS_DIR=$GD_PATH"/"$DIR_PATH

fi

echo $GLOCAL_DIR
echo $GTRANS_DIR

mkdir -p $GLOCAL_DIR
cd $GLOCAL_DIR
ln -s $Q_DIR"/"* $GLOCAL_DIR
ls -lrt $GLOCAL_DIR
echo $GLOCAL_DIR

echo $GDRIVE_SCRIPT push -quiet -ignore-name-clashes
$GDRIVE_SCRIPT push -quiet -ignore-name-clashes

echo ""
echo "Files moved from "$GLOCAL_DIR "...."
echo "to "$GTRANS_DIR
echo ""

$GDRIVE_SCRIPT ls 
