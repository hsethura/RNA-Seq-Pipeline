#!/bin/sh

source /idi/moc_ec/MOC/scripts/bash_header


ID=$1


### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

### determining paths and headers 
paths_and_headers $ID $@

ALL_REFS_OPT=`extract_option -all_refs - 1 $@`
GREF=`extract_option -gref Y 1 $@`

### set paths to directories, files, and scripts from config file
### default config file is /idi/moc_ec/MOC/config_files/Universal_config.yaml

CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/PC_config.yaml" 1 $@`

STOCK_LOG_PATH=`config_read $CONFIG_FILE stock_log_path`
STOCK_PLATE_MAP_PATH=`config_read $CONFIG_FILE stock_plate_map_path`
SAMPLE_PLATE_MAP_PATH=`config_read $CONFIG_FILE sample_plate_map_path`
STOCK_LOG_LOCAL=`config_read $CONFIG_FILE stock_log_local`
LINKED_STOCK_PATH=`config_read $CONFIG_FILE linked_stock_maps_path`
LINKED_OUT_PATH=`config_read $CONFIG_FILE linked_map_out_path`

GPD_PATH=`config_read $CONFIG_FILE gdriveDEV_path`
GTRACK_PATH=`config_read $CONFIG_FILE gdrivetrack_path`
GLOCAL_PATH=`config_read $CONFIG_FILE gdrivelocal_path`
GDRIVE_SCRIPT=`config_read $CONFIG_FILE gdrive_script`
TEMP_DIR=`config_read $CONFIG_FILE Temp_path`

mkdir -p $STOCK_LOG_LOCAL

TEMP_FILE=$TEMP_DIR"/"temp.txt

link ()
{
		LOG_FILE=$1
		PLATE_MAP=$2
		LINKED_FILE=$3
		OPTION=$4
		TEMP_FILE=$5
		
		echo "Linking log file $LOG_FILE and plate map $PLATE_MAP to create $LINKED_FILE..."
		
		cat $PLATE_MAP | sed 1d  > $TEMP_FILE
		
		if [ -f $LOG_FILE ] && [ -f $PLATE_MAP ];then

			while read line 
			do
				ROW=`echo $line | awk '{print $5}' | cut -d"," -f3 | sed 's/,/ /g'`
				COL=`echo $line | awk '{print $5}' | cut -d"," -f4 | sed 's/,/ /g'`
				BC=`echo $line | awk '{print $5}' | cut -d"," -f5 | sed 's/,/ /g'`
				
				if [ $OPTION == "Stock" ];then
					POSIT=$ROW$COL
					SPEC=`awk -F"\t" -v POSIT=$POSIT '{if(POSIT == $3) print $4}' $LOG_FILE`
					echo $NAME"	"$ROW"	"$COL"	"$BC"	"$SPEC | sed 's/ /_/g' >> $LINKED_FILE
				fi
				if [ $OPTION == "Plate" ];then
					SPEC=`cat $SPEC_LINKED_FILE | grep -w $BC | awk '{print $NF}'`
					echo $ROW"	"$COL"	"$BC"	"$SPEC | sed 's/ /_/g' >> $LINKED_FILE
				fi
			done < $TEMP_FILE
		else
			echo $LOG_FILE "or" $PLATE_MAP" not found!"

		fi

}



### Making local gdrive dir and move into it
GLOCAL_DIR=$GLOCAL_PATH"/"$STOCK_LOG_PATH
mkdir -p $GLOCAL_DIR
cd $GLOCAL_DIR
echo $GLOCAL_DIR


### Pulling all stock log desktop files from gdrive directory $STOCK_LOG_PATH 
echo ""
echo "Pulling all stock log desktop files from gdrive directory $STOCK_LOG_PATH"
$GDRIVE_SCRIPT pull -quiet

cd $GLOCAL_DIR
echo $GLOCAL_DIR


### Pulling all stock log sheets files from gdrive directory $STOCK_LOG_PATH using gids in desktop files
echo "Pulling all stock log sheets files from gdrive directory $STOCK_LOG_PATH using gids in desktop files"
echo ""

ALL_STOCKS=`ls -lrt | awk '{print $9}' | grep desktop`
SUFF="_log.txt"

echo $ALL_STOCKS

for STOCK in $ALL_STOCKS
do
	NAME=`echo $STOCK | rev | cut -d"/" -f1 | rev | sed 's/.desktop//g'`
	ID=`cat $STOCK | grep URL | cut -d"/" -f6`
	
	echo $SCRIPTS_DIR"GS_import.py" -s $ID -t \"Stock\" -p $NAME --Key_dir $STOCK_LOG_LOCAL -S $SUFF
 	$SCRIPTS_DIR"GS_import.py" -s $ID -t "Stock" -p $NAME --Key_dir $STOCK_LOG_LOCAL -S $SUFF
done

### Linking all barcodes in stock plate maps in $STOCK_PLATE_MAP_PATH to species names in cognate stock log files 
echo ""
echo "Linking all barcodes in stock plate maps in $STOCK_PLATE_MAP_PATH to species names in cognate stock log files" 
echo ""

mkdir -p $LINKED_STOCK_PATH
SPEC_LINKED_FILE=$LINKED_STOCK_PATH"/BMB_spec_linked.txt"
echo "" | sed 1d > $SPEC_LINKED_FILE

ALL_STOCK_PLATES=`ls -lrt $STOCK_PLATE_MAP_PATH* | awk '{print $9}' | grep "stock"`

ls -lrt $ALL_STOCK_PLATES

for STOCK_PLATE in $ALL_STOCK_PLATES
do
	NAME=`echo $STOCK_PLATE | rev | cut -d"/" -f1 | rev | sed 's/.desktop//g'`
	STOCK_LOG_FILE=$STOCK_LOG_LOCAL"/"$NAME$SUFF
	echo "link $STOCK_LOG_FILE $STOCK_PLATE $SPEC_LINKED_FILE Stock $TEMP_FILE"
	link $STOCK_LOG_FILE $STOCK_PLATE $SPEC_LINKED_FILE "Stock" $TEMP_FILE
	
done
  
ls -lrt $SPEC_LINKED_FILE


### Making output directory for linked sample files in local gdrive
  GLOCAL_OUT_DIR=$GLOCAL_PATH"/"$LINKED_OUT_PATH
mkdir $GLOCAL_OUT_DIR


### Making sample plate maps linked to species and moving them to $GLOCAL_OUT_DIR
echo ""
echo "Making sample plate maps linked to species and moving them to $GLOCAL_OUT_DIR\n"
echo ""
 
cd $GLOCAL_OUT_DIR
pwd
###### ISSUES ####

ALL_SAMP_PLATES=`ls -lrt $SAMPLE_PLATE_MAP_PATH* | awk '{print $9}'`

for SAMP_PLATE in $ALL_SAMP_PLATES
do
	
	NAME=`echo $SAMP_PLATE | rev | cut -d"/" -f1 | rev | sed 's/.desktop//g'`
	OUT_FILE=$GLOCAL_OUT_DIR"/"$NAME"_linked.tsv"
	echo "Row	Column	Barcode	Strain" > $OUT_FILE
	
	echo ""
	echo "Linking barcodes to species names in $SAMP_PLATE using $SPEC_LINKED_FILE and writing it to file $OUT_FILE "
	echo ""
	
	echo "link $SPEC_LINKED_FILE $SAMP_PLATE $OUT_FILE Plate $TEMP_FILE"
	link $SPEC_LINKED_FILE $SAMP_PLATE $OUT_FILE "Plate" $TEMP_FILE
	
	cat $OUT_FILE > $TEMP_FILE
	echo "Plate_ID:	"$NAME > $OUT_FILE
	echo  "" >> $OUT_FILE
	echo "Row	Column	Barcode	Strain" >> $OUT_FILE
	cat $TEMP_FILE >> $OUT_FILE
	ls -lrt $OUT_FILE
done

$GDRIVE_SCRIPT push -quiet -ignore-name-clashes

### change permissions for Results and temp dirs
change_perms $TEMP_DIR $STOCK_LOG_PATH


ALL_BC_STOCK=`cat $SPEC_LINKED_FILE | awk '{print $4}' | grep -v Barcode | grep -v EMPTY `
ALL_BC_USED=`cat $GLOCAL_OUT_DIR/*tsv | awk '{print $3}' | grep -v Barcode | grep -v EMPTY `
ALL_BC_REMAIN=`echo $ALL_BC_USED $ALL_BC_USED $ALL_BC_STOCK | tr ' ' '\n' | sort | uniq -c | sort | awk '{if($1==1) print $2}'`
echo $ALL_BC_REMAIN

for BC in $ALL_BC_REMAIN
do
	cat $SPEC_LINKED_FILE | grep $BC

done


ls -lrt $GLOCAL_OUT_DIR

ls -lrt $SPEC_LINKED_FILE
