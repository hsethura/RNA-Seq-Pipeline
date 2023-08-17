#!/bin/sh

FPKM_FILE=$1

SCRIPTS_PATH="/broad/IDP-Dx_storage/MOC/scripts/"
OUT_PATH=`dirname $FPKM_FILE`
OUT_DIR=$OUT_PATH"/FPKM_histos/"
TEMP_DIR="/broad/hptmp/MOC/"


mkdir -p $OUT_DIR
mkdir -p $TEMP_DIR

NF=`cat $FPKM_FILE | sed 1d | head -1 | awk '{print NF}'`

i=2

while [ $i -le $NF ]
do

	
	
	ID=`cat $FPKM_FILE | head -1 | awk -v i=$i '{print $i}'`
	TEMP_OUT=$TEMP_DIR$ID"_temp_coords.txt"
	OUTPUT_FILE=$OUT_DIR$ID"_fpkmHisto.pdf"

	cat $FPKM_FILE | sed 1d | awk -v i=$i '{print $i}' > $TEMP_OUT



	Rscript $SCRIPTS_PATH"ReadLen_histo.R" $TEMP_OUT $OUTPUT_FILE $ID

	i=`expr $i + 1`
	
	
	echo $OUTPUT_FILE
	
done

exit


