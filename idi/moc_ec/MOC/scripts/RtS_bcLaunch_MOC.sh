#!/bin/sh


PROJID=$1

SCRIPT_PATH="/broad/IDP-Dx_storage/MOC/scripts/"
OUTPUT_PATH="/broad/IDP-Dx_storage/RawSeqData/MOC/"
OUT_DIR="/broad/IDP-Dx_storage/MOC/RawSeqData/"$PROJID"/"
TEMP_DIR="/broad/hptmp/Dx/"$PROJID"/"




f=`ls $OUT_DIR/* | grep "1.fastq" | grep -v barcode`

echo $f 

for m in $f
do
	ROOT=`basename $m  | sed -e 's/.gz//g' -e 's/.fastq//g'`
	echo "qsub -l h_vmem=8g  -b y sh /broad/IDP-Dx_storage/BactRAP/scripts/RtS_bcSplit.sh "$PROJID" "$ROOT" "$m" -MEM 2000"
done