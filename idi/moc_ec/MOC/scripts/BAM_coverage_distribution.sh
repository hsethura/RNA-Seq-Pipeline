#!/bin/sh

source "/broad/software/scripts/useuse"
reuse -q use igvtools_2.3


TDF_FILE=$1
GFF_FILE=$2


INC=10

TEMP_DIR="/broad/hptmp/Dx/"

ROOT=`basename $TDF_FILE .tdf`


DEPTH_FILE=$TEMP_DIR$ROOT"_depth.txt"
HISTO_FILE=$TEMP_DIR$ROOT"_histo.txt"
GENE_FILE=$TEMP_DIR$ROOT"_gene.txt"
COV_DIST_FILE=$TEMP_DIR$ROOT"_covDist.txt"

#igvtools tdftobedgraph $TDF_FILE $DEPTH_FILE

cat $GFF_FILE | grep -v "#" | awk '{split($9, ar, ";"); print ar[1], $4, $5, $7}' > $GENE_FILE
SIZE=`cat $DEPTH_FILE | awk '{print $3}' | tail -1`


echo "/broad/IDP-Dx_storage/MOC/scripts/COV_dist2"

/broad/IDP-Dx_storage/MOC/scripts/COV_dist2 $GENE_FILE $DEPTH_FILE $SIZE 10 > $COV_DIST_FILE

echo $COV_DIST_FILE

typeset -i NB
typeset -i i

NB=`cat $COV_DIST_FILE | awk '{print NF}' | sort | uniq`

i=4

echo "" 

while [ $i -le $NB ]
do	
	TOTAL_BIN=`awk  -v i=$i -v total=0 '{
			
					total=total+($i*$2)
					print total/NR
	
		}' $COV_DIST_FILE | tail -1`
	
	printf "%s\t" $TOTAL_BIN
	i=`expr $i + 1`

done

echo ""