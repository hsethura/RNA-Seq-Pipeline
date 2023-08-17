#!/bin/sh

BC_db="/broad/IDP-Dx_storage/MOC/files/rts_bcs.txt"
DIST_FILE=$1

ALL_BC=`cat $DIST_FILE | awk '{print $1}'`

echo $ALL_BC

for BC in $ALL_BC
do

	ID=`cat $BC_db | grep $BC | awk '{print $1}'`
	STATS=`cat $DIST_FILE | grep $BC `
	echo $STATS, $ID
done | sort -k5
