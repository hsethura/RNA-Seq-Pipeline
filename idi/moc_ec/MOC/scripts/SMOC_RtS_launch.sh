#!/bin/sh

SMOC_ID=$1
shift

ALL_MOCS=`sh /broad/IDP-Dx_storage/MOC/scripts/MOC_DB_find.sh "SMOC_ID,"$SMOC_ID | grep -w "MOC ID:" | awk -F"\t" '{ print $2}' | sed 's/,/ /g'`

for MOC in $ALL_MOCS
do
	echo sh /broad/IDP-Dx_storage/MOC/scripts/RtS_Analysis_launch.sh $MOC $@
	
done
