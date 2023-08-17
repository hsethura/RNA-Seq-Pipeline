#!/bin/sh


DIR=$1
WEB_DIR=$2
URL="http://idpweb.broadinstitute.org/data/MOC/"

BAMS=`ls -lrt $DIR*bam | awk '{print $9}'`
cd $DIR

for BAM in $BAMS
do
	ROOT=`basename $BAM`
	echo $URL$ROOT
	qsub -l h_vmem=10g -V -b Y java -jar /broad/IDP-Dx_storage/BactRAP/scripts/bam2tdf-14/bam2tdf.jar $BAM
	mv $BAM".tdf" $WEB_DIR
done

