#!/bin/sh


ORIGIN_DIR=$1
DEST_PATH=$2


DIR_NAME=`basename $ORIGIN_DIR`
DEST_DIR=$DEST_PATH"/"$DIR_NAME

mkdir -p $DEST_DIR

echo "qsub -l h_vmem=4g -l os=RedHat7 -l h_rt=4:00:00 -b Y mv $ORIGIN_DIR/* $DEST_DIR"