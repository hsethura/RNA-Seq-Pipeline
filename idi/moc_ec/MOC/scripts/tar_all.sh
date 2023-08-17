#!/bin/sh


DIR=$1
TAG=$2

files=`ls $DIR | grep $TAG | grep -v tar`

cd $DIR

for file in $files
do
	echo "tar czf $m.tar.gz $file"
	qsub -b Y tar czf $DIR$file".tar.gz" $DIR$file

done