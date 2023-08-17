#!/bin/sh

for m in $@
do

	echo "Changing permissions in $m..."
	echo "chmod -R 777 $m"
	chmod -R 777 $m

done