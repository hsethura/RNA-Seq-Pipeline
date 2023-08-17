#!/bin/sh

DIR=$1

ALL_FILES=`ls -lrt $DIR | grep fasta | awk '{print $9}'`


for FILE in $ALL_FILES
do
	NAME=`basename $FILE | cut -d"_" -f1,2`
	#echo $NAME
	cat $FILE | awk -v NAME=$NAME '{
										if($1 ~ />/)
											print $1"|"NAME
										else
											print $1
										
									}'
done