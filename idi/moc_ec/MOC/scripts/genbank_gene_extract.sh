#!/bin/sh

GENE_FILE=$1
GENE=$2

TOP_DIR="/idi/moc_ec/MOC/genbank/"
ALL_DIRS=`ls $TOP_DIR`

echo $ALL_DIRS


echo "" | sed 1d > $GENE_FILE

for DIR in $ALL_DIRS
do
	echo $DIR
	NUM_FILES=`ls -lrt $TOP_DIR$DIR"/"*"/"*"rna"* | wc -l`
	if [ $NUM_FILES -gt 0 ];then
		FILES=`ls $TOP_DIR$DIR"/"*"/"*"rna"*`
		zcat $FILES  | awk -v flag="N" -v SPEC=$DIR -v GENE=$GENE '{
											if($0 ~ GENE && $1 ~ />/)
											{	
												if(flag=="Y")
													print ""
												flag="Y"
												#$1=$1"|"SPEC
												print $0
											}
											if($0 !~ GENE && $1 ~ />/)
											{	
												if(flag=="Y")
													print ""
												flag="N"
											}
											if($1 !~ />/ && flag=="Y")
												printf "%s", $0
													
													
													
										}' | sed 's/lcl/lcl|'$DIR'/g' >> $GENE_FILE
	fi
done

ls -lrt $GENE_FILE