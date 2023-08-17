#!/bin/sh


IN_DIR=$1
TRIM_DIR=$2
TRIM_1=$3
TRIM_2=$4

 

trim_fastq ()
{

	INPUT_FILE=$1
	TRIM_FILE=$2
	TRIM_LEN=$3
	
	zcat $INPUT_FILE | awk '{
								# qual
								if(NR%4==0)
								{
									y=substr($0,'$TRIM_LEN',10000)
									print y
								}
								
								# name
								if(NR%4==1)
									print $0
								
								# seq
								if(NR%4==2)
								{
									y=substr($0,'$TRIM_LEN',10000)
									print y
								}
								# +
								if(NR%4==3)
									print $0
							
							}' > $TRIM_FILE
							
							"Gzipping $TRIM_FILE..."

							gzip $TRIM_FILE
							

}

ALL_R1=`ls -lrt $IN_DIR/*R1*fastq* | grep -v "ndeterm" | awk '{print $9}'`
ALL_R2=`ls -lrt $IN_DIR/*R2*fastq* | grep -v "ndeterm" | awk '{print $9}'`

echo $ALL_R1
echo $ALL_R2

mkdir -p $TRIM_DIR

for FILE in $ALL_R1
do
	TRIM=$TRIM_1
	IN_FILE_NAME=`basename $FILE`
	TRIM_FILE_NAME=`basename $FILE .gz`
	
	INFILE=$IN_DIR"/"$IN_FILE_NAME
	TRIM_FILE=$TRIM_DIR"/"$TRIM_FILE_NAME
	
	echo trim_fastq $INFILE $TRIM_FILE $TRIM
	trim_fastq $INFILE $TRIM_FILE $TRIM
done

for FILE in $ALL_R2
do
	TRIM=$TRIM_2
	IN_FILE_NAME=`basename $FILE`
	TRIM_FILE_NAME=`basename $FILE .gz`
	
	INFILE=$IN_DIR"/"$IN_FILE_NAME
	TRIM_FILE=$TRIM_DIR"/"$TRIM_FILE_NAME
	
	echo trim_fastq $INFILE $TRIM_FILE $TRIM
	trim_fastq $INFILE $TRIM_FILE $TRIM
done
