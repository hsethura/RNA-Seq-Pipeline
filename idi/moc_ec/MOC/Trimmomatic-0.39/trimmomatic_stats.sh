#!/bin/sh

TRIMLOG_FILE=$1
TRIMMOMATIC_LOG_FILE=$2

TMP_FILE="${TRIMLOG_FILE/trimlog.txt/trimlog_tmp.txt}"

OP_FILE="${TRIMLOG_FILE/trimlog.txt/trimlog_stats.txt}"
echo "READ_LENGTH COUNT" > $OP_FILE

awk '($3 == 0 || $4 != 0 || $6 != 0) {print $3}' $TRIMLOG_FILE | sort | uniq -c | awk '{print $2, $1}' > $TMP_FILE
cat $TMP_FILE >> $OP_FILE

grep "Input Read Pairs" $TRIMMOMATIC_LOG_FILE >> $OP_FILE

awk '{ sum += $2 } END { print "Total reads trimmed by Trimmomatic:", sum }' $TMP_FILE >> $OP_FILE

rm $TMP_FILE
