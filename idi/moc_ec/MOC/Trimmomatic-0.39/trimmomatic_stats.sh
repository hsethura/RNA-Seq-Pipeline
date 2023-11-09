#!/bin/sh

# This code snippet parses trimmomatic log file and generates the distribution of read lengths after trimming by trimmomatic. 
# It also provides some stats such as %reads trimmed by Trimmomatic etc.

TRIMMOMATIC_LOG_FILE=$1 # trimmomatic log file
TRIMMOMATIC_SUMMARY_FILE=$2 # trimmomatic output file 

TMP_FILE="${TRIMMOMATIC_LOG_FILE/trimlog.txt/trimlog_tmp.txt}"

# OP_FILE="${TRIMMOMATIC_LOG_FILE/trimlog.txt/trimlog_stats.txt}"
OP_FILE="${TRIMMOMATIC_LOG_FILE/trimlog.txt/trimlog_stats.csv}"

# echo "READ_LENGTH COUNT" > $OP_FILE

# # $3 == 0 (check if read length after trimmomatic is zero)
# # $4 != 0 || $6 != 0 (checks if the read was trimmed either at the beginning or at the end)
# # sort | uniq -c | awk '{print $2, $1}' (sorts and prints the distribution of read lenghts)
# awk '($3 == 0 || $4 != 0 || $6 != 0) {print $3}' $TRIMMOMATIC_LOG_FILE | sort | uniq -c | awk '{print $2, $1}' > $TMP_FILE
# cat $TMP_FILE >> $OP_FILE

# grep "Input Read Pairs" $TRIMMOMATIC_SUMMARY_FILE >> $OP_FILE

# total_reads="$(wc -l $TRIMMOMATIC_LOG_FILE)"

# # calculates the total reads trimmed by trimmomatic by adding up the counts from the read distribution
# modified_reads="$(awk '{ sum += $2 } END { print sum }' $TMP_FILE)"

# echo "Total reads trimmed by Trimmomatic: $modified_reads" >> $OP_FILE
# awk -v total_reads="$total_reads" -v modified_reads="$modified_reads" 'BEGIN { printf "\% reads trimmed by Trimmomatic: %.2f", (100*modified_reads) / total_reads }' >> $OP_FILE
# rm $TMP_FILE

trimmed_read_pairs=$(awk '($3 == 0 || $4 != 0 || $6 != 0) {print $1}' $TRIMMOMATIC_LOG_FILE | uniq | wc -l)
total_read_pairs=$(awk '{print $1}' $TRIMMOMATIC_LOG_FILE | uniq | wc -l)
percent_trimmed_read_pairs=$(awk -v num="$trimmed_read_pairs" -v den="$total_read_pairs" 'BEGIN { printf "%.2f", (100*num) / den }')

empty_read_pairs=$(awk '($3 == 0) {print $1}' $TRIMMOMATIC_LOG_FILE | sort | uniq -c | awk '$1 == 2' | wc -l)
percent_empty_read_pairs=$(awk -v num="$empty_read_pairs" -v den="$total_read_pairs" 'BEGIN { printf "%.2f", (100*num) / den }')

trimmed_read1s=$(awk '$2 ~ /^1/' $TRIMMOMATIC_LOG_FILE | awk '($3 == 0 || $4 != 0 || $6 != 0) {print $1}' | wc -l)
total_read1s=$(awk '$2 ~ /^1/' $TRIMMOMATIC_LOG_FILE | wc -l)
percent_trimmed_read1s=$(awk -v num="$trimmed_read1s" -v den="$total_read1s" 'BEGIN { printf "%.2f", (100*num) / den }')

empty_read1s=$(awk '$2 ~ /^1/' $TRIMMOMATIC_LOG_FILE  | awk '($3 == 0) {print $1}' | wc -l)
percent_empty_read1s=$(awk -v num="$empty_read1s" -v den="$total_read1s" 'BEGIN { printf "%.2f", (100*num) / den }')

trimmed_read2s=$(awk '$2 ~ /^2/' $TRIMMOMATIC_LOG_FILE  | awk '($3 == 0 || $4 != 0 || $6 != 0) {print $1}' | wc -l)
total_read2s=$(awk '$2 ~ /^2/' $TRIMMOMATIC_LOG_FILE | wc -l)
percent_trimmed_read2s=$(awk -v num="$trimmed_read2s" -v den="$total_read2s" 'BEGIN { printf "%.2f", (100*num) / den }')

empty_read2s=$(awk '$2 ~ /^2/' $TRIMMOMATIC_LOG_FILE  | awk '($3 == 0) {print $1}' | wc -l)
percent_empty_read2s=$(awk -v num="$empty_read2s" -v den="$total_read2s" 'BEGIN { printf "%.2f", (100*num) / den }')

echo "total_read1s,trimmed_read1s,percent_trimmed_read1s,empty_read1s,percent_empty_read1s,total_read2s,trimmed_read2s,percent_trimmed_read2s,empty_read2s,percent_empty_read2s,total_read_pairs,trimmed_read_pairs,percent_trimmed_read_pairs,empty_read_pairs,percent_empty_read_pairs" > $OP_FILE
echo "$total_read1s,$trimmed_read1s,$percent_trimmed_read1s,$empty_read1s,$percent_empty_read1s,$total_read2s,$trimmed_read2s,$percent_trimmed_read2s,$empty_read2s,$percent_empty_read2s,$total_read_pairs,$trimmed_read_pairs,$percent_trimmed_read_pairs,$empty_read_pairs,$percent_empty_read_pairs" >> $OP_FILE