#!/bin/bash

read1_ip=$1
read2_ip=$2
read1_op=$3
read2_op=$4

# Create uncompressed file names
read1_uncomp_op="${read1_op/1.fastq.gz/1.fastq}"
read2_uncomp_op="${read2_op/2.fastq.gz/2.fastq}"

# Print statements to inform the user of progress
echo "Processing read1 file: $read1_ip"
echo "Processing read2 file: $read2_ip"

echo "Uncompressed output will be written to: $read1_uncomp_op and $read2_uncomp_op"
echo "Compressed output will be written to: $read1_op and $read2_op"

# Run the command using the file variables
paste <(zcat "$read1_ip" | paste - - - -) <(zcat "$read2_ip" | paste - - - -) | \
awk -F'\t' -v out1="$read1_uncomp_op" -v out2="$read2_uncomp_op" '{
    umi=substr($2, 1, 8); 
    gsub(/ /, "_" umi " ", $1);  # Replace space with "umi " in read1 ID
    gsub(/ /, "_" umi " ", $5);  # Replace space with "umi " in read2 ID
    print $1"\n"$2"\n"$3"\n"$4 >> out1; 
    print $5"\n"$6"\n"$7"\n"$8 >> out2
}'

# Print statements after writing the files
echo "Finished writing updated FASTQ files."

# Compress the output files
gzip -f "$read1_uncomp_op"
gzip -f "$read2_uncomp_op"

# Print statements after compression
echo "Compressed output files: $read1_op and $read2_op"

echo "Operation completed successfully."