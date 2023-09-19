# Usage: source generate_read_insertion_plots.sh <fastq.gz file>
# The script takes a fastq.gz (Read 2), identifies reads that contain adapter sequenes, and plots the insertion location of adapter sequences
# Adapters sequences should be provided in the adapters.fasta file
# If new adapter sequences are added, then a blast database needs to be created using `makeblastdb -type nucl adapters.fasta`

if [ $# -ne 1 ]; then
    echo "Error: No arguments provided"
    echo "Usage: generate_read_insertion_plots.sh <fastq.gz file>"
else
    use BLAST+

    fastq_gz_path="$1"
    fastq_path="${fastq_gz_path%.gz}"
    fasta_path="${fastq_path%.fastq}.fasta"

    # Decompressing
    echo "Decompressing: $fastq_gz_path to $fastq_path"
    gzip -d -c $fastq_gz_path > $fastq_path

    # converting fastq to fasta (obtained from chatgpt)
    echo "Converting fastq sequences to fasta: $fastq_path to $fasta_path"
    awk 'NR%4==1{printf ">%s\n", substr($0,2)} NR%4==2{print}' $fastq_path > $fasta_path

    # Maps each query (read) to the two adapters, the seed size is controlled by word size. Word size of 10 means there will be atleast 10 NT with perfect match between adapter and read
    # Output - csv format
    # csv columns - 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    blast_op_path="${fasta_path%.fasta}_blast_op_ws_10.csv"
    cmd="blastn -db adapters.fasta -query $fasta_path -out $blast_op_path -outfmt 10 -ungapped  -word_size 10 -strand plus"
    echo 'Using blastn to map reads to adapter sequences'
    echo "Command: $cmd"
    $cmd

    # Generating plots in Python
    cmd="python process_blast_output.py -csv $blast_op_path -fasta $fasta_path"
    echo "Generating plots in python"
    echo "Command: $cmd"
    $cmd

    # Remove intermediate files
    echo "Removing interemediate files"
    rm $fastq_path
    rm $fasta_path
fi