#!/bin/bash
#$ -l h_vmem=4G -pe smp 4 -binding linear:4
#$ -l h_rt=240:00:00

#Whitelist the Bead Barcode (BB) file to check for diversity
#It will output counts vs barcode plots.
#In the future maybe we can add a plot to show the distribution of the number of reads per BB, which can help QC
#The lavender sequence from the lucid chart is: GACGCTGCCGACGA
#Having this sequence in front of the BB read means the sequencing was done with illumina's Reverse Complement Workflow.
#This script will subsample $num_reads from the index read ($bb_file)
#Then it will produce the whitelist and we can check the diversity of outputted barcodes
#Will need to have umi_tools installed.

#User inputs
# outdir="/idi/cgtb/jbagnall/bacdrop/SCR-0001.2_Klebs_novaseq_230821/test_files/" #This is the desired output directory
# read1_file="/idi/cgtb/jbagnall/bacdrop/SCR-0001.2_Klebs_novaseq_230821/CGAGGCTG.unmapped.1.fastq.gz" #This is the read with the UMI (8bp) at the beginning
# bb_file="/idi/cgtb/jbagnall/bacdrop/SCR-0001.2_Klebs_novaseq_230821/CGAGGCTG.unmapped.index_1.fastq.gz" #This is the index read with the supposed bead barcode (BB, 16bp)

outdir=$1
read1_file=$2
bb_file=$3

num_reads=100000 #Number of reads to subset to

# mkdir -p "$outdir"

#Automatically generated names for whitelist output files
bb_name="$(basename -- $bb_file .fastq.gz)" #strips the directory and extension
wl_plot_prefix="${outdir}${bb_name}.wl_plot"
wl_log="${outdir}${bb_name}.wl.log"
wl_out="${outdir}${bb_name}.wl.txt"

echo "bb_name: $bb_name"
echo "wl_plot_prefix: $wl_plot_prefix"
echo "wl_log: $wl_log"
echo "wl_out: $wl_out"


umi_tools whitelist -I "$read1_file" \
                    --read2-in "$bb_file" \
                    --extract-method=string \
                    --bc-pattern=NNNNNNNN \
                    --bc-pattern2=CCCCCCCCCCCCCCCC \
                    --plot-prefix="$wl_plot_prefix" \
                    --method=umis \
                    --knee-method=distance \
                    --error-correct-threshold=2 \
                    --subset-reads="$num_reads" \
                    -L "$wl_log" > "$wl_out"
