#!/usr/bin/env Rscript

suppressMessages(library(argparse))

parser <- ArgumentParser(description='Generate histograms for length read collapse information')
parser$add_argument("-i", "--infile_str", dest = "infile_str", 
    type = "character", required = TRUE, help = "Path to the input file")
parser$add_argument("-o", "--outdir", dest = "outdir",
    type = "character", required = TRUE, help = "Path to the output directory")
parser$add_argument("-p", "--prefix", dest = "prefix",
    type = "character", required = TRUE, help = "Prefix names used in the output files")

args <- parser$parse_args()

# Collect the arguments
print(args)
infile_str <- args$infile_str
outdir <- args$outdir
prefix <- args$prefix

outfile_str <- paste0(outdir, "/", prefix, "_coll_len_hist.pdf")

# Read the data to the input file

info = file.info(infile_str)
lsize <- info$size
if (lsize == 0) {
    stop("The file is empty")
}

ltab <- read.table(infile_str, stringsAsFactors = FALSE)
ltab2 <- ltab$V1
xlab_str <- "Number of reads for a single collapse"
ylab_str <- "Collapse frequency"
min_str <- paste0("min: ", min(ltab2))
max_str <- paste0("max: ", max(ltab2))
median_str <- paste0("median: ", median(ltab2))
lbreaks = 0:101
ltab3 <- ltab2
ltab3[ltab2 > 100] <- 101
ltab_nol <- ltab2[!(ltab2 > 100)] 
single_count <- sum(ltab2 == 1)
single_count_p1 <- (single_count * 100)/ length(ltab2)
single_count_p <- format(round(single_count_p1, 2), nsmall = 2)
single_count_nol_p1 <- (single_count * 100)/length(ltab_nol)
single_count_nol_p <- format(round(single_count_nol_p1, 2), nsmall =2)
print(single_count_p)
print(single_count_nol_p)
#lbreaks = c(lbreaks1, max(ltab3))
total_reads_before <- sum(ltab2)
total_reads_after <- length(ltab2)
compress1 <- total_reads_before /total_reads_after
compress <- format(round(compress1, 2), nsmall = 2)
single_umi_p_str <- paste0("single_umi: ", single_count_p, "%")
single_umi_nol_p_str <- paste0("single_umi_tr_100: ", single_count_nol_p, "%")
total_str <- paste0("reads_before: ", total_reads_before, ", reads_after: ", total_reads_after, "\navg compression: ", compress, ", ", single_umi_p_str, "\n", single_umi_nol_p_str)
main_str = paste0(prefix, "\n", min_str, ", ", median_str, ", ", max_str, "\n", total_str)
xlim_range <- seq(1,101,20)
print(outfile_str)
pdf(outfile_str)
hist(ltab3, breaks = lbreaks, xlab = xlab_str, ylab = ylab_str, main = main_str, cex.main=0.8, axes = FALSE)
axis(side=1, at=xlim_range)
axis(side = 2)
dev.off()
