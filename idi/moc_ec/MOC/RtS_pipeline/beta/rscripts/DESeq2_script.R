#!/usr/bin/env Rscript

suppressMessages(library(argparse))
suppressMessages(library(DESeq2))
suppressMessages(library(data.table))

pos_to_sid <- function(sample_id_vec, pos_lst) {
    regex_str <- "[,|;]\\s*"
    pos_int_vec <- as.integer(unlist(strsplit(pos_lst, regex_str)))
    sid_vec <- sample_id_vec[pos_int_vec]
    return (sid_vec)

}

get_sample_id_vec <- function(key_path) {

    sample_id <- 'Sample_ID'
    ltab <- read.csv(key_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    sample_id_vec1 <- ltab[, sample_id]
    sample_id_vec2 <- sample_id_vec1[sample_id_vec1 != ""]
    # replace the empty spaces with underscore
    sample_id_vec <- gsub("\\s+", "_", sample_id_vec2)
    return(sample_id_vec)
}

format_sid_str <- function(sid_lst) {
    regex_str <- "[,|;]\\s*"
    sid_vec1 <- unlist(strsplit(sid_lst, regex_str))
    sid_vec = gsub("\\s+", "_", sid_vec1)
    return (sid_vec)
}


# Configure the parser
parser <- ArgumentParser(description='Process input for DESeq2')
parser$add_argument("-k", "--key_path", dest="key_path", type="character", 
    required = TRUE, help="Path to the key file")
parser$add_argument("-r", "--readcount_path", dest="count_path", 
    type = "character", required = TRUE, help = "Path to the readcount file")
parser$add_argument("-p", "--prefix", dest = "prefix", 
    type = "character", required = TRUE, help = "Prefix of names used in the output files")
parser$add_argument("-C", "--group-control", dest = "c_lst",
    type ="character", required = TRUE, help = "The list of names/positions for the control group")
parser$add_argument("-N", "--group-noncontrol", dest = "nc_lst",
    type ="character", required = TRUE, help = "The list of names/positions for the non-control group")
parser$add_argument("-t", "--input-type", dest = "input_type", 
    type = "character", default = "sample_id", help = "Optional/ How to specify the samples (sid/pos)" )
parser$add_argument("-o", "--outdir", dest = "outdir", 
    type = "character", default = ".", help = "The output directory for results")

args <- parser$parse_args()

# Collect the arguments

key_path = args$key_path
count_path = args$count_path
c_lst <- args$c_lst
nc_lst <- args$nc_lst
input_type <- args$input_type
outdir <- args$outdir
prefix <- args$prefix

# if the input type is in terms of position, open the key file and get the sample_id s 

c_sid_vec = NULL
nc_sid_vec = NULL

if (input_type == 'pos') {
    sample_id_vec <- get_sample_id_vec(key_path)
    c_sid_vec <- pos_to_sid(sample_id_vec, c_lst)
    nc_sid_vec <- pos_to_sid(sample_id_vec, nc_lst)
} else {
    c_sid_vec <- format_sid_str(c_lst)
    nc_sid_vec <- format_sid_str(nc_lst)
}

c_sid_str = paste(c_sid_vec, collapse = ", ")
nc_sid_str = paste(nc_sid_vec, collapse = ", ")

# Print the arguments for sanity check
print(paste0("c_lst: ", c_lst))
print(paste0("nc_lst: ",  nc_lst))
print(paste0("input_type: ", input_type))
print(paste0("Control samples: ", c_sid_str))
print(paste0("Non-control samples: ", nc_sid_str))

# Now that we know the sample ids, we can get the read counts for those samples 
# and can take DESeq2 from there.

count_tab <- read.csv(count_path, sep = "\t", header = TRUE, row.names = 1, 
    stringsAsFactors = FALSE, check.names = FALSE)

c_reads = count_tab[, c_sid_vec]
nc_reads = count_tab[, nc_sid_vec]

countData <- cbind(c_reads, nc_reads)
condVec <- c(rep("control", length(c_sid_vec)), rep("noncontrol", length(nc_sid_vec)))
colData <- data.frame(condVec)
colnames(colData) <- c("condition")
rownames(colData) <- colnames(countData)

colData

dds <- DESeqDataSetFromMatrix(countData = countData,
    colData = colData,
    design = ~ condition)
#dds

#dds <- dds[ rowSums(counts(dds)) > 1, ]
dds$condition <- relevel(dds$condition, ref="control")
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
resOrdered2 = setDT(data.frame(resOrdered), keep.rownames = TRUE)[]
colnames(resOrdered2)[1] <- "Gene_id"

dir.create(file.path(outdir))
outfile = paste0(outdir, "/", prefix, ".tsv")
write.table(resOrdered2, outfile, sep = "\t", row.names = FALSE)


ma_pdf_file = paste0(outdir, "/", prefix, "_MA.pdf")
pdf(ma_pdf_file)
plotMA(res, main="DESeq2", ylim=c(-2,2))
dev.off()
