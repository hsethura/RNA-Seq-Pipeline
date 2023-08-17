#!/usr/bin/env Rscript

library(docopt)
library(edgeR)

'edgeR pipeline for comparing between two groups.

Usage: edgeR_main.R [ -i infile -n <noncontrol> -c <control> -t <type> -o <outdir> --tag <tag>]

options:
  -i <infile> --infile <infile>
  -c <control> --control <control>
  -n <non-control> --nc <non-controlol>
  -t <type> --type <type>
  -o <outdir> --out <outdir>
  --tag <tag>' -> doc
  

opts <- docopt(doc)

infile_str <- opts$infile
ctrl_str <- opts$control
nctrl_str <- opts$nc
type_str <- tolower(opts$type)
outdir <- opts$out
tag_str <- opts$tag

print(paste0("infile: ", infile_str))
print(paste0("control: ", ctrl_str))
print(paste0("non_control: ", nctrl_str))
print(paste0("type: ", type_str))
print(paste0("out: ", outdir))
print(paste0("tag: ", tag_str))

# Load the data from the input file

get_splitted_parts <- function(lstr, reg_str) {
    part_lst <- strsplit(lstr, reg_str)
    parts <- part_lst[[1]]
    return (parts)
}

get_status_vec <- function(obj) {
    lstatus <- as.vector(decideTestsDGE(obj))
    lstatus[lstatus == -1] <- "Down"
    lstatus[lstatus == 1] <- "Up"
    lstatus[lstatus == 0] <- "non-DE"
    return(lstatus)
}

ltab <- read.table(infile_str, sep = "\t", header = TRUE, row.names = 1, 
    stringsAsFactors = FALSE, check.names=FALSE)

reg_str <- "\\s*[,|;]\\s*"
ctrl_parts <- get_splitted_parts(ctrl_str, reg_str)
nctrl_parts <- get_splitted_parts(nctrl_str, reg_str)

ctrl_data <- ltab[, ctrl_parts]
nctrl_data <- ltab[, nctrl_parts]

counts_in_use <- data.frame(ctrl_data, nctrl_data, check.names = FALSE)
ctrl_len <- length(ctrl_parts)
nctrl_len <- length(nctrl_parts)
group_arr <- c(rep("C", ctrl_len), rep("NC", nctrl_len))
group <- factor(group_arr)

gene_count <- length(rownames(counts_in_use))

y0 <- DGEList(counts = counts_in_use, group = group)
y1 <- calcNormFactors(y0)
design <- model.matrix(~0 + group)
y2 <- estimateDisp(y1, design, robust=TRUE)
NCvsC <- makeContrasts(groupNC - groupC, levels=design)

dir.create(outdir)
edgeR_tag <- "EdgeR_"
outfile_comp <- paste0(outdir, "/", edgeR_tag, tag_str, "_", type_str, "_comp.tsv")
pdf_file <- paste0(outdir, "/", edgeR_tag, tag_str, "_", type_str, "_comp.pdf")

if (type_str == "qlf") {
    fit <- glmQLFit(y2, design, robust=TRUE)
    qlf <- glmQLFTest(fit, contrast = NCvsC)

    lstatus <- get_status_vec(qlf)
    # Generate MA plot
    pdf(pdf_file)
    plotMD(qlf, status = lstatus, main = tag_str, values=c("Down","Up"), col=c("blue","red"))
    abline(h=c(-1, 1), col="blue")
    dev.off()

    # Generate final table
    final_tab <- topTags(qlf, n = gene_count)
    write.table(final_tab, outfile_comp, sep = "\t")

} else if (type_str == "lr") {

    fit <- glmFit(y2, design)
    lrt <- glmLRT(fit, contrast = NCvsC)

    lstatus <- get_status_vec(lrt)
    pdf(pdf_file)
    plotMD(lrt, status = lstatus, main = tag_str, values=c("Down", "Up"), col=c("blue","red"))
    abline(h=c(-1, 1), col="blue")
    dev.off()

    final_tab <- topTags(lrt, n = gene_count)
    write.table(final_tab, outfile_comp, sep = "\t")
} else {
    print(paste0("Unknown test type: ", type_str))
}



