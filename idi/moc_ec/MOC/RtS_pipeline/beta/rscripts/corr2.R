#!/usr/bin/env Rscript

# arg1 = FPKMO file
# arg2 = CORR file
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
FPKMO= read.table(file = args[1], header=T, row.names = 1, sep = "\t", quote = "", comment.char = "")
sorted_FPKMO <- FPKMO[ , order(names(FPKMO))]

log_FPKMO <- log10(sorted_FPKMO+1)

corr <- cor(log_FPKMO)

write.table(corr, file = args[2],  sep = "\t", quote = F, row.names = T) 

#my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
#pdf("HM.pdf", width=4, height=4)
#heatmap(as.matrix(corr, col=my_palette, scale="row", key=T, keysize=1.5,density.info="none", trace="none",cexCol=0.9, labRow=NA, Colv=NA))
#dev.off()" 



