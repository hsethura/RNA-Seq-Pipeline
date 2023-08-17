
args <- commandArgs(trailingOnly = TRUE)

gDat <- read.delim(args[1])

val <- gDat[,1]
log.val <- log(val)

pdf(args[2]) 
 
hist(log.val, col="lightblue", xlab="Log10 FPKMO", main=args[3], breaks=20)