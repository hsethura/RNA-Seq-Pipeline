#!/usr/bin/env Rscript


# Shane Roesemann's Amazing R-script.  And handsome to boot!
# "All the scripts that's fit to print()"

# arg 1 = metrics file
# arg 2 = corr file

install.packages("RColorBrewer", repos = "http://cran.us.r-project.org")
install.packages("ggplot2", repos = "http://cran.us.r-project.org")
install.packages("reshape2", repos = "http://cran.us.r-project.org")

library(RColorBrewer)
library(ggplot2)
library(reshape2)

##########################################
# Prepare R for taking in cmd arguments   
##########################################

args <- commandArgs(trailingOnly = T)
METRICS <- as.data.frame(read.table(file = args[1], header=T, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE))
CORR <- as.data.frame(read.table(file = args[2], header=T, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE))


##########################################
# This program's output will be a .pdf
##########################################

pdf(file = "shanetest.pdf", onefile = T, paper = "us")

##########################################
#Barplot of % of Reads Counted
##########################################
par(mfrow = c(1,1))

bp_reads <- barplot(METRICS$pcnt_counted, 
                    main = "Percent of Total reads Counted by Feature Counts", 
                    ylab = "Percent", 
                    ylim = c(0,100),
                    space = 1
                    )
end_point = 0.5 + nrow(METRICS) + nrow(METRICS)-1
#labels
text(seq(1.5,end_point,by=2), par("usr")[3]-1.5,  
     srt = 60, adj= 1, xpd = TRUE,
     labels = as.vector(METRICS$Sample), 
     cex=0.65)

##########################################
#Barplot of % Sense
##########################################

bp_sense <- barplot(METRICS$pcnt_sense, 
                    main = "Percent Aligned to Genes on the Sense Strand", 
                    ylab = "Percent", 
                    ylim = c(0,100),
                    space = 1
)
#labels
text(seq(1.5,end_point,by=2), par("usr")[3]-1.5, 
     srt = 60, adj= 1, xpd = TRUE,
     labels = paste(METRICS$Sample), cex=0.65)


##########################################
#Stacked plot of total counts + CDS counts
##########################################


#First, create a mini df with the appropriate metrics (NON_CDS, CDS)
Non_CDS <- as.vector(METRICS$total_counted - METRICS$CDS_total_counts)
CDS <- as.vector(METRICS$CDS_total_counts)
stacked_df <- as.data.frame(CDS)
stacked_df$Non_CDS <- Non_CDS

options(scipen = 5)
#Then, make bargraph
bp_total_counts <- barplot(as.matrix(t(stacked_df[,1:2])), 
                           main = "Total Reads Counted per Sample",
                           space = 1, col = c("grey70", "grey50"), las = 1)

# this adds the tilted sample labels on x-axis
text(seq(1.5,end_point,by=2), par("usr")[3]-1.5, 
     srt = 60, adj= 1, xpd = TRUE,
     labels = paste(METRICS$Sample), cex=0.65)

abline(1000000, 0, col = "red")
legend(x = "topleft", legend = colnames(stacked_df), 
       col = c("grey70", "grey50"), cex = 1, pch = 15, text.width = 3)



##########################################
#Stacked bar grpah of read type per sample
##########################################

#make a vectors of data you NEED
CDS <- as.vector(METRICS$CDS_pcnt_of_counted)
rRNA_pct_counted <- as.vector(METRICS$rRNA_pcnt_of_counted)
tRNA_pct_counted <- as.vector(METRICS$tRNA_pcnt_of_counted)
IGR <- as.vector(METRICS$IGR_pcnt_of_counted)
misc_rRNA <- as.vector(METRICS$misc_RNA_pcnt_of_counted)
AS_pct <- as.vector(100 - METRICS$pcnt_sense)

#put the vectors into a df
stackedbp2 <- as.data.frame(CDS)
stackedbp2$rRNA <- rRNA_pct_counted
stackedbp2$tRNA <- tRNA_pct_counted
stackedbp2$IGR <- IGR
stackedbp2$misc_rRNA <- misc_rRNA #this was added by shane - Jonathan did not ask for this
stackedbp2$AS_pct <- AS_pct #this was added by shane 


#plot
color1 <- rev(brewer.pal(n = 6,
                         name = "Blues"))

bp_stacked <- barplot(as.matrix(t(stackedbp2[,1:6])), ylim = c(0,100), 
                      col = color1, main = "Composition of Reads per Sample", space = 1)

legend(x = "bottomleft", legend = colnames(stackedbp2), 
       cex = 1, pch = 15, text.width = 3, col= color1)

text(seq(1.5,end_point,by=2), par("usr")[3]-1.5, 
     srt = 60, adj= 1, xpd = TRUE,
     labels = paste(METRICS$Sample), cex=0.65)

##########################################
#Correlation Heatmap
##########################################

#The dataset must be coerced into a numeric matrix before heat map construction.
corr<- apply(as.matrix(CORR[2:10,2:10]), 1, as.numeric)

#Heatmap construction
h_map <- heatmap(corr, symm = T, 
                labRow = METRICS$Sample[3:12],
                labCol = METRICS$Sample[3:12],
                main = "Correlation Heatmap")


##########################################
#End .pdf
##########################################

dev.off()







