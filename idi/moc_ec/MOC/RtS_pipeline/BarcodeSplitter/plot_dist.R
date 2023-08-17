args = commandArgs(TRUE)
infile <- args[1]
outfile <- args[2]
name_str <- args[3]

library(reshape2)
library(ggplot2)
res <- readLines(infile)
lastBar <- ""
lastZero <- 0.0
lastOne <- 0.0
lastPercent <- 0.0

totalP <- 0

# Build a dataframe where to store all the data

freqMat <- data.frame(nrow = 96, ncol = 3)
barVec <- c()
count <- 1
com_reg_str <- "([0-9]+\\.?[0-9]*e?\\-?[0-9]*)"
regex_zero <- paste0("^.*\\s+", com_reg_str, "%$")
regex_one <- regex_zero
regex_percent <- paste0("^.*\\s+", com_reg_str, "%\\)$")

for (lstr in res) {
	if (grepl("^Barcode", lstr)) {
		lastBar <- strsplit(lstr, " ")[[1]][2]
		
	} else if (grepl("^Zero base mismatch", lstr)) {
		lastZero <- as.numeric(gsub(regex_zero, "\\1", lstr))
	} else if (grepl("^One base mismatch", lstr)) {
		lastOne <- as.numeric(gsub(regex_one, "\\1", lstr))
	} else if (grepl("^Total read for this barcode", lstr)) {
		lastPercent <- as.numeric(gsub(regex_percent, "\\1", lstr))
		#cat("Debug: ", lastBar, " ", lastZero, " ", lastOne, " ", lastPercent, "\n")
		lastZeroT <- (lastZero * lastPercent)/100.0
		lastOneT <- (lastOne * lastPercent)/100.0
		if (is.na(lastZeroT) || is.na(lastOneT))
			next
		#cat (count, " ", lastBar, " ", lastZeroT, " ", lastOneT, "\n")
		totalP <- totalP + lastZeroT + lastOneT
		barVec <- c(barVec, lastBar)
		freqMat[count, 1] <- lastBar
		freqMat[count, 2] <- lastZeroT
		freqMat[count, 3] <- lastOneT
		count <- count + 1
	}
}

freqMat1 <- freqMat[complete.cases(freqMat),]

#rownames(freqMat) <- barVec
#plot_str <- strsplit(arg1, "\\.")[[1]][1]
plot_str <- name_str
colnames(freqMat1) <- c("barCode", "ZeroMM", "OneMM")
freqMat2 <- melt(freqMat1, id = "barCode")
cnames <- c("BarCode", "Mismatch", "PercentT")
colnames(freqMat2) <- cnames
# Now start working with ggplot2
lplot <- ggplot(data=freqMat2, aes(x=reorder(BarCode, -PercentT), y=PercentT, 
	fill=Mismatch)) + xlab("Barcode") +
	ylab("Percentage of total reads") +
	ggtitle(plot_str) + geom_bar(stat="identity") + 
	theme(axis.text.x = element_text(size = 4, angle = 90, hjust = 1), 
	legend.position = c(0.9, 0.9))

#plotname <- paste0(plot_str, ".pdf")
ggsave(outfile, lplot)
#cat(totalP, "\n")
