load("../predpow.rda")

predpow = sort(predpow, decreasing=TRUE)
hmark = grepl("^H3K", names(predpow))
chromod = names(predpow) %in% c("C5", "S9")

COL = rep('grey70', length(predpow))
COL[hmark] = 'orange2'
COL[chromod] = 'deeppink2'

pdf("Fig4a.pdf", useDingbats=FALSE)
barplot(predpow, ylab="Predictive power",
   col=COL, border=NA, axisnames=FALSE)
dev.off()

# Start writing to file.
sink(file("info_Fig4a.txt"))

cat("Sorted list of features.\n")
cat(paste(names(predpow), round(predpow, 3), sep="\t", collapse="\n"))

sink(NULL)
