# Read in the insertion data.
a = read.delim("../allprom_nochromP.txt", as.is=TRUE)
a = subset(a, prom != "p0")

# Use HMM to segment the genome in up/down domains.
library(HMMt)
b = bridge(a[,c(2,4,4,5)])
fit = BaumWelchT(b$x, b$series.length)

state = fit@ViterbiPath[b$nonvirtuals] - 1


# Transform 0/1 data in genomic domains.
source("../helpers/domainify.R")

dom1 = domainify(data.frame(a[,c(2,4,4)], state))
dom0 = domainify(data.frame(a[,c(2,4,4)], 1-state))
# Make GRange objects.
library(GenomicRanges)

gdom1 = GRanges(Rle(dom1$block),
   IRanges(start=dom1[,2], end=dom1[,3]))
gdom0 = GRanges(Rle(dom0$block),
   IRanges(start=dom0[,2], end=dom0[,3]))

save(gdom1, file="gdom1.rda")
save(gdom0, file="gdom0.rda")


# Start writing to file.
sink(file("info_domains.txt"))

cat("Mean expression in domains:\n")
tapply(INDEX=state, X=a$nexp, mean)

cat("Number of (down) domains:\n")
nrow(dom0)

cat("Mean size of up domains:\n")
mean(dom1$end-dom1$start)

cat("Size distribution of up domains:\n")
quantile(dom1$end-dom1$start)

cat("Mean size of down domains:\n")
mean(dom0$end-dom0$start)

cat("Size distribution of down domains:\n")
quantile(dom0$end-dom0$start)

# End writing to file.
sink(NULL)
