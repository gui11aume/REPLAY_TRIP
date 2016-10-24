# Read in the insertion data.
a = read.delim("allprom_nochromP.txt", as.is=TRUE)
a = subset(a, prom != "p0")

# Use HMM to segment the genome in up/down domains.
library(HMMt)
b = bridge(a[,c(2,4,4,5)])
fit = BaumWelchBin(b$x, b$series.length)

state = fit@ViterbiPath[b$nonvirtuals] - 1

# Mean expression in domains
tapply(INDEX=state, X=a$nexp, mean)


# Transform 0/1 data in genomic domains.
source("helpers/domainify.R")

dom1 = domainify(data.frame(a[,c(2,4,4)], state))
dom0 = domainify(data.frame(a[,c(2,4,4)], 1-state))

# Number of domains.
nrow(dom0)

# Size distribution of the domains.
mean(dom1$end-dom1$start)
quantile(dom1$end-dom1$start)

# Size distribution of inter-domains.
mean(dom0$end-dom0$start)
quantile(dom0$end-dom0$start)


# Make GRange objects.
library(GenomicRanges)

gdom1 = GRanges(Rle(dom1$block),
   IRanges(start=dom1[,2], end=dom1[,3]))
gdom0 = GRanges(Rle(dom0$block),
   IRanges(start=dom0[,2], end=dom0[,3]))
save(gdom1, file="gdom1.rda")
save(gdom0, file="gdom0.rda")

# Read genes for overlap.
genes = read.table("dmel_r5.57_FB2014_03/dmel-all-genes-r5.57.gff",
   sep="\t")
ggenes = GRanges(Rle(genes$V1), IRanges(start=genes$V4, end=genes$V5),
   strand=Rle(genes$V7))
gproms = promoters(ggenes, upstream=1, downstream=0)

# Get overlap between promoters and up domains.
ov = as.matrix(findOverlaps(gdom1, gproms))
mean(tapply(INDEX=ov[,1], X=ov[,1], length))
ov = as.matrix(findOverlaps(gdom0, gproms))
mean(tapply(INDEX=ov[,1], X=ov[,1], length))
