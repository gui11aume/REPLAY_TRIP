require(GenomicRanges)

# Read in preprocessed HiC data.
diags = read.table("20diags.txt.gz")
# Read in paradoxical domains.
prdx = read.delim("prdx_domains.txt")
# Read in color domains.
coldoms = read.delim("color_domains.txt")
RED = subset(coldoms, col=="RED")

prop5 = rowSums(diags[,5:6]) / rowSums(diags[,4:ncol(diags)])
save(prop5, file="/data/prop5.rda")

gdiags = GRanges(Rle(diags$V1),
   IRanges(start=diags$V2, end=diags$V3))
gprdx = GRanges(Rle(prdx$block),
   IRanges(start=prdx$start, end=prdx$end))
gRED = GRanges(Rle(RED$chr),
   IRanges(start=RED$start, end=RED$end))

ov1 = which(countOverlaps(gdiags, gprdx) > 0)
ov2 = which(countOverlaps(gdiags, gRED) > 0)

t.test(prop5[ov1], prop5[ov2[!ov2 %in% ov1]])
wilcox.test(prop5[ov1], prop5[ov2[!ov2 %in% ov1]])

# Do the linear regression "manually".
X = cbind(rep(1, 21), log(seq(1000, 41000, by=2000)))
Y = as.matrix(diags[,4:(ncol(diags)-1)] /
   rowSums(diags[,4:ncol(diags)]))
dual = (solve(t(X) %*% X) %*% t(X))[2,]
# Remove NAs so that computation can proceed.
decay = as.vector(Y %*% dual)
save(decay, file="/data/decay.rda")

t.test(decay[ov1], decay[ov2[!ov2 %in% ov1]])
wilcox.test(decay[ov1], decay[ov2[!ov2 %in% ov1]])

load("gdom0.rda")
load("gdom1.rda")

ova = which(countOverlaps(gdiags, gdom0) > 0)
ovb = which(countOverlaps(gdiags, gdom1) > 0)

t.test(decay[ova], decay[ovb])
wilcox.test(decay[ova], decay[ovb])

source("helpers/violin-plot.R")

pdf("violin_short_range.pdf", width=4)
violin.plot(list(prop5[ov1], prop5[ov2[!ov2 %in% ov1]],
   prop5[ova], prop5[ovb]), x.pos=c(1,2,3.5,4.5), hbars=FALSE)
violin.plot(list(decay[ov1], decay[ov2[!ov2 %in% ov1]],
   decay[ova], decay[ovb]), x.pos=c(1,2,3.5,4.5), hbars=FALSE)
dev.off()

pdf("anticorrelation.pdf")
plot(prop5, decay, pch='.', xlim=c(0,.3), ylim=c(-.05,0),
  col="#00000080", xlab="Proportion of contacts within 5 kb",
  ylab="Rate of decay")
dev.off()

prop = diags[,4:(ncol(diags)-1)] / rowSums(diags[,4:ncol(diags)])
y1 = colMeans(prop[ova,], na.rm=TRUE)
y2 = colMeans(prop[ovb,], na.rm=TRUE)

pdf("decay.pdf", useDingbats=FALSE)
plot(log10(seq(1000, 41000, by=2000)), log10(y1),
     xlab="Genomic distance [log10(bp)]",
     ylab="Average contact proportion [log10]", type="n")
rect(xleft=3, xright=4.6, ybottom=-2.15, ytop=-1, 
     col="grey95", border=NA)
# Horizontal lines.
segments(x0=rep(3.0,4), x1=rep(4.6,4),
   y0=seq(-2.2,-1, by=.4), y1=seq(-2.2,-1, by=.4),
   col="white", lwd=2)
# Vertical lines
segments(x0=c(3,3.5,4,4.5), x1=c(3,3.5,4,4.5),
   y0=rep(-2.2,4), y1=rep(-1,4),
   col="white", lwd=2)
points(log10(seq(1000, 41000, by=2000)), log10(y1))
points(log10(seq(1000, 41000, by=2000)), log10(y2), pch=19)
dev.off()

#a = subset(read.delim("allprom.txt"), prom != "p0")
#ga = GRanges(Rle(a$chr), IRanges(start=a$pos, width=1))

#GSE = read.delim("GSE36175_norm_aggregated_tiling_arrays.txt.gz")
## Remove "chr" for compatibility.
#GSE$chromosome = sub("^chr", "", GSE$chromosome)
#
#gGSE = GRanges(Rle(GSE$chromosome),
#   IRanges(start=GSE$start, end=GSE$end))
#
#ov = as.matrix(findOverlaps(ga, gdiags))
#x = prop5[ov[,2]]
#y = a$nexp[ov[,1]]
#lm(y~x)
#anova(lm(y~x))
#
#ov = as.matrix(findOverlaps(ga, gGSE))
#x = GSE$LAM[ov[,2]]
#y = a$nexp[ov[,1]]
#lm(y~x)
#anova(lm(y~x))
