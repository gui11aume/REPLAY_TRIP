require(GenomicRanges)

# Read in preprocessed HiC data.
diags = read.table("../20diags.txt.gz")
gdiags = GRanges(Rle(diags$V1),
   IRanges(start=diags$V2, end=diags$V3))

# Do the linear regression "manually".
X = cbind(rep(1, 21), log(seq(1000, 41000, by=2000)))
Y = as.matrix(diags[,4:(ncol(diags)-1)] /
   rowSums(diags[,4:ncol(diags)]))

dual = (solve(t(X) %*% X) %*% t(X))[2,]
decay = as.vector(Y %*% dual)

#save(decay, file="/data/decay.rda")

load("gdom0.rda")
load("gdom1.rda")

ova = which(countOverlaps(gdiags, gdom0) > 0)
ovb = which(countOverlaps(gdiags, gdom1) > 0)

prop = diags[,4:(ncol(diags)-1)] / rowSums(diags[,4:ncol(diags)])
y1 = colMeans(prop[ova,], na.rm=TRUE)
y2 = colMeans(prop[ovb,], na.rm=TRUE)

pdf("Fig3e.pdf", useDingbats=FALSE)
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
