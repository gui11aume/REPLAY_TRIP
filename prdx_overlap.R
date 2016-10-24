require(GenomicRanges)

load("dmel_r5.57_FB2014_03/genes_r5.57.rda")
load("dmel_r5.57_FB2014_03/act_genes_r5.57.rda")
load("dmel_r5.57_FB2014_03/inact_genes_r5.57.rda")
prdx = read.delim("prdx_domains.txt")

gprdx = GRanges(Rle(prdx$block),
             IRanges(start=prdx$start, end=prdx$end))

ov = as.matrix(findOverlaps(gprdx, genes_r5.57))
length(unique(ov[,2]))

IDs = unique(sub("^ID=([^;]*);.*", "\\1", mcols(genes_r5.57[ov[,2]])$attr))
write(IDs, file="prdx_idx.txt", ncol=1)

# Show the size of genes covered by paradoxical domains.
distr = width(genes_r5.57[ov[,2]])
quantile(distr)
quantile(width(act_genes_r5.57))
quantile(width(inact_genes_r5.57))

# Resample paradoxical domains 1000 times.

#if (file.exists("tmp.rda")) {
#   load("tmp.rda")
#} else {
#   set.seed(123)
#   sizes = list()
#   for (iter in 1:1000) {
#      blocklist = list()
#
#      for (blockname in unique(prdx$block)) {
#         block = subset(prdx, block == blockname)
#         nblock = nrow(block)
#         idx = sample(nblock, size=1)
#         if (idx > 1) {
#            shift1 = block$start[idx] - block$start[1]
#            shift2 = block$end[nblock] - block$end[idx-1]
#            block$start[idx:nblock] = block$start[idx:nblock] - shift1
#            block$end[idx:nblock] = block$end[idx:nblock] - shift1
#            block$start[1:(idx-1)] = block$start[1:(idx-1)] + shift2
#            block$end[1:(idx-1)] = block$end[1:(idx-1)] + shift2
#            block = rbind(block[c(idx:nblock,1:(idx-1)),])
#         }
#         blocklist[[blockname]] = block
#      } 
#
#      rsmpl = Reduce(rbind, blocklist)
#      grsmpl = GRanges(Rle(rsmpl$block),
#         IRanges(start=rsmpl$start, end=rsmpl$end))
#
#      ov = as.matrix(findOverlaps(grsmpl, genes_r5.57))
#      sizes[[iter]] = width(genes_r5.57[ov[,2]])
#
#   }
#   save(sizes, file="tmp.rda")
#}
#
#allmed = sapply(sizes, median, na.rm=TRUE)
#mean(median(distr, na.rm=TRUE) > allmed)

# Show the significance after resampling.
#pdf("/data/hist.pdf")
#hist(allmed, breaks=20, xlim=c(2000,12000))
#abline(v=median(distr, na.rm=TRUE), col=3)
#dev.off()

# Make a stripchart to compare.
pdf("stripsize.pdf", width=4, height=8, useDingbats=FALSE)

COLa = c("#EFB642", "#C14730", "#1E2A44")
COLb = c("#EFB64280", "#C1473080", "#1E2A4480")

x1 = log10(distr)
x2 = log10(sample(width(act_genes_r5.57), size=length(distr)))
x3 = log10(sample(width(inact_genes_r5.57), size=length(distr)))

plot(c(.5,3.5), c(1.7,5.6), type="n", bty="n", ylab="", xlab="",
     xaxt="n", yaxt="n")
rect(xleft=.4, ybottom=1.55, xright=3.6, ytop=5.75,
   col="grey95", border=NA)
segments(x0=rep(.4,4), x1=rep(3.6,4), y0=2:5, y1=2:5,
   lwd=2, col="white")
boxplot(list(x1,x2,x3), outline=FALSE, names=FALSE,
   whisklty=2, whisklwd=2, staplelty=1, staplelwd=2, bty="n",
   notch=TRUE, boxlwd=2, ylab="Gene size [log10(bp)]",
   col=COLb, cex.axis=1.3, cex.lab=1.3, add=TRUE)
points(1+runif(length(distr), -.25,.25), x1,
    pch=21, cex=.6, col=COLa[1], bg=COLb[1])
points(2+runif(length(distr), -.25,.25), x2,
    pch=21, cex=.6, col=COLa[2], bg=COLb[2])
points(3+runif(length(distr), -.25,.25), x3,
    pch=21, cex=.6, col=COLa[3], bg=COLb[3])

dev.off()
