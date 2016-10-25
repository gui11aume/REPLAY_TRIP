require(GenomicRanges)

load("../dmel_r5.57_FB2014_03/act_exons_r5.57.rda")
load("../dmel_r5.57_FB2014_03/inact_exons_r5.57.rda")

ov = function(...) { suppressWarnings(countOverlaps(...)) }

a = read.delim("../allprom_nochromP.txt", comment.char="#")

matexons   = matrix(NA, nrow=4, ncol=5)

promoters = c("p0","pI", "pII", "pIII", "pIV")

for (i in 1:5) {

   # Take only current promoter and create a GenomicRanges object.
   p = subset(a, prom == promoters[i])
   pGR = GRanges(p$chr, IRanges(start=p$pos, width=1), p$strand)

   # Extract insertions in active exons and
   # introns, and compute the median.
   x1 = median(p$nexp[ov(pGR, act_exons_r5.57) > 0], na.rm=TRUE)
   x2 = median(p$nexp[ov(pGR, inact_exons_r5.57) > 0], na.rm=TRUE)

   # Reverse the strand and repeat.
   p$strand = ifelse(p$strand == "+", "-", "+")
   pGR = GRanges(p$chr, IRanges(start=p$pos, width=1), p$strand)

   x3 = median(p$nexp[ov(pGR, act_exons_r5.57) > 0], na.rm=TRUE)
   x4 = median(p$nexp[ov(pGR, inact_exons_r5.57) > 0], na.rm=TRUE)

   matexons[,i] = c(x1,x3,x2,x4)

}

COL = colorRampPalette(c("seagreen3","royalblue4","purple4"))(6)

pdf("SuppFigX.pdf", useDingbats=FALSE, height=5, width=7)
barplot(matexons, beside=TRUE, names.arg=promoters,
  ylab="Median expression (log2)", col=rep(COL, each=4))
legend(x="topright", bty="n", cex=.8,
  legend=c("Active forward", "Active reverse",
  "Inactive forward", "Inactive reverse"))
dev.off()
