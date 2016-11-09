# Read in the insertion data.
a = read.delim("../allprom_nochromP.txt", as.is=TRUE)
a = subset(a, prom != "p0")

# Subset promoter II vs all the rest on chromosome 2L.
aII = subset(a, prom == 'pII' & chr == '2L' & pos < 2e7 & pos > 1e7)
anII = subset(a, prom != 'pII' & chr == '2L' & pos < 2e7 & pos > 1e7)

require(GenomicRanges)
load("gdom0.rda")
gaII = GRanges(Rle(aII$chr), IRanges(start=aII$pos, width=1))
ganII = GRanges(Rle(anII$chr), IRanges(start=anII$pos, width=1))

ovII = countOverlaps(gaII, gdom0) > 0
ovnII = countOverlaps(ganII, gdom0) > 0
#dom2L = subset(as.data.frame(gdom0),
#   seqnames=="2L" & start < 2e7 & end > 1e7)

# These are seagreen3 and purple4 with transparency.
colors = c("#43CD8070", "#551A8B70")

pdf("Fig3b.pdf", height=6, width=10)
plot(aII$pos/1e6, aII$nexp,
     xlim=c(10,20),
     type="h", bty="n", yaxt="n", ylab="",
     xlab="Position on chr2L (Mb)",
     col=colors[1+ovII], ylim=c(-18,6))
axis(side=2, at=c(-5,0,5))
# Overlay plots on the same canvas.
par(new=TRUE)
plot(anII$pos/1e6, anII$nexp,
     xlim=c(10,20),
     type="h", bty="n", xaxt="n", yaxt="n", xlab="", ylab="",
     col=colors[1+ovnII], ylim=c(-6,18))
axis(side=2, at=c(-5,0,5))
#rect(xleft=dom2L$start/1e6, xright=dom2L$end/1e6,
#     ybottom=-6, ytop=18, border=NA, col="#551A8B20")
dev.off()
