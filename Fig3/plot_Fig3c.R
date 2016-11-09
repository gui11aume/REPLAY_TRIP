source("../helpers/violin-plot.R")
require(GenomicRanges)

a = read.delim("../allprom.txt", comment.char="#")

load("gdom0.rda")
load("gdom1.rda")

allvio = list()
for (thisprom in c("p0","pI","pII","pIII","pIV")) {
   p = subset(a, prom == thisprom)
   gp = GRanges(Rle(p$chr), IRanges(start=p$pos, width=1))
   dom0 = subset(p, countOverlaps(gp, gdom0) > 0)
   dom1 = subset(p, countOverlaps(gp, gdom1) > 0)
   allvio = c(allvio, list(dom0$nexp, dom1$nexp))
}
COL = colorRampPalette(c("seagreen3", "royalblue4","purple4"))(6)
pdf("Fig3c.pdf", height=4, width=4)
violin.plot(allvio, col=rep(COL, each=2),
   x.pos=c(1,2, 4,5, 7,8, 10,11, 13,14))
dev.off()
