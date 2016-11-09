source("../helpers/violin-plot.R")

a = read.delim("../allprom_nochromP.txt", comment.char="#")

allvio = list()
for (thisprom in c("p0","pI","pII","pIII","pIV")) {
   p = subset(a, prom == thisprom)
   Het = subset(p, grepl("Het|U", p$chr))
   noHet = subset(p, !grepl("Het|U", p$chr))
   allvio = c(allvio, list(noHet$nexp, Het$nexp))
}
COL = colorRampPalette(c("seagreen3", "royalblue4","purple4"))(6)
pdf("Fig3a.pdf", height=4, width=4)
violin.plot(allvio, col=rep(COL, each=2),
   x.pos=c(1,2, 4,5, 7,8, 10,11, 13,14))
dev.off()
