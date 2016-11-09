source("../helpers/violin-plot.R")
require(GenomicRanges)

# Read in the insertion data.
a = read.delim("../allprom.txt", as.is=TRUE)
a = subset(a, prom != "p0")

# Read in color domains.
coldoms = read.delim("../color_domains.txt")
BLACK = subset(coldoms, col=="BLACK")

ga = GRanges(Rle(a$chr), IRanges(start=a$pos, width=1))
gcol = GRanges(Rle(coldoms$chr),
   IRanges(start=coldoms$start, end=coldoms$end))

ov = as.matrix(findOverlaps(ga, gcol))
l = tapply(X=a$nexp[ov[,1]], INDEX=coldoms$col[ov[,2]], c)

ord = order(c("BLACK","RED", "YELLOW", "BLUE", "GREEN"))
colors=c("black","red", "gold2", "blue", "seagreen3")

pdf("Fig3d.pdf", width=4, height=4)
violin.plot(l[ord], col=colors)
dev.off()

