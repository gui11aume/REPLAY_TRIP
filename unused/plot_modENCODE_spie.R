require(GenomicRanges)

S9 = read.table("modENCODE_9states.txt")
prdx = read.delim("prdx_domains.txt")

gS9 = GRanges(Rle(S9$V2), IRanges(start=S9$V3, end=S9$V4))
gprdx = GRanges(Rle(prdx$block),
   IRanges(start=prdx$start, end=prdx$end))

inter = intersect(gS9, gprdx)
ov = as.matrix(findOverlaps(inter, gS9))
coverage = width(ranges(inter)[ov[,1],])
states = S9[ov[,2],1]


source("helpers/spie.R")
expt = tapply(INDEX=S9$V1, X=S9$V4-S9$V3, sum)
obs = tapply(INDEX=states, X=coverage, sum)

pdf("modENCODE_spie.pdf")
par(mar=c(0,0,0,0))
spiechart(observed=obs, expected=expt, col=1:9)
dev.off()
