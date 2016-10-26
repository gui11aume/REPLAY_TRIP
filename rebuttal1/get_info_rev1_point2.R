library(GenomicRanges)

load("../dmel_r5.57_FB2014_03/exons_r5.57.rda")
load("../dmel_r5.57_FB2014_03/introns_r5.57.rda")

ov = function(...) { suppressWarnings(countOverlaps(...)) }
dT = function(...) { suppressWarnings(distanceToNearest(...)) }

# Create a 'GRanges' object from insertion data.
a = subset(read.delim("../allprom_nochromP.txt"), prom != "p0")
ga = GRanges(Rle(a$chr), IRanges(a$pos, width=1))

ovexons = ov(ga, exons_r5.57, ignore.strand=TRUE)
ovintrons = ov(ga, introns_r5.57, ignore.strand=TRUE)

# Remove insertions in introns and exons.
a = subset(a, !ovexons & !ovintrons)
ga = GRanges(Rle(a$chr), strand=Rle(a$strand), IRanges(a$pos, width=1))

# Extract active TSS from the positions of active genes.
load("../dmel_r5.57_FB2014_03/act_genes_r5.57.rda")
gTSS = promoters(act_genes_r5.57, up=0, down=1)

dplus = mcols(dT(ga, gTSS))$distance

# Reverse insertions in order to find closest hit on the
# opposite strand.
strand(ga) = ifelse(strand(ga) == "+", "-", "+")

dminus = mcols(dT(ga, gTSS))$distance

aplus = a[which(dplus < dminus),]
aminus = a[which(dminus < dplus),]

med = median(abs(dplus-dminus), na.rm=TRUE)
pmed = median(pmin(dplus, dminus), na.rm=TRUE)

# Start writing to file.
sink(file("info_rev1_point2.txt"))


cat(paste("Closer to promoter in same orientation: ",
   nrow(aplus), "\n", sep=""))
cat(paste("Closer to promoter in opposite orientation: ",
   nrow(aminus), "\n", sep=""))
cat(paste("Median absolute difference of distance: ",
   med, "\n", sep=""))
cat(paste("Median distance to nearest: ",
   pmed, "\n", sep=""))

wilcox.test(aplus$nexp, aminus$nexp)

# End writing to file.
sink(NULL)
