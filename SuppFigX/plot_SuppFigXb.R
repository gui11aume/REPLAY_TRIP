require(GenomicRanges)

load("../dmel_r5.57_FB2014_03/exons_r5.57.rda")
load("../dmel_r5.57_FB2014_03/introns_r5.57.rda")

ov = function(...) { suppressWarnings(countOverlaps(...)) }

# Create a 'GRanges' object from insertion data.
a = subset(read.delim("../allprom.txt"), prom != "p0")
ga = GRanges(Rle(a$chr), IRanges(a$pos, width=1))

ovexons = ov(ga, exons_r5.57, ignore.strand=TRUE)
ovintrons = ov(ga, introns_r5.57, ignore.strand=TRUE)

# Remove insertions in introns and exons.
a = subset(a, !ovexons & !ovintrons)
ga = GRanges(Rle(a$chr), IRanges(a$pos, width=1), nexp=a$nexp)

###########     Begin copy/paste     #########

predpow_out_genes = rep(NA, 118)

# Get the segmentation in 5 colors and 9 states.
C5 = read.delim("../color_domains.txt")
gC5 = GRanges(Rle(C5$chr), IRanges(start=C5$start, end=C5$end),
   col=C5$col)

S9 = read.table("../modENCODE_9states.txt")
gS9 = GRanges(Rle(S9$V2), IRanges(start=S9$V3, end=S9$V4),
   state=as.integer(S9$V1))

# Compute the predictive power of these models directly
# i.e. without binning by 2000 because it would break
# the domains.

ov = findOverlaps(ga, gC5)
y = mcols(ga[queryHits(ov)])$nexp
x = mcols(gC5[subjectHits(ov)])$col

modC5 = lm(y ~ x)
predpow_out_genes[117] = var(modC5$fitted.values) /
   (var(modC5$residuals + modC5$fitted.values))

ov = findOverlaps(ga, gS9)
y = mcols(ga[queryHits(ov)])$nexp
x = as.factor(mcols(gS9[subjectHits(ov)])$state)

modS9 = lm(y ~ x)
predpow_out_genes[118] = var(modS9$fitted.values) /
   (var(modS9$residuals + modS9$fitted.values))


load("/data/potEF1.rda")
load("/data/potGSE.rda")
load("/data/potprom.rda")
load("/data/poterm.rda")

makepos = function(x) seq(from=1, by=2000, length.out=x)
lEF1 = sapply(pot.EF1, length)
lGSE = sapply(pot.GSE, length)
lprom = sapply(pot.prom, length)
lterm = sapply(pot.term, length)


EF1 = data.frame(chr=rep(names(pot.EF1), times=lEF1),
         pos=unlist(sapply(lEF1, makepos)),
         value=log(0.01+unlist(pot.EF1)))
GSE = data.frame(chr=rep(names(pot.GSE), times=lGSE),
         pos=unlist(sapply(lGSE, makepos)),
         value=log(0.01+unlist(pot.GSE)))
prom = data.frame(chr=rep(names(pot.prom), times=lprom),
         pos=unlist(sapply(lprom, makepos)),
         value=log(0.01+unlist(pot.prom)))
term = data.frame(chr=rep(names(pot.term), times=lterm),
         pos=unlist(sapply(lterm, makepos)),
         value=log(0.01+unlist(pot.term)))

gEF1 = GRanges(Rle(EF1$chr), IRanges(start=EF1$pos, width=2000))
gGSE = GRanges(Rle(GSE$chr), IRanges(start=GSE$pos, width=2000))
gprom = GRanges(Rle(prom$chr), IRanges(start=prom$pos, width=2000))
gterm = GRanges(Rle(term$chr), IRanges(start=term$pos, width=2000))

ov = as.matrix(findOverlaps(ga, gEF1))
val = tapply(INDEX=ov[,2], X=a$nexp[ov[,1]], mean)
y = rep(NA, nrow(EF1))
y[as.integer(names(val))] = val

x1 = EF1$value;  z1 = x1^3
x2 = GSE$value;  z2 = x2^3
x3 = prom$value; z3 = x3^3
x4 = term$value; z4 = x4^3

mod1 = lm(y ~ x1+z1)
predpow_out_genes[1] = var(mod1$fitted.values) /
   (var(mod1$residuals + mod1$fitted.values))
mod2 = lm(y ~ x2+z2)
predpow_out_genes[2] = var(mod2$fitted.values) /
   (var(mod2$residuals + mod2$fitted.values))
mod3 = lm(y ~ x3+z3)
predpow_out_genes[3] = var(mod3$fitted.values) /
   (var(mod3$residuals + mod3$fitted.values))
mod4 = lm(y ~ x4+z4)
predpow_out_genes[4] = var(mod4$fitted.values) /
   (var(mod4$residuals + mod4$fitted.values))

GSE = read.delim("../GSE36175_norm_aggregated_tiling_arrays.txt.gz")
GSE$chromosome = sub("^chr", "", GSE$chromosome)
gGSE = GRanges(Rle(GSE$chromosome),
         IRanges(start=GSE$start, end=GSE$end))

variables = list()
ov = as.matrix(findOverlaps(gGSE, gEF1))
for (i in 1:112) {
   val = tapply(INDEX=ov[,2], X=GSE[ov[,1],i+4], mean)
   x = rep(NA, nrow(EF1))
   x[as.integer(names(val))] = val; z = x^3
   mod = lm(y ~ x+z)
   predpow_out_genes[i+4] = var(mod$fitted.values) /
      (var(mod$residuals + mod$fitted.values))
   nm = colnames(GSE)[i+4]
   variables[[paste(nm, 1, sep="_")]] = x
   variables[[paste(nm, 2, sep="_")]] = x^3
}

###########     End copy/paste     #########


load("../predpow.rda")
predpow_out_genes = predpow_out_genes[order(predpow, decreasing=TRUE)]

pdf("SuppFigXb.pdf", useDingbats=FALSE, height=5, width=7)
barplot(predpow_out_genes, ylab="Predictive power", axisnames=FALSE)
title(xlab="Features in the same order as in Figure 4a", line=.5)
dev.off()

