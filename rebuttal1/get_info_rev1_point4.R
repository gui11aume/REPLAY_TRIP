library(GenomicRanges)

# Create a 'GRange' object for genes, without strand information.
genes = read.table("../dmel_r5.57_FB2014_03/data/Kc_exp_color_MC.txt")
ggenes = GRanges(Rle(genes$V1), IRanges(start=genes$V2, end=genes$V3),
   expr=log2(.01 + genes$V5))

ov = function(...) { suppressWarnings(findOverlaps(...)) }

# Create a 'GRanges' object from insertion data.
a = subset(read.delim("../allprom.txt"), prom != "p0")
ga = GRanges(Rle(a$chr), IRanges(a$pos, width=1), nexp=a$nexp)

# First approach: get the flanking genes and possibly the gene(s)
# where the reporter is inserted, compute their average and use
# it as a predictor.

# Get genes overlapping with the insertions.
ovg = ov(ga, ggenes, ignore.strand=TRUE)

# Get the gene before and the gene after.
before = follow(ga, ggenes, ignore.strand=TRUE)
after = precede(ga, ggenes, ignore.strand=TRUE)

X = c(
   mcols(ggenes)$expr[subjectHits(ovg)],
   mcols(ggenes)$expr[before],
   mcols(ggenes)$expr[after]
)

INDEX = c(
   queryHits(ovg),
   1:length(ga),
   1:length(ga)
)

# Nearest neighbor mean.
nnmean = tapply(X=X, INDEX=INDEX, FUN=mean, na.rm=TRUE)
aid = as.integer(names(nnmean))

y1 = nnmean
x1 = mcols(ga[aid])$nexp; z1 = x1^2

mod1 = lm(y1 ~ x1+z1)
predpow1 = var(mod1$fitted.values) /
   (var(mod1$residuals + mod1$fitted.values))


# Second approach, take a window of fixed size around the
# insertion and compute the average of all genes that intersect
# the window.
ga20kb = resize(ga, width=20000)

# Get genes overlapping with the windows
ov20kb = ov(ga20kb, ggenes, ignore.strand=TRUE)

X = mcols(ggenes)$expr[subjectHits(ov20kb)]
INDEX = queryHits(ov20kb)

# Nearest neighbor mean.
nnmean = tapply(X=X, INDEX=INDEX, FUN=mean, na.rm=TRUE)
aid = as.integer(names(nnmean))

y2 = nnmean
x2 = mcols(ga20kb[aid])$nexp; z2 = x2^2

mod2 = lm(y2 ~ x2+z2)
predpow2 = var(mod2$fitted.values) /
   (var(mod2$residuals + mod2$fitted.values))

# Start writing to file.
sink(file("info_rev1_point4.txt"))

cat("Predictive power of the first approach (using two\n")
cat("flanking genes plus eventual insertion gene to predict\n")
cat("the expression of the reporter):\n")

predpow1

cat("\n")
cat("Predictive power of the second approach (using windows\n")
cat("to predict the expression of the reporter):\n")

predpow2

# End writing to file.
sink(NULL)
