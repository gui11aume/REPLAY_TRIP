require(GenomicRanges)

# Read in color domains.
coldoms = read.delim("color_domains.txt")
RED = subset(coldoms, col=="RED")
YELLOW = subset(coldoms, col=="YELLOW")

# Read in paradoxical domains.
prdx = read.delim("prdx_domains.txt")

# Read in enhancers.
EF1 = read.table("/data/eEF1delta_S2_rep1.peaks.txt")
GSE = read.table("/data/GSE57876_pnr_S2.peaks.txt")

# Remove "chr" for compatibility.
EF1$V1 = sub("^chr", "", EF1$V1)
GSE$V1 = sub("^chr", "", GSE$V1)

gprdx = GRanges(Rle(prdx$block),
   IRanges(start=prdx$start, end=prdx$end))
gcol = GRanges(Rle(coldoms$chr),
   IRanges(start=coldoms$start, end=coldoms$end))
gRED = GRanges(Rle(RED$chr),
   IRanges(start=RED$start, end=RED$end))
gYELLOW = GRanges(Rle(YELLOW$chr),
   IRanges(start=YELLOW$start, end=YELLOW$end))

gEF1 = GRanges(Rle(EF1$V1), IRanges(start=EF1$V2, width=1))
gGSE = GRanges(Rle(GSE$V1), IRanges(start=GSE$V2, width=1))

a = 1e6 * (sum(countOverlaps(gEF1, gprdx) > 0) +
   sum(countOverlaps(gGSE, gprdx) > 0)) / sum(sum(coverage(gprdx)))

b = 1e6 * (sum(countOverlaps(gEF1, gRED) > 0) +
   sum(countOverlaps(gGSE, gRED) > 0)) / sum(sum(coverage(gRED)))

1e6 * (sum(countOverlaps(gEF1, gYELLOW) > 0) +
   sum(countOverlaps(gGSE, gYELLOW) > 0)) / sum(sum(coverage(gYELLOW)))

c = 1e6 * (sum(countOverlaps(gEF1, gcol) > 0) +
   sum(countOverlaps(gGSE, gcol) > 0)) / sum(sum(coverage(gcol)))

# Data for statistical tests.
print(sum(countOverlaps(gEF1, gprdx) > 0) +
   sum(countOverlaps(gGSE, gprdx) > 0))
print(sum(sum(coverage(gprdx))))

print(sum(countOverlaps(gEF1, gRED) > 0) +
   sum(countOverlaps(gGSE, gRED) > 0))
print(sum(sum(coverage(gRED))))


ov1 = countOverlaps(gEF1, gprdx) > 0
ov2 = countOverlaps(gGSE, gprdx) > 0
d = (sum(EF1$V3[ov1]) + sum(GSE$V3[ov2])) / (sum(ov1) + sum(ov2))

ov1 = countOverlaps(gEF1, gRED) > 0
ov2 = countOverlaps(gGSE, gRED) > 0
e = (sum(EF1$V3[ov1]) + sum(GSE$V3[ov2])) / (sum(ov1) + sum(ov2))

ov1 = countOverlaps(gEF1, gcol) > 0
ov2 = countOverlaps(gGSE, gcol) > 0
f = (sum(EF1$V3[ov1]) + sum(GSE$V3[ov2])) / (sum(ov1) + sum(ov2))

pdf("barplot_enhancers.pdf")
par(mfrow=c(1,2))
barplot(c(a,b,c))
barplot(c(d,e,f))
dev.off()
