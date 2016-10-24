require(GenomicRanges)

# Read in data.
chromP = read.delim("GSE36175_norm_aggregated_tiling_arrays.txt.gz")
coldoms = read.delim("color_domains.txt")
prdx = read.delim("prdx_domains.txt")

# Keep only RED domains.
reddoms = subset(coldoms, col == "RED")
# Remove "chr" from chromosome names
chromP$chromosome = sub("^chr", "", chromP$chromosome)

gchromP = GRanges(Rle(chromP$chromosome),
            IRanges(start=chromP$start, end=chromP$end))
# Shorten the intervals of GATC fragments to make sure that
# they do not overlap with each other (can bias the overlaps).
gchromP = narrow(gchromP, start=3, end=-3)
gred = GRanges(Rle(reddoms$chr),
            IRanges(start=reddoms$start, end=reddoms$end))
gprdx = GRanges(Rle(prdx$block),
            IRanges(start=prdx$start, end=prdx$end))

ov1 = countOverlaps(gchromP, gred) > 0
ov2 = countOverlaps(gchromP, gprdx) > 0

pdf("prdx_vs_RED.pdf", useDingbats=FALSE)
x = colMeans(chromP[ov1,5:ncol(chromP)], na.rm=TRUE)
y = colMeans(chromP[ov2,5:ncol(chromP)], na.rm=TRUE)
plot(x,y, type="n",
     xlab = "RED domains", ylab= "Paradoxical domains")
rect(xleft=-1.5, xright=2, ybottom=-1.5, ytop=2, col="grey96",
     border=NA)
segments(x0=-1.5, x1=2, y0=c(-1,0,1), y1=c(-1,0,1),
         lwd=2, col="white")
segments(x0=c(-1,0,1,2), x1=c(-1,0,1,2), y0=-1.5, y1=2,
         lwd=2, col="white")
points(x, y, pch=19, cex=.8)
dev.off()
