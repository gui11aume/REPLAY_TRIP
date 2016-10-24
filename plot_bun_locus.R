require(GenomicRanges)

# Coordinates of the bun gene.
ID = "FBgn0259176"
chrom    = "2L"
bunstart = 12455540L
bunend   = 12546611L

genes = subset(
   read.table("dmel_r5.57_FB2014_03/dmel-all-genes-r5.57.gff",
      sep="\t", comm="", quote=""), V1 == "2L")

bun = subset(genes, grepl("ID=FBgn0259176;", V9))

# 22 exons of the bun gene.
exons = subset(read.table("dmel_r5.57_FB2014_03/dmel-all-exons-r5.57.gff"),
   grepl("^Name=bun:", V9))

left  = bun$V4 - 0.15*(bun$V5 - bun$V4)
right = bun$V5 + 0.1*(bun$V5 - bun$V4)

# CAGE tag counts on 2L by 300 bp windows.
load("CAGE2L.rda")
inside = (1:length(CAGE2Lplus)*300 > left) &
   (1:length(CAGE2Lplus)*300 < right)

# Read in the insertion data.
a = read.delim("allprom_nochromP.txt", as.is=TRUE)
a = subset(a, prom != "p0" & chr == "2L" & pos > left & pos < right)

# Read in enhancers.
EF1 = subset(read.table("/data/eEF1delta_S2_rep1.peaks.txt"),
   V1 == "chr2L" & V2 > left & V2 < right)
GSE = subset(read.table("/data/GSE57876_pnr_S2.peaks.txt"),
   V1 == "chr2L" & V2 > left & V2 < right)

# Read in chromatin proteins
chromP = subset(
   read.delim("GSE36175_norm_aggregated_tiling_arrays.txt.gz"),
   chromosome == "chr2L" & start > left & end < right)



pdf("bun_locus.pdf")
plot(x=c(left, right) / 1e6, y=c(0,20), type="n", xaxt="n",
     xlab="Position on chr2L (Mb)", ylab="", yaxt="n", bty="n")
segments(x0=bun$V4 / 1e6, x1=bun$V5 /1e6, y0=0.5, y1=0.5)
rect(xleft=exons$V4 / 1e6, xright=exons$V5/ 1e6, ytop=1,
     ybottom=0, col="black")
axis(1, at=seq(from=12.44, to=12.56, by=.02))

par(new=TRUE)
plot(x=c(left, right) / 1e6, y=c(-4,16), type="n",
     xlab="", ylab="", yaxt="n", xaxt="n", bty="n")
points(a$pos/1e6, a$nexp / 3, type="h")
segments(x0=left / 1e6, x1=right /1e6, y0=0, y1=0, col="grey50")

par(new=TRUE)
plot(x=c(left, right) / 1e6, y=c(-8,12), type="n",
     xlab="", ylab="", yaxt="n", xaxt="n", bty="n")
points((1:length(CAGE2Lplus)*300)[inside]/1e6,
   CAGE2Lplus[inside]/1000, type="h")
points((1:length(CAGE2Lminus)*300)[inside]/1e6,
   -CAGE2Lminus[inside]/1000, type="h")
segments(x0=left / 1e6, x1=right /1e6, y0=0, y1=0, col="grey50")

par(new=TRUE)
plot(x=c(left, right) / 1e6, y=c(-10.5,10.5), type="n",
     xlab="", ylab="", yaxt="n", xaxt="n", bty="n")
points(EF1$V2/1e6, EF1$V3 / 12, type="h")
segments(x0=left / 1e6, x1=right /1e6, y0=0, y1=0, col="grey50")

par(new=TRUE)
plot(x=c(left, right) / 1e6, y=c(-11.5,9.5), type="n",
     xlab="", ylab="", yaxt="n", xaxt="n", bty="n")
points(GSE$V2/1e6, GSE$V3 / 12, type="h")
segments(x0=left / 1e6, x1=right /1e6, y0=0, y1=0, col="grey50")

par(new=TRUE)
plot(x=c(left, right) / 1e6, y=c(-14,6), type="n",
     xlab="", ylab="", yaxt="n", xaxt="n", bty="n")
points(chromP$start/1e6, chromP[["H1"]]/5, type="h")
segments(x0=left / 1e6, x1=right /1e6, y0=0, y1=0, col="grey50")

par(new=TRUE)
plot(x=c(left, right) / 1e6, y=c(-15,5), type="n",
     xlab="", ylab="", yaxt="n", xaxt="n", bty="n")
points(chromP$start/1e6, chromP[["TBP"]]/5, type="h")
segments(x0=left / 1e6, x1=right /1e6, y0=0, y1=0, col="grey50")

par(new=TRUE)
plot(x=c(left, right) / 1e6, y=c(-16,4), type="n",
     xlab="", ylab="", yaxt="n", xaxt="n", bty="n")
points(chromP$start/1e6, chromP[["RPII18"]]/5, type="h")
segments(x0=left / 1e6, x1=right /1e6, y0=0, y1=0, col="grey50")

par(new=TRUE)
plot(x=c(left, right) / 1e6, y=c(-17,3), type="n",
     xlab="", ylab="", yaxt="n", xaxt="n", bty="n")
points(chromP$start/1e6, chromP[["BRM"]]/5, type="h")
segments(x0=left / 1e6, x1=right /1e6, y0=0, y1=0, col="grey50")

par(new=TRUE)
plot(x=c(left, right) / 1e6, y=c(-18,2), type="n",
     xlab="", ylab="", yaxt="n", xaxt="n", bty="n")
points(chromP$start/1e6, chromP[["H3K4ME2"]]/5, type="h")
segments(x0=left / 1e6, x1=right /1e6, y0=0, y1=0, col="grey50")

dev.off()
