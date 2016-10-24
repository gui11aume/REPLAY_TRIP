require(GenomicRanges)

load("dmel_r5.57_FB2014_03/act_promoters_r5.57.rda")

chroms = c("2L", "2R", "3L", "3R", "4", "X",
   "2LHet", "2RHet", "3LHet", "3RHet", "YHet")

YELLOW = subset(read.delim("color_domains.txt"), col=="YELLOW")
gY = GRanges(Rle(YELLOW$chr), IRanges(start=YELLOW$start, end=YELLOW$end))
ov = countOverlaps(act_promoters_r5.57, gY) > 0
gpromY = act_promoters_r5.57[ov]

pot.prom = list()
pot.promY = list()

for (chrom in chroms) {
   fname = paste("/data/", chrom, ".mat.gz", sep="")
   # Takes long...
   HiC = read.table(fname)
   g. = GRanges(Rle(chrom),
        IRanges(start=seq(1, 1+(nrow(HiC)-1)*2000, by=2000),
                      end=seq(2000, nrow(HiC)*2000, by=2000)))
   HiC.n = as.matrix(HiC / rowSums(HiC))

   ov = countOverlaps(g., act_promoters_r5.57)
   x = HiC.n %*% ov
   pot.prom[[chrom]] = x
   ov2 = countOverlaps(g., gpromY)
   x = HiC.n %*% ov2
   pot.promY[[chrom]] = x

}

save(pot.prom, file="potprom.rda")
save(pot.promY, file="potpromY.rda")
