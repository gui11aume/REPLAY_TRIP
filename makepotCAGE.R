require(GenomicRanges)

a = subset(read.delim("allprom_nochromP.txt"), prom != "p0")
ga = GRanges(Rle(a$chr), IRanges(start=a$pos, width=1))

CAGE = read.table("/data/mappgedGAGE.txt")

chroms = c("2L", "2R", "3L", "3R", "4", "X",
   "2LHet", "2RHet", "3LHet", "3RHet", "YHet")

pot.CAGE = list()

for (chrom in chroms) {
   fname = paste("/data/", chrom, ".mat.gz", sep="")
   # Takes long...
   HiC = read.table(fname)
   g. = GRanges(Rle(chrom),
        IRanges(start=seq(1, 1+(nrow(HiC)-1)*2000, by=2000),
                      end=seq(2000, nrow(HiC)*2000, by=2000)))
   HiC.n = as.matrix(HiC / rowSums(HiC))

   CAGE. = subset(CAGE, V1 == chrom)
   if (nrow(CAGE.) > 0) {
      gCAGE. = GRanges(Rle(CAGE.$V1), IRanges(start=CAGE.$V2, width=1))
      ov = countOverlaps(g., gCAGE.)
      pot.CAGE[[chrom]] = HiC.n %*% ov
   }

}

save(pot.CAGE, file="potCAGE.rda")
