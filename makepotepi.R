require(GenomicRanges)

# Load the H3K4me1 profile from Corces lab
EPI = read.table('/data/GSM890120_H3K4me1_Kc.bed')

#chroms = c("2L", "2R", "3L", "3R", "4", "X",
#   "2LHet", "2RHet", "3LHet", "3RHet", "YHet")

chroms = '2L'
pot.EPI = list()

for (chrom in chroms) {
   fname = paste("/data/", chrom, ".mat.gz", sep="")
   # Takes long...
   HiC = read.table(fname)
   g. = GRanges(Rle(chrom),
        IRanges(start=seq(1, 1+(nrow(HiC)-1)*2000, by=2000),
                      end=seq(2000, nrow(HiC)*2000, by=2000)))
   HiC.n = as.matrix(HiC / rowSums(HiC))

   EPI. = subset(EPI, V1 == chrom)
   if (nrow(EPI.) > 0) {
      # We use the middle point of the peak to assign it to 1 single window
      gEPI. = GRanges(Rle(EPI.$V1), IRanges(start=EPI.$V2, end=EPI.$V3))
#      ov = as.matrix(findOverlaps(g., gEF1.))
#      val = tapply(INDEX=ov[,1], X=EF1.$V3[ov[,2]], sum)
#      pot = rep(0, nrow(HiC.n))
#      pot[as.integer(names(val))] = val
#      x = HiC.n %*% pot
      ov = countOverlaps(g., gEPI.)
      x = HiC.n %*% ov
      pot.EPI[[chrom]] = x
   }
}

save(pot.EPI, file="potEPI.rda")
