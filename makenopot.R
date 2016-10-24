require(GenomicRanges)

a = subset(read.delim("allprom_nochromP.txt"), prom != "p0")
ga = GRanges(Rle(a$chr), IRanges(start=a$pos, width=1))

EF1 = read.table("/data/eEF1delta_S2_rep1.peaks.txt")
GSE = read.table("/data/GSE57876_pnr_S2.peaks.txt")

# Remove "chr".
EF1$V1 = sub("^chr", "", EF1$V1)
GSE$V1 = sub("^chr", "", GSE$V1)

chroms = c("2L", "2R", "3L", "3R", "4", "X",
   "2LHet", "2RHet", "3LHet", "3RHet", "YHet")

pot.EF1 = list()
pot.GSE = list()

for (chrom in chroms) {
   fname = paste("/data/", chrom, ".mat.gz", sep="")
   write(paste("reading", chrom), file=stderr())
   # Takes long...
   HiC = as.matrix(read.table(fname))

   # Uniformize the HiC map.
   m = nrow(HiC)
   require(Matrix)
   toe = toeplitz(1:m)
   means = tapply(X=HiC, INDEX=toe, mean, na.rm=TRUE)
   tmp = matrix(c(means,0), nrow=m, ncol=m)
   tmp = lower.tri(tmp)*tmp
   tmp = tmp + t(tmp)
   diag(tmp) = means[1]
   HiC = tmp

   g. = GRanges(Rle(chrom),
        IRanges(start=seq(1, 1+(nrow(HiC)-1)*2000, by=2000),
                      end=seq(2000, nrow(HiC)*2000, by=2000)))
   HiC.n = as.matrix(HiC / rowSums(HiC))

   EF1. = subset(EF1, V1 == chrom)
   if (nrow(EF1.) > 0) {
      gEF1. = GRanges(Rle(EF1.$V1), IRanges(start=EF1.$V2, width=1))
      ov = as.matrix(findOverlaps(g., gEF1.))
      val = tapply(INDEX=ov[,1], X=EF1.$V3[ov[,2]], sum)
      pot = rep(0, nrow(HiC.n))
      pot[as.integer(names(val))] = val
      x = HiC.n %*% pot
      pot.EF1[[chrom]] = x
   }

   GSE. = subset(GSE, V1 == chrom)
   if (nrow(GSE.) > 1) {
      gGSE. = GRanges(Rle(GSE.$V1), IRanges(start=GSE.$V2, width=1))
      ov = as.matrix(findOverlaps(g., gGSE.))
      val = tapply(INDEX=ov[,1], X=GSE.$V3[ov[,2]], sum)
      pot = rep(0, nrow(HiC.n))
      pot[as.integer(names(val))] = val
      x = HiC.n %*% pot
      pot.GSE[[chrom]] = x
   }

}

save(pot.EF1, file="nopotEF1.rda")
save(pot.GSE, file="nopotGSE.rda")
