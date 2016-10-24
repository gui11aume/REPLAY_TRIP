require(GenomicRanges)

GSE = read.delim("GSE36175_norm_aggregated_tiling_arrays.txt.gz")

# Remove "chr".
GSE$chromosome = sub("^chr", "", GSE$chromosome)

chroms = c("2L", "2R", "3L", "3R", "4", "X",
   "2LHet", "2RHet", "3LHet", "3RHet", "YHet")

pot.all = list()

for (chrom in chroms) {
   fname = paste("/data/", chrom, ".mat.gz", sep="")
   # Takes long...
   HiC = read.table(fname)
   g. = GRanges(Rle(chrom),
        IRanges(start=seq(1, 1+(nrow(HiC)-1)*2000, by=2000),
                      end=seq(2000, nrow(HiC)*2000, by=2000)))

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

   HiC.n = as.matrix(HiC / rowSums(HiC))

   GSE. = subset(GSE, chromosome == chrom)
   if (nrow(GSE.) > 1) {
      gGSE. = GRanges(Rle(GSE.$chromosome),
            IRanges(start=(GSE.$start+GSE.$end)/2, width=1))
      ov = as.matrix(findOverlaps(g., gGSE.))
      val = matrix(0, nrow=nrow(HiC), ncol=ncol(GSE)-4)
      for (i in 1:ncol(val)) {
         tmp = tapply(X=GSE.[,i+4][ov[,2]], INDEX=ov[,1], mean)
         val[as.integer(names(tmp)),i] = tmp 
      }
      val[is.na(val)] = 0
      pot.all[[chrom]] = HiC.n %*% val 
   }   

}

save(pot.all, file="potall.rda")
