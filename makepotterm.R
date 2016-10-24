require(GenomicRanges)

load("dmel_r5.57_FB2014_03/act_genes_r5.57.rda")
gterm = resize(act_genes_r5.57, fix="end", width=1)

chroms = c("2L", "2R", "3L", "3R", "4", "X",
   "2LHet", "2RHet", "3LHet", "3RHet", "YHet")

pot.term = list()

for (chrom in chroms) {
   fname = paste("/data/", chrom, "_norm.mat.gz", sep="")
   # Takes long...
   HiC = read.table(fname)
   g. = GRanges(Rle(chrom),
        IRanges(start=seq(1, 1+(nrow(HiC)-1)*2000, by=2000),
                      end=seq(2000, nrow(HiC)*2000, by=2000)))
   HiC.n = as.matrix(HiC / rowSums(HiC))

   ov = countOverlaps(g., gterm)
   pot.term[[chrom]] = HiC.n %*% ov

}

save(pot.term, file="/data/poterm.rda")
