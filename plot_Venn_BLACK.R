# Read in color domains.
coldoms = read.delim("color_domains.txt")
BLACK = subset(coldoms, col=="BLACK")

library(GenomicRanges)
gBLACK = GRanges(Rle(BLACK$chr),
   IRanges(start=BLACK$start, end=BLACK$end))

# Venn diagram.
load("gdom0.rda")
AB = sum(sum(coverage(intersect(gdom0, gBLACK)))) / 1e6
A = sum(sum(coverage(gdom0))) / 1e6
B = sum(sum(coverage(gBLACK))) / 1e6

source("helpers/Venn.R")
pdf("Venn.pdf")
Venn(A, B, AB)
dev.off()

# Resample 1000 times with circular permutation.
orig = tapply(INDEX=coldoms$chr, X=1:nrow(coldoms),
   FUN=function(x) coldoms[x,])

set.seed(123)

overlap = rep(NA, 1000)
for (iter in 1:1000) {
   colchrom = orig
   for (i in 1:length(colchrom)) {
      block = colchrom[[i]]
      nblock = nrow(block)
      idx = sample(nblock, size=1)
      if (idx > 1) {
         shift1 = block$start[idx] - block$start[1]
         shift2 = block$end[nblock] - block$end[idx-1]
         block$start[idx:nblock] = block$start[idx:nblock] - shift1
         block$end[idx:nblock] = block$end[idx:nblock] - shift1
         block$start[1:(idx-1)] = block$start[1:(idx-1)] + shift2
         block$end[1:(idx-1)] = block$end[1:(idx-1)] + shift2
         block = rbind(block[c(idx:nblock,1:(idx-1)),])
      }
      colchrom[[i]] = block
   } 

   resampled = Reduce(rbind, colchrom)
   BLACK = subset(resampled, col=="BLACK")

   gBLACK = GRanges(Rle(BLACK$chr),
      IRanges(start=BLACK$start, end=BLACK$end))

   # The resampling takes some time, mostly because of 'intersect()'.
   overlap[iter] = sum(sum(coverage(intersect(gdom0, gBLACK)))) / 1e6

}

max(overlap)
mean(overlap >= AB)
