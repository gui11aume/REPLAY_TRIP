library(GenomicRanges)

load("../dmel_r5.57_FB2014_03/act_genes_r5.57.rda")
write.table(as.data.frame(act_genes_r5.57)[,1:5],
   file="act_genes_r5.57.txt", sep="\t", quote=FALSE, row.names=FALSE)
