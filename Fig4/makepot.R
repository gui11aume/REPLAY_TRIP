require(GenomicRanges)

load("../dmel_r5.57_FB2014_03/act_genes_r5.57.rda")
gprom = promoters(act_genes_r5.57, up=1, down=0)
gterm = resize(act_genes_r5.57, fix="end", width=1)

EF1 = read.table("/data/eEF1delta_S2_rep1.peaks.txt")
GSE = read.table("/data/GSE57876_pnr_S2.peaks.txt")
# Remove "chr".
EF1$V1 = sub("^chr", "", EF1$V1)
GSE$V1 = sub("^chr", "", GSE$V1)

gEF1 = GRanges(Rle(EF1$V1), IRanges(start=EF1$V2, width=1))
gGSE = GRanges(Rle(GSE$V1), IRanges(start=GSE$V2, width=1))

H3K27ac = subset(read.table("SRR442091_H3K27Ac_discretized.tsv"), V4 == 1)
H3K4me1 = subset(read.table("SRR442090_H3K4me1_discretized.tsv"), V4 == 1)

gH3K27ac = GRanges(H3K27ac$V1, IRanges(start=H3K27ac$V2, end=H3K27ac$V3))
gH3K4me1 = GRanges(H3K4me1$V1, IRanges(start=H3K4me1$V2, end=H3K4me1$V3))

ov = countOverlaps(gH3K27ac, gH3K4me1)
genh = gH3K27ac[ov > 0]

chroms = c("2L", "2R", "3L", "3R", "4", "X")

pot.prom = list()
pot.term = list()
pot.EF1  = list()
pot.GSE  = list()
pot.enh  = list()

for (chrom in chroms) {

   cat(paste(chrom, "\n", sep=""), file=stderr())
   fname = paste("/data/", chrom, "_norm.mat.gz", sep="")
   # Takes long...
   HiC = read.table(fname)
   g. = GRanges(Rle(chrom),
        IRanges(start=seq(1, 1+(nrow(HiC)-1)*2000, by=2000),
                      end=seq(2000, nrow(HiC)*2000, by=2000)))

   # Normalize in order to get a probability distribution.
   HiC.n = as.matrix(HiC / rowSums(HiC))

   ov = countOverlaps(g., gprom)
   pot.prom[[chrom]] = HiC.n %*% ov

   ov = countOverlaps(g., gterm)
   pot.term[[chrom]] = HiC.n %*% ov

   ov = countOverlaps(g., gEF1)
   pot.EF1[[chrom]] = HiC.n %*% ov

   ov = countOverlaps(g., gGSE)
   pot.GSE[[chrom]] = HiC.n %*% ov

   ov = countOverlaps(g., genh)
   pot.enh[[chrom]] = HiC.n %*% ov

}

save(pot.term, file="potterm.rda")
save(pot.prom, file="potprom.rda")
save(pot.EF1,  file="potEF1.rda")
save(pot.GSE,  file="potGSE.rda")
save(pot.enh,  file="potenh.rda")
