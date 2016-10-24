require(GenomicRanges)

# Keep complete cases and remove control promoter.
a = read.delim("allprom.txt")
a = subset(a, prom != "p0")
a = subset(a, complete.cases(a))

ga = GRanges(Rle(a$chr), IRanges(start=a$pos, width=1))
 
load("gdom0.rda")

ov0 = countOverlaps(ga, gdom0) > 0
a0 = as.matrix(a[ov0,8:ncol(a)])

pca = prcomp(a0[,8:ncol(a0)], scale.=TRUE)

pdf("PCA.pdf", width=6, height=6)
plot(pca$x, pch=19, cex=.3, col="#551A8B20")
dev.off()

# The outliers are approximately below the line (-10,-1).
outliers = pca$x[,2] + pca$x[,1] < -10

chromout = colMeans(a0[which(outliers),8:ncol(a0)])
chromin  = colMeans(a0[which(!outliers),8:ncol(a0)])
names    = colnames(a0)[8:ncol(a0)]

x = cbind(chromin, chromout)
h = hclust(dist(x))
x = x[h$order,]
colnames(x) = NULL

pdf("outliers_heatmap.pdf")
heatmap(x, Rowv=NA, Colv=NA, scale="none", cexRow=.6)
dev.off()

# Here show the min and the max for scaling purposes.
min(x) # -0.9166672
max(x) # 1.516331

# Use a binary HMM to find paradoxical domains.
library(HMMt)
source("helpers/BaumWelchBin.R")

prdxcal = rep(0L, nrow(a))
prdxcal[which(ov0)[outliers]] = 1L

mean(prdxcal)

b = bridge(data.frame(a[,c(2,4,4)], prdxcal))

fit = BaumWelchBin(b$x, b$series.length, m=2)

state = fit$ViterbiPath[b$nonvirtuals]-1

source("helpers/domainify.R")
prdx = domainify(data.frame(a[,c(2,4,4)], state))
write.table(prdx, file="prdx_domains.txt", sep="\t",
   quote=FALSE, row.names=FALSE)

dim(prdx)
