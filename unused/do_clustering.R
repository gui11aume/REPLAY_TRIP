library(GenomicRanges)
a = read.delim("allprom.txt")
a = subset(a, prom != "p0")
a = subset(a, complete.cases(a))
#ga = GRanges(Rle(a$chr), IRanges(start=a$pos, width=1))
 
#load("gdom0.rda")
#load("gdom1.rda")

#ov0 = countOverlaps(ga, gdom0) > 0
#ov1 = countOverlaps(ga, gdom1) > 0

#a0 = as.matrix(a[ov0,8:ncol(a)])
#a1 = as.matrix(a[ov1,8:ncol(a)])

#if (!file.exists("h0.rda")) {
#   h0 = hclust(dist(a0))
#   save(h0, file="h0.rda")
#} else {
#   load("h0.rda")
#}
#
#if (!file.exists("h1.rda")) {
#   h1 = hclust(dist(a1))
#   save(h1, file="h1.rda")
#} else {
#   load("h1.rda")
#}
#
#g = hclust(dist(t(rbind(a0,a1))))
#
#bw=colorRampPalette(c("purple3", "purple3", "white", "seagreen3", "seagreen3"))(512)
#
#n0 = nrow(a0)
#n1 = nrow(a1)
#
#png("/data/heatmap.png", width=(n0+n1)/16, height=(n0+n1)/16)
#image(t(rbind(a1[h1$order,g$order], a0[h0$order,g$order])),
#     col=bw, xaxt="n", yaxt="n")
#abline(h=n1/(n0+n1), lwd=10)
#dev.off()

#a = rbind(data.frame(x=0,a0),data.frame(x=1,a1))
#a = a[sample(nrow(a)),]

#pca = prcomp(a[,-1], scale.=TRUE)
pca = prcomp(a[,8:ncol(a)], scale.=TRUE)
pdf("/data/PCA.pdf")
plot(pca$x, pch=19, cex=.3,
   col=c("#551A8B20","#43CD8020")[1 + (a$nexp > 0)])
plot(pca$x[,3:2], pch=19, cex=.3,
   col=c("#551A8B20","#43CD8020")[1 + (a$nexp > 0)])
dev.off()

#colnames(a0)[g$order]
