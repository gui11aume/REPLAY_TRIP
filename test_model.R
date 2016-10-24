require(GenomicRanges)

#diags = read.table("20diags.txt.gz")
# Read in TRIP data.
a = subset(read.delim("allprom.txt"), prom != "p0")
ga = GRanges(Rle(a$chr), IRanges(start=a$pos, width=1))

#prop5 = rowSums(diags[,5:6]) / rowSums(diags[,4:ncol(diags)])
#save(prop5, file="/data/prop5.rda")

#gdiags = GRanges(Rle(diags$V1),
#   IRanges(start=diags$V2, end=diags$V3))

#ov = as.matrix(findOverlaps(ga, gdiags))
#x = prop5[ov[,2]]
#y = a$nexp[ov[,1]]
#lm(y~x)
#anova(lm(y~x))

load("/data/potEF1.rda")
load("/data/potGSE.rda")
load("/data/potprom.rda")
load("/data/poterm.rda")
#load("potCAGE.rda")
#load("potall.rda")
#load("nopotall.rda")
#load("nopotprom.rda")
#load("potpromY.rda")

makepos = function(x) seq(from=1, by=2000, length.out=x)
lEF1 = sapply(pot.EF1, length)
lGSE = sapply(pot.GSE, length)
lprom = sapply(pot.prom, length)
lterm = sapply(pot.term, length)
#lCAGE = sapply(pot.CAGE, length)
#lall = sapply(pot.all, nrow)
#lnoall = sapply(nopot.all, nrow)
#lnoprom = sapply(nopot.prom, length)
#lpromY = sapply(pot.promY, length)


EF1 = data.frame(chr=rep(names(pot.EF1), times=lEF1),
         pos=unlist(sapply(lEF1, makepos)),
         value=log(0.01+unlist(pot.EF1)))
GSE = data.frame(chr=rep(names(pot.GSE), times=lGSE),
         pos=unlist(sapply(lGSE, makepos)),
         value=log(0.01+unlist(pot.GSE)))
prom = data.frame(chr=rep(names(pot.prom), times=lprom),
         pos=unlist(sapply(lprom, makepos)),
         value=log(0.01+unlist(pot.prom)))
term = data.frame(chr=rep(names(pot.term), times=lterm),
         pos=unlist(sapply(lterm, makepos)),
         value=log(0.01+unlist(pot.term)))
#CAGE = data.frame(chr=rep(names(pot.CAGE), times=lCAGE),
#         pos=unlist(sapply(lCAGE, makepos)), value=log(unlist(pot.CAGE)))
#aall = data.frame(chr=rep(names(pot.all), times=lall),
#         pos=unlist(sapply(lall, makepos)) , Reduce(rbind, pot.all))
#noaall = data.frame(chr=rep(names(nopot.all), times=lnoall),
#         pos=unlist(sapply(lnoall, makepos)) , Reduce(rbind, nopot.all))
#noprom = data.frame(chr=rep(names(nopot.prom), times=lnoprom),
#         pos=unlist(sapply(lnoprom, makepos)), value=unlist(nopot.prom))
#promY = data.frame(chr=rep(names(pot.promY), times=lpromY),
#         pos=unlist(sapply(lpromY, makepos)), value=unlist(pot.promY))

gEF1 = GRanges(Rle(EF1$chr), IRanges(start=EF1$pos, width=2000))
gGSE = GRanges(Rle(GSE$chr), IRanges(start=GSE$pos, width=2000))
gprom = GRanges(Rle(prom$chr), IRanges(start=prom$pos, width=2000))
gterm = GRanges(Rle(term$chr), IRanges(start=term$pos, width=2000))
#gCAGE = GRanges(Rle(CAGE$chr), IRanges(start=CAGE$pos, width=2000))
#gall = GRanges(Rle(aall$chr), IRanges(start=aall$pos, width=2000))
#gnoall = GRanges(Rle(noaall$chr), IRanges(start=noaall$pos, width=2000))
#gnoprom = GRanges(Rle(noprom$chr), IRanges(start=noprom$pos, width=2000))
#gpromY = GRanges(Rle(promY$chr), IRanges(start=promY$pos, width=2000))


predpow = rep(NA, 116)

ov = as.matrix(findOverlaps(ga, gEF1))
val = tapply(INDEX=ov[,2], X=a$nexp[ov[,1]], mean)
y = rep(NA, nrow(EF1))
y[as.integer(names(val))] = val

#x1 = EF1$value[ov[,2]]
#z1 = EF1$value[ov[,2]]^2

x1 = EF1$value
x2 = GSE$value
x3 = prom$value
x4 = term$value

mod1 = lm(y ~ x1+x1^2)
predpow[1] = var(mod1$fitted.values) /
   (var(mod1$residuals + mod1$fitted.values))
mod2 = lm(y ~ x2+x2^3)
predpow[2] = var(mod2$fitted.values) /
   (var(mod2$residuals + mod2$fitted.values))
mod3 = lm(y ~ x3+x3^3)
predpow[3] = var(mod3$fitted.values) /
   (var(mod3$residuals + mod3$fitted.values))
mod4 = lm(y ~ x4+x4^3)
predpow[4] = var(mod4$fitted.values) /
   (var(mod4$residuals + mod4$fitted.values))

GSE = read.delim("GSE36175_norm_aggregated_tiling_arrays.txt.gz")
GSE$chromosome = sub("^chr", "", GSE$chromosome)
gGSE = GRanges(Rle(GSE$chromosome),
         IRanges(start=GSE$start, end=GSE$end))

variables = list()
ov = as.matrix(findOverlaps(gGSE, gEF1))
for (i in 1:112) {
   val = tapply(INDEX=ov[,2], X=GSE[ov[,1],i+4], mean)
   x = rep(NA, nrow(EF1))
   x[as.integer(names(val))] = val
   mod = lm(y ~ x+x^3)
   predpow[i+4] = var(mod$fitted.values) /
      (var(mod$residuals + mod$fitted.values))
   nm = colnames(GSE)[i+4]
   variables[[paste(nm, 1, sep="_")]] = x
   variables[[paste(nm, 2, sep="_")]] = x^3
}

save(predpow, file="predpow.rda")

dat = data.frame(y, as.data.frame(variables))
mod = lm(dat)
var(mod$fitted.values) /
      (var(mod$residuals + mod$fitted.values))

dat = data.frame(y, as.data.frame(variables), x4,x4^3)
mod = lm(dat)
var(mod$fitted.values) /
      (var(mod$residuals + mod$fitted.values))

write.table(dat, file="/data/allpred.txt",
   sep="\t", quote=FALSE, row.names=FALSE)

library(glmnet)
dat = subset(dat, complete.cases(dat))
y = dat[,1]
x = scale(as.matrix(dat[,-1]), scale=TRUE)
cvfit = cv.glmnet(x=x, y=y)
save(cvfit, file="glmnet_fit.rda")

pred = predict(cvfit, newx=x, s="lambda.1se")
cor(pred, y)^2
