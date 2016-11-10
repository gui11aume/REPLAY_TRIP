# Load the data
setwd('/data/rep_timing_raw')

# Early replicating
early.rep1 <- read.table('SRR586004_RT.tsv')
names(early.rep1) <- c('chr','reads')
early.rep2 <- read.table('SRR586008_RT.tsv')
names(early.rep2) <- c('chr','reads')
early <- data.frame(chr=early.rep1$chr, 
                    ereads=rowMeans(cbind(early.rep1$reads,early.rep2$reads)))

# MidEarly replicating
midearly.rep1 <- read.table('SRR586005_RT.tsv')
names(midearly.rep1) <- c('chr','reads')
midearly.rep2 <- read.table('SRR586009_RT.tsv')
names(midearly.rep2) <- c('chr','reads')
midearly <- data.frame(chr=midearly.rep1$chr, 
                    mereads=rowMeans(cbind(midearly.rep1$reads,midearly.rep2$reads)))

# Midlate replicating
midlate.rep1 <- read.table('SRR586006_RT.tsv')
names(midlate.rep1) <- c('chr','reads')
midlate.rep2 <- read.table('SRR586010_RT.tsv')
names(midlate.rep2) <- c('chr','reads')
midlate <- data.frame(chr=midlate.rep1$chr, 
                    mlreads=rowMeans(cbind(midlate.rep1$reads,midlate.rep2$reads)))

# Late replicating
late.rep1 <- read.table('SRR586007_RT.tsv')
names(late.rep1) <- c('chr','reads')
late.rep2 <- read.table('SRR586011_RT.tsv')
names(late.rep2) <- c('chr','reads')
late <- data.frame(chr=late.rep1$chr, 
                    lreads=rowMeans(cbind(late.rep1$reads,late.rep2$reads)))

# All times
allt <- data.frame(early$chr,early$ereads,midearly$mereads,midlate$mlreads,late$lreads)
names(allt)[1] <- 'chr'

# Normalize for sequencing coverage to compare among times (permilion proportion)
alltn <- data.frame(allt[,1],apply(allt[2:5],2, function(col) {col * 1e6/sum(col)}))
names(alltn) <- c('chr','early','midearly','midlatea','late')

# 2 tracks are enough Early and Late
alltn$soon <- alltn$early + alltn$midearly
alltn$later <- alltn$midlate + alltn$late

alltn <- alltn[,c(1,6,7)]

write.table(alltn,file = '/data/rep_timing_tracks.bed',
            col.names=F,row.names=F, quote = F, sep = '\t')


# Prepare chr2L
twoL <- subset(alltn, alltn$chr == '2L')
twoL$pos <- 1:length(twoL$chr)

# Up and Down domains in chr2L
library(GenomicRanges)
load('/home/REPLAY_TRIP/gdom0.rda') # Down domains
load('/home/REPLAY_TRIP/gdom1.rda') # Up domains
ups <- as.data.frame(gdom0)
ups <- ups[,1:3]
downs <- as.data.frame(gdom1)
downs <- downs[,1:3]
ups2L <- subset(ups, ups$seqnames == '2L')
downs2L <- subset(downs, downs$seqnames == '2L')


# Load Matrix chr2L
mtx <- read.table('/data/2L.mat.gz')

# Plot HiC + rep timing + up and downs
png(file='/data/rep_timing_chr2L.png',width = 1000, height = 1000)
layout(matrix(c(1, 2, 3), 3, 1), heights=c(2,2,8))
par(mar=c(0.25, 0.5, 0.25, 0.5))

## Top pannel
plot(twoL$pos,twoL$soon,type = 'h',
     col='seagreen',ylim = c(-120,120),
     lwd=1.5,xlab=NA,ylab=NA,axes=F,
     xaxt='n', yaxt='n', xaxs='i', yaxs='i')
lines(twoL$pos,-twoL$later,type = 'h',
      col='plum3')

## Second pannel
plot(0, 
     ylim=c(0, 1), xlab=NA, ylab=NA,
     xaxt='n', yaxt='n', xaxs='i',
     yaxs='i',axes=F)
rect(xleft=ups2L$start,xrigth=ups2L$end,
     ytop=10,ybottom=0,col='seagreen')
rect(xleft=downs2L$start,xrigth=downs2L$end,
     ytop=0,ybottom=-10, col='plum3')

## Third pannel
par(mar=c(0.5, 0.5, 0, 0.5))
grad = colorRampPalette(c('white','black'))(12)
image(log(1 + mat2L), col=grad)





dev.off()

