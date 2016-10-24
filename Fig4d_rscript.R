# Load GRanges library
library(GenomicRanges)

#Load  kc167 endogenous expression data set:
endogenous<- read.delim('/rscripts/Kc_exp_color.bed')

#Make sure that genes in minus strand have the correct start:
positions<- rep(NA,times=nrow(endogenous))
positions[endogenous$strand=='+']<- endogenous$start[endogenous$strand=='+']
positions[endogenous$strand=='-']<- endogenous$end[endogenous$strand=='-']

#load buds and overlap both data sets:
buds<-read.table('/rscripts/BUDs.txt', stringsAsFactors=F, header=T)

endoGR<- GRanges(seqname= endogenous$chr, IRanges(start= positions, end= positions+1))
budsGR<- GRanges(seqnames= buds$chrom, ranges=IRanges(start= buds$start, end=buds$end))
ov<-findOverlaps(endoGR, budsGR, select='first')
endogenous$bud<- buds$puff[ov]


#Make correct format for plotting:
mat<-table(endogenous$color, endogenous$bud)

par(mfrow=c(3,1))
pie(table(endogenous$bud), col=c('slategrey', 'plum3'), labels="", border=F)
pie(mat[,2], col=names(mat[,1]), labels="", border=F)
pie(mat[,1], col=names(mat[,1]), labels="", border=F)

