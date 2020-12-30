# Created on 09/07/13 by jdutheil
# Get block counts

loghist<-function(x, ...) {
  par(xaxt="n", mar=c(4,4,1,1)+0.1)
  hist(log(x, base=10), xlab="Block length", main=NA, xlim=c(0, 4), freq=TRUE, cex.axis=1.25, cex.lab=1.25, ...)
  par(xaxt="s")
  axis(1, at=0:4, labels=10^(0:4), cex.axis=1.25, cex.lab=1.25)
}

blocks.start<-read.table("chr1.statistics.csv", header=TRUE)
nrow(blocks.start)
sum(blocks.start$BlockLength)
loghist(blocks.start$BlockLength, col="cadetblue")

blocks.min100<-read.table("chr1.min100.statistics.csv", header=TRUE)
nrow(blocks.min100)
sum(blocks.min100$BlockLength)
loghist(blocks.min100$BlockLength, col="cornsilk3")

blocks.5sp<-read.table("chr1.5ways.statistics.csv", header=TRUE)
nrow(blocks.5sp)
sum(blocks.5sp$BlockLength)
loghist(blocks.5sp$BlockLength, col="cornsilk3")

blocks.merged<-read.table("chr1.merged.statistics.csv", header=TRUE)
nrow(blocks.merged)
sum(blocks.merged$BlockLength)
loghist(blocks.merged$BlockLength, col="cornsilk3")

blocks.nogap<-read.table("chr1.nogaps.statistics.csv", header=TRUE)
nrow(blocks.nogap)
sum(blocks.nogap$BlockLength)
loghist(blocks.nogap$BlockLength, col="cornsilk3")

blocks.window<-read.table("chr1.window1k.statistics.csv", header=TRUE)
nrow(blocks.window)
sum(blocks.window$BlockLength)
loghist(blocks.window$BlockLength, breaks=c(0,2.9,3.1,10), col="cornsilk3")

library(ape)

trees<-read.tree("chr1.window1k.dnd")
trees2<-lapply(trees, function(t) { t$edge.length<-rep(1, length(t$edge.length)); return(t) })
matrices<-lapply(trees2, cophenetic)
test1<-function(m) {
  rownames(m)<-sapply(strsplit(rownames(m), "\\."), function(x) x[1])
  colnames(m)<-sapply(strsplit(colnames(m), "\\."), function(x) x[1])
  return(m["hg19", "gorGor1"] == min(m[m>0]) & m["hg19", "panTro2"] < m["hg19", "ponAbe2"] & m["hg19", "ponAbe2"] < m["hg19", "rheMac2"])
}
sum(sapply(matrices, test1)) #613
test2<-function(m) {
  rownames(m)<-sapply(strsplit(rownames(m), "\\."), function(x) x[1])
  colnames(m)<-sapply(strsplit(colnames(m), "\\."), function(x) x[1])
  return(m["panTro2", "gorGor1"] == min(m[m>0]) & m["panTro2", "hg19"] < m["panTro2", "ponAbe2"] & m["panTro2", "ponAbe2"] < m["panTro2", "rheMac2"])
}
sum(sapply(matrices, test2)) #547

library(phangorn)
trees[sapply(matrices, test2)]
t<-trees[sapply(matrices, test2)]
plot(midpoint(ladderize(t[[8]])), no.margin=TRUE, edge.width=5, cex=1.5, x.lim=.1)

# Profiling
library(chron)
prof<-read.table("syrupy_20131115092051.ps.log", header=TRUE, stringsAsFactors=FALSE)
prof$TIME<-chron(times=prof$TIME)
prof$TIME<-prof$TIME-prof$TIME[1]
layout(matrix(1:2, nrow=2))
plot(CPU~TIME, prof, type="l", main="Processor usage (%)", ylim=c(0,100))
plot(RSS~TIME, prof, type="l", main="Memory usage (kB)")

dev.print(pdf, file="Profiling.pdf", width=8, height=12)
dev.print(svg, file="Profiling.svg", width=8, height=12)
