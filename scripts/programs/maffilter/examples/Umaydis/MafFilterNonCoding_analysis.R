# Created on July 11, 2013 by jdutheil
# Analize output of maffilter.

d<-read.table("Umaydis.noncoding.statistics.csv", header=TRUE)
hist(log(d$BlockLength))
shapiro.test(log(sample(d$BlockLength,5000)))

#Try to fit an exponential distribution (and a gamma one):
library(MASS)
e<-fitdistr(d$BlockLength, "exponential")
g<-fitdistr(d$BlockLength, "gamma")
l<-fitdistr(d$BlockLength, "log-normal")
hist(d$BlockLength, prob=TRUE, nclass=300, xlim=c(0,3000))
curve(dexp(x, rate=e$estimate), col="blue", add=TRUE, lwd=2)
curve(dgamma(x, shape=g$estimate[1], rate=g$estimate[2]), col="red", add=TRUE, lwd=2)
curve(dlnorm(x, meanlog=l$estimate[1], sdlog=l$estimate[2]), col="red", add=TRUE, lwd=2)
# The log-normal fit is good

1/e$estimate #1069
1/g$estimate[2] #825
exp(l$estimate[1]) #694

sum(d$BlockLength>=300)

# GC content
d$GC<-(d$Count.G + d$Count.C) / (d$Count.A + d$Count.C + d$Count.G + d$Count.T)

gc<-d$GC[d$BlockLength>=300]
hist(gc, nclass=50, prob=TRUE)
n<-fitdistr(gc, "normal")
h<-fitdistr(gc, "cauchy")
k<-fitdistr(gc, "log-normal")
l<-fitdistr(gc, "logistic")
curve(dnorm(x, mean=n$estimate[1], sd=n$estimate[2]), col="blue", add=TRUE, lwd=2)
curve(dcauchy(x, location=n$estimate[1], scale=n$estimate[2]), col="red", add=TRUE, lwd=2)
curve(dlnorm(x, meanlog=k$estimate[1], sdlog=k$estimate[2]), col="green", add=TRUE, lwd=2)
curve(dlogis(x, location=l$estimate[1], scale=l$estimate[2]), col="orange", add=TRUE, lwd=2)

l$estimate

# Check per chromosome:
chrsizes<-c(chr01=2476500, chr02=1879391, chr03=1633472, chr04=885077, chr05=1393418, chr06=1031380, chr07=957188,
            chr08=813246, chr09=733964, chr10=692355, chr11=690617, chr12=650985, chr13=606072, chr14=611467,
            chr15=575418, chr16=552767, chr17=576627, chr18=560726, chr19=571809, chr20=523884, chr21=470506, chr22=403592, chr23=344927)
d$ChrSize<-chrsizes[d$Chr]
boxplot(GC~ChrSize, subset=BlockLength>=300, d)
cor.test(~GC+ChrSize,d, subset=BlockLength>=300, method="kendal")

# Figure for article:
library(gplots)

la<-matrix(1:4, nrow=2, byrow=FALSE)
la[2,2]<-3
layout(la)
hist(d$BlockLength, prob=TRUE, nclass=100, xlim=c(0,5000), main="Distribution of block length", xlab="Length of intergenic region", col=grey(0.7), ylim=c(0, 1e-3))
#curve(dexp(x, rate=e$estimate), col="black", add=TRUE, lwd=2)
l<-fitdistr(d$BlockLength, "log-normal")
curve(dlnorm(x, meanlog=l$estimate[1], sdlog=l$estimate[2]), col="black", add=TRUE, lwd=2)

hist(gc, nclass=50, prob=TRUE, main="Distribution of GC content", xlab="GC content in block", col=grey(0.7))
l<-fitdistr(gc, "logistic")
curve(dlogis(x, location=l$estimate[1], scale=l$estimate[2]), col="black", add=TRUE, lwd=2)

lst<-with(subset(d, BlockLength>=300), split(GC, ChrSize))

#plotmeans(GC~ChrSize, subset=BlockLength>=300, d)
av<-sapply(lst, mean)
dn<-sapply(lst, quantile, prob=0.05)
up<-sapply(lst, quantile ,prob=0.95)
plotCI(as.numeric(names(lst)), av, uiw=up-av, liw=av-dn, xlab="Chromosome size", ylab="GC content", main="GC content per chromosome", pch=19, col=grey(0.6))
m<-lm(GC~ChrSize, subset=BlockLength>=300, d)
abline(m, lwd=2)

dev.print(pdf, "Umaydis_intergenic_GC.pdf", width=10, height=6)
dev.print(svg, "Umaydis_intergenic_GC.svg", width=10, height=6)

### Why a logistic distribution ? ###

# If GC is homogeneous along a given region,
# The probability of being GC at a given position is a bernoulli variable of parameter p,
# The GC content of a segment of lengh n follow a binomial distribution of parameters p and n.

gc<-function(p, n) rbinom(1, prob=p, size=n) / n

sim<-function(rep, p , n) {
  x<-numeric(rep)
  for (i in 1:rep) {
    x[i]<-gc(p, n)
  }
  fl<-fitdistr(x, "logistic")
  fn<-fitdistr(x, "normal")
  flf<-function(x) dlogis(x, location=fl$estimate[1], scale=fl$estimate[2])
  fnf<-function(x) dnorm(x, mean=fl$estimate[1], sd=fl$estimate[2])
  mx<-max(flf(p), fnf(p))

  hist(x, prob=TRUE, main=paste("Segment length =", n), ylim=c(0, mx), xlim=c(0.4, 0.6), col=grey(0.7), border=grey(0.7), nclass=30)
  curve(flf, col="black", add=TRUE, lwd=2)
  curve(fnf, col="black", add=TRUE, lwd=2, lty="dashed")
  legend(0.4, mx, legend=c("Logistic", "Normal"), lty=c("solid", "dashed"), lwd=2, col="black", box.lty=0)
  return(x)
}

simRnd<-function(rep, p, meanLog, sdLog) {
  x<-numeric(rep)
  for (i in 1:rep) {
    x[i]<-gc(p, round(rlnorm(1, meanlog=meanLog, sdlog=sdLog)))
  }
  fl<-fitdistr(x, "logistic")
  fn<-fitdistr(x, "normal")
  flf<-function(x) dlogis(x, location=fl$estimate[1], scale=fl$estimate[2])
  fnf<-function(x) dnorm(x, mean=fl$estimate[1], sd=fl$estimate[2])
  mx<-max(flf(p), fnf(p))

  hist(x, prob=TRUE, main="Segment length = LogNormal", ylim=c(0, mx), xlim=c(0.4, 0.6), col=grey(0.7), border=grey(0.7), nclass=30)
  curve(flf, col="black", add=TRUE, lwd=2)
  curve(fnf, col="black", add=TRUE, lwd=2, lty="dashed")
  legend(0.4, mx, legend=c("Logistic", "Normal"), lty=c("solid", "dashed"), lwd=2, col="black", box.lty=0)
  return(x)
}

layout(matrix(1:4, nrow=2, byrow=TRUE))
sim(20000, 0.5, 100)
sim(20000, 0.5, 500)
sim(20000, 0.5, 1000)
simRnd(20000, 0.5, 6.54, 0.96)

dev.print(pdf, file="SimulateGC50.pdf", width=12, height=12)

simRnd2<-function(rep, p, meanLog, sdLog, sdGC) {
  x<-numeric(rep)
  for (i in 1:rep) {
    #g<-rnorm(1, mean=p, sd=sdGC)
    g<-sample(c(0.3, 0.7), 1)
    if (g<0) g<-0
    if (g>1) g<-1
    x[i]<-gc(g, round(rlnorm(1, meanlog=meanLog, sdlog=sdLog)))
  }
  fl<-fitdistr(x, "logistic")
  fn<-fitdistr(x, "normal")
  flf<-function(x) dlogis(x, location=fl$estimate[1], scale=fl$estimate[2])
  fnf<-function(x) dnorm(x, mean=fl$estimate[1], sd=fl$estimate[2])
  mx<-max(flf(p), fnf(p))

  hist(x, prob=TRUE, main=paste("Sd GC =", sdGC), ylim=c(0, mx), xlim=c(0.2, 0.8), col=grey(0.7), border=grey(0.7), nclass=30)
  curve(flf, col="black", add=TRUE, lwd=2)
  curve(fnf, col="black", add=TRUE, lwd=2, lty="dashed")
  legend(0.2, mx, legend=c("Logistic", "Normal"), lty=c("solid", "dashed"), lwd=2, col="black", box.lty=0)
  return(x)
}

layout(matrix(1:4, nrow=2, byrow=TRUE))
simRnd2(20000, 0.5, 6.54, 0.96, 0.01)
simRnd2(20000, 0.5, 6.54, 0.96, 0.1)
simRnd2(20000, 0.5, 6.54, 0.96, 0.5)
simRnd2(20000, 0.5, 6.54, 0.96, 1)


# Profiling
library(chron)
prof<-read.table("syrupy_20131115102802.ps.log", header=TRUE, stringsAsFactors=FALSE)
prof$TIME<-chron(times=prof$TIME)
prof$TIME<-prof$TIME-prof$TIME[1]
layout(matrix(1:2, nrow=2))
plot(CPU~TIME, prof, type="l", main="Processor usage (%)", ylim=c(0,100))
plot(RSS~TIME, prof, type="l", main="Memory usage (kB)")
max(prof$RSS)
dev.print(pdf, file="Profiling.pdf", width=8, height=12)
dev.print(svg, file="Profiling.svg", width=8, height=12)

