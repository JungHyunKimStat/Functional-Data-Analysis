subway <- read.csv("subway1.csv",header=T)
matplot(subway,type='l',xlab="시간",ylab="인원")
library(fda)
subway.basis <- create.bspline.basis(c(1,21),nbasis=23,norder=4,breaks=(1:21))

subway<- as.matrix(subway)

loglam <- -5:0
nlam = length(loglam)
gcvsave = rep(NA,nlam)
for (i in 1:nlam) {
  #i=1
  lambda   = 10^loglam[i]
  fdParobj = fdPar(subway.basis, int2Lfd(2), lambda)
  smoothlist = smooth.basis(1:21, subway,fdParobj,dfscale=1.2)
  gcvsave[i] = sum(smoothlist$gcv)
}
windows()
plot(x=loglam,y=gcvsave,type="l",xlab="Log(lam)",ylab="gcv")
lamb <- 10^(loglam[which.min(gcvsave)])

subwayd <- smooth.basisPar(argvals=1:21,
                           y=subway, fdobj=subway.basis,
                           Lfdobj=int2Lfd(2), lambda=lamb)
subway.fd = subwayd$fd
# open a new graph window in R in a Windows laptop.
# run quartz() if you are using a Mac computer
windows()
plot(subway.fd,xlab="시간",ylab="인원")

#  a principal components analysis
nharm = 3
#nonparametric fpca
subway.pcalist.bs = pca.fd(subway.fd, nharm=3)
#variance proportation explained
subway.pcalist.bs$varprop

sum(subway.pcalist.bs$varprop)
# Varimax rotation 해야??

coef.bs <- subway.pcalist.bs$harmonics$coefs
timepts <- seq(1, 21, length.out=200)
randfd1.bs = fd(coef.bs[,1],subway.basis)
pc1.bs = eval.fd(timepts,randfd1.bs)
randfd2.bs = fd(coef.bs[,2],subway.basis)
pc2.bs = eval.fd(timepts,randfd2.bs)
randfd3.bs = fd(coef.bs[,3],subway.basis)
pc3.bs = eval.fd(timepts,randfd3.bs)

windows()
par(mfrow=c(1,3),lwd=3,cex.lab=1.2,oma=c(0,1,0,0))
plot(x=timepts,y=pc1.bs,type="l",lty=1,xlab="시간",ylab="FPC 1")
plot(x=timepts,y=pc2.bs,type="l",lty=1,xlab="시간",ylab="FPC 2")
plot(x=timepts,y=pc3.bs,type="l",lty=1,xlab="시간",ylab="FPC 3")
# 해석하기가 어려워...