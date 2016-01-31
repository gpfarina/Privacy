aP=bP=0
sens=0
eps=0
distance=NULL
db<-NULL
bitPrec<-1000
library(Rmpfr)
betaPrec<-function(a,b, prec){
    z<-mpfr(a, round(prec/2))
    w<-mpfr(b, round(prec/2))
    d<-beta(z,w)
    return(d)
}
hellinger <- function (a1,b1,a2,b2){
  #    d<-sqrt(1.0 - betaPrec((a1+a2)/2,(b1+b2)/2, bitPrec)/sqrt(betaPrec(a1,b1,bitPrec)*betaPrec(a2,b2,bitPrec)))
         d<-sqrt(1.0 - beta((a1+a2)/2,(b1+b2)/2)/sqrt(beta(a1,b1)*beta(a2,b2)))
         return (d)
}
kl <- function(a1,b1,a2,b2){
    d<-log(beta(a2,b2)/beta(a1,b1)) - (a2-a1)*digamma(a1) - (b2-b1)*digamma(b1) + (a2-a1+b2-b1)*digamma(a1+b1)
    return(d)
}
klSymmetric <- function(a1,b1,a2,b2){
    d<-kl(a1,b1,a2,b2)+kl(a1,b2,a1,b1)
    return (d)
}
totalVariation <-function(a1,b1,a2,b2){
    integ <- function(x) {
#        first<-mpfr(x^(a1-1), round(bitPrec/2))
 #       second<-mpfr((1-x)^(b1-1),round(bitPrec/2))
  #      third<-mpfr(x^(a2-1), round(bitPrec/2))
   #     fourth<-mpfr((1-x)^(b2-1),round(bitPrec/2))
    #    fifth<-mpfr((first*second)/betaPrec(a1,b1,round(bitPrec/2)), round(bitPrec/2))
     #   sixth<-mpfr((third*fourth)/betaPrec(a2,b2,round(bitPrec/2)),round(bitPrec/2))
      #  final<-mpfr(fifth-sixth,round(bitPrec/2))
        #return(as.numeric(abs(final)))
        return(abs(
        ((x^(a1-1)*(1-x)^(b1-1))/beta(a1,b1))
        -
        ((x^(a2-1)*(1-x)^(b2-1))/beta(a2,b2))))
    }
    d<-integrate(integ, lower = 0, upper = 1)
    return ((d[1]$value)/2)
}
rangeBetas <- function(prior, n){
      r<-NULL
      for (a in 0:n){
      	  r <- append(r, list(c(prior[1]+a,prior[2]+n-a)))
      }
 return (r)
}
rangeBetasFull <- function(prior, n){
      r<-NULL
      for(i in 0:n){
           r<-append(r, rangeBetas(prior, i))
      }
      return(r)
}
trues <-function(db){
      sum(db==TRUE)
}
betaRet<-function(prior, db){
        a<-trues(db)
	return(c(prior[1]+a, prior[2]+length(db)-a))	
}
setEps <-function(ppeps=0.1){eps<<-ppeps}
setPrior <-function(pa=1,pb=1){ aP<<- pa
    bP<<-pb}
setDist<-function(ppdist){distance<<-ppdist}
setSens<-function(ppsens){sens<<-ppsens}
setDb<-function(ppdb){db<<-ppdb}
genDb<-function(theta,n){
   d<- sample(c(TRUE,FALSE),n, replace=T, prob=c(theta,1-theta))
}
setParameters<- function(paP1=1, pbP1=1, peps1=0.5, pdist1=hellinger, psens1=hellinger(2,1,1,2),
                         pdb1=c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,TRUE)){
    setEps(peps1)
    setPrior(paP1, pbP1)
    setDist(pdist1)
    setSens(psens1)
    setDb(pdb1)
}
