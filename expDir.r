rm(list=ls())
library(SimplicialCubature)
library(stats)
library(gtools)
source("expBeta.r")
## The Simplicial Package only works for Dirichilet distribution with dimension >=2 
#For the beta distribution: we use the previously written code

#x is the first k-1 components, the last component is determined by 1-sum(x)
#a is k dimensional parameter.
#e.g. k=3, x = c(0.2,0.3) and a = c(1,2,3)
dirchPdf<-function(x,a){
  if(length(c(a)) - length(c(x)) != 1){stop("Error; dimension of a and x dont match")}
  if(sum(x)>1 || sum(x)==1){ return(0)}
  if(sum(x)<1){ 
    new.x = c(x,1-sum(x));
    ans = ddirichlet(new.x, a );
    return(ans)
  }
}


#TotalVariationD: computes the total variation distance between two dirichlet distributions (with the same number of arguments)
#a: vector of real numbers representing the parameters of the first dirichlet distribution
#b: vector of real numbers representing the parameters of the second dirichlet distribution

totalVariationD<-function(a,b){
  d<- -1
  if (length(a)!=length(b)){
          stop("Error; dirichlets dimensions are not matching")
   }
  else{
      if(length(a)==2){return(totalVariation(a[1],a[2],b[1],b[2]))}
      else if(length(a)>2){
          f<-function(x){return(abs(dirchPdf(x,a)-dirchPdf(x,b)))}
          d<-adaptIntegrateSimplex(f, CanonicalSimplex(length(a)-1))}
  }
  return(d$integral/2)
}
#betam: is the generalized beta function
#a: vector of parameters 
betam<-function(a){
#   return(Reduce(function(x,y){x*y},sapply(a,gamma),1 )/gamma(Reduce(function(x,y){x+y},a,0)))
   return(Reduce(function(x,y){x*y},new("mpfr",unlist(sapply(a,function(x){gamma(as(x,"mpfr"))}))), as(1,"mpfr") )/gamma(Reduce(function(x,y){x+y},a,as(0, "mpfr"))))
}
hellingerD <- function (a,b){
      d<-NULL
       if(identical(a, b)) {
           show("A distance of 0 has been reached")
           d<-as.numeric(as(0, "mpfr"))
           }
       else{
                d<-as.numeric(sqrt(1.0 - betam((a+b)/2)/sqrt(betam(a)*betam(b))))
   	}
       return (d)
}
#this function generates all the possible parameters of a dirichlet with k categories and n observations
#n : positive integer, is the number of observations
#k:  integer greater or equal than 2, is the number of categories
f<-function(n,k){
    if(k==1){return(matrix(n,1,1))}
    else{
        out<-matrix(0,0,k)
        for(i in 0: n){
            rec<-f(n-i,k-1)
            for(j in 1:nrow(rec)){
                out<-rbind(out,c(i,rec[j,]))
            }
        }
        return(out)
    }
}
#given a prior this gives you all the possible dirichlet distributions given the size of the database
#prior: vector of real non negative numbers representing the prior dirichlet distribution
#n: size of the database
rangeDir<-function(prior, n){
    return(t(t(f(n,length(prior)))+prior))
}
#we learn the posterior , incrementing the counts, and compute the distribution of the exponential mechanism
#priors: vector of real non negative numbers representing  the prior dirichlet distribution
#db: database with values coming from "length(priors)" possible categories
computePost<-function(priors, db){
    prior<-priors
    for(i in 0:length(prior)-1){
        prior[i+1]<-prior[i+1]+length(db[db==i])
    }
    return(prior)
}

createDistr<-function(priors, db, eps,distance,sens){
    M<-rangeDir(priors, length(db))
    posts<-computePost(priors, db)
    dists<-apply(M, 1,function(x){distance(x,posts)})
    X<-new("mpfr",unlist(sapply(dists, function(x){exp((-eps*x)/(2*sens))})))
    return(X/sum(X))
}

#sampling function
expMechDir<-function(priors, db, eps, dist, sens){
    distr<-createDistr(priors,db, eps, dist, sens)
    return(sample(apply(rangeDir(priors, length(db)),1,function(x){return(list(x))}) ,1,prob=distr))
}
#adds noise to the parameters
laplaceNoiseD<-function(sens, eps, prior, db){
    beta<- computePost(prior, db)
    return(sapply(beta, function(b){rlaplace(1, b, (sens/eps))}))
}
#make constraints for the LP solver
makeConstraints <- function(n){
    constr<- NULL
    for(i in 1: n ){
        cstr1<- rep(0, 2 * n)
        cstr2<- rep(0, 2 * n)
        cstr1[i]<- -1
        cstr2[i]<- -1
        cstr1[n+i] <- -1
        cstr2[n+i] <- +1
        constr<- rbind(constr, cstr1, cstr2)
    }
    for(i in 1: n ){
        cstr<- rep(0,  2 * n)
        cstr[n+i] <- 1
        constr<-rbind(constr, cstr)
    }
    f1<-c(rep(0, n), rep(1, n))
    constr<-rbind(constr, f1, f1)
    return(matrix(constr, nrow=3*n+2))
}
#adds noise but returns consistent values
laplaceNoiseLPD <-function(sens, eps, prior, db, value, int){
    z<-rep(0, 2*length(value))
    for(i in 1:(2*length(value))){
            z[i]<-  ((-1)**(i))*(value[ceiling(i/2)])
    }
    tt.obj  <-  c(rep(1, length(value)), rep(0, length(value)))
    tt.con <-  makeConstraints(length(prior))
    tt.dir  <-   c( rep("<=", 2*length(prior)), rep(">=", length(prior)) , ">=", "<="  )
    tt.rhs  <-  c(z, prior, length(db)+sum(prior), length(db)+sum(prior))
    if(int){
        d<-lp ("min", tt.obj, tt.con, tt.dir, tt.rhs,int.vec=seq(length(prior)+1, 2*length(prior))) 
    }
    else{
        d<-lp ("min", tt.obj, tt.con, tt.dir, tt.rhs)
    }
    return(d$solution[c(seq(length(prior)+1, 2*length(prior)))])
}
#generates a db of n elements coming from a size categories with proability theta
genDbD<- function(n, size, theta){
    return(sample(seq(0, size-1), n, prob=theta, replace=TRUE) )
}
#returns a beta which minimizes the hellinger distance from the one with noisy oarameters
#the returned one has parameters which are consistent
laplaceNoisePostHellingerD <- function (sens, eps, prior, db, noisy){
    beta<- computePost(prior, db)
    r<-rangeDir(prior, length(db))
    v<-rep(0, nrow(r))
    for (t in 1 :nrow(r)){
     	if(sum(noisy<0)>0)
                v[t]<- 1
        else
            	v[t]<- hellingerD(noisy,r[t,])
    }
    v<-as(unlist(v),"mpfr")
    return(r[which.min(v),])
}

lpvsHellPostD<-function(eps, prior, n, size, times, theta){
    db       <- genDbD(n, size, theta)
    real    <- computePost(prior, db)
    histLp  <- NULL
    histH   <-  NULL
    r<-rangeDir(prior, length(db))
    for(i in 1: times){
        hell<-laplaceNoisePostHellingerD(2, eps, prior,db)
        lp<-laplaceNoiseLPD(2, eps, prior, db)
        histLp[i]<-hellingerD(lp, real)
        histH[i]<-hellingerD(hell, real)
    }
    D1 <- hist(histLp, plot=FALSE,  breaks=seq(0.0,1.0, by=0.1))$counts
    D2 <- hist(histH,  plot=FALSE,  breaks=seq(0.0,1.0, by=0.1))$counts
    dat <- rbind(D1, D2)
    barplot(dat, col=c("red", "black"),beside=TRUE, space=c(0, 0.1), names.arg=seq(0,0.9, by=0.1),
            main=sprintf("L1-Norm Vs Hellinger minimization Post Processing, eps:%.3f, samples:%d prior=%s, size=%d, theta=%s", eps, times,
               sprintf("(%s)", paste(prior, collapse=" ")),
               size,
               sprintf("(%s)", paste(theta, collapse=" "))
            ))
    legend("topright", 
       legend = c("L1-Norm", "Hellinger"), 
       fill = c("red", "black"))
    return(0)
}
#convert to strings the name of the function
dist2Str<-function(distance){
if (identical(distance, totalVariationD)) return ("TV")
	else if(identical(distance, hellingerD)) return ("Hell")
		#else if(identical(distance, klD)) return ("KL")
		#	else if(identical(distance, klSymmetricD))return("KLSym")
}



#n: number of times we will add noise, so in our histogram will consist of n points in 10 bins (0.0, by 0.1, 1.0)
#teps: epsilon parameter
#theta: vector of probabilities used to generate the database (it has to come from the len-simplex
#len: number of categorie
#size: size of the db
#prior: prior distribution: has to have len parameters
#distance: the distance to be used, usually hellingerD
#sensH: sensitivity used for the exponential mechanism, or better the sensitivity of the statistic metric
#sensL1: sensitivity of the l1 norm over the parameters: it should be 2
plotAccuraciesD<-function(n, teps, theta, len, size, prior, distance, sensH, sensL1){
  hist0<-1:n
  hist1<-1:n
  hist2<-1:n
  hist3<-1:n
  data<-matrix(0, n, 2)
  db<-genDbD(size, len, theta)
  real<-computePost(prior, db)
  r<-rangeDir(prior, len)
#  distr<-createDistr(prior, db, teps, distance, sensH) # we are not using expmech so we don't create it for now
  failed<-0
  for(i in 1:n){
      data[i,]<-laplaceNoiseD(sensL1, teps, prior, db)
  }
  for(i in 1:n){
     # show(i)
      	if(sum((data[i,])<0)>0){
                hist0[i]<-1
                hist2[i]<-1
                hist3[i]<-1
                failed<-failed+1
            }
        else{
            hist0[i]<-as.numeric(distance(data[i,],real))
        #    v3<-laplaceNoisePostHellingerD(sensL1, teps, prior, db, data[i,])
       #     hist3[i]<-as.numeric(distance(v3, real))
            #distr<-as.numeric(distr)
         #   v2<- sample(x=as.list(data.frame(t(r))), 1, replace = T, prob=as.vector(t(distr))) #expmech
        #    hist2[i]<-as.numeric(distance(v2[[1]], real))
        }
        v1<-laplaceNoiseLPD(sensL1, teps, prior, db, data[i,], TRUE)
        hist1[i]<-as.numeric(distance(v1,real))
      }
  par(mfrow = c(2,2))
  histPercent(hist0,  n, failed, main =   "Laplace Noise no Post", xlab=dist2Str(distance))
  histPercent(hist1,  n, 0, main =   "Laplace Noise Mult LPminL1Norm", xlab=dist2Str(distance))
  histPercent(hist2,  n, failed, main =   "ExpMech on Mult no Post"   , xlab=dist2Str(distance))
  histPercent(hist3,  n, failed, main =   "Laplace Noise on Mult LPminHell"   , xlab=dist2Str(distance))
  string<-sprintf("(%s)", paste(theta, collapse=" "))
  string3<-sprintf("(%s)", paste(prior, collapse=" "))
  string2<-sprintf("eps=%.3f, samples=%d, theta=%s, size db=%d, prior=%s, dist=hellinger, sens L1=%.2f, sens Exp=%.5f",teps, n, string, len,
                    string3, sensL1, sensH)
  mtext(outer=TRUE, string2 , line=-1.2)
  return(failed)
}
histPercent <- function(x, n,failed,...) {
   H <- hist(x, plot = FALSE, breaks=seq(0.0,1.0, by=0.1))
   H$density <- with(H, 100 * density*diff(breaks)[1])
   labs <- paste(round(H$density), "%", sep="")
   plot(H, freq=FALSE, labels = labs, col="gray", ylim=c(0, 1.10*max(H$density)),...)
   usr <- par( "usr" )
   text( usr[ 2 ], usr[ 4 ], sprintf("failed: %.1f %%", (failed/n)*100), adj = c( 1, 1 ), col = "red" )
}

#lpsolve come svegliere fra varie soluzioni
#tutte le soluzioni su lpMin e poi pescare...perché ho molti minimi
#stampare le beta nei grafici
#n: number of times we will add noise, so in our histogram will consist of n points in 10 bins (0.0, by 0.1, 1.0)
#teps: epsilon parameter
#theta: vector of probabilities used to generate the database (it has to come from the len-simplex
#len: number of categorie
#size: size of the db
#prior: prior distribution: has to have len parameters
plotAccuraciesNoVsLp<-function(n, teps, theta, len, prior){
  hist0<-1:n
  hist1<-1:n
  hist2<-1:n
  data<-matrix(0, n, 2)
  db<-genDbD(len, length(prior), theta)
  real<-computePost(prior, db)
  failed<-0
  for(i in 1:n){
        data[i,]<-laplaceNoiseD(2, teps, prior, db)
        v<-laplaceNoiseLPD(2, teps, prior, db, data[i,], FALSE)
        s<-laplaceNoiseLPD(2, teps, prior, db, data[i,], TRUE)
        hist0[i]<-abs(data[i,1]-real[1])+abs(data[i,2]-real[2])
        hist1[i]<-abs(v[1]-real[1])+abs(v[2]-real[2])
        hist2[i]<-abs(s[1]-real[1])+abs(s[2]-real[2])
      #  show(i)
      }
  par(mfrow = c(2,2))
  histPercent1(hist0,  n, failed, main =   "Laplace Noise no Post", xlab=dist2Str(distance))
  histPercent1(hist1,  n, 0, main =   "Laplace Noise Mult LPminL1Norm No Int Constraints", xlab=dist2Str(distance))
  histPercent1(hist2,  n, 0, main =   "Laplace Noise Mult LPminL1Norm Int Constraints", xlab=dist2Str(distance))
  string1<-sprintf("(%s)", paste(theta, collapse=" "))
  string2<-sprintf("(%s)", paste(prior, collapse=" "))
  string<- sprintf("eps=%.3f, samples=%d, theta=%s, size db=%d, prior=%s",teps, n, string1, len,
                    string2)
  mtext(outer=TRUE, string , line=-1.2)
  return(0)
}
histPercent1 <- function(x, n,failed,...) {
   H <- hist(x, plot=FALSE, breaks=seq(0,max(x)+5, by=5))
   H$density <- with(H, 100 * density*diff(breaks)[1])
   labs <- paste(H$density, "%", sep="")
   plot(H, freq=FALSE, labels = labs, col="gray", ylim=c(0, 1.10*max(H$density)),...)
}

compAverage<-function(simplexSamples, n, teps, len, prior){
    diff0<-matrix(0, simplexSamples+1, 2)
    diff1<-matrix(0, simplexSamples+1, 2)
    t=1
    for(i in seq(0, 1, by = (1/simplexSamples))){
          theta<-c(i,1-i)
          db<-genDbD(len, length(prior), theta)
          real<-computePost(prior, db)
          data<-matrix(0, n, 2)
          hist0<-1:n
          hist1<-1:n
          hist2<-1:n
          for(j in 1:n){
                data[j,]<-laplaceNoiseD(2, teps, prior, db)
                v<-laplaceNoiseLPD(2, teps, prior, db, data[j,], FALSE)
                s<-laplaceNoiseLPD(2, teps, prior, db, data[j,], TRUE)
                  hist0[j]<-abs(data[j,1]-real[1])+abs(data[j,2]-real[2])
                  hist1[j]<-abs(v[1]-real[1])+abs(v[2]-real[2])
                  hist2[j]<-abs(s[1]-real[1])+abs(s[2]-real[2])
                  #hist0[j]<-sqrt((data[j,1]-real[1])**2+(data[j,2]-real[2])**2)
                  #hist1[j]<-sqrt((v[1]-real[1])**2+(v[2]-real[2])**2)
                  #hist2[j]<-sqrt((s[1]-real[1])**2+(s[2]-real[2])**2)
           }
          H0<-hist(hist0, plot=FALSE, breaks=seq(0,max(hist0)+5, by=5))
          H1<-hist(hist1, plot=FALSE, breaks=seq(0,max(hist1)+5, by=5))
          H2<-hist(hist2, plot=FALSE, breaks=seq(0,max(hist2)+5, by=5))
          par(mfrow = c(2,2))
          diff0[t,]<-c(H1$counts[1], H0$counts[1])
          diff1[t,]<-c(H2$counts[1], H1$counts[1])
          t<-t+1
          histPercent1(hist0,  n, 0, main =   "Laplace Noise no Post", xlab=dist2Str(distance))
          histPercent1(hist1,  n, 0, main =   "Laplace Noise Mult LPminL1Norm No Int Constraints", xlab=dist2Str(distance))
          histPercent1(hist2,  n, 0, main =   "Laplace Noise Mult LPminL1Norm Int Constraints", xlab=dist2Str(distance))
          string1<-sprintf("(%s)", paste(theta, collapse=" "))
          string2<-sprintf("(%s)", paste(prior, collapse=" "))
          string<- sprintf("eps=%.3f, samples=%d, theta=%s, size db=%d, prior=%s",teps, n, string1, len,
                           string2)
          mtext(outer=TRUE, string , line=-1.2)
    }
    return(cbind(diff0, as.logical(diff0[,1]>=diff0[,2]), diff1, as.logical(diff1[,1]>=diff1[,2])))
}
#TODO:
#geondi: how to generalize to a stastical distance and how to sample
comp1<-function(prior, eps, db, times,d){
    real<-computePost(prior,db)
    v1<-0
    v2<-0
    v3<-0
    for(i in 1:times){

        	nL<-laplaceNoiseD(2, eps, prior, db)
        	nLPNI<-laplaceNoiseLPD(2, eps, prior, db, nL, FALSE)
        	nLPI<-laplaceNoiseLPD(2, eps, prior, db, nL, TRUE)
                if(d){
                    if(!(is.nan(hellingerD(real, nL)))){
                        v1<-v1+(hellingerD(real, nL))/times
                    }
                    else{
                        v1<-v1+1/times
                    }
                    v2<-v2+(hellingerD(real, nLPNI))/times
                    v3<-v3+(hellingerD(real, nLPI))/times
            }
               else{
                v1<-v1+(((nL[1]-real[1])**2)+(nL[2]-real[2])**2)/times
                v2<-v2+(((nLPNI[1]-real[1])**2)+(nLPNI[2]-real[2])**2)/times
                v3<-v3+(((nLPI[1]-real[1])**2)+(nLPI[2]-real[2])**2)/times
            }
            }
    return(c(v1,v2,v3))
}

MSE0<-function(prior, db, intervals, times,d){
    i<-1
    v<-matrix(0, intervals, 3)
#    for(eps in seq(0.2/intervals, 0.2, by = (0.2/intervals))){
        for(eps in seq(1/intervals, 1, by = (1/intervals))){
        v[i,]<-comp1(prior, eps, db, times,d)
        i<-i+1
    }
    return(v)
}


intervals<-50
d<-TRUE
times<- 50
sizedb<-100
theta<-c(0.5,0.5)
res<-MSE0(c(10,10), genDbD(sizedb, 2, theta), intervals, times, d)
eps <- seq(1/intervals, 1, by = (1/intervals))
plot(eps, res[,3],col='green')
lines(eps, res[,3],col='green')
par(new=T)
lines(eps, res[,2],col='red')
par(new=T)
lines(eps, res[,1],col='black')
legend('topright',c("POST-INT-POS-SUM","POST-POS-SUM","NO-POST") , 
   lty=1, col=c('green', 'red', 'black'), bty='n', cex=.75)
