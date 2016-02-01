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
       d<-sqrt(1.0 - betam((a+b)/2)/sqrt(betam(a)*betam(b)))
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
laplaceNoiseD<-function(sens, eps, prior, db){
    beta<- computePost(prior, db)
    return(sapply(beta, function(b){rlaplace(1, b, (sens/eps))}))
}
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
laplaceNoiseLPD <-function(sens, eps, prior, db, value){
    z<-rep(0, 2*length(value))
    for(i in 1:(2*length(value))){
            z[i]<-  ((-1)**(i))*(value[ceiling(i/2)])
    }
    tt.obj  <-  c(rep(1, length(value)), rep(0, length(value)))
    tt.con <-  makeConstraints(length(prior))
    tt.dir  <-   c( rep("<=", 2*length(prior)), rep(">=", length(prior)) , ">=", "<="  )
    tt.rhs  <-  c(z, prior, length(db)+sum(prior), length(db)+sum(prior))
 #   d<-lp ("min", tt.obj, tt.con, tt.dir, tt.rhs,int.vec=seq(length(prior)+1, 2*length(prior))) #removed the integer constraints 
    d<-lp ("min", tt.obj, tt.con, tt.dir, tt.rhs)
    return(d$solution[c(seq(length(prior)+1, 2*length(prior)))])
}

genDbD<- function(n, size, theta){
    return(sample(seq(0, size-1), n, prob=theta, replace=TRUE) )
}
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
dist2Str<-function(distance){
if (identical(distance, totalVariationD)) return ("TV")
	else if(identical(distance, hellingerD)) return ("Hell")
		else if(identical(distance, klD)) return ("KL")
			else if(identical(distance, klSymmetricD))return("KLSym")
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
  db<-genDbD(len, size, theta)
  real<-computePost(prior, db)
  r<-rangeDir(prior, len)
#  distr<-createDistr(prior, db, teps, distance, sensH) # we are not using expmech so we don't create it for now
  failed<-0
  for(i in 1:n){
      data[i,]<-laplaceNoiseD(sensL1, teps, prior, db)
  }
  for(i in 1:n){
      show(i)
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
        v1<-laplaceNoiseLPD(sensL1, teps, prior, db, data[i,])
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
