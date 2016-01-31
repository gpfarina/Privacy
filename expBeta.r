source("utils.r")
library(VGAM)
library(rgl)
library(lpSolve)
laplaceNoise<-function(sens, eps, prior, db){
    beta<-betaRet(prior, db)
    return( c(rlaplace(1, beta[1], sens/(eps/2)), rlaplace(1, beta[2], sens/(eps/2))))
}
laplaceNoiseLP <-function(sens, eps, prior, db, val){
    beta<-betaRet(prior, db)
    xnn <- rlaplace(1, beta[1], sens/(eps/2)) #real numbers
    ynn <- rlaplace(1, beta[2], sens/(eps/2)) #real numbers
    f.obj<-c(1,1,0,0)
    if(val){
        f.con<-matrix(c(
            		-1,0,-1,0,
                        -1,0,1,0,
                        0,-1,0,-1,
                        0,-1,0,+1,
                        0,0,1,0,
            		0,0,0,1,
                        0,0,1,1),nrow=7, byrow=TRUE)
        f.dir<-c("<=","<=","<=","<=",">=",">=","<=")
        f.rhs<-c(-xnn,xnn,-ynn,ynn,prior[1],prior[2],length(db)+prior[1]+prior[2])
    }
    else{
        f.con<-matrix(c(-1,0,-1,0,
                        	  -1,0,1,0,
                                  0,-1,0,-1,
                                  0,-1,0,+1,
                                  0,0,1,0,
                                  0,0,0,1,
                                  0,0,1,1,
                                  0,0,1,1),nrow=8, byrow=TRUE)
        f.dir<-c("<=","<=","<=","<=",">=",">=",">=","<=")
        f.rhs<-c(-xnn,xnn,-ynn,ynn,prior[1],prior[2],length(db)+prior[1]+prior[2], length(db)+prior[1]+prior[2])
    }
    d<-lp ("min", f.obj, f.con, f.dir, f.rhs,int.vec=c(3,4))
    return(c(d$solution[3],d$solution[4]))#pair of integers
}
laplaceNoisePostHellinger <- function (sens, eps, prior, db){
    beta<- betaRet(prior, db)
    r<-rangeBetas(prior, length(db))
    xnn <- rlaplace(1, beta[1], sens/(eps/2)) #real numbers
    ynn <- rlaplace(1, beta[2], sens/(eps/2)) #real numbers
    v<-NULL
    for (t in 1 :length(r)){
     	if(xnn<=0 || ynn <=0)
                v[t]<- 1
        else
            	v[t]<-hellinger(xnn,ynn,r[[t]][1],r[[t]][2])
    }
    return(r[which.min(v)])
}
lpvsHellPost<-function(sens, eps, prior, size, times, theta){
    db       <- genDb(theta, size)
    beta    <- betaRet(prior, db)
    histLp  <- NULL
    histH   <-  NULL
    for(i in 1: times){
        hell<-laplaceNoisePostHellinger(1, eps, prior,db)
        lp<-laplaceNoiseLP(1, eps, prior, db, FALSE)
        histLp[i]<-hellinger(lp[1],lp[2], beta[1], beta[2])
        histH[i]<-hellinger(hell[[1]][1], hell[[1]][2],beta[1], beta[2])
    }
    D1 <- hist(histLp, plot=FALSE,  breaks=seq(0.0,1.0, by=0.1))$counts
    D2 <- hist(histH,  plot=FALSE,  breaks=seq(0.0,1.0, by=0.1))$counts
    dat <- rbind(D1, D2)
    barplot(dat, col=c("red", "black"),beside=TRUE, space=c(0, 0.1), names.arg=seq(0,0.9, by=0.1),
            main=sprintf("L1-Norm Vs Hellinger minimization Post Processing, eps:%f, prior=Beta(%d,%d), size=%d, theta=%f", eps, prior[1], prior[2], size, theta)
            )
    legend("topright", 
       legend = c("L1-Norm", "Hellinger"), 
       fill = c("red", "black"))
    return(0)
}
computeDenum<-function(prior, eps, sens, distance, db, val){
    	if(val){
             r<-rangeBetasFull(prior, length(db))
         }
         else{
             r<-rangeBetas(prior, length(db))
         }
	k<-0.0
        nTrues<-trues(db)
	nFalses<-length(db)-nTrues
	for (a in r){
	    k=k+exp(-(eps/(2*sens))*distance(nTrues+prior[1],nFalses+prior[2],a[1],a[2]))
	}
	return(k)
}
distrExpMech <- function(prior, eps, sens, distance, db,val){
	 probDist<-c()
	 nTrues=trues(db)
	 nFalses=length(db)-nTrues
	 realBeta<-betaRet(prior, db)
	 denum<-computeDenum(prior, eps, sens, distance, db, val)
         if(val){
             r<-rangeBetasFull(prior, length(db))
         }
         else{
             r<-rangeBetas(prior, length(db))
         }
         #probDist=mpfrArray(x=c(0), precBits=bitPrec, dim=length(r)) #to have high precision for probabilities maybe we don't need it
	 for (i in 1:length(r)){
	     probDist[i]<- as.numeric(exp(-(eps/(2*sens))*distance(nTrues+prior[1],nFalses+prior[2],r[[i]][1],r[[i]][2]))/denum)
	 }
    return(probDist)
}
expMechSamples <- function(prior, eps, sens, distance, db, times,val) {
        if(val){
             r<-rangeBetasFull(prior, length(db))
         }
         else{
             r<-rangeBetas(prior, length(db))
         }
        sample(x=r, times, replace = T, prob = distrExpMech(prior, eps, sens, distance, db, val))
}
expMech<-function(prior, eps, sens, distance, db,val){
   expMechSamples(prior, eps, sens, distance, db, 1,val)
}

dist2Str<-function(distance){
if (identical(distance, totalVariation)) return ("TV")
	else if(identical(distance, hellinger)) return ("Hell")
		else if(identical(distance, kl)) return ("KL")
			else if(identical(distance, klSymmetric))return("KLSym")
}
plotAccuracies<-function(n, teps, theta, size, prior, distance, sens, val=FALSE, dist3d=FALSE){
  hist1<-1:n
  hist2<-1:n
  db<-genDb(theta,size)
  real<-betaRet(prior, db)
  if(val){
          r<-rangeBetasFull(prior, length(db))
   }
  else{
          r<-rangeBetas(prior, length(db))
   }
  distr<-distrExpMech(prior, teps*2, sens, distance, db, val)
  for(i in 1:n){
      v<-laplaceNoiseLP(1, teps, prior, db,val)
      hist1[i]<-as.numeric(distance(v[1],v[2],real[1],real[2]))
      v1<-sample(x=r, 1, replace = T, prob=distr)
      hist2[i]<-as.numeric(distance(v1[[1]][1],v1[[1]][2],real[1],real[2])) #as.numeric to convert to standard numbers....because so far I don't know how to sample
      #from distributions where the probabilities are MPFR numbers....
      }
  par(mfrow = c(2,1))
  hist(hist1,  breaks=seq(0.0,1.0, by=0.1), main = "Laplace Noise on a and b", xlab=dist2Str(distance))
  hist(hist2,   breaks=seq(0.0,1.0, by=0.1), main = "ExpMech on Beta(a,b)",xlab=dist2Str(distance))
  title(line=-1,outer=TRUE,sprintf("eps=%f, samples=%d, theta=%f, size db=%d, prior=(%d,%d),dist=hellinger",teps, n, theta, size, prior[1],prior[2]))
  if(dist3d){
      fun<-function(x,y){distance(real[1],real[2], x, y)}
      f<- Vectorize(fun)
      persp3d(f, xlim=c(prior[1],length(db)+prior[1]), ylim=c(prior[2],length(db)+prior[2]), xlab="a", ylab="b")
  }
  return(distr)
}
