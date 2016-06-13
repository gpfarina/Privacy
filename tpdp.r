require("Rmpfr")
hellinger <- function (a1,b1,a2,b2){sqrt(1.0 - beta((a1+a2)/2,(b1+b2)/2)/sqrt(beta(a1,b1)*beta(a2,b2)))}
l1sens <- function(a1,b1,a2,b2){abs(a1-a2)+abs(b1-b2)}
rangeBetas <- function(prior, n){
      r<-NULL
      for (a in 0:n){
      	  r <- append(r, list(c(prior[1]+a,prior[2]+n-a)))
      }
     return (r)
}
trues <-function(db){
      sum(db==TRUE)
}
betaPost<-function(prior, db){
        a<-trues(db)
        return(c(prior[1]+a, prior[2]+length(db)-a))	
    }
genDb<-function(theta,n){
   d<- sample(c(TRUE,FALSE),n, replace=T, prob=c(theta,1-theta))
}
computeDenum<-function(prior, eps, sens, db){
    r<-rangeBetas(prior, length(db))
  	k<-0.0
    nTrues<-trues(db)
	nFalses<-length(db)-nTrues
	for (a in r){
	    k=k+exp(-(eps/(2*sens))*hellinger(nTrues+prior[1],nFalses+prior[2],a[1],a[2]))
	}
	return(k)
}
distrExpMech <- function(prior, eps, sens, distance, db){
	 probDist<-c()
	 nTrues=trues(db)
	 nFalses=length(db)-nTrues
	 realBeta<-betaPost(prior, db)
	 denum<-computeDenum(prior, eps, sens, db)
     r<-rangeBetas(prior, length(db))
	 for (i in 1:length(r)){
	     probDist[i]<- exp(-(eps/(2*sens))*distance(nTrues+prior[1],nFalses+prior[2],r[[i]][1],r[[i]][2]))/denum
	 }
    return(probDist)
 }
dist2Str<-function(distance){
	if(identical(distance, hellinger)) return ("Hell")
	else return("")
}



repsidednoisymax<-function(prior, db, distance,eps,sens){
    r<-rangeBetas(prior, length(db))
    real<-betaPost(prior, db)
    t<-1:length(r)
    for (i in 1:length(r)){
        t[i]<- (1-(distance(r[[i]][1],r[[i]][2], real[1],real[2])))+rexp(1, r=eps/(sens))
    }
    return(r[which.max(t)])
}

plotAccuracies<-function(theta, size, n, eps, prior, distance){
  hist1<-1:n
  hist2<-1:n
  hist3<-1:n
  db<-genDb(theta,size)
  real<-betaPost(prior, db)
  r<-rangeBetas(prior, length(db))
  distr1<-distrExpMech(prior, eps, 2, distance, db)
  distr2<-distrExpMech(prior, eps, distance(prior[1]+length(db),prior[2], prior[1]+length(db)-1,prior[2]+1), hellinger, db)
  for(i in 1:n){
      v1<-sample(x=r, 1, replace = T, prob=distr1)
      hist1[i]<-distance(v1[[1]][1],v1[[1]][2],real[1],real[2])
      v2<-sample(x=r, 1, replace = T, prob=distr2)
      hist2[i]<-distance(v2[[1]][1],v2[[1]][2],real[1],real[2])
      t<-repsidednoisymax(prior,db,distance, eps, distance(prior[1]+length(db),prior[2], prior[1]+length(db)-1,prior[2]+1))
      hist3[i]<-distance(t[[1]][1],t[[1]][2],real[1],real[2])
      }
  par(mfrow = c(2,2))
  hist(hist1,   breaks=seq(0.0,1.0, by=0.1), main = "ExpMech on Beta(a,b)/L1 Sens",  xlab=dist2Str(distance))
  hist(hist2,   breaks=seq(0.0,1.0, by=0.1), main = "ExpMech on Beta(a,b)/Hell Sens",        xlab=dist2Str(distance))
  hist(hist3,   breaks=seq(0.0,1.0, by=0.1), main = "Rep1SidedNoisyMax on Beta(a,b)/Hell Sens",        xlab=dist2Str(distance))
  title(line=-1,outer=TRUE,sprintf("eps=%3.2f, n=%d, theta=%3.2f, sizeDb=%d, prior=(%d,%d)",eps, n, theta, size, prior[1],prior[2]))
}






