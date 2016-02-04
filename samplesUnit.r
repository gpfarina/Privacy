library(hitandrun)
source("expDir.r")
#
#samples from the 2-simplex a pair (x, y) and then
#using that pair generates a database of sizedb elements
#using the prior computes the posterior given the observations in sizedb
#now for m times adds noise to the posterior using LP and LP with linear post processing
#then computes the distances z and z' from the real posterior and averages it
#then draws a poitn in (x,y,z) and (x,y,z')

sampleFromSimplex<-function(m,nPoints, sizedb,prior, eps){
    avgs<-rep(0, nPoints)
    avgs2<-rep(0, nPoints)
    cat<-length(prior)
    prob<-matrix(rep(0, nPoints*cat), ncol=cat)
    for(i in 1:nPoints){
        prob[i,]<- simplex.sample(cat, 1, sort=FALSE)$samples[1,]
        db<-genDbD(sizedb, cat, prob[i,])
        postReal <- computePost(prior, db)
        for(l in 1:m){
            noisy2 <-laplaceNoise(2, eps, prior, db) 
            noisy<-laplaceNoiseLPD(2, eps, prior, db, noisy2)
            avgs[i]<- avgs[i]+as.numeric((hellingerD(postReal, noisy)))/m
            if(sum(noisy2<0)<=0){
                avgs2[i]<- avgs2[i]+as.numeric((hellingerD(postReal, noisy2)))/m
            }
        }
    }
    data <- rbind(cbind(prob,avgs),cbind(prob,avgs2))
    plot3d(x=data[,1], y=data[,2], z=data[,3], xlab="prob", ylab="prob", zlab="avg hell-distance", col=c("red"))
      legend3d("topright", 
      legend = c("Laplace Noise LP",
           		 "Laplace Noise No Post"), 
       fill =        c("red", "black"),
           )
  #  bgplot3d({
  #plot.new()
  #title(main = sprintf("Lap-Noise, w/post processing LP, eps:%3.2f,sizedb:%d, prior:%s, nPoints:%d ", eps, sizedb, sprintf("(%s)", paste(prior, collapse=" ")), nPoints), line = 3.0)
#})
    return(data)
}


#nPoints:number of points sampled from the 2-simples.
#m: number of times we add noise to the posterior
#sizedb: number of rows in the db
#prior: prior beta distribution
#eps: epsilon parameter
#distX: 0<=distX<=1, this is the mximum distance from the real posterior we as limit
#v: if true LPvsHell if false LPvsNoPost
#This function samples nPoints form the simplex and then for every of this points
#generates m times a db with length(prior) categories and of size sizedb.
#after this computes the real posterior given the database and
#then adds noise to the posterior using laplace.
#It then performs postprocessing (LP and Hell) if the
#posterior noisy distributions are closer than distX to the real increases teh relative counters.

#v: if TRUE compares  LP(black) vs Hell(red)
# if FALSE compares   LP(black) vs no-Post(red)


exp<-function(nPoints, m, sizedb, prior, eps, distX, v){
    cat<-length(prior)
    prob<-matrix(rep(0, nPoints*cat), ncol=cat)
    less1<-rep(0, nPoints)
    less2<-rep(0, nPoints)
    less3<-rep(0, nPoints)
    for(i in 1:nPoints){
        prob[i,]<- simplex.sample(cat, 1, sort=FALSE)$samples[1,]
        db<-genDbD(sizedb, cat, prob[i,])
        postReal <- computePost(prior, db)
         for(l in 1:m){
            noisy  <-laplaceNoise(2, eps, prior, db)
            noisy1<-laplaceNoiseLPD(2, eps, prior, db, noisy)
            if(as.numeric(hellingerD(noisy1, postReal ))<=distX){less1[i]<-less1[i]+1}
             if(v){
                noisy2<-laplaceNoisePostHellingerD(2, eps, prior, db, noisy)
                if(as.numeric(hellingerD(noisy2, postReal ))<=distX){less3[i]<-less3[i]+1}
            }
             else{
              	if(!(sum(noisy<0)>0)){
                    if(as.numeric(hellingerD(noisy  , postReal ))<=distX){
                        less2[i]<-less2[i]+1
                     }
                }
                else
                    if(sum(noisy<0)>0){show(sprintf("This ended up with negative parameters: %s, ", sprintf("(%s)", paste(noisy, collapse=" "))))}
            }
        }
        show(sprintf("---> %d-th loop", i))
    }
    less1<-(less1/m)*100
    less2<-(less2/m)*100
    less3<-(less3/m)*100
    if(v){
        data<-rbind(cbind(prob, less1),  cbind(prob, less3))
        str<-"Laplace Noise Hell min"
    }
    else{
        data<-rbind(cbind(prob, less1), cbind(prob, less2))
        str<-"Laplace Noise  no Post "
    }
    show(data)
    r<-plot3d(x=data[,1], y=data[,2], z=data[,3], xlab="prob", ylab="prob", zlab=sprintf("#<=%.1f",distX), col=c(rep("black",nPoints), rep("red",nPoints)))
    legend3d("topright", 
    legend = c("Laplace Noise LP",
           		 str), 
       fill =        c("black", "red"),
           )

    return(r)
}


exp2<-function(nPoints, m, sizedb, prior, eps, distX, v){
    cat<-length(prior)
    prob<-matrix(rep(0, nPoints*cat), ncol=cat)
    lp<-rep(0, nPoints)
    nopost<-rep(0, nPoints)
    hell<-rep(0, nPoints)
    for(i in 1:nPoints){
        prob[i,]<- simplex.sample(cat, 1, sort=FALSE)$samples[1,]
        db<-genDbD(sizedb, cat, prob[i,])
        postReal <- computePost(prior, db)
         for(l in 1:m){
            noisy  <-laplaceNoise(2, eps, prior, db)
            noisy1<-laplaceNoiseLPD(2, eps, prior, db, noisy)
            if(is.nan(hellingerD(noisy1, postReal ))){
                show(noisy1)
                show(postReal)
                lp[i]<-lp[i]+1
            }
            else
                if(hellingerD(noisy1, postReal )<=distX){
                    lp[i]<-lp[i]+1
            }
            if(v){
                noisy2<-laplaceNoisePostHellingerD(2, eps, prior, db, noisy)
                if(hellingerD(noisy2, postReal )<=distX){hell[i]<-hell[i]+1}
            }
            else{
              	if(!(sum(noisy<0)>0)){
                        if(hellingerD(noisy  , postReal )<=distX){
                            nopost[i]<-nopost[i]+1
                     }
                }
                else
                    if(sum(noisy<0)>0){show(sprintf("This ended up with negative parameters: %s, ", sprintf("(%s)", paste(noisy, collapse=" "))))}
            }
        }
        show(sprintf("---> %d-th loop", i))
    }
    if(v){plot(prob[,1], hell-lp)}
    else {plot(prob[,1], nopost-lp)}
}


