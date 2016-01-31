source("/home/gpietro/Dropbox/PrivacyWithStatisticalMetrics/Code/utils.r")
test<-function(a,b,n){
    print(sprintf("Prior: a:%d, b:%d, n:%d",a,b, n))
    for (i in 0:n){
        print(sprintf("a:%d, b:%d, h:%f, tv:%f, kl:%f, klSym:%f",
                      a+i,
                      b+n-i,
                      hellinger(a+i+1,b+n-i,a+i,b+n-i+1),
                      totalVariation(a+i+1,b+n-i,a+i,b+n-i+1),
                      kl(a+i+1,b+n-i,a+i,b+n-i+1),
                      klSymmetric(a+i+1,b+n-i,a+i,b+n-i+1) ))
    }
}
