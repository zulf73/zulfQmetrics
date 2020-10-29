#  A deeper look at marital satisfaction
# data set is tcm

# First look at single variable predictions
uvstats <- list()
for ( nfeat in 2:31){
  name_var <- names(tcm)[nfeat]
  print(paste("Univariate Feature ",name_var ))
  uvmod <- lm( target ~ tcm[,nfeat], data=tcm)
  S<-summary(uvmod)
  uvstats<-append(uvstats, list(var=name_var,
                                stat = S$r.squared))
  
}
 # build triangular

mod1 <- lm( target ~ tcm[,3] + tcm[,18])

for ( k in 2:31){ 
  mod<-lm( target ~ tcm[,3]+tcm[,18]+tcm[,22]+
             tcm[,12] +tcm[,2] + tcm[,27] + tcm[,25]+
             tcm[,13]
           + tcm[,k], data=tcm); 
  print(paste(k,summary(mod)$r.squared))
}