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
  print(paste(k,names(tcm)[k],summary(mod)$r.squared))
}

# Regression without N
for ( k in 2:31){ 
  mod<-lm( target ~ tcm[,18]+tcm[,22]+
             tcm[,12] + tcm[,27] + tcm[,25]+
             tcm[,13] + tcm[,17] + tcm[,19] + 
             tcm[,20] + tcm[,31]
           + tcm[,k], data=xtcm); 
  print(paste(k,names(tcm)[k],summary(mod)$r.squared))
}
#> names(xtcm)[var_seq2]
#[1] "Ideas"             "Altruism"         
#[3] "ExcitementSeeking" "Order"            
#[5] "TenderMindedness"  "PositiveEmotions" 
#[7] "Actions"           "Values"           
#[9] "Trust"             "Deliberation"    



###################################################
#  Univariate aggregated models
# Neuroticism
indices<-1:6

univariate_aggregated_rsq<-function( indices ){
  agg_prof_diff <- data.frame( x = exp(-rowMeans(abs_centroid[,indices+1])),
                               y=centroid_and_marsat[,"target"])
  uvagg<-lm(y ~ x, data=agg_prof_diff)
  summary(uvagg)$r.squared
}

rsquared_neoac<-c( 0.791, 0.823, 0.845,0.579,0.694)

# Pretty good univariate prediction with R2=0.905
univariate_aggregated_rsq( c(2,3,7,8,17,27))
# For aggregate univariate predictions does not 
# seem much better is possible


