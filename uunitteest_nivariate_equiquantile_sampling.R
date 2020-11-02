# UNIT TEST FOR EQUIDISTRIBUTED QUANTILE SAMPLING
m<- 0.3
s<- 0.45
g<- -0.1
l <- 1.2
a <- 1.3
testUnivariateDist <- ghyp(mu=m,sigma=s, gamma=g,lambda=l,alpha.bar=a)


indexInOrderedVec<-function( val, vec ){
  idx <- 1
  w <-which( vec > val)
  if (length(w)>0){
    idx <- head(w , 1)
  } else {
    idx <- length(vec)
  }
  idx
}

univariateQuantiles<-function( ncuts, distrib ){
  probs <- 1/ncuts*seq(1,ncuts)
  cuts <- qghyp(probs,distrib)
  cuts
}

testBinSamples <-function( nsamples, ncuts, dist ){
  samples <- rghyp( nsamples, dist )
  bins <- rep(0, ncuts)
  cuts <- univariateQuantiles( ncuts, dist)
  for ( x in samples ){
    idx <- indexInOrderedVec(x, cuts)
    bins[ idx ] <- bins[ idx ] + 1
  }
  bins
}

# plot utility
pplot<-function(x){
  ggplot(data.frame(z=x),aes(x=1:length(x),y=z))+
    geom_line()
}
