# Determine generalized hyperbolic parameters of Jung 
# transformed big five using big five parameters
# fitted parameters of big five
library(ghyp)
library(Matrix)
library(matlib)
lambda <- 1.000000
alpha.bar <-  1.788706

mu <-c( 3.053392,
        3.021512, 
        3.183587,
        3.154377,
        3.308383)

sigmaBaseVec <- c( 0.1314515387, 0.0003453978, 0.03400488, 0.02127298, 0.01147807,
                   0.0003453978, 0.5063454311, 0.05195018, 0.06297821, 0.02006035,
                   0.0340048782, 0.0519501804, 0.14269251, 0.03386428, 0.02592150,
                   0.0212729771, 0.0629782131, 0.03386428, 0.16951071, 0.04414651,
                   0.0114780685, 0.0200603541, 0.02592150, 0.04414651, 0.16969360 )

sigmaBase <- matrix( sigmaBaseVec, nrow=5)
gamma <- c( -0.030213436, 0.002608754, -0.028674928, -0.031492846, -0.043415295)

bigfiveDist <- ghyp(lambda=1, alpha.bar=alpha.bar, mu=mu, sigma=sigmaBase, gamma=gamma)


# y = t(A)x + e
# then A = cov(xy)/cov(xx)
linear_regr_transform<-function( corrxy,x ){
  cx <- cov(x)
  icx <- inv(cx)
  sx <- sqrt(diag(cx))
  sy <- rep(1, dim(corrxy)[1]) # rows
  covxy <-  corrxy
  for (j in 1:dim(corrxy)[1]){
    for (k in 1:dim(corrxy)[2]){
      covxy[j,k] <- corrxy[j,k]*sx[j]*sy[k]
    }
  }
  #A <- covxy %*% icx
  A <- corrxy %*% icx
  y <-  x %*% t(A)
  y
}

###################################################
# Zulf's jungTransfrom function
#

jung_bigfive<-base::matrix( nrow=4, ncol=5)

# from McRae/Costa Reinterpreting MBTI paper 1989
#
jung_bigfive[1,]<-c( 0.16,-0.74, 0.03, -0.03, 0.08)
jung_bigfive[2,]<-c( -0.06,0.1, 0.72, 0.04, -0.15)
jung_bigfive[3,]<-c( 0.06,0.19, 0.02, 0.44, -0.15)
jung_bigfive[4,]<-c( 0.11,0.15, 0.30, -0.06, -0.49)

vector_by_gram_schmidt <- function( J0 ){
  z <- rep(0.2,5)
  jx <- J0[1, ]/norm( J0[1,], type="2")
  z <- z - dot( z, jx)*jx
  jx <- J0[2, ]/norm( J0[2,], type="2")
  z <- z - dot( z, jx)*jx
  jx <- J0[3, ]/norm( J0[3,], type="2")
  z <- z - dot( z, jx)*jx
  jx <- J0[4, ]/norm( J0[4,], type="2")
  z <- z - dot( z, jx)*jx
  z <- 0.1*z /sum(abs(z))
  z  
}


jungTransform<-function( ocean_matrix ){
  J0 <- jung_bigfive
  J <- base::matrix( 0, nrow = 5 , ncol = 5 )
  J[ 1:4, ] <- J0
  J[ 5, ] <- vector_by_gram_schmidt( J0 )
  #J <- J + diag(5)
  #print(dim(J))
  x <- matrix(ocean_matrix, ncol=5)
  print(head(x))
  #print( J %*% x)
  #print(dim(x))
  
  #T <- (0.2 *  x) %*% J
  T <- 0.2*linear_regr_transform( t(J), x)
  #print(dim(T))
  #maxT <- max(max(T))  
  #minT <-min(min(T))
  #out <- (T-minT)/(maxT-minT)
  out <-T
  print(head(out))
  mjv<-colMeans(out) - rep(0.5,5)
  out <- apply( out, 1, function(z) z - mjv )
  #print(head(out))
  out
}

vf <- rghyp( 100000, prob.fem )
jtvf <- jungTransform(vf)

print("before fitting")
# fit generalized hyperbolic distribution
ghfit<-fit.ghypmv( t(jtvf), lambda=1.3, alpha.bar = 1, reltol=5e-7 )

print("fitting complete")
# obtain diagonal of covariance produce univariate
# ghyps
# use quantiles to get cut points
# design a class to hold all the cuts
# with method to give index of the personality 
# type.  let's keep it 2^P
#
# pseudocode
# pt <- PersonalityTypes( P, ghfit )
# vf <- rghyp( 1000, bigfiveDist)
# vm <- rghyp( 1000, bigfiveDist)
# typesf <- pt$personalityType( vf )
# typesm <- pt$personalityType( vm )
# countMatrix <- matrix( rep(0,nm*nf),nrow=nm)
# for (p in 1:length(typesm)){
#    for (q in 1:length(typesf)){
#       dist <- metric( vm[p], vf[q])
#       if( dist > tail_thresh ) {
#           countMatrix[p,q] = countMatrix[p,q]+1
#        }
#    }
# }
# Use S4 OOP
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

#library(methods)

# class method
#
#setClass( "PersonalityTypes",
#  slots = c( levels = "numeric",
#          ghdist = "ghyp",
#          cuts = "matrix"))

#setMethod("initCuts",
#          "PersonalityTypes",
#            function( ){
#              ncuts <- 2^levels
#              C <- sigma(ghdist)
#              M <- mu(ghdist)
#             a <- alpha(ghdist)
#              l <- lambda(ghdist)
#              for (j in 1:5){
#                m <- M[j]
#                s <- C[j,j]
#                g <- gamma[j]
#                dist <- ghyp(l,m,g,l)
#                cuts[j,]<-univariateQuantiles( ncuts, dist)
#              }
#            }
#          )

#setMethod( "personalityType", 
#           "PersonalityTypes",
#              return(ans)
#           })

alpha.bar.get<-function(sgh){
  tail(sgh$trace.pars$alpha.bar,1)
}

gamma.get<-function(sgh){
  tail(sgh$trace.pars$gamma,1)
}

lambda.get<-function(sgh){
  tail(sgh$trace.pars$lambda,1)
}


getPersonalityType <- function( bfv, level, ghdist )
{
  sgh <- ghyp.fit.info(ghdist)
  # determine cuts
  ncuts <- 2^level
  cuts <- matrix(rep(0, 5*ncuts),nrow=5)
  C <- vcov(ghdist)
  M <- mean(ghdist)
  a <- alpha.bar.get(sgh)
  l <- lambda.get(sgh)
  g <- gamma.get(sgh)
  for (j in 1:5){
    m <- M[j]
    s <- sqrt(C[j,j])
    g <- gamma[j]
    dist <- ghyp(mu=m,sigma=s, gamma=g,lambda=l,alpha.bar=a)
    dcuts <-univariateQuantiles( ncuts, dist)
    dcuts[is.infinite(dcuts)]<-5.0
    cuts[j,]<-dcuts
  }
  
  # use the cuts to determine the indices

  r <- length(cuts[1,])
  N <- dim(bfv)[1]
  out <- rep( 0, N)
  for ( vr in 1:N){
    idx <- rep(0,5)
    ans <- 0
    x <- bfv[vr, ]
    for (j in 1:5){
      idx[j]<-indexInOrderedVec( x[j], cuts[j,])
      ans <- ans + r^(j-1) * (idx[j]-1)
    }
    out[vr] <- ans
  }
  out
}

metric <- function( a, b ){
  Q(a,b) + Q2(a,b)
}

do_not_run <- function(){
nsample <- 1000000
vf <- matrix(rghyp( nsample, bigfiveDist),ncol=5)
vm <- matrix(rghyp( nsample, bigfiveDist), ncol=5)

jvf <- t(jungTransform( vf ))
jvm <- t(jungTransform( vm ))
typesf <- getPersonalityType( jvf, level = 2, ghfit )
typesm <- getPersonalityType( jvm, level = 2, ghfit )
nm <- 2*4^5
nf <- 2*4^5
countMatrix <- Matrix( nrow=nm, ncol=nf, data=0, sparse=T)

tail_thresh <- -0.5
for (p in 1:nsample){
  dist <- metric( vm[p,], vf[p,])
  if( dist > tail_thresh ) {
    tp <- typesf[p]
    tq <- typesm[p]
    if ((!is.na(tp)) && (!is.na(tq)) 
        && tp <= nm && tq <= nm ){
      #print(paste(tp,tq,dist))
      countMatrix[tp,tq] = countMatrix[tp,tq]+1
    }
  }
  
  if (p %% 200 == 0){
    density <- sum(countMatrix>0)/(nf*nm)
    print(paste(p, 'density=',density))
  }
}

}

scale_count_matrix <- function( A, target_width ){
  n<-dim(A)[1]
  if ( n != dim(A)[2]){
    return
  }
  if ( n %% target_width != 0 ){
    return
  }
  B <- matrix( data=0, ncol=target_width, 
              nrow = target_width )
  K <- n/target_width
  for (p in 1:target_width){
    for (q in 1:target_width){
    # determine the top left and bottom right
      x0 <- K*p
      x1 <- (K+1)*p -1
      y0 <- K*q
      y1 <- (K+1)*q -1
      new_count <- sum(A[x0:x1,y0:y1])
      B[p,q] <- new_count
    }
  }
  B
}
