# correlation of v_f v_m
.libPaths("C:\\Users\\Zulfi\\userLibs")
library(ghyp)
Bvals = c( 25, -8, -4, -3, -9,
           3, 10, 11, 5, -3,
           5, 10, 14, -3, -7,
           -5, 3, 1, 10, -5,
           1, -4, -11, 4, 11)

B = matrix( Bvals/100.0, nrow=5 )

# fitted parameters of big five
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
sigmaBase <- sigmaBase + diag(5)*0.30

gamma <- c( -0.030213436, 0.002608754, -0.028674928, -0.031492846, -0.043415295)

# make a 10 x 10 sigma matrix that includes
# the correlations
Sigma2 <- cbind( rbind(sigmaBase, B ), rbind( B, sigmaBase))
gamma2 <- cbind(gamma,gamma)
mu2 <- cbind(mu,mu)

prob.fem <- ghyp(lambda=1, alpha.bar=alpha.bar, mu=mu, sigma=sigmaBase, gamma=gamma)
prob.mal <- ghyp(lambda=1, alpha.bar=alpha.bar, mu=mu, sigma=sigmaBase, gamma=gamma)
prob.comb <- ghyp(lambda=1, alpha.bar=alpha.bar, mu=mu2, sigma=Sigma2, gamma=gamma2)


Q<- function( vf, vm ){
pf <- pghyp( vf, object=prob.fem)
pm <- pghyp( vm, object=prob.mal)
pc <- pghyp( c(vf,vm), object=prob.comb)

ans <- log( pc/(pm*pf) + 1e-5)
if (is.nan(ans) || is.infinite(ans)){
  ans <- log(1e-5)
}
ans
}


sampleRomanticDistance<-function( ntimes ){
  ans <- c()
  for (k in 1:ntimes){
    # now sim from fem, mal and print Q
    vf <- rghyp( 1, prob.fem )
    vm <- rghyp( 1, prob.mal )
    dist <- Q(vf,vm)
    ans <-append(ans, dist)
  }
  ans
}


GenTraitCorr <- c(
10,13,6,45,-34,
35,-5,-9,60,28,
-4,-4,6,4,71,
-10,2,3,66,8,
40,0,16,49,8,
3,-9,45,46,4,
-65,-9,-3,46,17,
-22,2,20,66,-9,
-15,0,23,65,-2,
-5,-5,27,62,-6,
-5,5,25,68,-3,
-15,6,13,63,8,
12,-7,7,21,70,
-2,4,10,77,4,
-8,3,21,73,0,
9,-2,-12,44,8,
1,-4,18,71,-19,
8,2,70,19,8,
13,22,43,16,40,
-2,0,65,-4,6,
57,16,-2,9,22,
11,23,43,6,46,
7,26,31,22,48,
0,60,14,-14,28,
11,49,29,0,56,
68,4,14,32,-9,
30,50,25,-4,41,
34,11,29,18,60,
-2,14,50,11,14,
48,26,14,0,21,
20,5,12,66,-2,
-43,53,-5,4,6,
3,-22,44,25,-11,
37,-13,-2,53,27,
15,51,11,8,10,
-8,-12,45,37,2,
-24,53,-12,8,1,
3,-19,42,38,2,
7,50,-6,4,-12,
2,-18,33,60,1,
29,43,-4,-18,11,
40,0,16,49,8,
13,39,4,1,48,
21,-11,1,61,15,
5,42,-4,9,-4)

B0 <- matrix( GenTraitCorr/100.0, ncol=5)

Q2 <- function( vf, vm ) {
  vd <- matrix( vf-vm, ncol=1)
  #print(norm(vd,"1"))
  w <- B0 %*% (vd)
  nw <- length(w)
  #print(nw)
  ans <-  exp(- norm(w,type="1")/length(w))
  ans
}

sampleQ2<-function( ntimes ){
  ans <- c()
  for (k in 1:ntimes){
    # now sim from fem, mal and print Q
    vf <- rghyp( 1, prob.fem )
    vm <- rghyp( 1, prob.mal )
    dist <- Q2(vf,vm)
    ans <-append(ans, dist)
  }
  ans
}

.libPaths("C:\\Users\\Zulfi\\userLibs")

# male/fem complainer/target
mct <- data.matrix(read.table( "mctraits.csv", sep=","))
fct <- data.matrix(read.table( "fctraits.csv", sep=","))
mtt <- data.matrix(read.table("mttraits.csv", sep=","))
ftt <- data.matrix(read.table("fttraits.csv", sep=","))


# unrealistic sampling
#vm <- as.matrix(rnorm(5, 2.5, 1.2))
#vf <- as.matrix(rnorm(5, 2.5, 1.2))



#scen = hpp.scenario(rate = 30, num.events = 300, num.sims = 100)


simRel <- function( lambda, totalTime = 20, xm, xf)
{
  curTime <- 0
  curRomanticCapital <- 200
  romanticCapital <- c(curRomanticCapital)
  eventTimes <- c( 0 )
  while (curTime < totalTime ){
    newTime <- rexp( n=1,rate = lambda )
    curTime <- curTime + newTime
    
    for (k in 1:15){
      signalMVal <- runif(n=1)
      signalFVal <- runif(n=1)
      if ( signalMVal > xm[k]  && signalFVal < xf[k] ) {
        eventTimes <- append( curTime, eventTimes )
        penalty <- abs(xm[k]) + abs(xf[k])
        penalty<-penalty*10
        #print(penalty)
        curRomanticCapital <- curRomanticCapital - penalty 
        romanticCapital <- append( curRomanticCapital, romanticCapital)
      }
    }
    #print( paste( curTime, curRomanticCapital))
    
  }
  list( x=eventTimes, y=romanticCapital)
}

#sR <- simRel( lambda=15,totalTime=10, xm=xm, xf=xf)
#plot(sR$x,sR$y, pch=',')


Q3<-function( vm, vf){
  xm <- abs((0.01/10.) *(mct + mtt) %*% vm)
  xf <- abs((0.01/10.) * (fct + ftt) %*% vf)
  sR <- simRel( lambda=100,totalTime=10, xm=xm, xf=xf)
  lateTimes<-sR$x[sR$y<0]
  ans<-lateTimes[1]
  print(lateTimes)
  ans
}


sampleQ3<-function( ntimes ){
  ans <- c()
  for (k in 1:ntimes){
    # now sim from fem, mal and print Q
    vf <- matrix(rghyp( 1, prob.fem ), ncol=1)
    vm <- matrix(rghyp( 1, prob.mal ),ncol=1)
    dist <- Q3(vf,vm)
    ans <-append(ans, dist)
  }
  ans
}
