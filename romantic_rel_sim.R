.libPaths("C:\\Users\\Zulfi\\userLibs")

# male/fem complainer/target
mct <- data.matrix(read.table( "mctraits.csv", sep=","))
fct <- data.matrix(read.table( "fctraits.csv", sep=","))
mtt <- data.matrix(read.table("mttraits.csv", sep=","))
ftt <- data.matrix(read.table("fttraits.csv", sep=","))


# unrealistic sampling
vm <- as.matrix(rnorm(5, 2.5, 1.2))
vf <- as.matrix(rnorm(5, 2.5, 1.2))
xm <- abs((0.01/10.) *(mct + mtt) %*% vm)
xf <- abs((0.01/10.) * (fct + ftt) %*% vf)


#scen = hpp.scenario(rate = 30, num.events = 300, num.sims = 100)


simRel <- function( lambda, totalTime = 20, xm, xf)
{
  curTime <- 0
  curRomanticCapital <- 1000
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
        print(penalty)
        curRomanticCapital <- curRomanticCapital - penalty 
        romanticCapital <- append( curRomanticCapital, romanticCapital)
        }
    }
    print( paste( curTime, curRomanticCapital))
    
  }
  list( x=eventTimes, y=romanticCapital)
}

sR <- simRel( lambda=15,totalTime=10, xm=xm, xf=xf)
plot(sR$x,sR$y, pch=',')