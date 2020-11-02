# MPQ Big Five Corr
library(ghyp)
library(glmnet)
library(geometry)

construct_big_five_dist <- function(){
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
  bigfiveDist  
}


facet_transform<-function( v, post=""){
  facet_names<-c(
    "Anxiety",
    "AngryHostility",
    "Depression",
    "Selfconsciousness",
    "Impulsiveness",
    "Vulnerability",
    "Warmth",
    "Gregariousnness",
    "Assertiveness",
    "Activity",
    "ExcitementSeeking",
    "PositiveEmotions",
    "Fantasy",
    "Aesthetics",
    "Feeling",
    "Actions",
    "Ideas",
    "Values",
    "Trust",
    "Straightforwardness",
    "Altruism",
    "Compliance",
    "Modesty",
    "TenderMindedness",
    "Competence",
    "Order",
    "Dutifulness",
    "AchievementStriving",
    "SelfDiscipline",
    "Deliberation"
  )
  
  facet_state<-list(
    list("is not anxious", "is anxious"),
    list( "is not angry and hostile","is angry and hostile"),
    list( "is not depressed","is depressed"),
    list( "is not self-consciousness","is self-conscious"),
    list( "is not impulsive","is impulsive"),
    list( "is not vulnerable","is vulnerable"),
    list( "is cold", "is warm"),
    list( "is quiet", "is talkative"),
    list( "is not assertive","is assertive"),
    list( "is passive","is active"),
    list( "does not seek excitement"," seeks excitement"),
    list( "is not filled with positive emotions",
          "is filled with positive emotions"),
    list( "is not imaginative","is imaginative"),
    list( "is not beauty oriented", "is beauty oriented"),
    list( "is not feeling driven","is feeling driven"),
    list( "does not like to act","likes to act"),
    list( "is not filled with ideas",
          "is filled with ideas"),
    list( "cares for values","does not care for values"),
    list( "is untrusting","is trusting"),
    list( "is not straightforward","is straightforward"),
    list( "is mean","is generous"),
    list( "is rebellious", "is compliant"),
    list( "is immodest","is modest"),
    list( "is tough minded","is tender-minded"),
    list( "is not competence oriented","is competence oriented"),
    list( "is disorganised","is organised"),
    list( "is unreliable", "is reliable"),
    list( "is a slacker","is achievement striving"),
    list( "has no self-discipline","has self-discipline"),
    list( "is unfocused","is deliberate and focused")
  )
  
  facet_data <- c( 
    -10, 0,	-5,	0,	71,
    -1,	1,	1,	-47,	51,
    -3,	-18,	-9,	4,	65,
    -3,	-9,	-28,	17,	59,
    5,	-22,	32,	-8,	37,
    -16,	-28,	-8,	14,	59,
    8,	4,	66,	38,	-11,
    -9,	-8,	57,	14,	-9,
    29,	21,	38,	-46,	-24,
    2,	39,	42,	-29,	-4,
    11,	-10,	38,	-39,	-7,
    19,	11,	53,	13,	-10,
    47,	-24,	7,	-12,	22,
    55,	5,	4,	22,	9,
    42,	17,	41,	4,	27,
    45,	2,	13,	10,	-36,
    67,	11,	-6,	-9,	-9,
    51,	4,	13,	12,	-8,
    12,	63,	18,	9,	-31,
    5,	63,	-8,	22,	1,
    12,	49,	43,	19,	0,
    -9,	71,	-18,	-9,	-18,
    5,	67,	-10,	9,	9,
    24,	53,	27,	11,	10,
    20,	-13,	15,	47,	-28,
    -27,	-7,	-5,	64,	5,
    8,	19,	8,	53,	-4,
    6,	-27,	15,	56,	-5,
    -6,	-2,	1,	63,	-19,
    -3,	8,	-27,	38,	-7
  )
  
  facet_mtx <-matrix( facet_data, ncol=5)
  
  ans <- t((facet_mtx/100.) %*% matrix(v,nrow=5))
  colnames(ans) <- sapply(facet_names, function(x) 
    paste(x,post,collapse=""))
    
  ans
}


# We want to sample from big five and then use
# Richard Robins, Avshalom Caspi and Terrie E. Moffitt, 
# "Two Personalities, One Relationship: 
# Both Partners' Personality Traits Shape the Quality
# of Their Relationship" Journal of Personality and Social Psychology, 2000, 
# Vol. 79 (2) p. 251-259

robins_marital_sat<-function( vf, vm ){
  
  mpq_mtx_vals <- c( 
    31, 10, 31, 12, -39,
    35, 4, 42, -26, -16,
    32, 42, 7, -7, -9,
    10, 3, 61, 27, 3,
    -2, -2, -8, -12, 73,
    -7, -17, 0, -50, 11,
    9, -13, -10, -21, 27,
    -24, 52, -11, 3, -8,
    -35, 15, 4, 15, 23,
    -28, 23, 10, 14, -6,
    40, 11, 4, -13, 9 )
  
  mpq_mtx <- matrix( data=mpq_mtx_vals, ncol=5)

  A <-mpq_mtx/100.
  
  xf <- matrix( vf, nrow=5)
  yf <-A %*% xf

  xm <- matrix( vm, nrow=5)
  ym <-A %*% xm
  
  male_sat_vec <- c(4,-7,-3,4,-22,-4,
                           -14,11,-2,-5)
  fem_sat_vec <- c(11,7,0,16,-7,-18,
                          -6,16,-6,11)
  male_sat <- geometry::dot(yf[1:10],male_sat_vec)/100.
  fem_sat <- geometry::dot(ym[1:10],fem_sat_vec)/100.
  ans<- 0.5*(male_sat + fem_sat)
  ans
}

create_marsat_dataset <- function( nsamples ){
  
  big_five_dist <- construct_big_five_dist()
  data_mtx <- matrix( data=0, ncol= 5+5+30+1, 
                      nrow=nsamples)
  bff <- rghyp(nsamples, big_five_dist )
  bfm <- rghyp(nsamples, big_five_dist )

  facets_f <- facet_transform( bff )
  facets_m <- facet_transform( bfm )
  
  for (j in 1:nsamples){
    data_mtx[j,1:5]<-bff[j,]
    data_mtx[j,6:10]<-bfm[j,]
    data_mtx[j,11:40] <-
        exp(-abs(facets_f[j,] -facets_m[j,]))
    data_mtx[j,41] <-robins_marital_sat(bff[j,],bfm[j,])
  }
  data_mtx
}

df <- data.frame( create_marsat_dataset(100) )
names(df)<-c("o1","c1","e1","a1", "n1",
             "o2","c2","e2","a2", "n2",
             facet_names, "target")
mod <- lm(  target ~ .-o1-c1-e1-a1-n1-o2-c2-e2-a2-n2 -target, 
            data=df )
summary(mod)



