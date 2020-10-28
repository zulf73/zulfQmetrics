# We want to consider 2048 types
# in terms of the McRae Facets, six per
# OCEAN variable then consider holistic
# narrative self in terms of these values

library(matlib)
library(data.table)
library(dplyr)




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


facet_profile<-function( v, t=1){
  #w0<-linear_regr_transform( facet_mtx/100., matrix(v,ncol=1))
  w0 <- (facet_mtx/100.) %*% matrix(v,ncol=1)
  w <- t(w0)
  out<-list(type=t)
  for (j in 1:30){
    key <- facet_names[j]
    value <- w[j]
    out[[key]] <- value
  }
  out
}

center_and_rescale<-function( df ) {
  dft <- df
  nc <- ncol(dft)
  for (k in 2:nc){
    v <- unlist(df[,k ])
    v <- (v - min(v))/(max(v)-min(v))
    meanv <- mean(v)
    print(paste('meanv',meanv))
    v <- v - meanv + 0.5
    dft[,k] <- v
  }
  dft
}

binarise<-function( df ) {
  dft <- df
  nc <- ncol(dft)
  for (k in 2:nc){
    v <- unlist(df[,k ])
    v[v>0.5]<-1
    v[v <= 0.5]<-0
    dft[,k] <- v
  }
  dft
}

histTotalHigh<-function( facets_types ){
  totalOneCount <- rep(0, 30)
  ftmRowSums<-rowSums(ftm_bin[,2:31])
  for (r in 1:781){
    v <- ftmRowSums[r]
    totalOneCount[ v ] <- totalOneCount[v] + 1
  }
  totalOneCount
}

narrative_from_binary<-function(binary_code, name="Tony "){
  out <- ""
  n <- length(binary_code)
  for (j in 1:n){
    code <- binary_code[j]
    if (code>1){
      name <- paste("Type",code,sep="")
    }
    if (code==0 || code==1){
      if (code==1){
         out <- paste(out, name," ", facet_state[[j]][[2]],". ",sep="", collapse="")
      }
    }
  }
  out
}



nsample <- 100000
vf <- matrix(rghyp( nsample, bigfiveDist),ncol=5)

jvf <- t(jungTransform( vf ))
typesf <- getPersonalityType( jvf, level = 2, ghfit )
nf <- 2*4^5

proflist <-list()
for (p in 1:nsample){
  tp <- typesf[p]
  nlist <- facet_profile( vf[p,], t=tp)
  proflist[[p]] <- nlist
}



#-------- Running code ----
df <- rbindlist( proflist)
facet_type_means <- df %>% 
  group_by(type) %>% 
  summarise(across(everything(),mean))
ftm <- center_and_rescale(facet_type_means)
#ggplot( data=ftm, aes(x=TenderMindedness)) +geom_histogram()
bin_ftm<-binarise(ftm)