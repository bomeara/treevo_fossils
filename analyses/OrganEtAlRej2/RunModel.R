library(TreEvo)

options(warn=1)
#options(error = utils::recover)
source("/Users/bomeara/Documents/MyDocuments/GitClones/treevo_fossils/data/OrganEtAl/ProcessOrganData.R")


#assume generation time of 20 years. Tree is in MY time units.
TreeYears=1e6
generation.time=500000
#timeStep<-1000/TreeYears
timeStep<-generation.time/TreeYears
continuous.rate.guess <- geiger::fitContinuous(phy, trait)$opt$sigsq #units are variance / MY
sd.per.gen.guess <- sd(rnorm(1e5,0,sqrt(continuous.rate.guess*timeStep)))



intrinsicFn=boundaryMinIntrinsic
extrinsicFn=nullExtrinsic
startingPriorsFns=c("uniform")
startingPriorsValues=matrix(range(trait),nrow=2,byrow=FALSE) #assume that the min value is the root state
intrinsicPriorsFns=c("exponential","uniform") 
intrinsicPriorsValues=matrix(c(
  rep(1/sd.per.gen.guess , 2),
  c(0, min(trait))), nrow=2, byrow=FALSE)
extrinsicPriorsFns=c("fixed")
extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE)

abcTolerance<-0.01


#dealing with zero length branches
phy$edge.length[which(phy$edge.length<0.01)] <- 0.01 



system.time(results <- doRun_rej(
  phy = phy,
  traits = trait,
  intrinsicFn=intrinsicFn,
  extrinsicFn=extrinsicFn,
  startingPriorsFns=startingPriorsFns,
  startingPriorsValues=startingPriorsValues,
  intrinsicPriorsFns=intrinsicPriorsFns,
  intrinsicPriorsValues=intrinsicPriorsValues,
  extrinsicPriorsFns=extrinsicPriorsFns,
  extrinsicPriorsValues=extrinsicPriorsValues,
  TreeYears=TreeYears,
  standardDevFactor=0.2,
  StartSims=10,
 jobName='OrganEtAlRej',
  multicore=FALSE,
  coreLimit=1,
  savesims=TRUE
))


save(list=ls(), file="OrganEtAlResults.rda")
