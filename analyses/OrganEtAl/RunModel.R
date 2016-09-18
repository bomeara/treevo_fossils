source("/Users/bomeara/Documents/MyDocuments/GitClones/treevo_fossils/data/OrganEtAl/ProcessOrganData.R")

changingRateIntrinsic <-function(params, states, timefrompresent) {
  #params[1] is sd at present,
  #params[2] is multiplier for parabolic change
  #overall model is sd_now = params[1] * (1 + params[2]*(timefrompresent)^2)
  newdisplacement <- rnorm(n=length(states), mean=0, sd=params[1] * (1 + params[2]*(timefrompresent^2)))
  return(newdisplacement)
}


#assume generation time of 20 years. Tree is in MY time units.
TreeYears=1e6
timeStep<-20/TreeYears

continuous.rate.guess <- geiger::fitContinuous(phy, trait)$opt$sigsq #units are variance / MY
sd.per.gen.guess <- sd(rnorm(1e5,0,sqrt(continuous.rate.guess*timeStep)))

max.time.squared <- max(phytools::nodeHeights(phy))^2

#Assume rate at beginning is not too much greater / less than rate at end: say two standard dev include 10fold diff

sd.param2 <- (c(10,-10)/2)/max.time.squared

intrinsicFn=changingRateIntrinsic
extrinsicFn=nullExtrinsic
startingPriorsFns=c("uniform")
startingPriorsValues=matrix(range(trait),nrow=2,byrow=FALSE) #assume that the min value is the root state
intrinsicPriorsFns=c("exponential","normal") #do fixed for param2 b/c dont have a lot of data to use
intrinsicPriorsValues=matrix(c(
  rep(1/sd.per.gen.guess , 2),
  sd.param2), nrow=2, byrow=FALSE)
extrinsicPriorsFns=c("fixed")
extrinsicPriorsValues=matrix(c(0, 0), nrow=2, byrow=FALSE)

abcTolerance<-0.01

results <- doRun_prc(
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
  plot=FALSE,
  StartSims=100,
  epsilonProportion=0.1,
  epsilonMultiplier=0.7,
  nStepsPRC=5,
  numParticles=2000,
  jobName='OrganEtAl',
  stopRule=FALSE,
  multicore=TRUE,
  coreLimit=3
)

save(list=ls(), file="OrganEtAlResults.rda")
