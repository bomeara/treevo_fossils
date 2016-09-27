library(TreEvo)
library(geiger)
library(ape)
library(phytools)
setwd("/Users/bomeara/Documents/MyDocuments/GitClones/treevo_fossils/analyses/OrganEtAlLocalRuns")
options(warn=1)
#options(error = utils::recover)
#source("/Users/bomeara/Documents/MyDocuments/GitClones/treevo_fossils/data/OrganEtAl/ProcessOrganData.R")
#traits <- trait
load("InputData.rda")
#assume generation time of 20 years. Tree is in MY time units.
TreeYears=1e6
generation.time=10000
#timeStep<-1000/TreeYears
timeStep<-generation.time/TreeYears
continuous.rate.guess <- geiger::fitContinuous(phy, traits)$opt$sigsq #units are variance / MY
sd.per.gen.guess <- sd(rnorm(1e5,0,sqrt(continuous.rate.guess*timeStep)))

trait <- traits

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

StartSims=50000
nrepSim<-StartSims
multicore=FALSE
checkpointFile=NULL
niter.goal=5


#dealing with zero length branches
phy$edge.length[which(phy$edge.length<0.01)] <- 0.01 

timeStep<-generation.time/TreeYears
	#splits<-getSimulationSplits(phy) #initialize this info
	taxon.df <- getTaxonDFWithPossibleExtinction(phy)
startingValuesGuess=c()
 intrinsicStatesGuess=c()
extrinsicStatesGuess=c()

	#figure out number of free params
	numberParametersTotal<-dim(startingPriorsValues)[2] +  dim(intrinsicPriorsValues)[2] + dim(extrinsicPriorsValues)[2]
	numberParametersFree<-numberParametersTotal
	numberParametersStarting<-0
	numberParametersIntrinsic<-0
	numberParametersExtrinsic<-0
	freevariables<-matrix(data=NA, nrow=2, ncol=0)
	titlevector<-c()
	freevector<-c()

	namesForPriorMatrix<-c()
	PriorMatrix<-matrix(c(startingPriorsFns, intrinsicPriorsFns, extrinsicPriorsFns), nrow=1, ncol=numberParametersTotal)
	for (a in 1:dim(startingPriorsValues)[2]) {
		namesForPriorMatrix<-c(paste("StartingStates", a, sep=""))
	}
	for (b in 1:dim(intrinsicPriorsValues)[2]) {
		namesForPriorMatrix<-append(namesForPriorMatrix, paste("IntrinsicValue", b, sep=""))
	}
	for (c in 1:dim(extrinsicPriorsValues)[2]) {
		namesForPriorMatrix <-append(namesForPriorMatrix, paste("ExtrinsicValue", c, sep=""))
	}
	PriorMatrix<-rbind(PriorMatrix, cbind(startingPriorsValues, intrinsicPriorsValues, extrinsicPriorsValues))
	colnames(PriorMatrix)<-namesForPriorMatrix
	rownames(PriorMatrix)<-c("shape", "value1", "value2")

	for (i in 1:dim(startingPriorsValues)[2]) {
		priorFn<-match.arg(arg=startingPriorsFns[i],choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE)
		if (priorFn=="fixed") {
			numberParametersFree<-numberParametersFree-1
			freevector<-c(freevector, FALSE)
		}
		else {
			numberParametersStarting<-numberParametersStarting+1
			freevariables<-cbind(freevariables, startingPriorsValues[, i])
			titlevector <-c(titlevector, paste("Starting", numberParametersStarting))
			freevector<-c(freevector, TRUE)
		}
	}
	for (i in 1:dim(intrinsicPriorsValues)[2]) {
		priorFn<-match.arg(arg=intrinsicPriorsFns[i],choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE)
		if (priorFn=="fixed") {
			numberParametersFree<-numberParametersFree-1
			freevector<-c(freevector, FALSE)
		}
		else {
			numberParametersIntrinsic<-numberParametersIntrinsic+1
			freevariables<-cbind(freevariables, intrinsicPriorsValues[, i])
			titlevector <-c(titlevector, paste("Intrinsic", numberParametersIntrinsic))
			freevector<-c(freevector, TRUE)
		}
	}

	for (i in 1:dim(extrinsicPriorsValues)[2]) {
		priorFn<-match.arg(arg=extrinsicPriorsFns[i],choices=c("fixed", "uniform", "normal", "lognormal", "gamma", "exponential"),several.ok=FALSE)
		if (priorFn=="fixed") {
			numberParametersFree<-numberParametersFree-1
			freevector<-c(freevector, FALSE)
		}
		else {
			numberParametersExtrinsic<-numberParametersExtrinsic+1
			freevariables<-cbind(freevariables, extrinsicPriorsValues[, i])
			titlevector <-c(titlevector, paste("Extrinsic", numberParametersExtrinsic))
			freevector<-c(freevector, TRUE)
		}
	}

	#initialize guesses, if needed
	if (length(startingValuesGuess)==0) { #if no user guesses, try pulling a value from the prior
		startingValuesGuess<-rep(NA,length(startingPriorsFns))
		for (i in 1:length(startingPriorsFns)) {
			startingValuesGuess[i]<-pullFromPrior(startingPriorsValues[,i],startingPriorsFns[i])
		}
	}
	if (length(intrinsicStatesGuess)==0) { #if no user guesses, try pulling a value from the prior
		intrinsicStatesGuess<-rep(NA,length(intrinsicPriorsFns))
		for (i in 1:length(intrinsicPriorsFns)) {
			intrinsicStatesGuess[i]<-pullFromPrior(intrinsicPriorsValues[,i],intrinsicPriorsFns[i])
		}
	}
	if (length(extrinsicStatesGuess)==0) { #if no user guesses, try pulling a value from the prior
		extrinsicStatesGuess<-rep(NA,length(extrinsicPriorsFns))
		for (i in 1:length(extrinsicPriorsFns)) {
			extrinsicStatesGuess[i]<-pullFromPrior(extrinsicPriorsValues[,i],extrinsicPriorsFns[i])
		}
	}

	brown<-fitContinuous(phy=phy, dat=traits, model="BM", ncores=1, control=list(niter=100)) #it actually runs faster without checking for cores. And we parallelize elsewhere
	lambda<-fitContinuous(phy=phy, dat=traits, model="lambda", ncores=1, control=list(niter=100))
	delta<-fitContinuous(phy=phy, dat=traits, model="delta", ncores=1, control=list(niter=100))
	ou<-fitContinuous(phy=phy, dat=traits, model="OU", ncores=1, control=list(niter=100))
	white<-fitContinuous(phy=phy, dat=traits, model="white", ncores=1, control=list(niter=100))

	cat("Setting number of starting points for Geiger optimization to")
	niter.brown.g <- round(max(10, min(niter.goal/TreEvo:::solnfreq(brown),100)))
	cat(paste("\n",niter.brown.g, "for Brownian motion"))
	niter.lambda.g <- round(max(10, min(niter.goal/TreEvo:::solnfreq(lambda),100)))
	cat(paste("\n",niter.lambda.g, "for lambda"))
	niter.delta.g <- round(max(10, min(niter.goal/TreEvo:::solnfreq(delta),100)))
	cat(paste("\n",niter.delta.g, "for delta"))
	niter.OU.g <- round(max(10, min(niter.goal/TreEvo:::solnfreq(ou),100)))
	cat(paste("\n",niter.OU.g, "for OU"))
	niter.white.g <- round(max(10, min(niter.goal/TreEvo:::solnfreq(white),100)))
	cat(paste("\n",niter.white.g, "for white noise"))

	load("CheckPoint.trueFreeValuesANDSummaryValues.Rsave")
	trueFreeValuesMatrix<-trueFreeValuesANDSummaryValues[,1:numberParametersFree]
	summaryValuesMatrix<-trueFreeValuesANDSummaryValues[,-1:-numberParametersFree]


	res<-PLSRejection(summaryValuesMatrix, trueFreeValuesMatrix, phy, traits, abcTolerance)
	input.data<-rbind(jobName, length(phy[[3]]), timeStep, StartSims, standardDevFactor, abcTolerance)
  print(res)
	rejectionResults<-vector("list")
	#names(rejectionResults)<-c("input.data", "PriorMatrix", "phy", "traits")
#save(trueFreeValuesMatrix,res, file="BarbsTestofDistanceCalc2.Rdata")
	rejectionResults$input.data<-traits
	rejectionResults$PriorMatrix<-PriorMatrix
	rejectionResults$phy<-phy
	rejectionResults$traits<-traits
	rejectionResults$trueFreeValuesANDSummaryValues<-trueFreeValuesANDSummaryValues
  rejectionResults$SimTime<-1
  rejectionResults$abcDistances<-res$abcDistances
	rejectionResults$particleDataFrame<-res$particleDataFrame
  rejectionResults$CredInt<-CredInt(res$particleDataFrame)
	rejectionResults$HPD<-HPD(res$particleDataFrame)
save(rejectionResults, file="RejectionResults.rda")
#plotPosteriors(rejectionResults$particleDataFrame, rejectionResults$PriorMatrix)
param.names <- c("Root state (pg of DNA)", "Rate of evolution (ln(pg)^2 / MY)", "Minimum (pg of DNA)")
pdf(file="OrganEmpiricalResults.pdf", width=15, height=4)
for (i in sequence(3)) {
	all.vals <- trueFreeValuesMatrix[,i]
	best.vals <- trueFreeValuesMatrix[which(rejectionResults$abcDistances < quantile(rejectionResults$abcDistances, abcTolerance)),i]
	if(i!=2) {
		all.vals <- exp(all.vals) 	
		best.vals <- exp(best.vals)
	}
	posteriorCurve=getUnivariatePosteriorCurve(best.vals, from=min(all.vals), to=max(all.vals))
	priorCurve=getUnivariatePosteriorCurve(all.vals, from=min(all.vals), to=max(all.vals))
	plotUnivariatePosteriorVsPrior(posteriorCurve, priorCurve, label=param.names[i])
}
dev.off()

param.names <- c("Root state (ln(pg) of DNA)", "Rate of evolution (ln(pg)^2 / MY)", "Minimum (ln(pg) of DNA)")
pdf(file="OrganEmpiricalResultsNonLog.pdf", width=15, height=4)
for (i in sequence(3)) {
	all.vals <- trueFreeValuesMatrix[,i]
	best.vals <- trueFreeValuesMatrix[which(rejectionResults$abcDistances < quantile(rejectionResults$abcDistances, abcTolerance)),i]
	#maybe do prior curve instead
	posteriorCurve=getUnivariatePosteriorCurve(best.vals, from=min(all.vals), to=max(all.vals))
	priorCurve=getUnivariatePosteriorCurve(all.vals, from=min(all.vals), to=max(all.vals))
	plotUnivariatePosteriorVsPrior(posteriorCurve, priorCurve, label=param.names[i])
}
dev.off()

traits.df <- data.frame(exp(trait))
colnames(traits.df) <- "Genome size (pg)"
ggplot(data=traits.df) + geom_bar(stat="identity")


pdf(file="OrganEtAlTree.pdf", width=10, height=7)
plot(contMap(phy, exp(trait)))
dev.off()
system("open *.pdf")