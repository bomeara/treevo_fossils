library(TreEvo)
library(ape)
library(phytools)
library(geiger)
original.dir <- getwd()
setwd("/Users/bomeara/Documents/MyDocuments/GitClones/treevo_fossils/data/OrganEtAl")
phy <- read.tree("tree4.tre")
data <- read.csv("supplementaltable8.csv", stringsAsFactor=FALSE, strip.white=TRUE)
data <- data[!is.na(data$Ln.Genome.Size),]
trait <- data$Ln.Genome.Size
names(trait) <- data$Taxa
for (i in sequence(Ntip(phy))) {
	phy$tip.label[i] <- data$Taxa[which(as.character(data$Codes.Tree.4)==phy$tip.label[i])]	
}

#plot(contMap(phy, trait))

dinos <- getMRCA(phy, c("Scutellosaurus lawleri", "Gallus gallus"))
phy <- extract.clade(phy, dinos)
pruned <- treedata(phy, trait, warnings=FALSE)
phy <- pruned$phy
trait <- pruned$data[,1]
names(trait) <- rownames(pruned$data)
rm(dinos, data, pruned, i)
setwd(original.dir)