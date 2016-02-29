source("getConsensusOvarianSubtypes.R")
source("getRandomForestConsensusOvarianSubtypes.R")
## This file is produced from classificationAcrossDatasets.Rnw
load("esets.with.survival.RData")

# rescale per gene
esets.with.survival.scaled <- lapply(esets.with.survival, function(eset) {
  exprs(eset) <- t(scale(t(exprs(eset))))
  return(eset)
})

dataset.names <- names(esets.with.survival.scaled)

classification.vals.pam <- list()
classification.vals.rf <- list()

for(dataset.name in dataset.names) {
  left.out.dataset <- esets.with.survival.scaled[[dataset.name]]
  training.dataset.names <- dataset.names[dataset.names != dataset.name]
  
  consensus.classifier.output.pam <- getConsensusOvarianSubtypes(left.out.dataset, .dataset.names.to.keep = training.dataset.names)
  consensus.classifier.output.rf <- getRandomForestConsensusOvarianSubtypes(left.out.dataset, .dataset.names.to.keep = training.dataset.names)
  
  classification.vals.pam[[dataset.name]] <- consensus.classifier.output.pam$Annotated.eset$Ovarian.subtypes
  classification.vals.rf[[dataset.name]] <- consensus.classifier.output.rf$Annotated.eset$Ovarian.subtypes.rf
}

# Make some tables

print(table(unlist(lapply(esets.with.survival.scaled, function(eset) eset$Helland.subtypes)), unlist(classification.vals)))
print(table(unlist(lapply(esets.with.survival.scaled, function(eset) eset$Verhaak.subtypes)), unlist(classification.vals)))
print(table(unlist(lapply(esets.with.survival.scaled, function(eset) eset$Konecny.subtypes)), unlist(classification.vals)))
