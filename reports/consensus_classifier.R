source("getConsensusOvarianSubtypes.R")
## This file is produced from classificationAcrossDatasets.Rnw
#load("esets.with.survival.RData")
load("esets.not.rescaled.classified.RData")
# rescale per gene
#   esets.with.survival.scaled <- lapply(esets.with.survival, function(eset) {
#     exprs(eset) <- t(scale(t(exprs(eset))))
#     return(eset)
#   })

esets.scaled <- lapply(esets.not.rescaled.classified, function(eset) {
  exprs(eset) <- t(scale(t(exprs(eset))))
  return(eset)
})

# dataset.names <- names(esets.with.survival.scaled)
dataset.names <- names(esets.scaled)
                       
classification.vals <- list()
sample_ids <- list()

for(dataset.name in dataset.names) {
  # left.out.dataset <- esets.with.survival.scaled[[dataset.name]]
  left.out.dataset <- esets.scaled[[dataset.name]]
  training.dataset.names <- dataset.names[dataset.names != dataset.name]
  
  consensus.classifier.output <- getConsensusOvarianSubtypes(left.out.dataset, .dataset.names.to.keep = training.dataset.names)
  
  classification.vals[[dataset.name]] <- consensus.classifier.output$Annotated.eset$Ovarian.subtypes %>% as.vector
  sample_ids[[dataset.name]] <- paste(dataset.name,".",consensus.classifier.output$Annotated.eset %>% exprs %>% colnames, sep="")
}

# Make some tables

print(table(unlist(lapply(esets.scaled, function(eset) eset$Helland.subtypes)), unlist(classification.vals)))
print(table(unlist(lapply(esets.scaled, function(eset) eset$Verhaak.subtypes)), unlist(classification.vals)))
print(table(unlist(lapply(esets.scaled, function(eset) eset$Konecny.subtypes)), unlist(classification.vals)))

# Check if classification.vals matches the prediction for the training set- highly pure subtypes

# Check if classification.vals matches the prediction for the pure subtypes 
compare_subtypes <- function (sample_id)
{
pure_subtype = intersection_pooled.subtypes[sample_id,"Verhaak"]
predicted = 
  unlist(classification.vals)[((sample_ids %>% unlist) == sample_id) %>% which]
predicted_subtype = substr(predicted,start = 1,stop = 3)
pure_subtype == predicted_subtype
}

## the following gives the amounts of consensus between the pure subtypes and the purest subtypes respectively

lapply(rownames(intersection_pooled.subtypes), compare_subtypes) %>% unlist %>% table
lapply(rownames(Filtered_intersection_pooled.subtypes), compare_subtypes) %>% unlist %>% table