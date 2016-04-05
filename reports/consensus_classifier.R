library(Biobase)
library(pamr)
library(randomForest)
library(survcomp)

load("Filtered_intersection_pooled.subtypes.RData")

source("getConsensusOvarianSubtypes.R")
source("getRandomForestConsensusOvarianSubtypes.R")
source("~/repos/MetaGx/R/datasetMerging.R")
source("~/repos/MetaGx/R/stripWhiteSpace.R")
## This file is produced from classificationAcrossDatasets.Rnw
#load("esets.with.survival.RData")
load("esets.not.rescaled.classified.RData")

esets.scaled <- lapply(esets.not.rescaled.classified, function(eset) {
  exprs(eset) <- t(scale(t(exprs(eset))))
  return(eset)
})

dataset.names <- names(esets.scaled)

classification.vals.pam <- list()
classification.vals.rf <- list()
classification.vals.pam.probs <- list()
classification.vals.rf.prob <- list()
sample_ids <- list()

for(dataset.name in dataset.names) {
  # left.out.dataset <- esets.with.survival.scaled[[dataset.name]]
  left.out.dataset <- esets.scaled[[dataset.name]]
  training.dataset.names <- dataset.names[dataset.names != dataset.name]
  
  consensus.classifier.output.pam <- getConsensusOvarianSubtypes(left.out.dataset, .dataset.names.to.keep = training.dataset.names)
  consensus.classifier.output.rf <- getRandomForestConsensusOvarianSubtypes(left.out.dataset, .dataset.names.to.keep = training.dataset.names)
  
  sample_ids[[dataset.name]] <- paste(dataset.name,".",consensus.classifier.output.pam$Annotated.eset %>% exprs %>% colnames, sep="")
  classification.vals.pam[[dataset.name]] <- consensus.classifier.output.pam$Annotated.eset$Ovarian.subtypes
  classification.vals.pam.probs[[dataset.name]] <- consensus.classifier.output.pam$posterior.probs 
  classification.vals.rf[[dataset.name]] <- consensus.classifier.output.rf$Annotated.eset$Ovarian.subtypes.rf
  classification.vals.rf.prob[[dataset.name]] <- consensus.classifier.output.rf$rf.probs
}

# check if the pam classifier and random forest classifier result in the same subtypes
table(unlist(classification.vals.pam) == unlist(classification.vals.rf))

# Get the margins, given the list of probabilities of datasets
get_margins <- function(dataset.probs) {
  lapply(1:nrow(dataset.probs), function (i)
  sort(cl[i,],partial=4)[4] - sort(cl[i,],partial=3)[3]) %>% unlist
} 

margins.rf = lapply(classification.vals.rf.prob, get_margins) %>% unlist
names(margins.rf) <- unlist(sample_ids)

## Filter the predictions
# rule 1 : Extract those samples whose probabilities are above 0.5
get_top_prob <- function (dataset.number) {
  lapply(1:length(classification.vals.rf[[dataset.number]]), function (i)
    post.probs[i,classification.vals.rf[[dataset.number]][[i]]]) %>% unlist
}

top_prob.rf <- lapply(1:length(classification.vals.rf), get_top_prob) %>% unlist
names(top_prob.rf) <- names(margins.rf)

top_prob_above0.5 <- top_prob.rf > 0.5
table(top_prob_above0.5)

# rule 2 : Extract those samples whose margins are above prob_cutoff
prob_cutoff <- 0.05
top_prob_margin_cutoff <- margins.rf > prob_cutoff
table(top_prob_margin_cutoff)

(top_prob_above0.5 & top_prob_margin_cutoff) %>% table

## The filtered dataset 
predicted.rf <- unlist(classification.vals.rf)[top_prob_above0.5 & top_prob_margin_cutoff]
names(top_prob.rf) <- names(margins.rf)
# Make some tables
helland <- (unlist(lapply(esets.scaled, function(eset) eset$Helland.subtypes)))[names(predicted.rf)]
verhaak <- (unlist(lapply(esets.scaled, function(eset) eset$Verhaak.subtypes)))[names(predicted.rf)]
konecny <- (unlist(lapply(esets.scaled, function(eset) eset$Konecny.subtypes)))[names(predicted.rf)]

print(table(helland, predicted.rf))
print(table(verhaak, predicted.rf))
print(table(konecny, predicted.rf))

# Check if classification.vals matches the prediction for the training set- highly pure subtypes

# Check if classification.vals matches the prediction for the pure subtypes 
compare_subtypes <- function (sample_id)
{
pure_subtype = intersection_pooled.subtypes[sample_id,"Verhaak"]
predicted = 
  predicted.rf[((sample_ids %>% unlist) == sample_id) %>% which]
predicted_subtype = substr(predicted,start = 1,stop = 3)
pure_subtype == predicted_subtype
}

## the following gives the amounts of consensus between the pure subtypes and the purest subtypes respectively
print("Subtypes common to pure subtypes: \n")
lapply(rownames(intersection_pooled.subtypes), compare_subtypes) %>% unlist %>% table
print("Subtypes common to purest subtypes: \n")
lapply(rownames(Filtered_intersection_pooled.subtypes), compare_subtypes) %>% unlist %>% table


########## Venn Diagrams 
# overlaps = get_overlaps_margin[[margin]][[subtype_order]]
# #overlaps.parts = lapply(overlaps$parts, delete_leading_dataset_string)
# IMR = unlist(classification.vals) %in% "IMR_consensus" %>% which
# overlaps$parts[[1]] %in% unlist(sample_ids)[IMR]
########## Survival Curves

## below is the survival plot AFTER removing the patients with low margin 
survival_combined.df <- data.frame(sample.names = names(predicted.rf), prediction = predicted.rf)

survival_combined_pooled.df <- pooled.subtypes[names(predicted.rf),]
survival_combined_pooled.df$groups <- predicted.rf
# patients with survival data
survival_combined_pooled.df = na.omit(survival_combined_pooled.df)
# since we name the subtypes of the filtered set to that of Verhaak

pval <- summary(coxph(surv.obj ~ groups + strata(data.source), survival_combined_pooled.df))$sctest["pvalue"]

hr.out <- survcomp::hazard.ratio(x=survival_combined_pooled.df$groups, surv.time=survival_combined_pooled.df$years_to_death, surv.event=survival_combined_pooled.df$vital_status, strat=survival_combined_pooled.df$data.source)

if(length(hr.out$hazard.ratio) == 1) {
  text <- paste0(sprintf("HR: %.3f (%.3f-%.3f)\n", hr.out$hazard.ratio, hr.out$lower, hr.out$upper), sprintf("Logrank p = %.1E", pval))
} else {
  for(i in 1:length(hr.out$hazard.ratio)) {
    text <- paste0(text, sprintf("HR %s: %.3f (%.3f-%.3f)\n", levels(survival_filtered.df$groups)[i+1], hr.out$hazard.ratio[i], hr.out$lower[i], hr.out$upper[i]))
  }
  text <- paste0(text, sprintf("Logrank p = %.1E", pval))
}
cols <- 1:4

title <- "Combined Classifier"

km.coxph.plot(surv.obj ~ groups, survival_combined_pooled.df, x.label="Time (years)", y.label = "Overall Survival", main.title="", show.n.risk = FALSE, n.risk.step=2, leg.text = levels(survival_combined_pooled.df$groups), leg.pos="topright", leg.inset=0, leg.bty="n", n.risk.cex=0.85, cex=0.4, o.text="", .col=cols, cex.lab=1.5)
title(title, cex.main=2)
text(0,0.05, text, cex=0.85, pos=4)
