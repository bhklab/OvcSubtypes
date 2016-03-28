source("getConsensusOvarianSubtypes.R")
source("getRandomForestConsensusOvarianSubtypes.R")
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

dataset.names <- names(esets.with.survival.scaled)

classification.vals.pam <- list()
classification.vals.rf <- list()

for(dataset.name in dataset.names) {
  # left.out.dataset <- esets.with.survival.scaled[[dataset.name]]
  left.out.dataset <- esets.scaled[[dataset.name]]
  training.dataset.names <- dataset.names[dataset.names != dataset.name]
  
  consensus.classifier.output.pam <- getConsensusOvarianSubtypes(left.out.dataset, .dataset.names.to.keep = training.dataset.names)
  consensus.classifier.output.rf <- getRandomForestConsensusOvarianSubtypes(left.out.dataset, .dataset.names.to.keep = training.dataset.names)
  
  classification.vals.pam[[dataset.name]] <- consensus.classifier.output.pam$Annotated.eset$Ovarian.subtypes
  classification.vals.rf[[dataset.name]] <- consensus.classifier.output.rf$Annotated.eset$Ovarian.subtypes.rf
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
survival_combined.df <- data.frame(sample.names = unlist(sample_ids), prediction = unlist(classification.vals))

survival_combined_pooled.df <- pooled.subtypes[unlist(sample_ids),]
survival_combined_pooled.df$groups <- unlist(classification.vals)
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
