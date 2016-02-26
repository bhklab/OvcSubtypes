## if threshold.auto, use pamr.adaptthresh to select a threshold; otherwise, use the threshold parameter

getConsensusOvarianSubtypes <- function(eset, .dataset.names.to.keep=names(esets.with.survival.scaled), threshold.auto=TRUE, threshold=1) {
  
  ### Load training data
  print("Loading training data")
  ## This file is produced from classificationAcrossDatasets.Rnw
  load("esets.with.survival.RData")
  
  # rescale per gene
  esets.with.survival.scaled <- lapply(esets.with.survival, function(eset) {
    exprs(eset) <- t(scale(t(exprs(eset))))
    return(eset)
  })
  
  dataset.names.to.keep <- .dataset.names.to.keep
  
  esets.with.survival.scaled <- esets.with.survival.scaled[dataset.names.to.keep]
  
  esets.survival.scaled.merged <- MetaGx::datasetMerging(esets.with.survival.scaled, method = "intersect")
  
  subtype.correspondances <- data.frame(Konecny=c("C1_immL", "C2_diffL", "C3_profL", "C4_mescL"),
                                        Verhaak=c("IMR", "DIF", "PRO", "MES"),
                                        Helland=c("C2", "C4", "C5", "C1"))
  
  # Only keep cases that are concordant for all three classifiers
  
  cases.to.keep <- match(esets.survival.scaled.merged$Konecny.subtypes, subtype.correspondances$Konecny) ==
    match(esets.survival.scaled.merged$Verhaak.subtypes, subtype.correspondances$Verhaak) &
    match(esets.survival.scaled.merged$Verhaak.subtypes, subtype.correspondances$Verhaak) ==
    match(esets.survival.scaled.merged$Helland.subtypes, subtype.correspondances$Helland)
  
  training.dataset <- esets.survival.scaled.merged[,cases.to.keep]
  
  ### Once we are happy with the normalization / removal of discordant cases, this eset should be a package data file.
  
  train.labels <- training.dataset$Verhaak.subtypes
  levels(train.labels) <- paste0(levels(train.labels), "_consensus")
 
  ### The function pamr.adaptthresh seems to have a scoping issue, so a workaround is to save 
  # required data in the global environment
  
  consensus.global.obj <- list(training.dataset=training.dataset, train.labels=train.labels, eset=eset)
  
  assign("consensus.global.obj", consensus.global.obj, .GlobalEnv)
   
  pamr.model <- pamr.train(data=list(x=exprs(consensus.global.obj$training.dataset), y=consensus.global.obj$train.labels),
                           gene.subset = intersect(rownames(exprs(consensus.global.obj$training.dataset)), rownames(exprs(consensus.global.obj$eset))))
  #pamr.model <- pamr.train(data=list(x=exprs(training.dataset), y=train.labels),
  #                         gene.subset = intersect(rownames(exprs(training.dataset)), rownames(exprs(eset))))
  
  if(threshold.auto == TRUE) {
    threshold <- pamr.adaptthresh(pamr.model)
  } else {
    threshold = threshold
  }
  
  assign("pamr.model", pamr.model, .GlobalEnv)
  
  my.predictions <- pamr.predict(pamr.model, newx=exprs(eset)[intersect(rownames(exprs(training.dataset)), rownames(exprs(eset))),], threshold = threshold)
  my.predictions.posterior <- pamr.predict(pamr.model, newx=exprs(eset)[intersect(rownames(exprs(training.dataset)), rownames(exprs(eset))),], threshold = threshold, type="posterior")
  
  eset$Ovarian.subtypes <- my.predictions
  
  return(list(Annotated.eset=eset, pamr.model=pamr.model, posterior.probs=my.predictions.posterior))
}