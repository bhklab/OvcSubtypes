
getRandomForestConsensusOvarianSubtypes <- function(eset, .dataset.names.to.keep=names(esets.with.survival.scaled)) {
  
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
  
  # Convert pData factor columns to strings
  #esets.with.survival.scaled <- lapply(esets.with.survival.scaled, function(eset) {
  #  ind <- sapply(pData(eset), is.factor)
  #  pData(eset)[ind] <- lapply(pData(eset)[ind], as.character)
  #  return(eset)
  #})
  
  esets.survival.scaled.merged <- MetaGx::datasetMerging(esets.with.survival.scaled, method = "intersect", standardization = "none")
  
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
  
  intersecting.gene.names <- as.character(intersect(rownames(exprs(training.dataset)), rownames(exprs(eset))))
  
  print("Training Random Forest...")
  rf.model <- randomForest(x=t(exprs(training.dataset)[intersecting.gene.names,]), y=train.labels)
  
  my.predictions <- predict(rf.model, newdata = t(exprs(eset)))
  my.predictions.probs <- predict(rf.model, newdata = t(exprs(eset)), type = 'prob')
  
  eset$Ovarian.subtypes.rf <- my.predictions
  
  return(list(Annotated.eset=eset, rf.probs=my.predictions.probs))
}