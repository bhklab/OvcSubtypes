
getGenePairConsensusOvarianSubtypes <- function(eset, .dataset.names.to.keep=names(esets.not.rescaled.classified), purest.subtypes = FALSE) {
  
  ### Load training data
  print("Loading training data")
  ## This file is produced from classificationAcrossDatasets.Rnw
  load("esets.not.rescaled.classified.RData")
  
#  library("org.Hs.eg.db")
#  konecny.supplementary.data <- read.xls(system.file("extdata", "jnci_JNCI_14_0249_s05.xls", package="MetaGx"), sheet=4)
#  konecny.entrez.ids <- konecny.supplementary.data$EntrezGeneID
#  
#  supplementary.type.1 <- read.xls(system.file("extdata", "journal.pone.0018064.s015.XLS", package="MetaGx"), sheet=1)
#  supplementary.type.2 <- read.xls(system.file("extdata", "journal.pone.0018064.s015.XLS", package="MetaGx"), sheet=2)
#  supplementary.type.4 <- read.xls(system.file("extdata", "journal.pone.0018064.s015.XLS", package="MetaGx"), sheet=3)
#  supplementary.type.5 <- read.xls(system.file("extdata", "journal.pone.0018064.s015.XLS", package="MetaGx"), sheet=4)
#  supplementary.tables <- list(C1=supplementary.type.1, C2=supplementary.type.2, C4=supplementary.type.4, C5=supplementary.type.5)
#  entrez.id.logFC.list <- lapply(supplementary.tables, function(x) {
#    ## Use the supplementary table's listed probe id and gene name to determine the Entrez ID
#    # If there is only one EntrezID that maps to a probe in hgu133plus2.db, use that Entrez ID.
#    # If there are multiple EntrezIDs that map to a probe, then use the EntrezID (if any) that corresponds to the provided gene symbol.
#    current.mapping <- suppressWarnings(AnnotationDbi::select(hgu133plus2.db, as.character(x$ID), c("ENTREZID", "SYMBOL")))
#    current.mapping <- current.mapping[ !is.na(current.mapping$ENTREZID), ]
#    colnames(x)[1:2] <- c("PROBEID", "SYMBOL")
#    mappings.with.unique.probeid <- current.mapping[ !(current.mapping$PROBEID %in% current.mapping$PROBEID[duplicated(current.mapping$PROBEID)]),]
#    mappings.with.duplicate.probeid <- current.mapping[ current.mapping$PROBEID %in% current.mapping$PROBEID[duplicated(current.mapping$PROBEID)],]
#    mappings.with.duplicate.probeid <- merge(x, mappings.with.duplicate.probeid, by=c("PROBEID", "SYMBOL"))[, c("PROBEID", "ENTREZID", "SYMBOL")]
#    mappings.with.duplicate.probeid <- unique(mappings.with.duplicate.probeid)
#    current.mapping <- rbind(mappings.with.unique.probeid, mappings.with.duplicate.probeid)
#    to.return <- merge(x, current.mapping, by="PROBEID")[, c("ENTREZID", "PROBEID", "logFC")]
#    return(to.return)
#    })
#  helland.entrez.ids <- do.call(c,sapply(entrez.id.logFC.list, function(x) x$ENTREZID))
#  names(helland.entrez.ids) <- NULL
#  
#  supplementary.data.sheet7 <- read.xls(system.file("extdata", "JCI65833sd1.xls", package="MetaGx"), sheet=7, skip=1)
#  
#  verhaak.gene.symbols <- lapply(levels(supplementary.data.sheet7$CLASS), function(y) as.character(supplementary.data.sheet7[supplementary.data.sheet7$CLASS==y,1]))
#  verhaak.gene.symbols <- do.call(c, verhaak.gene.symbols)
#  
#  do.call(c, as.list(org.Hs.egALIAS2EG)[verhaak.gene.symbols])
#  
#  verhaak.entrez.ids <- do.call(c, as.list(org.Hs.egALIAS2EG)[verhaak.gene.symbols])
#  
#  entrez.id.union <- unique(do.call(c, list(helland.entrez.ids, konecny.entrez.ids, verhaak.entrez.ids)))
#  save(verhaak.entrez.ids, file="verhaak.entrez.ids.RData")
#  save(entrez.id.union, file="entrez.id.union.RData")
  load("entrez.id.union.RData")
  load("verhaak.entrez.ids.RData")
  
  dataset.names.to.keep <- .dataset.names.to.keep
  
  esets.not.rescaled.classified <- esets.not.rescaled.classified[dataset.names.to.keep]
  
  esets.merged <- MetaGx::datasetMerging(esets.not.rescaled.classified, method = "intersect", standardization = "none")
  
  esets.merged <- esets.merged[fData(esets.merged)$EntrezGene.ID %in% verhaak.entrez.ids,]
  
  subtype.correspondances <- data.frame(Konecny=c("C1_immL", "C2_diffL", "C3_profL", "C4_mescL"),
                                        Verhaak=c("IMR", "DIF", "PRO", "MES"),
                                        Helland=c("C2", "C4", "C5", "C1"))
  
  # Only keep cases that are concordant for all three classifiers
  
  cases.to.keep <- 
    if (purest.subtypes) { 
      print("Purest subtypes")
      purest_subtypes <- rownames(Filtered_intersection_pooled.subtypes)
      (exprs(esets.merged) %>% colnames) %in% delete_leading_dataset_string(purest_subtypes)
    } else {
      match(esets.merged$Konecny.subtypes, subtype.correspondances$Konecny) ==
        match(esets.merged$Verhaak.subtypes, subtype.correspondances$Verhaak) &
        match(esets.merged$Verhaak.subtypes, subtype.correspondances$Verhaak) ==
        match(esets.merged$Helland.subtypes, subtype.correspondances$Helland)
    }  
  training.dataset <- esets.merged[,cases.to.keep]
  
  ### Once we are happy with the normalization / removal of discordant cases, this eset should be a package data file.
  
  train.labels <- training.dataset$Verhaak.subtypes
  levels(train.labels) <- paste0(levels(train.labels), "_consensus")
  
  intersecting.entrez.ids <- as.character(intersect(fData(training.dataset)$EntrezGene.ID, fData(eset)$EntrezGene.ID))
  
  print("Training Random Forest...")
  
  train.expression.matrix <- t(exprs(training.dataset)[match(intersecting.entrez.ids, fData(training.dataset)$EntrezGene.ID),])
  
  train.pairwise.matrix <-
    apply(combn(1:length(intersecting.entrez.ids),2), 2, function(pair) train.expression.matrix[,pair[1]] > train.expression.matrix[,pair[2]])
  
  rf.model <- randomForest(x=train.pairwise.matrix, y=train.labels)
  
  test.expression.matrix <- t(exprs(eset)[match(intersecting.entrez.ids, fData(eset)$EntrezGene.ID),])
  
  test.pairwise.matrix <-
    apply(combn(1:length(intersecting.entrez.ids),2), 2, function(pair) test.expression.matrix[,pair[1]] > test.expression.matrix[,pair[2]])
  
  my.predictions <- predict(rf.model, newdata = test.pairwise.matrix)
  my.predictions.probs <- predict(rf.model, newdata = test.pairwise.matrix, type = 'prob')
  
  eset$Ovarian.subtypes.rf <- my.predictions
  
  return(list(Annotated.eset=eset, rf.probs=my.predictions.probs))
}