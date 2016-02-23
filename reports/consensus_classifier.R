library(MetaGx)
library(MetaGxOvarian)

load("esets.with.survival.RData")

# rescale per gene
esets.with.survival.scaled <- lapply(esets.with.survival, function(eset) {
  exprs(eset) <- t(scale(t(exprs(eset))))
  return(eset)
})

esets.survival.scaled.merged <- MetaGx::datasetMerging(esets.with.survival.scaled)

subtype.correspondances <- data.frame(Konecny=c("C1_immL", "C2_diffL", "C3_profL", "C4_mescL"),
                                      Verhaak=c("IMR", "DIF", "PRO", "MES"),
                                      Helland=c("C2", "C4", "C5", "C1"))

cases.to.keep <- match(esets.survival.scaled.merged$Konecny.subtypes, subtype.correspondances$Konecny) ==
  match(esets.survival.scaled.merged$Verhaak.subtypes, subtype.correspondances$Verhaak) &
  match(esets.survival.scaled.merged$Verhaak.subtypes, subtype.correspondances$Verhaak) ==
  match(esets.survival.scaled.merged$Helland.subtypes, subtype.correspondances$Helland)

# leave-one-dataset-out thing
classification.vals <- list()
for(dataset.name in names(esets.with.survival.scaled)) {
  left.out.dataset <- esets.with.survival.scaled[[dataset.name]]
  training.datasets <- esets.with.survival.scaled[names(esets.with.survival.scaled) != dataset.name]
  training.datasets.merged <- MetaGx::datasetMerging(training.datasets, method="intersect")
  
  pamr.model <- pamr.train(data=list(x=exprs(training.datasets.merged), y=training.datasets.merged$Verhaak.subtypes),
                           gene.subset = intersect(rownames(exprs(training.datasets.merged)), rownames(exprs(left.out.dataset))))
  
  my.predictions <- pamr.predict(pamr.model, newx=exprs(left.out.dataset)[intersect(rownames(exprs(training.datasets.merged)), rownames(exprs(left.out.dataset))),], threshold = 1)
  
  classification.vals[[dataset.name]] <- my.predictions
}