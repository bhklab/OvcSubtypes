### module load R/3.3.0

# .libPaths("/mnt/work1/users/bhklab/Rlib")

### install CRAN and Bioconductor packages
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("devtools", "gdata", "knitr", "HiDimDA", "survival", "reshape2", "genefu", "annotate", "hgu133plus2.db", "stringr", "survMisc", "xtable", "gridExtra", "Biobase", "GSVA", "sparsediscrim", "survcomp", "NMF", "ggplot2", "e1071", "randomForest", "clue", "GEOquery", "logging", "metafor"))

### build and install MetaGx
# library(devtools)
# devtools::install_github(repo="bhklab/MetaGx")
# library(MetaGx)

### Load compendium of ovarian cancer datasets
# install.packages("/mnt/work1/users/bhklab/Data/MetaGxData-private/MetaGxOvarian_0.99.0.tar.gz", repos=NULL)

### clone the github repository
# system("git clone https://github.com/bhklab/OvcSubtypes.git")

########################
### load libraries
########################

library(gdata)
library(survival)
library(survcomp)
library(reshape2)
library(genefu)
library(annotate)
library(HiDimDA)
library(hgu133plus2.db)
library(stringr)
library(survMisc) 
library(xtable)
ibrary(grid)
library(gridExtra)
library(Biobase)
library(GSVA)
library(sparsediscrim)
library(MetaGxOvarian)
library(ggplot2)
library(e1071)
library(randomForest)
library(knitr)
library(xtable)
library(Biobase)
library(randomForest)
library(metafor)
library(dplyr)
library(VennDiagram)
library(GGally)


library(MetaGx)
library(MetaGxOvarian)


########################
### knit files
########################

# knitr::knit(file.path("reports", "reproduceResults.Rnw"))
knitr::knit(file.path("reproduceResults.Rnw"))
knitr::knit(file.path("classificationAcrossDatasets.Rnw"))
knitr::knit(file.path("robustness.Rnw"))
