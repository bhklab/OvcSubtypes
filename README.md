# Analysis and report generation for ovarian subtyping project


**_Instructions:_**

First, install MetaGxOvarian (contact Greg Chen or Deena Gendoo for the latest version).
Install MetaGx with the commands

```
git clone https://github.com/bhklab/MetaGx.git
R CMD BUILD MetaGx
R CMD INSTALL MetaGx
```

Clone this repository:

```
git clone https://github.com/bhklab/OvcSubtypes.git
```

In the "reports" directory, there are three main knitr files:
- reproduceResults.Rnw
- classificationAcrossDatasets.Rnw
- robustness.Rnw
Note that classificationAcrossDatasets.Rnw produces the file "esets.not.rescaled.RData" which contains samples of high-grade serious ovarian cancer, with genes **not** z-score rescaled by gene. In order to ensure consistency in all analyses, this file is used by robustness.Rnw and batch.cluster.all.R.

The three knitr files can be run to produce pdf files within an R session with the command:
```
library(knitr)
knit("reproduceResults.Rnw")
knit("classificationAcrossDatasets.Rnw")
knit("robustness.Rnw")
```

Note that several dependencies may need to be installed. Packages on CRAN can be installed with, for example,
```
install.packages("xtable")
```
Packages on Bioconductor can be installed with, for example,
```
source("https://bioconductor.org/biocLite.R")
biocLite("survcomp")
```

Within robustness.Rnw, there are two main components:
- a reproduction of clustering algorithms and comparison to the cluster labels given in original supplementary texts (producing heatmaps demonstrating concordance between our implementation and the original results)
- evaluation of robustness using prediction strength

To build this file without the prediction strength, comment out the line:
```
\Sexpr{knit_child('robustness_validation.Rnw')}
```

In the prediction strength analysis, we re-cluster each dataset with each algorithm 100 times. This takes a substantial amount of computation, so this is parallelized on an SGE cluster with the files "submit_batch_clustering.sh" and "batch.cluster.all.R". Ensuring that "esets.not.rescaled.RData" exists, the command
```
./submit_batch_clustering.sh
```
submits an array job to the cluster, creating a directory with clustering output. This directory is read in "robustness_validation.Rnw".
