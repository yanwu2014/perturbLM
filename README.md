# perturbLM
Linear model for analyzing genetic perturbations with scRNA-seq readouts

## Install instructions
* Open RStudio or R session
* Install devtools if not installed: `install.packages(devtools)`
* Install liger for GSEA enrichment: `devtools::install_github("JEFworks/liger")`
* Install qvalue for FDR correction: `source("https://bioconductor.org/biocLite.R"); biocLite("qvalue");`
* Install swne for helper functions: `devtools::install_github("yanwu2014/swne")`
* Install perturbLM: `devtools::install_github("yanwu2014/perturbLM")`
