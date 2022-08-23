# perturbLM
A lightweight R package that uses linear models to analyze the effects of chemical/genetic perturbations, 
conditions, and disease states on single gene expression. Inspired by the MIMOSCA python package (https://github.com/asncd/MIMOSCA)

## Install instructions
```
if (!require("remotes")) {
  install.packages("remotes")
}
remotes::install_github("yanwu2014/perturbLM")
```

## Tutorials
A [tutorial](https://yanwu2014.github.io/perturbLM/demos/elasticnet_demo.nb.html) on how to fit an 
ElasticNet model on single cell perturbational response data
and evaluate model performance on a per perturbation basis. This enables interpretation of 
perturbation effects on the transcriptome as well as other covariates (such as cell type).

Coming soon:
1. A tutorial on how to use linear models to infer non-linear interactions between perturbations
2. A tutorial on how to infer cell type specific perturbation effects
