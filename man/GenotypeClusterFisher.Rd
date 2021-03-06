% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment.R
\name{GenotypeClusterFisher}
\alias{GenotypeClusterFisher}
\title{Use Fisher's exact test to calculate enrichment for each genotype in each cluster}
\usage{
GenotypeClusterFisher(genotype.cluster.counts)
}
\arguments{
\item{genotype.cluster.counts}{Genotype x cluster counts matrix}
}
\value{
Genotype x cluster p-values matrix
}
\description{
Use Fisher's exact test to calculate enrichment for each genotype in each cluster
}
