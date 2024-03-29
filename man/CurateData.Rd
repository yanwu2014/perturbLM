% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/curation.R
\name{CurateData}
\alias{CurateData}
\title{Curate perturbation dataset (in Seurat object format)}
\usage{
CurateData(
  obj,
  guide,
  gene,
  ctrl.gene,
  var.genes,
  dataset.name,
  out.dir,
  dose.val = NULL,
  min.cells = 25,
  batch.col = NULL,
  cell.type.col = NULL,
  nfeature.col = "nFeature_RNA",
  size.factor.col = "nCount_RNA",
  mito.col = "percent.mt",
  nfolds = 5,
  seed = 23,
  export = T
)
}
\arguments{
\item{obj}{Input Seurat object to curate}

\item{guide}{Dictionary or vector defining guide RNA perturbations}

\item{gene}{Dictionary or vector defining gene targets}

\item{ctrl.gene}{Control gene target (i.e. NTC, AAVS)}

\item{var.genes}{Variable genes}

\item{dataset.name}{Name of curated dataset}

\item{out.dir}{Directory to write curated dataset}

\item{dose.val}{Perturbation dose (set to 1 if NULL)}

\item{min.cells}{Minimum cells per gene perturbation}

\item{batch.col}{Metadata column containing technical/biological batches}

\item{cell.type.col}{Metadata column containing cell types}

\item{nfeature.col}{Metadata column containing number of genes expressed}

\item{size.factor.col}{Metadata column containing number of total counts}

\item{mito.col}{Metadata column containing fraction of mitochondrial reads}

\item{nfolds}{Number of folds to use for test/train splitting}

\item{seed}{Random seed for defining folds}

\item{export}{If TRUE, write to file (default: TRUE)}
}
\description{
Curate perturbation dataset (in Seurat object format)
}
