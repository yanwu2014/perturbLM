#### Functions for data curation

#' Efficiently load tabular datasets using data.table
#'
#' @param fi File path
#' @param sep Delimiter character (default tab)
#' @import data.table
#' @export
#'
FastRead <- function(fi, sep = "\t") {
  con <- file(fi, "r")
  cells <- readLines(con, n = 1)
  close(con)

  cells <- strsplit(cells, split = sep)[[1]]
  cells <- cells[!is.na(cells)]
  cells <- cells[cells != ""]
  cells <- gsub("\"", "", cells)

  m <- fread(fi, sep = sep, header = F, skip = 1)
  m <- as.matrix(m, rownames = 1)
  colnames(m) <- cells

  m
}

#' Curate perturbation dataset (in Seurat object format)
#'
#' @param obj Input Seurat object to curate
#' @param guide Dictionary or vector defining guide RNA perturbations
#' @param gene Dictionary or vector defining gene targets
#' @param ctrl.gene Control gene target (i.e. NTC, AAVS)
#' @param var.genes Variable genes
#' @param dataset.name Name of curated dataset
#' @param out.dir Directory to write curated dataset
#' @param dose.val Perturbation dose (set to 1 if NULL)
#' @param min.cells Minimum cells per gene perturbation
#' @param batch.col Metadata column containing technical/biological batches
#' @param cell.type.col Metadata column containing cell types
#' @param nfeature.col Metadata column containing number of genes expressed
#' @param size.factor.col Metadata column containing number of total counts
#' @param mito.col Metadata column containing fraction of mitochondrial reads
#' @param nfolds Number of folds to use for test/train splitting
#' @param seed Random seed for defining folds
#' @param export If TRUE, write to file (default: TRUE)
#'
#' @export
#'
CurateData <- function(obj,
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
                       export = T) {
  require(Seurat)
  require(SeuratDisk)

  DefaultAssay(obj) <- "RNA"

  if (is.list(condition)) {
    condition <- FlattenGroups(condition)
  }

  if (!is.vector(condition)) {
    stop("Condition must be a list or vector")
  }

  conditions.keep <- names(table(condition)[table(condition) >= min.cells])
  condition <- condition[condition %in% conditions.keep]

  cells.use <- intersect(colnames(obj), names(condition))
  if (length(cells.use) == 0) {
    stop("colnames(obj) must overlap with names(condition)")
  }

  obj <- obj[,cells.use]
  obj[["condition"]] <- condition[cells.use]

  if (is.null(dose.val)) {
    dose.val <- 1
  }
  obj$dose_val <- dose.val
  obj$control <- obj$condition == ctrl
  if (sum(obj$control) == 0) {
    stop("No groups named ctrl in condition")
  }

  if (!is.null(cell.type.col)) {
    obj[["cell_type"]] <- obj[[cell.type.col]]
  }

  if (!is.null(batch.col)) {
    obj[["batch"]] <- obj[[batch.col]]
  }

  if (mito.col != "percent.mt") {
    obj[["percent.mt"]] <- obj[[mito.col]]
  }

  if (nfeature.col != "nFeature_RNA") {
    obj[["nFeature_RNA"]] <- obj[[nfeature.col]]
  }

  if (size.factor.col != "nCounts_RNA") {
    obj[["nCounts_RNA"]] <- obj[[nfeature.col]]
  }

  obj@assays$RNA@var.features <- var.genes
  obj@assays$RNA@meta.features$hvg <- rownames(obj) %in% var.genes

  ## Get cross validation folds for data
  folds <- perturbLM::SplitFoldsByGroup(obj$condition, nfolds = nfolds, seed = seed)
  for (i in 1:length(folds)) {
    obj[[paste0("split", i)]] <- folds[[i]][colnames(obj)]
  }

  ## Export to file
  if (export) {
    out.seurat.file <- paste0(out.dir, "/", dataset.name, ".rds")
    saveRDS(obj, file = out.seurat.file)

    out.h5seurat.file <- paste0(out.dir, "/", dataset.name, ".h5Seurat")
    SaveH5Seurat(obj, filename = out.h5seurat.file, overwrite = T)
    Convert(out.h5seurat.file, dest = "h5ad", overwrite = T)
    unlink(out.h5seurat.file)
  }

  obj
}
