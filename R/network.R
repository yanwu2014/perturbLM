
#' Build Shared Nearest Neighbors matrix. From Seurat (https://github.com/satijalab/seurat)
#'
#' @param data.use Cells x features matrix to construct SNN with
#' @param k.param Number of neighbors
#' @param k.scale Neighborhood size to consider
#' @param prune.SNN Min fraction of shared nearest neighbors to consider
#' @return Sparse shared nearest neighbors matrix
#'
#' @importFrom FNN get.knn
#' @import Seurat
#' @export
#'
BuildSNN <- function(data.use, k.param = 10, k.scale = 10, prune.SNN = 1/15) {
  n.cells <- nrow(data.use)

  my.knn <- get.knn(data = data.use, k = min(k.scale * k.param, n.cells - 1))
  nn.ranked <- cbind(1:n.cells, my.knn$nn.index[, 1:(k.param - 1)])
  nn.large <- my.knn$nn.index

  w <- Seurat:::CalcSNNSparse(cell.names = rownames(data.use), k.param = k.param, nn.large = nn.large,
                              nn.ranked = nn.ranked, prune.SNN = prune.SNN, print.output = T)

  return(w)
}


#' Calculate effect of perturbations on gene modules
#'
#' @param coefs.matrix Regression coefficients matrix
#' @param gene.clusters.list List of gene clusters
#' @param gene.module.mapping Mapping gene modules to descriptions
#' @param min.coef Minimum cofficient magnitude for a gene to be considered
#' @return Matrix of perturbation effects on each gene modules
#'
#' @export
#'
CalcGeneModuleEffect <- function(coefs.matrix, gene.clusters.list, gene.module.mapping = NULL,
                                 min.coef = 0.025) {
  coefs.matrix <- coefs.matrix[apply(coefs.matrix, 1, function(x) any(abs(x) > min.coef)),]

  if (!is.null(gene.module.mapping)) {
    gene.modules <- gene.module.mapping$Description
    names(gene.modules) <- gene.module.mapping$Module
    names(gene.clusters.list) <- sapply(names(gene.clusters.list), function(x) gene.modules[[x]])
  }
  gene.clusters <- FlattenGenotypeList(gene.clusters.list)

  coefs.matrix <- coefs.matrix[rownames(coefs.matrix) %in% names(gene.clusters),]
  gene.clusters <- gene.clusters[rownames(coefs.matrix)]

  tfs.modules.matrix <- apply(coefs.matrix, 2, function(x) {
    tapply(x, gene.clusters, mean)
  })
  tfs.modules.matrix <- tfs.modules.matrix[,order(colnames(tfs.modules.matrix))]
  tfs.modules.matrix
}
