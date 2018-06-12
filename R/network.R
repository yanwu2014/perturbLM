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
    gene.clusters.list <- gene.clusters.list[names(gene.clusters.list) %in% names(gene.modules)]
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
