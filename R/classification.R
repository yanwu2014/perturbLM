

#' Map cells to a reference dataset using a kNN classifier
#'
#' @param query.data Query matrix
#' @param ref.data Reference matrix
#' @param ref.groups Reference cell types
#' @param k Number of neighbors
#'
#' @return Matrix of cell type probabilities
#'
#' @import FNN
#' @export
#'
MapCellsKNN <- function(query.data, ref.data, ref.groups, k = 30) {
  knn.res <- FNN::get.knnx(t(ref.data), t(query.data), k = k)

  ref.groups <- as.character(ref.groups)
  unique.grps <- unique(ref.groups)

  grp.frac <- apply(knn.res$nn.index, 1, function(x) {
    grp.counts <- vector(mode = "integer", length = length(unique.grps))
    names(grp.counts) <- unique.grps

    grp.tbl <- table(ref.groups[x])
    grp.counts[names(grp.tbl)] <- grp.tbl
    grp.counts/k
  })
  colnames(grp.frac) <- colnames(query.data)

  return(grp.frac)
}
MapCellsKNN <- compiler::cmpfun(MapCellsKNN)


#' Map mouse genes to human genes
#'
#' @param mouse.genes Vector of mouse genes
#' @param mouse.human.gene.file File with mouse to human gene mappings
#'
#' @return vector of human genes with mouse genes as names
#' @export
#'
MouseHumanMapping <- function(mouse.genes, mouse.human.gene.file) {
  mouse.gene.mapping.df <- read.table(mouse.human.gene.file, sep = "\t", header = T, stringsAsFactors = F)
  mouse.gene.mapping.df <- mouse.gene.mapping.df[!duplicated(mouse.gene.mapping.df$Gene.name),]

  mouse2human <- mouse.gene.mapping.df$Human.gene.name
  names(mouse2human) <- mouse.gene.mapping.df$Gene.name

  mouse.genes <- mouse.genes[mouse.genes %in% names(mouse2human)]
  return(mouse2human[mouse.genes])
}
MouseHumanMapping <- compiler::cmpfun(MouseHumanMapping)
