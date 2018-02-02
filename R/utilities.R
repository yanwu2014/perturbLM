
#' Read data from either a tab separated values file or a 10X genomics sparse matrix directory
#'
#' @param matrix.dir Input matrix file or directory
#' @return Matrix
#' @export
#'
ReadData <- function(matrix.dir) {
  if (dir.exists(matrix.dir)) {
    counts <- as.matrix(Seurat::Read10X(matrix.dir))
  } else {
    counts <- read.table(matrix.dir, sep = "\t", header = T, row.names = 1)
  }
  colnames(counts) <- make.names(colnames(counts))
  return(counts)
}


#' Filter rows and columns of counts matrix
#'
#' @param counts Input counts matrix
#' @param min.cells.frac Min number of cells to keep a gene
#' @param min.genes Min genes expressed to keep a cell
#' @param min.expr Minimum level to count as expressed
#' @param max.transcripts Max number of UMIs per cell
#'
#' @return Filtered counts matrix
#' @import Matrix
#' @export
#'
FilterData <- function(counts, min.cells.frac, trim, min.genes = 500, min.expr = 0, max.transcripts = 50000) {
  min.cells <- round(ncol(counts)*min.cells.frac)
  counts <- counts[ , Matrix::colSums(counts) < max.transcripts]
  counts <- counts[ , Matrix::colSums(counts > min.expr) > min.genes]
  counts <- counts[Matrix::rowSums(counts > min.expr) > min.cells, ]
  if (trim > 0) {
    counts <- t(.winsorize_matrix(t(counts), trim = trim))
  }
  return(counts)
}


## Trim input matrices. From PAGODA2
.winsorize_matrix <- function(mat, trim) {
  if(trim  >  0.5) { trim <- trim/ncol(mat)  }
  wm <- winsorizeMatrix(mat, trim)
  rownames(wm) <- rownames(mat)
  colnames(wm) <- colnames(mat)
  return(wm)
}


#' Convert dataframe to matrix, specifying all column names
#'
#' @param df Dataframe
#' @param output.name Column of df to use as matrix values
#' @param row.col Column of df to use as matrix rows
#' @param col.col Column of df to use as matrix columns
#'
#' @return Matrix
#' @importFrom reshape2 acast
#' @export
#'
UnflattenDataframe <- function(df, output.name, row.col = 'Gene', col.col = 'Group') {
  df <- df[c(row.col, col.col, output.name)]
  colnames(df) <- c('row', 'column', output.name)
  mat.out <- reshape2::acast(df, row~column, value.var = output.name)
  return(mat.out)
}
