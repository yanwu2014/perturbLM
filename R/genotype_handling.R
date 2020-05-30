#### Functions for handling genotype/phenotype dictionaries


#' Convert genotypes list to design matrix
#'
#' @param genotypes.list Named list of cell vectors for each genotype
#' @param max.genotypes Maximum number of genotypes per cell
#' @param min.cells Minumum number of cells per genotype, including combo genotypes
#'
#' @return Sparse binary matrix where rows are cells, and columns are genotypes. 1 means the cell received the genotype.
#'
#' @import Matrix
#' @import hash
#' @export
#'
DesignMatrixGenotypes <- function(genotypes.list, max.genotypes = 2, min.cells = 10) {
  cell.names <- unique(unlist(genotypes.list, F, F))
  single.genotypes <- names(genotypes.list)

  single.mat <- Matrix::Matrix(0, length(cell.names), length(single.genotypes),
                               dimnames = list(cell.names, single.genotypes))
  for (i in 1:ncol(single.mat)) {
    single.mat[genotypes.list[[i]], i] <- 1
  }
  single.mat <- single.mat[Matrix::rowSums(single.mat) > 0, ]

  combo.rows <- single.mat[Matrix::rowSums(single.mat) > 1,]
  combo.genotypes <- unique(apply(combo.rows, 1, function(x) paste(names(x[x == 1]), collapse = ":", sep = "")))
  if (max.genotypes > 1 && length(combo.genotypes) > 0) {
    combo.genotypes.list <- hash::hash(keys = combo.genotypes)
    for (i in 1:nrow(combo.rows)) {
      mat.row <- combo.rows[i,]
      genotype <- paste(names(mat.row[mat.row == 1]), collapse = ":", sep = "")
      cell <- rownames(combo.rows)[[i]]
      combo.genotypes.list[[genotype]] <- c(combo.genotypes.list[[genotype]], cell)
    }
    combo.genotypes.list <- as.list(combo.genotypes.list)
    combo.genotypes.list[['keys']] <- NULL
    combo.genotypes.list <- combo.genotypes.list[combo.genotypes]

    combo.mat <- Matrix(0, nrow(single.mat), length(combo.genotypes), dimnames = list(rownames(single.mat), combo.genotypes))
    for (i in 1:ncol(combo.mat)) {
      combo.mat[combo.genotypes.list[[i]],i] <- 1
    }
    design.mat <- cbind(single.mat, combo.mat)

  } else {
    design.mat <- single.mat
  }

  design.mat <- design.mat[,Matrix::colSums(design.mat) > min.cells]
  return(as(design.mat, "dgCMatrix"))
}
DesignMatrixGenotypes <- compiler::cmpfun(DesignMatrixGenotypes)



#' Pad design matrix with zeros for cells without a called genotype
#'
#' @param X Design matrix
#' @param row.names Rownames (must contain at least all rows of X)
#'
#' @return Design matrix, with all row.names not in X set to all zeros
#' @export
#'
PadDesignMatrix <- function(X, row.names) {
  Y <- matrix(0, length(row.names), ncol(X), dimnames = list(row.names, colnames(X)))
  Y[rownames(X),] <- X
  return(Y)
}



#' Filters out cells that belong to control genotype and other genotypes
#'
#' @param design.mat Design matrix
#' @param ctrl Control genotype
#'
#' @return Filtered design matrix
#' @export
#'
CleanDesignCtrl <- function(design.mat, ctrl) {
  ctrl.iy <- which(colnames(design.mat) == ctrl)
  design.mat <- Matrix(t(apply(design.mat, 1, function(x) {
    if (x[[ctrl.iy]] == 1 && sum(x) > 1) { x[[ctrl.iy]] <- 0 }
    return(x)
  })))
  design.mat <- design.mat[,!(grepl(":", colnames(design.mat)) & grepl(ctrl, colnames(design.mat)))]
  return(design.mat)
}
CleanDesignCtrl <- compiler::cmpfun(CleanDesignCtrl)



#' Counts the number of cells in each group overlap combination
#'
#' @param group.list.1 Named list of cell vectors for each genotype
#' @param group.list.2 Named list of cell vectors for each cluster
#' @return Matrix of cell counts that overlap the two sets of groups
#' @export
#'
GroupOverlapCounts <- function(group.list.1, group.list.2) {
  df <- matrix(0, length(group.list.1), length(group.list.2))
  rownames(df) <- names(group.list.1)
  colnames(df) <- names(group.list.2)

  for(i in 1:length(group.list.1)) {
    group.1.cells <- group.list.1[[i]]
    for(j in 1:length(group.list.2)) {
      group.2.cells <- group.list.2[[j]]
      df[i,j] <- length(intersect(group.1.cells, group.2.cells))
    }
  }
  return(df)
}
