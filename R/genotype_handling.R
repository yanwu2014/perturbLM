#### Functions for handling genotype/phenotype dictionaries

#' Convert genotypes dictionary from list to named vector. All cells with more than one genotype are removed
#'
#' @param genotypes.list Named list of cell vectors for each genotype
#' @return Named vector with the values as genotypes and the names as cell names.
#' @export
#'
FlattenGenotypeList <- function(genotypes.list) {
  cell.names <- unlist(genotypes.list, F, F)
  if (length(cell.names) > length(unique(cell.names))) {
    print("Warning: removing all cells with more than one genotype")
    genotypes.list <- .get_single_genotypes(genotypes.list)
    cell.names <- unlist(genotypes.list, F, F)
  }

  genotypes <- rep(NA, length(cell.names))
  names(genotypes) <- cell.names
  for (g in names(genotypes.list)) {
    genotypes[genotypes.list[[g]]] <- g
  }

  if (any(is.na(genotypes))) { stop("Some unassigned cells"); }
  return(genotypes)
}
FlattenGenotypeList <- compiler::cmpfun(FlattenGenotypeList)


.get_single_genotypes <- function(genotypes.list, min.cells = 1) {
  genotypes.matrix <- DesignMatrixGenotypes(genotypes.list, max.genotypes = 1, min.cells = min.cells)
  genotypes.list <- lapply(colnames(genotypes.matrix), function(g) {
    rownames(genotypes.matrix)[which(genotypes.matrix[,g] == 1)]
  })
  names(genotypes.list) <- colnames(genotypes.matrix)
  return(genotypes.list)
}
.get_single_genotypes <- compiler::cmpfun(.get_single_genotypes)



#' Convert genotypes from named vector to list
#'
#' @param cell.genotypes Named vector of cell genotypes
#' @param min.cells Minimum cells for a genotype to be included in the list
#' @return Named list of cell vectors for each genotype
#'
#' @export
#'
UnflattenCellGenotypes <- function(cell.genotypes, min.cells = 1) {
  genotypes.list <- c()
  genotypes <- unique(cell.genotypes)
  for (g in genotypes) {
    g.cells <- names(cell.genotypes[cell.genotypes == g])
    if (length(g.cells) > min.cells) {
      genotypes.list[[g]] <- g.cells
    }
  }
  return(genotypes.list)
}
UnflattenCellGenotypes <- compiler::cmpfun(UnflattenCellGenotypes)


#' Read in genotypes dictionary from csv file. Returns an R list.
#'
#' @param pheno.dict.file Genotype dictionary file in csv format
#' @param sep.char Value separating columns in genotype dictionary file
#' @return Named list of cell vectors for each genotype
#'
#' @export
#'
ReadGenotypes <- function(pheno.dict.file, sep.char = ",") {
  geno.data <- read.table(pheno.dict.file, sep = sep.char, header = F, stringsAsFactors = F)
  genotypes <- geno.data[[1]]
  genotypes.list <- lapply(rownames(geno.data), function(i) sapply(strsplit(geno.data[i,2], split = ',')[[1]], trimws))
  genotypes.list <- lapply(genotypes.list, function(x) make.names(x))
  names(genotypes.list) <- genotypes
  return(genotypes.list)
}
ReadGenotypes <- compiler::cmpfun(ReadGenotypes)


#' Write genotypes to csv file
#'
#' @param genotypes.list Named list of cell vectors for each genotype
#' @param out.file Output file name
#' @return Writes genotypes.list to out.file
#' @export
#'
WriteGenotypes <- function(genotypes.list, out.file) {
  geno.data <- sapply(genotypes.list, function(x) paste('\"', paste(x, collapse = ", "), '\"', sep = ""))
  geno.data <- sapply(names(geno.data), function(x) paste(x, geno.data[[x]], sep = ","))
  fileConn <- file(out.file, "w")
  writeLines(geno.data, fileConn)
  close(fileConn)
}
WriteGenotypes <- compiler::cmpfun(WriteGenotypes)


#' Convert genotypes list to design matrix
#'
#' @param genotypes.list Named list of cell vectors for each genotype
#' @param max.genotypes Maximum number of genotypes per cell
#' @param min.cells Minumum number of cells per genotype
#' @param drop.cells Drop cells in genotypes that are filtered out
#'
#' @return Sparse binary matrix where rows are cells, and columns are genotypes. 1 means the cell received the genotype.
#'
#' @import Matrix
#' @import hash
#' @export
#'
DesignMatrixGenotypes <- function(genotypes.list, max.genotypes = 2, min.cells = 5, drop.cells = T) {
  cell.names <- unique(unlist(genotypes.list, F, F))
  single.genotypes <- names(genotypes.list)

  single.mat <- Matrix::Matrix(0, length(cell.names), length(single.genotypes), dimnames = list(cell.names, single.genotypes))
  for (i in 1:ncol(single.mat)) {
    single.mat[genotypes.list[[i]], i] <- 1
  }
  if (drop.cells) single.mat <- single.mat[Matrix::rowSums(single.mat) <= max.genotypes, ]
  n_cells_single <- Matrix::colSums(single.mat)
  single.mat <- single.mat[ ,names(n_cells_single[n_cells_single >= min.cells])]
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

    # Filter out combos that have too few cells
    n_cells_combo <- Matrix::colSums(combo.mat)
    combos_keep <- names(n_cells_combo[n_cells_combo >= min.cells])

    # Remove filtered combo cells from design matrix
    if (length(combos_keep) > 0) {
      combo.mat <- combo.mat[ ,combos_keep]
      design.mat <- cbind(single.mat, combo.mat)
      filtered.cells <- (Matrix::rowSums(single.mat) > 1) & (Matrix::rowSums(combo.mat) == 0)

    } else {
      design.mat <- single.mat
      filtered.cells <- (Matrix::rowSums(single.mat) > 1)

    }
    if (drop.cells) design.mat <- design.mat[!filtered.cells,]
  } else {
    design.mat <- single.mat
  }
  return(Matrix(design.mat))
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
  design.mat <- Matrix(apply(design.mat, 1, function(x) {
    if (x[[ctrl.iy]] == 1 && sum(x) > 1) { x[[ctrl.iy]] <- 0 }
    return(x)
  }))
  return(Matrix::t(design.mat))
}
CleanDesignCtrl <- compiler::cmpfun(CleanDesignCtrl)


#' Counts the number of cells in each genotype x cluster combination
#'
#' @param genotypes.list Named list of cell vectors for each genotype
#' @param clusters.list Named list of cell vectors for each cluster
#' @return Matrix of genotypes and cluster counts
#' @export
#'
GenotypeClusterCounts <- function(genotypes.list, clusters.list) {
  df <- matrix(0, length(genotypes.list), length(clusters.list))
  rownames(df) <- names(genotypes.list)
  colnames(df) <- names(clusters.list)

  for(i in 1:length(genotypes.list)) {
    genotype.cells <- genotypes.list[[i]]
    for(j in 1:length(clusters.list)) {
      cluster.cells <- clusters.list[[j]]
      df[i,j] <- length(intersect(genotype.cells, cluster.cells))
    }
  }
  return(df)
}
GenotypeClusterCounts <- compiler::cmpfun(GenotypeClusterCounts)
