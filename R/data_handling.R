#### Functions for handling perturbation dictionaries and other input data

#' Extract a field from a delimited string
#'
#' @param string Input string
#' @param field Field to extract
#' @param delim Character delimiter (default: "_")
#'
#' @export
#'
ExtractField <- function (string, field = 1, delim = "_") {
  fields <- as.numeric(unlist(strsplit(x = as.character(x = field), split = ",")))
  if (length(fields) == 1) {
    return(strsplit(string, split = delim)[[1]][field])
  }
  return(paste(strsplit(string, split = delim)[[1]][fields], collapse = delim))
}


#' Convert list of sample groups into a flat vector.
#' Removes all samples belonging to multiple groups
#'
#' @param groups.list List of sample groups
#'
#' @return A character vector of the group each sample belongs to. The vector names are the sample names.
#'
#' @export
#'
FlattenGroups <- function(groups.list) {
  cell.names <- unlist(groups.list, F, F)
  if (length(cell.names) > length(unique(cell.names))) {
    print("Warning: removing all samples belonging to multiple groups")
    cell.tbl <- table(cell.names)
    cells.keep <- names(cell.tbl[cell.tbl == 1])
  } else {
    cells.keep <- cell.names
  }

  groups <- rep(NA, length(cells.keep))
  names(groups) <- cells.keep
  for (g in names(groups.list)) {
    g.cells <- intersect(groups.list[[g]], cells.keep)
    groups[g.cells] <- g
  }

  if (any(is.na(groups))) { stop("Some unassigned cells"); }
  return(groups)
}


#' Converts a flat sample groups character vector into a list format
#'
#' @param groups Character vector of sample groups
#'
#' @return List of groups: each list element contains the samples for that group
#'
#' @export
#'
UnflattenGroups <- function(groups) {
  groups.list <- lapply(unique(groups), function(g){
    names(groups[groups == g])
  })
  names(groups.list) <- unique(groups)
  return(groups.list)
}


#' Stratify groups by additional factor
#'
#' @param group.list Perturbation dictionary (in list format)
#' @param stratify Factor to stratify by
#' @param sep Character delimiter to add to the stratified perturbation names (default: "_")
#'
#' @return Perturbation dictionary stratified by factor
#' @export
StratifyGroupList <- function(group.list, stratify, sep = '_') {
  stratify <- as.factor(stratify)
  group.list.stratified <- c()
  for (g in names(group.list)) {
    g.cells <- group.list[[g]]
    for (l in levels(stratify)) {
      l.cells <- names(stratify[stratify == l])
      k <- paste(g, l, sep = '_')
      group.list.stratified[[k]] <- intersect(g.cells, l.cells)
    }
  }
  return(group.list.stratified)
}


#' Merge perturbations with the same pattern in names.
#' For example, we might want to collapse all guide RNAs into the gene they
#' target (i.e. TULP2-1, TULP2-2, TULP2-3 into TULP2)
#'
#' @param group.list Perturbation dictionary (in list format)
#' @param pattern Character pattern to collapse (i.e. TULP2)
#' @param new.group.name Name of new perturbation (default: pattern)

#' @return List of perts with all perturbations containing the pattern merged
#' @export
MergeGroups <- function(group.list, pattern, new.group.name = NULL) {
  if (is.null(new.group.name)) {
    new.group.name <- pattern
  }

  guides <- names(group.list)[grepl(pattern, names(group.list))]
  guides.cells <- unique(unlist(group.list[guides], F, F))

  group.list[guides] <- NULL
  group.list[[new.group.name]] <- guides.cells

  return(group.list)
}


#' Merge two perturbation dicts
#'
#' @param group.list.1 1st Perturbation dictionary to merge
#' @param group.list.2 1st Perturbation dictionary to merge
#'
#' @return Merged perturbation dictionary
#' @export
MergeGroupLists <- function(group.list.1, group.list.2) {
  common.groups <- intersect(names(group.list.1), names(group.list.2))
  group.list <- c()
  for (p in common.groups) {
    group.list[[p]] <- union(group.list.1[[p]], group.list.2[[p]])
  }
  unique.groups.1 <- names(group.list.1)[!names(group.list.1) %in% common.groups]
  group.list[unique.groups.1] <- group.list.1[unique.groups.1]

  unique.groups.2 <- names(group.list.2)[!names(group.list.2) %in% common.groups]
  group.list[unique.groups.2] <- group.list.2[unique.groups.2]

  return(group.list)
}


#' Convert genotypes list to design matrix.
#' DesignMatrixGenotypes is deprecated, use CreateDesignMatrix.
#'
#' @export
#'
DesignMatrixGenotypes <- function(genotypes.list, max.genotypes = 2, min.cells = 10) {
  .Deprecated("DesignMatrixGenotypes", package="perturbLM",
              msg = "DesignMatrixGenotypes is deprecated, use CreateDesignMatrix",
              old = as.character(sys.call(sys.parent()))[1L])
  CreateDesignMatrix(genotypes.list, max.genotypes, min.cells)
}


#' Convert perturbation dictionary to design matrix
#'
#' @param group.list Named list of cell vectors for each perturbation
#' @param max.groups Maximum number of perturbations per cell
#' @param min.cells Minumum number of cells per perturbation, including combo perturbations
#'
#' @return Sparse binary matrix where rows are cells, and columns are perturbation. 1 means the cell received the perturbation.
#'
#' @import Matrix
#' @import hash
#' @export
#'
CreateDesignMatrix <- function(group.list, max.groups = 3, min.cells = 10) {
  cell.names <- unique(unlist(group.list, F, F))
  single.groups <- names(group.list)

  single.mat <- Matrix::Matrix(0, length(cell.names), length(single.groups),
                               dimnames = list(cell.names, single.groups))
  for (i in 1:ncol(single.mat)) {
    single.mat[group.list[[i]], i] <- 1
  }
  single.mat <- single.mat[Matrix::rowSums(single.mat) > 0, ]

  combo.rows <- single.mat[Matrix::rowSums(single.mat) > 1,]
  combo.groups <- unique(apply(combo.rows, 1, function(x) {
    paste(names(x[x == 1]), collapse = ":", sep = "")
  }))
  if (max.groups > 1 && length(combo.groups) > 0) {
    combo.group.list <- hash::hash(keys = combo.groups)
    for (i in 1:nrow(combo.rows)) {
      mat.row <- combo.rows[i,]
      genotype <- paste(names(mat.row[mat.row == 1]), collapse = ":", sep = "")
      cell <- rownames(combo.rows)[[i]]
      combo.group.list[[genotype]] <- c(combo.group.list[[genotype]], cell)
    }
    combo.group.list <- as.list(combo.group.list)
    combo.group.list[['keys']] <- NULL
    combo.group.list <- combo.group.list[combo.groups]

    combo.mat <- Matrix(0, nrow(single.mat), length(combo.groups),
                        dimnames = list(rownames(single.mat), combo.groups))
    for (i in 1:ncol(combo.mat)) {
      combo.mat[combo.group.list[[i]],i] <- 1
    }
    design.mat <- cbind(single.mat, combo.mat)

  } else {
    design.mat <- single.mat
  }

  design.mat <- design.mat[,Matrix::colSums(design.mat) > min.cells]
  return(as(design.mat, "dgCMatrix"))
}


#' Filters out cells that belong to control genotype and other genotypes
#' CleanDesignCtrl is deprecated, use CleanControls instead
#'
#' @export
#'
CleanDesignCtrl <- function(design.mat, ctrl) {
  .Deprecated("CleanDesignCtrl", package="perturbLM",
              msg = "CleanDesignCtrl is deprecated, use CleanControls instead",
              old = as.character(sys.call(sys.parent()))[1L])

  ctrl.iy <- which(colnames(design.mat) == ctrl)
  design.mat <- Matrix(t(apply(design.mat, 1, function(x) {
    if (x[[ctrl.iy]] == 1 && sum(x) > 1) { x[[ctrl.iy]] <- 0 }
    return(x)
  })))
  design.mat <- design.mat[,!(grepl(":", colnames(design.mat)) & grepl(ctrl, colnames(design.mat)))]
  return(design.mat)
}
CleanDesignCtrl <- compiler::cmpfun(CleanDesignCtrl)


#' Remove all cells in the control perturbation(s) that also received a non-control perturbation
#'
#' @param group.list Perturbation dictionary
#' @param ctrl.groups Control perturbation(s)
#'
#' @return Filtered perturbation dictionary
#' @export
#'
CleanControls <- function(group.list, ctrl.groups) {
  pert.groups <- names(group.list)[!names(group.list) %in% ctrl.groups]
  pert.cells <- unique(unlist(group.list[pert.groups], F, F))
  for (g in ctrl.groups) {
    g.cells <- group.list[[g]]
    g.cells <- g.cells[!g.cells %in% pert.cells]
    group.list[[g]] <- g.cells
  }

  return(group.list)
}


#' Counts the number of cells in each group overlap combination
#'
#' @param group.list.1 Named list of cell vectors for each genotype
#' @param group.list.2 Named list of cell vectors for each cluster
#' @return Matrix of cell counts that overlap the two sets of groups
#' @export
#'
GroupOverlapCounts <- function(group.list.1, group.list.2) {
  .Deprecated("GroupClusterCounts", package="perturbLM",
              msg = "Use GroupClusterCounts instead",
              old = as.character(sys.call(sys.parent()))[1L])
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


#' Counts the number of cells in each perturbation/cluster
#'
#' @param group.list Perturbation dictionary (in list format)
#' @param clusters Cluster annotations
#'
#' @return Dataframe of cluster cell counts for each perturbation
#' @export
GroupClusterCounts <- function(group.list, clusters) {
  clusters <- factor(clusters)
  df <- matrix(0, length(group.list), nlevels(clusters))
  rownames(df) <- names(group.list)
  colnames(df) <- levels(clusters)

  for(i in 1:length(group.list)) {
    group.cells <- group.list[[i]]
    for(cl in levels(clusters)) {
      cl.cells <- names(clusters[clusters == cl])
      df[i,cl] <- length(intersect(group.cells, cl.cells))
    }
  }
  return(df)
}


#' Get cluster composition of each perturbation
#'
#' @param group.list Perturbation dictionary (in list format)
#' @param clusters Cluster annotations
#' @param samples Samples (i.e. mice, patients, biological replicates, default: NULL)
#' @param stabilize.var Apply variance stabilizing transformation (default: TRUE)
#' @param ctrl.group Control perturbation to compute logFC against (for paired samples only, default: NULL)
#'
#' @return Dataframe of cluster compositions or logFCs for each perturbation
#' @export
GroupClusterComp <- function(group.list,
                             clusters,
                             samples = NULL,
                             stabilize.var = TRUE,
                             ctrl.group = NULL)
{
  if (!is.null(samples)) {
    group.list <- StratifyGroupList(group.list, samples, sep = '_')
  }

  group_cluster_counts <- GroupClusterCounts(group.list, clusters)
  group_cluster_frac <- group_cluster_counts/rowSums(group_cluster_counts)

  if (stabilize.var) {
    group_cluster_frac <- asin(sqrt(group_cluster_frac))
  }

  if (!is.null(ctrl.group)) {
    if (is.null(samples)) {
      stop("Must specify samples to compute logFCs")
    }
    row_samples <- sapply(rownames(group_cluster_frac),
                          ExtractField,
                          delim = '_',
                          field = 2)
    logfc_list <- lapply(unique(row_samples), function(s) {
      ctrl.row <- paste0(ctrl.pert, '_', s)
      fracs <- group_cluster_frac[which(row_samples == s),]
      logfcs <- t(apply(fracs, 1, function(x) {
        log2((x + 0.05)/(fracs[ctrl.row,] + 0.05))
      }))
      logfcs <- logfcs[rownames(logfcs) != ctrl.row,]
      logfcs
    })
    group_cluster_frac <- do.call(rbind, logfc_list)
  }

  return(group_cluster_frac)
}


#' Splits data into folds stratified by a perturbation dictionary
#'
#' @param groups.list Perturbation dictionary
#' @param nfolds Number of folds
#' @param seed Set a random seed for reproducibility
#'
#' @return List of fold splits
#' @export
#'
SplitFoldsByGroup <- function(group.list, nfolds, seed = NULL) {
  require(caret)
  all.cells <- unique(unlist(group.list, F, F))
  folds.list <- lapply(group.list, function(cells) {
    set.seed(seed)
    lapply(caret::createFolds(cells, k = nfolds), function(ix) cells[ix])
  })

  lapply(1:nfolds, function(i) {
    test.cells <- unique(unlist(lapply(folds.list, function(cells.split) cells.split[[i]]), F, F))

    fold.split <- rep("train", length(all.cells))
    names(fold.split) <- all.cells

    fold.split[test.cells] <- "test"
    return(fold.split)
  })
}
