#### Linear models ####
#' @import glmnet
#' @import matrixStats
#'
#' @useDynLib perturbLM
#' @importFrom Rcpp sourceCpp
NULL


#' Extract coefficient matrix from a multigaussian glmnet object
#'
#' @param mfit Glmnet object (multigaussian)
#' @param best.lambda Optimal lambda regularization parameter
#'
#' @return Matrix of regression coefficients
#' @export
#'
GetCoefMatrix <- function(mfit, best.lambda) {
  cfs <- coef(mfit, s = best.lambda)
  cfs <- lapply(cfs, function(x) {y <- as.numeric(x); names(y) <- rownames(x); return(y)})
  cfs <- do.call('rbind', cfs)
  colnames(cfs) <- gsub('genotype', '', colnames(cfs))
  return(cfs)
}


#' Calculate regression coefficients and p-values via permutation testing
#'
#' @param design.matrix Design matrix for regression
#' @param metadata Additional covariates to take into account (will not be permuted)
#' @param y Linear model response
#' @param alpha Elasticnet ratio: (0 is fully L2, 1 fully is L1)
#' @param lambda Coefficient regularization parameter
#' @param family Regression family to use
#' @param n.rand Number of permutations for calculating coefficient significance
#' @param n.cores Number of cores to use
#' @param n.bins If binning genes, number of bins to use
#' @param use.quantiles If binning genes, whether or not to bin by quantile
#' @param output.name Column name of regression coefficient
#'
#' @return If binning genes, returns a dataframe with Gene, Perturbation, Regression Coefficient (cf), and P-value.
#'         Otherwise returns a list of regression coefficients and p-values
#' @import snow
#' @import lsr
#' @import qvalue
#' @importFrom data.table rbindlist
#' @export
#'
CalcGlmnetPvals <- function(design.matrix, metadata, y, alpha, lambda, family, n.rand, n.cores,
                            n.bins = 10, use.quantiles = T, output.name = "cf") {
  mfit <- glmnet::glmnet(Matrix(cbind(design.matrix, metadata)), y = y, family = family, alpha = alpha, lambda = lambda,
                         standardize = F)
  cfs <- GetCoefMatrix(mfit, lambda)

  cols.keep <- colnames(cfs)[2:ncol(cfs)]
  cols.keep <- cols.keep[!cols.keep %in% colnames(metadata)]

  cl <- snow::makeCluster(n.cores, type = "SOCK")
  snow:: clusterExport(cl, c("design.matrix", "y", "metadata", "alpha", "lambda", "family", "GetCoefMatrix"),
                       envir = environment())
  source.log <- snow::parLapply(cl, 1:n.cores, function(i) library(glmnet))
  cfs.rand <- snow::parLapply(cl, 1:n.rand, function(i) {
    design.matrix.permute <- design.matrix[sample(1:nrow(design.matrix)),]
    mfit <- glmnet::glmnet(Matrix(cbind(design.matrix.permute, metadata)), y = y, family = family, alpha = alpha, lambda = lambda,
                           standardize = F)
    GetCoefMatrix(mfit, lambda)
  })

  gene_vars <- matrixStats::colVars(y)
  gene_avgs <- colMeans(y)
  gene.covs <- data.frame(gene_avgs, gene_vars)

  cfs <- cfs[,cols.keep]
  cfs.rand <- lapply(cfs.rand, function(x) x[,cols.keep])

  null.coefs.df <- lapply(cfs.rand, .flatten_score_matrix, output.name = output.name,
                          gene.covs = gene.covs, group.covs = NULL)
  null.coefs.df <- rbindlist(null.coefs.df)
  null.coefs.df <- .get_multi_bins(null.coefs.df, colnames(gene.covs), n.bins, use.quantiles, bin.group = T)
  null.binned.dat <- .get_binned_list(null.coefs.df, output.name)
  rm(null.coefs.df)
  invisible(gc())

  coefs.df <- .flatten_score_matrix(cfs, output.name, gene.covs, NULL)
  coefs.df <- .get_multi_bins(coefs.df, colnames(gene.covs), n.bins, use.quantiles, bin.group = T)
  coefs.df <- .calc_emp_pvals(coefs.df, null.binned.dat, output.name = 'cf', n.cores = round(n.cores/4))
  coefs.df <- coefs.df[order(coefs.df$p_val),]

  coefs.df[c("gene_avgs", "gene_vars", "bin.index")] <- NULL
  return(coefs.df)
}

#### Permutation testing and P-value matrix manipulation ####

## Given a matrix of regression coefficients, and a 3D array of permuted coefficients
## calculate empirical p-values
.calc_pvals_cfs <- function(cfs, cfs.rand) {
  p.vals <- matrix(1, nrow(cfs), ncol(cfs))
  colnames(p.vals) <- colnames(cfs)
  rownames(p.vals) <- rownames(cfs)
  for (i in 1:nrow(cfs)) {
    for (j in 1:ncol(cfs)) {
      v <- sort(na.omit(cfs.rand[i,j,]))
      v.pos <- v[v > 0]
      v.neg <- v[v <= 0]
      b <- cfs[i,j]
      if (b > 0) {
        p.vals[i,j] <- calcPvalGreaterCpp(v.pos, b)
      } else {
        p.vals[i,j] <- calcPvalLessCpp(v.neg, b)
      }
    }
  }
  return(p.vals)
}


## Helper function for p-value calculation.
.flatten_score_matrix <- function(score.matrix, output.name, gene.covs = NULL, group.covs = NULL) {
  scores.df <- as.data.frame.table(score.matrix, stringsAsFactors = F, responseName = output.name)
  colnames(scores.df) <- c('Gene','Group', output.name)
  num.cfs <- nrow(scores.df)

  if (!is.null(gene.covs)) {
    for (cov in colnames(gene.covs)) {
      scores.df[[cov]] <- numeric(num.cfs)
    }
  }

  if (!is.null(group.covs)) {
    for (cov in colnames(group.covs)) {
      scores.df[[cov]] <- numeric(num.cfs)
    }
  }

  st <- 1
  en <- nrow(score.matrix)
  for (i in 1:ncol(score.matrix)) {
    for (cov in colnames(gene.covs)) {
      scores.df[st:en, cov] <- gene.covs[[cov]]
    }
    for (cov in colnames(group.covs)) {
      scores.df[st:en, cov] <- group.covs[i,cov]
    }

    st <- st + nrow(score.matrix)
    en <- en + nrow(score.matrix)
  }

  return(scores.df)
}


## Helper function for p-value calculation.
.get_multi_bins <- function(scores.df, all.covs, n.bins, use.quantiles = F, bin.direction = F, bin.group = F) {
  if (length(n.bins) == 1) {
    n.bins <- rep(n.bins, length(all.covs))
  }
  names(n.bins) <- all.covs

  bin.cov.list <- c()
  for (cov in all.covs) {
    bin.cov <- paste(cov, 'binned', sep = '_')
    if (use.quantiles) {
      bin.cov.list[[bin.cov]] <- factor(lsr::quantileCut(scores.df[[cov]], n = n.bins[[cov]], labels = F, include.lowest = T))
    } else {
      bin.cov.list[[bin.cov]] <- factor(cut(scores.df[[cov]], breaks = n.bins[[cov]], labels = F, include.lowest = T))
    }
  }

  if (bin.direction) {
    bin.cov.list[['Direction']] <- factor(scores.df$Direction)
  }

  if (bin.group) {
    bin.cov.list[['Group']] <- factor(scores.df$Group)
  }

  bin <- interaction(bin.cov.list, drop = T, sep = '_')
  bin.names <- levels(bin)

  bin.idx <- 1:length(bin.names)
  names(bin.idx) <- bin.names

  bin <- as.character(bin)

  scores.df$bin.index <- vapply(bin, function(x) bin.idx[[x]], 1)
  return(scores.df)
}


## Helper function for p-value calculation.
.get_binned_list <- function(scores.df, output.name) {
  bin.names <- sort(unique(scores.df$bin.index))

  binned.dat <- lapply(bin.names, function(x) {
    scores <- scores.df[scores.df$bin.index == x,];
    if (nrow(scores) > 0) {
      return(sort(omitNaCpp(scores[[output.name]]), decreasing = F))
    } else {
      return(c())
    }
  })
  return(binned.dat)
}


## Helper function for p-value calculation.
.fast_pvals <- function(r, binned.dat.up, binned.dat.dn) {
  bin.idx <- r[[1]]
  x <- r[[2]]
  if (x >= 0) {
    v <- binned.dat.up[[bin.idx]]
    p.val <- calcPvalGreaterCpp(v,x)
  } else {
    v <- binned.dat.dn[[bin.idx]]
    p.val <- calcPvalLessCpp(v,x)
  }
  return(p.val)
}


## Helper function for p-value calculation.
.calc_emp_pvals <- function(scores.df, binned.dat, output.name, n.cores = 1, direction = c("both", "lower"),
                           correct = "BH") {
  scores.mat <- as.matrix(scores.df[c('bin.index', output.name)])

  binned.dat.up <- lapply(binned.dat, function(v) v[v >= 0])
  binned.dat.dn <- lapply(binned.dat, function(v) v[v <= 0])
  rm(binned.dat)

  if (n.cores > 1) {
    cl <- makeCluster(n.cores, type = 'SOCK')
    snow::clusterExport(cl, c('binned.dat.up', 'binned.dat.dn'), envir = environment())
    print('Starting P-Value Calculations')
    coefs.pvals <- snow::parApply(cl, scores.mat, 1, .fast_pvals, binned.dat.up = binned.dat.up,
                                  binned.dat.dn = binned.dat.dn)
    snow::stopCluster(cl)
  } else {
    print('Starting P-Value Calculations')
    coefs.pvals <- apply(scores.mat, 1, .fast_pvals, binned.dat.up = binned.dat.up,
                         binned.dat.dn = binned.dat.dn)
  }

  scores.df$p_val <- coefs.pvals
  if (correct == "qvalue") {
    scores.df$FDR <- qvalue::qvalue(coefs.pvals)$qvalues
  } else {
    scores.df$FDR <- p.adjust(coefs.pvals, method = "BH")
  }
  scores.df <- scores.df[order(scores.df$p_val),]
  return(scores.df)
}
