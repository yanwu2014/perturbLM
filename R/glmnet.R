#### Linear models ####
#' @import glmnet
#' @import matrixStats
#' @import compiler
#'
#' @useDynLib perturbLM
#' @importFrom Rcpp sourceCpp
NULL


#' Calculate log-FC for each genotype against the control
#'
#' @param X Design matrix
#' @param Y Gene expression
#' @param ctrl Control genotype
#'
#' @return log-FC matrix
#' @export
#'
CalcLogFC <- function(X, Y, ctrl, eps = 1e-3) {
  mean.exprs <- apply(X, 2, function(ix) Matrix::colMeans(Y[ix == 1,]))
  log.fc <- apply(mean.exprs, 2, function(x) log2((x + eps)/(mean.exprs[,ctrl] + eps)))
  return(log.fc)
}


#' Extract coefficient matrix from a multigaussian glmnet object
#'
#' @param mfit Glmnet object (multigaussian)
#' @param best.lambda Optimal lambda regularization parameter
#'
#' @return Matrix of regression coefficients
#' @export
#'
GetCoefMatrix <- function(mfit, lambda) {
  cfs <- glmnet::coef.glmnet(mfit, s = lambda)
  cfs <- sapply(cfs, function(x) {
    y <- as.numeric(x)
    names(y) <- rownames(x)
    return(y)
  })
  colnames(cfs) <- gsub('genotype', '', colnames(cfs))
  return(cfs)
}
GetCoefMatrix <- compiler::cmpfun(GetCoefMatrix)


#' Calculate penalized linear regression with glmnet
#'
#' @param design.matrix Design matrix for regression
#' @param metadata Additional covariates to take into account (will not be permuted)
#' @param y Linear model response
#' @param alpha Elasticnet ratio: (0 is fully L2, 1 fully is L1)
#' @param lambda Coefficient regularization parameter
#' @param lambda.seq Regularization sequence to follow
#' @param family Regression family to use
#' @param ctrl Control variable in the design matrix (optional)
#'
#' @return Matrix of regression coefficients
#' @export
#'
CalcGlmnet <- function(design.matrix, metadata, y, alpha, lambda, lambda.seq, family, ctrl = NULL) {
  if (is.null(ctrl)) {
    x <- as.matrix(cbind(design.matrix, metadata))
    mfit <- glmnet::glmnet(x, y = y, family = family, alpha = alpha, lambda = lambda.seq,
                           standardize = F)
    cfs <- GetCoefMatrix(mfit, lambda)
    cfs <- cfs[,colnames(design.matrix)]
  } else {
    stopifnot(ctrl %in% colnames(design.matrix))
    genotypes <- colnames(design.matrix)[colnames(design.matrix) != ctrl]
    cfs <- vapply(genotypes, function(g) {
      ix <- Matrix::rowSums(design.matrix[,c(g,ctrl)]) > 0
      x <- as.matrix(cbind(design.matrix[ix, g], metadata[ix,]))
      mfit <- glmnet::glmnet(x, y = y[ix,], family = family, alpha = alpha, lambda = lambda.seq,
                             standardize = F)
      return(GetCoefMatrix(mfit, lambda)[,2])
    }, rep(0, ncol(y)))
    colnames(cfs) <- genotypes
  }

  return(cfs)
}
CalcGlmnet <- compiler::cmpfun(CalcGlmnet)


#' Calculate regression coefficients and p-values via permutation testing
#'
#' @param design.matrix Design matrix for regression
#' @param metadata Additional covariates to take into account (will not be permuted)
#' @param y Linear model response
#' @param alpha Elasticnet ratio: (0 is fully L2, 1 fully is L1)
#' @param lambda Coefficient regularization parameter
#' @param lambda.seq Regularization sequence to follow
#' @param family Regression family to use
#' @param ctrl Control variable to compare against (optional)
#' @param n.rand Number of permutations for calculating coefficient significance
#' @param n.cores Number of cores to use
#' @param n.bins If binning genes, number of bins to use
#' @param use.quantiles If binning genes, whether or not to bin by quantile
#' @param output.name Column name of regression coefficient
#'
#' @return Returns a dataframe with Gene, Perturbation, log-FC, p-value
#' @import snow
#' @import lsr
#' @import qvalue
#' @importFrom data.table rbindlist
#' @export
#'
CalcGlmnetPvals <- function(design.matrix, metadata, y, alpha, lambda, lambda.seq, family,
                            ctrl = NULL, n.rand = 20, n.cores = 16, n.bins = 10, use.quantiles = T,
                            output.name = "cf") {

  metadata <- as.matrix(metadata)
  cfs <- CalcGlmnet(design.matrix, metadata, y, alpha, lambda, lambda.seq, family, ctrl)

  cl <- snow::makeCluster(n.cores, type = "SOCK")
  snow::clusterExport(cl, c("design.matrix", "y", "metadata", "alpha", "lambda", "lambda.seq", "family",
                            "GetCoefMatrix"),
                      envir = environment())
  source.log <- snow::parLapply(cl, 1:n.cores, function(i) library(glmnet))
  cfs.rand <- snow::parLapply(cl, 1:n.rand, function(i) {
    design.matrix.permute <- design.matrix[sample(1:nrow(design.matrix)),]
    CalcGlmnet(design.matrix.permute, metadata, y, alpha, lambda, lambda.seq, family, ctrl)
  })
  snow::stopCluster(cl)

  gene.covs <- data.frame(gene_avgs = colMeans(y), gene_vars = matrixStats::colVars(y))
  null.coefs.df <- lapply(cfs.rand, .flatten_score_matrix, output.name = output.name,
                          gene.covs = gene.covs, group.covs = NULL)
  null.coefs.df <- rbindlist(null.coefs.df)
  null.coefs.df <- .get_multi_bins(null.coefs.df, colnames(gene.covs), n.bins, use.quantiles, bin.group = T)
  null.binned.dat <- .get_binned_list(null.coefs.df, output.name)
  rm(null.coefs.df); invisible(gc());

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
.calc_pvals_cfs <- compiler::cmpfun(.calc_pvals_cfs)


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
.flatten_score_matrix <- compiler::cmpfun(.flatten_score_matrix)


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
.get_multi_bins <- compiler::cmpfun(.get_multi_bins)


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
.get_binned_list <- compiler::cmpfun(.get_binned_list)


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
.fast_pvals <- compiler::cmpfun(.fast_pvals)


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
.calc_emp_pvals <- compiler::cmpfun(.calc_emp_pvals)

