#### Linear models ####
#' @import glmnet
#' @import matrixStats
#' @import compiler
#'
#' @useDynLib perturbLM
#' @importFrom Rcpp sourceCpp
NULL


#' Helper function for running cross validation on a sequence of lambda values
#'
#' @param x Design matrix + covariates matrix
#' @param y Expression response
#' @param groups Perturbation dictionary (in list format) or named vector
#' @param metric Metric to use for evaluation (pearson, spearman)
#' @param alpha Alpha value
#' @param family GLM family to use for elasticnet (default: mgaussian)
#' @param lambda lambda values to test
#' @param folds List of train/test splits
#'
#' @return Matrix of cross validation correlations per fold
#' @import glmnet
#'
cross_validate_lambda <- function(x,
                                  y,
                                  groups,
                                  alpha,
                                  metric,
                                  family,
                                  lambda,
                                  folds) {
  if (!is.list(groups)) {
    groups <- UnflattenGroups(groups)
  }

  fold_cor <- lapply(folds, function(split.ix) {
    train.ix <- split.ix == "train"
    test.ix <- split.ix == "test"

    y.test <- as.matrix(y[test.ix,])
    y.test <- y.test[,colSums(y.test) > 0]
    y.test <- apply(y.test, 2, scale)
    rownames(y.test) <- rownames(y)[test.ix]

    groups.sub <- lapply(groups, function(cells) intersect(cells, rownames(y.test)))
    groups.sub <- groups.sub[sapply(groups.sub, length) > 1]
    y.test.avg <- lapply(groups.sub, function(ix) colMeans(y.test[ix,]))
    y.test.avg <- do.call(rbind, y.test.avg)

    fit_train <- glmnet(x[train.ix,], y[train.ix, colnames(y.test)], family = family, alpha = alpha, lambda = lambda)

    sapply(lambda, function(s) {
      y.pred <- predict(fit_train, newx = x[test.ix, ], s = s)[,,1]
      y.pred <- apply(y.pred, 2, scale)
      rownames(y.pred) <- rownames(y)[test.ix]
      y.pred[is.na(y.pred) | is.nan(y.pred)] <- 0

      y.pred.avg <- lapply(groups.sub, function(ix) colMeans(y.pred[ix,]))
      y.pred.avg <- do.call(rbind, y.pred.avg)

      pred.cors <- sapply(1:nrow(y.pred.avg), function(i) cor(y.pred.avg[i,], y.test.avg[i,], method = metric))
      pred.cors[is.na(pred.cors)] <- 0

      mean(pred.cors)
    })
  })

  do.call(rbind, fold_cor)
}



#' Run cross validation for both alpha and lambda
#'
#' @param x Design matrix + covariates matrix
#' @param y Expression response
#' @param groups Perturbation dictionary (in list format) or named vector
#' @param alpha.seq Sequence of alpha values to test
#' @param plot Plot cross validation results (default: True)
#' @param n.cores Number of cores to use (default: 4)
#' @param family GLM family to use for elasticnet (default: mgaussian)
#' @param nlambda Number of lambda values to test (default: 10)
#' @param lambda.min.ratio Sets the minimum lambda value to test (default: 1e-6)
#' @param nfolds Number of folds to cross validate over (default: 5)
#'
#' @return List of results containing cross validation objects (cv.list),
#'         cross validation summary stats (cv.summary),
#'         and the optimal lambda/alpha combo (alpha.min, lambda.min)
#' @import snow
#' @import glmnet
#' @import ggplot2
#' @export
#'
RunCrossValidation <- function(x,
                               y,
                               groups,
                               metric = "pearson",
                               alpha.seq = seq(0, 1, by = 0.2),
                               plot = T,
                               n.cores = 4,
                               family = "mgaussian",
                               nlambda = 10,
                               lambda.min.ratio = 0.01,
                               nfolds = 4,
                               seed = NULL) {
  require(snow)

  valid_families <- c("mgaussian")
  if (!family %in% valid_families) {
    stop(paste0("Only ", paste(valid_families, collapse = ", "), " families are currently implemented"))
  }

  ## Get folds
  if (!is.list(groups)) {
    groups <- UnflattenGroups(groups)
  }
  groups <- lapply(groups, function(cells) intersect(cells, rownames(x)))
  folds <- SplitFoldsByGroup(groups, nfolds, seed = seed)

  cl <- snow::makeCluster(n.cores, type = "SOCK")
  snow::clusterExport(cl, c("x", "y", "groups", "metric", "family", "nlambda", "lambda.min.ratio", "folds"),
                      envir = environment())
  invisible(snow::parLapply(cl, 1:n.cores, function(i) {
    library(glmnet)
    library(Matrix)
  }))
  cv_list <- parLapply(cl, alpha.seq, function(a) {
    ## Get lambda sequence
    fit_full <- glmnet(x, y, family = family, alpha = a, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio)
    lambda <- fit_full$lambda

    fold_cor <- cross_validate_lambda(x = x,
                                      y = y,
                                      groups = groups,
                                      alpha = a,
                                      metric = metric,
                                      family = family,
                                      lambda = lambda,
                                      folds = folds)

    mean.cor <- apply(fold_cor, 2, function(x) mean(na.omit(x)))
    data.frame(lambda = lambda, log.lambda = log(lambda), alpha = a, R = mean.cor)
  })
  on.exit(snow::stopCluster(cl))


  cv_df <- do.call(rbind, cv_list)
  cv_df$alpha_fac <- as.factor(cv_df$alpha)
  rownames(cv_df) <- NULL

  if (plot) {
    print(ggplot(data=cv_df, aes(x=log.lambda, y=R, color=alpha_fac)) +
            geom_line() +
            geom_point() +
            theme_classic())
  }

  alpha.min <- cv_df[which.max(cv_df$R), "alpha"]
  lambda.min <- cv_df[which.max(cv_df$R), "lambda"]
  print(paste0("Best alpha: ", alpha.min))
  print(paste0("Best lambda: ", lambda.min))

  list(cv.summary = cv_df,
       alpha.min = alpha.min,
       lambda.min = lambda.min)
}




#' Run per-gene linear mixed models using glmmLasso
#'
#' @param model_df Dataframe containing covariates used in the model formulas
#' @param Y Gene expression matrix
#' @param fixed_form Model formula for fixed effects
#' @param rand_form Model formula for random effects
#' @param lambda Lasso regularization parameter (default: 0)
#' @param fam Distribution that describes gene expression (default: gaussian())
#' @param n_cores Number of cores to use for parallel processing (default: 8)
#' @param verbose Print error statements (default: False)
#'
#' @import glmmLasso
#' @return List of fitted linear mixed model objects (one per gene)
#' @export
#'
RunMixedLM <- function(model_df,
                       Y,
                       fixed_form,
                       rand_form,
                       lambda = 0,
                       fam = gaussian(),
                       n_cores = 8,
                       batch.size = 20,
                       verbose = F) {
  require(snow)

  model_df$condition <- droplevels(model_df$condition)
  gene.batch <- split(colnames(Y), ceiling(seq_along(colnames(Y))/batch.size))
  Y.batch <- lapply(gene.batch, function(g) Y[,g])
  rm(Y); invisible(gc());

  cl <- makeCluster(n_cores)
  clusterExport(cl, c("model_df", "fixed_form", "rand_form", "lambda", "fam"),
                envir = environment())
  parLapply(cl, 1:length(cl), function(i) {
    require(glmmLasso)
  })
  fit_list <- parLapply(cl, Y.batch, function(y) {
    return(apply(as.matrix(y), 2, function(x) {
      model_df$gene <- x
      tryCatch({
        fit <- glmmLasso(fix = fixed_form,
                         rnd = rand_form,
                         data = model_df,
                         lambda = lambda,
                         family = fam)
        fit$W <- fit$X <- fit$data <- fit$fitted.values <- fit$y <- NULL
        return(fit)
      }, error = function(e) {
        if (verbose) print(e)
        return(NULL)
      })
    }))
  })
  stopCluster(cl)
  names(fit_list) <- NULL
  fit_list <- unlist(fit_list, F, T)

  if (verbose) {
    print(paste0(sum(sapply(fit_list, is.null)), " genes had errors"))
  }

  fit_list <- fit_list[!sapply(fit_list, is.null)]
  return(fit_list)
}


#' Predict expression given a new model dataframe
#'
#' @param fit_list List of fitted glmmLasso objects
#' @param new_model_df New model dataframe to predict on (must contain same covariate columns as original dataframe)
#'
#' @return Matrix of predicted values
#' @export
#'
PredictMixedLM <- function(fit_list, new_model_df) {
  new_model_df$gene <- 0
  pred_list <- lapply(fit_list, function(fit) {
    as.numeric(predict(fit, newdata = new_model_df))
  })
  pred_mat <- do.call(rbind, pred_list)
  colnames(pred_mat) <- rownames(new_model_df)

  return(pred_mat)
}


#' Extract coefficient matrix from a multigaussian glmnet object
#'
#' @param mfit Glmnet object (multigaussian)
#'
#' @return Matrix of regression coefficients
#' @export
#'
GetCoefMatrix <- function(mfit) {
  cfs <- t(sapply(mfit$beta, function(x) {
    y <- as.numeric(x)
    names(y) <- rownames(x)
    return(y)
  }))
  colnames(cfs) <- gsub('genotype', '', colnames(cfs))
  return(cfs)
}
GetCoefMatrix <- compiler::cmpfun(GetCoefMatrix)


#' Calculate cluster enrichment with logistic regression
#'
#' @param design.mat Design matrix for regression
#' @param cov.mat Matrix of covariates to regress out
#' @param clusters Cluster assignments
#' @param alpha Elasticnet ratio: (0 is fully L2, 1 fully is L1)
#' @param lambda.use Shrinkage parameter
#' @param ctrl Control genotype to compare all other genotypes to (e.g. NTC, AAVS)
#'
#' @return Matrix of regression coefficients
#' @export
#'
CalcClusterEnrich <- function(design.mat, cov.mat, clusters, alpha, lambda.use, ctrl = NULL) {
  clusters <- as.factor(clusters)
  if(is.null(ctrl)) {
    design.mat.full <- as.matrix(cbind(design.mat, cov.mat))
    fit <- glmnet(design.mat.full, as.factor(clusters), family = "multinomial",
                  alpha = alpha, lambda = lambda.use)
    return(GetCoefMatrix(fit))
  } else {
    genotypes <- colnames(design.mat)[colnames(design.mat) != ctrl]
    cfs <- vapply(genotypes, function(g) {
      if (grepl(":", g)) {
        g1 <- strsplit(g, split = ":")[[1]][[1]]
        g2 <- strsplit(g, split = ":")[[1]][[2]]
        ix <- rowSums(design.mat[,c(g, g1, g2, ctrl)]) > 0
        design.mat.full <- as.matrix(cbind(design.mat[ix, c(g, g1, g2, ctrl)], cov.mat[ix,]))
        colnames(design.mat.full) <- c(g, g1, g2, ctrl, colnames(cov.mat))
        fit <- glmnet(design.mat.full, as.factor(clusters)[rownames(design.mat.full)], family = "multinomial",
                      alpha = alpha, lambda = lambda.use)
      } else {
        ix <- rowSums(design.mat[,c(g, ctrl)]) > 0
        design.mat.full <- as.matrix(cbind(design.mat[ix, c(g, ctrl)], cov.mat[ix,]))
        colnames(design.mat.full) <- c(g, ctrl, colnames(cov.mat))
        fit <- glmnet(design.mat.full, as.factor(clusters)[rownames(design.mat.full)], family = "multinomial",
                      alpha = alpha, lambda = lambda.use)
      }
      return(GetCoefMatrix(fit)[,g])
    }, rep(0, nlevels(clusters)))
  }
}
CalcClusterEnrich <- compiler::cmpfun(CalcClusterEnrich)



#' Calculate cluster enrichment with logistic regression
#'
#' @param design.mat Design matrix for regression
#' @param cov.mat Matrix of covariates to regress out
#' @param clusters Cluster assignments
#' @param alpha Elasticnet ratio: (0 is fully L2, 1 fully is L1)
#' @param lambda.use Shrinkage parameter
#' @param n.rand Number of permutations
#' @param n.cores Number of multiprocessing cores
#' @param ctrl Control genotype to compare all other genotypes to (e.g. NTC, AAVS)
#'
#' @return Matrix of regression coefficients
#' @import snow
#' @import glmnet
#' @import abind
#'
#' @export
#'
CalcClusterEnrichPvals <- function(design.mat, cov.mat, clusters, alpha, lambda.use,
                                   n.rand = 1000, n.cores = 8, ctrl = NULL) {
  design.mat <- as.matrix(design.mat)

  if (!is.null(cov.mat)) {
    cov.mat <- as.matrix(cov.mat);
    stopifnot(all(rownames(design.mat) == rownames(cov.mat)))
  }
  design.mat.full <- cbind(design.mat, cov.mat)

  cfs <- CalcClusterEnrich(design.mat, cov.mat, clusters, alpha, lambda.use, ctrl = ctrl)

  cl <- snow::makeCluster(n.cores, type = "SOCK")
  snow::clusterExport(cl, c("design.mat", "cov.mat", "clusters", "alpha", "lambda.use", "ctrl",
                            "GetCoefMatrix"), envir = environment())
  source.log <- snow::parLapply(cl, 1:n.cores, function(i) library(glmnet))
  cfs.rand <- snow::parLapply(cl, 1:n.rand, function(i) {
    design.mat.permute <- design.mat[sample(1:nrow(design.mat)),]
    rownames(design.mat.permute) <- rownames(design.mat)
    CalcClusterEnrich(design.mat.permute, cov.mat, clusters, alpha, lambda.use, ctrl = ctrl)
  })
  snow::stopCluster(cl)
  cfs.rand <- abind::abind(cfs.rand, along = 3)

  return(.calc_pvals_cfs(cfs, cfs.rand))
}
CalcClusterEnrichPvals <- compiler::cmpfun(CalcClusterEnrichPvals)



#' Calculate penalized linear regression with glmnet
#'
#' @param design.matrix Design matrix for regression
#' @param metadata Additional covariates to take into account (will not be permuted)
#' @param y Linear model response
#' @param alpha Elasticnet ratio: (0 is fully L2, 1 fully is L1)
#' @param lambda.use Coefficient regularization parameter
#' @param family Regression family to use
#' @param ctrl Control variable in the design matrix (optional)
#'
#' @return Matrix of regression coefficients
#' @export
#'
CalcGlmnet <- function(design.matrix, metadata, y, alpha, lambda.use, family, ctrl = NULL) {
  if (is.null(ctrl)) {
    x <- Matrix(cbind(design.matrix, metadata))
    mfit <- glmnet::glmnet(x, y = y, family = family, alpha = alpha, lambda = lambda.use,
                           standardize = F)
    cfs <- GetCoefMatrix(mfit)
    cfs <- cfs[,colnames(design.matrix)]
  } else {
    genotypes <- colnames(design.matrix)[colnames(design.matrix) != ctrl]
    cfs <- vapply(genotypes, function(g) {
      if (grepl(":", g)) {
        g1 <- strsplit(g, split = ":")[[1]][[1]]
        g2 <- strsplit(g, split = ":")[[1]][[2]]
        ix <- Matrix::rowSums(design.matrix[,c(g, g1, g2, ctrl)]) > 0
        x <- Matrix(cbind(design.matrix[ix, c(g, g1, g2)], metadata[ix,]))
        colnames(x) <- c(g, g1, g2, colnames(metadata))
        mfit <- glmnet::glmnet(x, y = y[ix,], family = family, alpha = alpha, lambda = lambda.use,
                               standardize = F)
      } else {
        ix <- Matrix::rowSums(design.matrix[,c(g, ctrl)]) > 0
        x <- Matrix(cbind(design.matrix[ix, g], metadata[ix,]))
        colnames(x) <- c(g, colnames(metadata))
        mfit <- glmnet::glmnet(x, y = y[ix,], family = family, alpha = alpha, lambda = lambda.use,
                               standardize = F)
      }
      return(GetCoefMatrix(mfit)[,g])
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
#' @param lambda.use Coefficient regularization parameter
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
CalcGlmnetPvals <- function(design.matrix, metadata, y, alpha, lambda.use, family,
                            ctrl = NULL, n.rand = 20, n.cores = 16, n.bins = 10,
                            use.quantiles = T, output.name = "cf") {

  metadata <- as.matrix(metadata)
  cfs <- CalcGlmnet(design.matrix, metadata, y, alpha, lambda.use, family, ctrl)

  cl <- snow::makeCluster(n.cores, type = "SOCK")
  snow::clusterExport(cl, c("design.matrix", "y", "metadata", "alpha", "lambda.use",
                            "family", "GetCoefMatrix"), envir = environment())
  source.log <- snow::parLapply(cl, 1:n.cores, function(i) library(glmnet))
  cfs.rand <- snow::parLapply(cl, 1:n.rand, function(i) {
    design.matrix.permute <- design.matrix[sample(1:nrow(design.matrix)),]
    CalcGlmnet(design.matrix.permute, metadata, y, alpha, lambda.use, family, ctrl)
  })
  snow::stopCluster(cl)

  y <- as.matrix(y)
  gene.covs <- data.frame(gene_avgs = colMeans(y), gene_vars = matrixStats::colVars(y))
  null.coefs.df <- lapply(cfs.rand, .flatten_score_matrix, output.name = output.name,
                          gene.covs = gene.covs, group.covs = NULL)
  null.coefs.df <- rbindlist(null.coefs.df)
  null.coefs.df <- .get_multi_bins(null.coefs.df, colnames(gene.covs), n.bins, use.quantiles, bin.group = T)
  null.binned.dat <- .get_binned_list(null.coefs.df, output.name)
  rm(null.coefs.df); invisible(gc());

  coefs.df <- .flatten_score_matrix(cfs, output.name, gene.covs, NULL)
  coefs.df <- .get_multi_bins(coefs.df, colnames(gene.covs), n.bins, use.quantiles, bin.group = T)
  coefs.df <- .calc_emp_pvals(coefs.df, null.binned.dat, output.name = output.name, n.cores = round(n.cores/2))
  coefs.df <- coefs.df[order(coefs.df$p_val),]

  coefs.df[c("gene_avgs", "gene_vars", "bin.index")] <- NULL
  return(coefs.df)
}
CalcGlmnetPvals <- compiler::cmpfun(CalcGlmnetPvals)



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
      b <- cfs[i,j]
      if (b > 0) {
        p.vals[i,j] <- calcPvalGreaterCpp(v, b)
      } else if (b < 0) {
        p.vals[i,j] <- -1*calcPvalLessCpp(v, b)
      } else {
        p.vals[i,j] <- 1
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
      bin.cov.list[[bin.cov]] <- factor(lsr::quantileCut(scores.df[[cov]], n = n.bins[[cov]],
                                                         labels = F, include.lowest = T))
    } else {
      bin.cov.list[[bin.cov]] <- factor(cut(scores.df[[cov]], breaks = n.bins[[cov]],
                                            labels = F, include.lowest = T))
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

