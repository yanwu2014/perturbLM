#### Linear models ####
#' @import glmnet
#' @import matrixStats
#' @import compiler
#'
#' @useDynLib perturbLM
#' @importFrom Rcpp sourceCpp
NULL

#' Set up design matrices and response given a Seurat object
#'
#' @param obj Seurat object
#' @param pert.col Metadata column containing perturbation info (required)
#' @param batch.col Metadata column containing batch info (default: NULL)
#' @param size.factor.col Metadata column containing library size info (default: NULL)
#' @param mito.col Metadata column containing mitochondrial fraction info (default: NULL)
#' @param response.col Metadata column containing the desired response.
#'                     If NULL, the normalized data in the Seurat object will be used as the response.
#' @param features.use If response.col is NULL, a subset of features to use for the response
#' @param test.prop Proportion of cells held out for model evaluation
#'
#' @return List of x (perturbation + covariate design matrix)
#'         y (response vector or matrix), and family (response family)
#' @export
#'
SetupLinearModel <- function(obj,
                             pert.col,
                             batch.col = NULL,
                             size.factor.col = NULL,
                             mito.col = NULL,
                             response.col = NULL,
                             features.use = NULL,
                             test.prop = 0.2,
                             seed = 12345) {

  if (!require("Seurat")) {
    stop("Install Seurat before running this function")
  }

  ## Set up design matrix
  perts <- obj@meta.data[[pert.col]]
  names(perts) <- colnames(obj)
  pert.mat <- CreateDesignMatrix(perts)

  if (is.null(batch.col)) {
    cov.mat <- Matrix(1, nrow = nrow(pert.mat), ncol = 1,
                      dimnames = c(rownames(pert.mat), c("batch")))
  }
  batch <- obj@meta.data[[batch.col]]
  names(batch) <- colnames(obj)
  cov.mat <- CreateDesignMatrix(as.factor(batch[rownames(pert.mat)]))

  if (!is.null(size.factor.col)) {
    size.factor <- scale(obj@meta.data[rownames(pert.mat), size.factor.col], center = F)
    cov.mat <- cbind(cov.mat, size.factor)
    colnames(cov.mat)[ncol(cov.mat)] <- "size.factor"
  }

  if (!is.null(mito.col)) {
    percent.mt <- scale(obj@meta.data[rownames(pert.mat), mito.col], center = F)
    cov.mat <- cbind(cov.mat, percent.mt)
    colnames(cov.mat)[ncol(cov.mat)] <- "percent.mt"
  }

  x <- cbind(pert.mat, cov.mat)

  ## Set up response
  if (is.null(response.col)) {
    print("Assuming response is normalized expression data ")
    y <- t(GetAssayData(obj, slot = "data"))[rownames(pert.mat),]

    if (!is.null(features.use)) {
      y <- y[,features.use]
    }

    fam <- "mgaussian"
  } else {
    y <- obj@meta.data[[response.col]]
    names(y) <- colnames(obj)
    y <- y[rownames(pert.mat)]

    if (is.factor(y) | is.character(y)) {
      print("Assuming response is categorical")
      y <- as.factor(y)
      fam <- "multinomial"
    } else if (is.numeric(y)) {
      print("Assuming response is numeric")
      fam <- "gaussian"
    } else {
      stop("Response must be either a factor, character, numeric, or NULL")
    }
  }

  ## Get train/test split
  if (test.prop < 0.01 | test.prop > 1/3) {
    stop("test.prop must be between 0.01 and 1/3")
  }
  nfolds <- round(1/test.prop)
  split <- SplitFoldsByGroup(perts, nfolds, seed = seed)[[1]]
  train.cells <- names(split[split == "train"])
  test.cells <- names(split[split == "test"])

  ## Return design matrices & response
  list(x = x,
       y = y,
       perts = perts,
       train.cells = train.cells,
       test.cells = test.cells,
       family = fam)
}


#' Helper function for running cross validation on a sequence of lambda values for the mgaussian family
#'
#' @param x Design matrix + covariates matrix
#' @param y Expression response
#' @param groups Perturbation dictionary (in list format) or named vector
#' @param metric Metric to use for evaluation (pearson, spearman)
#' @param alpha Alpha value
#' @param family GLM family to use for elasticnet (default: mgaussian)
#' @param lambda lambda values to test
#' @param folds List of train/test splits
#' @param seq.lambda.pred Predict expression at each lambda value sequentially to save memory
#' @return Matrix of cross validation correlations per fold
#' @import glmnet
#'
cross_validate_lambda <- function(x,
                                  y,
                                  groups,
                                  alpha,
                                  lambda,
                                  folds,
                                  seq.lambda.pred) {
  if (!is.list(groups)) {
    groups <- UnflattenGroups(groups)
  }

  fold_cor <- lapply(folds, function(split.ix) {
    train.ix <- split.ix == "train"
    test.ix <- split.ix == "test"

    y.test <- as.matrix(y[test.ix,])
    y.test <- y.test[,colSums(y.test) > 0]

    fit_train <- glmnet(x[train.ix,], y[train.ix, colnames(y.test)], family = "mgaussian", alpha = alpha, lambda = lambda, standardize = F)


    if (seq.lambda.pred) {
      lambda.cors <- sapply(lambda, function(s) {
        y.pred <- predict(fit_train, newx = x[test.ix, ], s = s)[,,1]
        return(mean((y.pred - y.test)^2))
      })
    } else {
      y.pred.list <- predict(fit_train, newx = x[test.ix,], s = lambda)
      y.pred.list <- lapply(1:dim(y.pred.list)[[3]], function(i) y.pred.list[,,i])
      lambda.cors <- sapply(y.pred.list, function(y.pred) {
        mean((y.pred - y.test)^2)
      })
      names(lambda.cors) <- lambda
    }

    lambda.cors
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
#' @param lambda.min.ratio Sets the minimum lambda value to test (default: 0.01)
#' @param nfolds Number of folds to cross validate over (default: 5)
#' @param seed Random seed for fold reproducibility
#' @param seq.lambda.pred Predict expression at each lambda value sequentially to save memory (default: F)
#'
#' @return List of results containing cross validation objects (cv.list),
#'         cross validation summary stats (cv.summary),
#'         and the optimal lambda/alpha combo (alpha, lambda)
#' @import snow
#' @import glmnet
#' @import ggplot2
#' @export
#'
RunCrossValidation <- function(x,
                               y,
                               groups,
                               alpha.seq = c(0.1, 0.4, 0.7, 0.9, 0.95, 1),
                               plot = T,
                               n.cores = 4,
                               family = "mgaussian",
                               nlambda = 10,
                               lambda.min.ratio = 0.01,
                               nfolds = 4,
                               seed = NULL,
                               seq.lambda.pred = F) {
  require(snow)

  if ((family == "mgaussian") && (length(dim(y)) != 2)) {
    stop("Response y must be a 2D matrix when fitting mgaussian model")
  }

  if ((family == "multinomial") && (!is.factor(y))) {
    stop("Response y must be a factor when fitting multinomial model")
  }

  ## Get folds
  if (!is.list(groups)) {
    groups <- UnflattenGroups(groups)
  }
  groups <- lapply(groups, function(cells) intersect(cells, rownames(x)))
  folds <- SplitFoldsByGroup(groups, nfolds, seed = seed)

  vars.export <- c("x", "y", "groups", "family", "nlambda", "lambda.min.ratio", "folds", "seq.lambda.pred")
  cl <- snow::makeCluster(n.cores, type = "SOCK")
  snow::clusterExport(cl, vars.export,
                      envir = environment())
  invisible(snow::parLapply(cl, 1:n.cores, function(i) {
    library(glmnet)
    library(Matrix)
  }))
  cv_list <- parLapply(cl, alpha.seq, function(a) {
    if (family == "mgaussian") {
      ## Get lambda sequence
      fit_full <- glmnet(x, y, family = family, alpha = a, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, standardize = F)
      lambda <- fit_full$lambda

      fold_cor <- cross_validate_lambda(x = x,
                                        y = y,
                                        groups = groups,
                                        alpha = a,
                                        lambda = lambda,
                                        folds = folds,
                                        seq.lambda.pred = seq.lambda.pred)

      mean.cor <- apply(fold_cor, 2, function(x) mean(na.omit(x)))
      return(data.frame(lambda = lambda, log.lambda = log(lambda), alpha = a, err = mean.cor))
    } else {
      cv <- cv.glmnet(x, y, family = family, alpha = a, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, standardize = F)
      return(data.frame(lambda = cv$lambda, log.lambda = log(cv$lambda), alpha = a, err = cv$cvm))
    }
  })
  on.exit(snow::stopCluster(cl))



  cv_df <- do.call(rbind, cv_list)
  cv_df$alpha_fac <- as.factor(cv_df$alpha)
  rownames(cv_df) <- NULL

  if (plot) {
    print(ggplot(data=cv_df, aes(x=log.lambda, y=err, color=alpha_fac)) +
            geom_line() +
            geom_point() +
            theme_classic())
  }
  alpha.use <- cv_df[which.min(cv_df$err), "alpha"]
  lambda.use <- cv_df[which.min(cv_df$err), "lambda"]

  print(paste0("Best alpha: ", alpha.use))
  print(paste0("Best lambda: ", lambda.use))

  list(cv.summary = cv_df,
       alpha = alpha.use,
       lambda = lambda.use)
}


#' Helper function for computing cross entropy
cross_entropy <- function(y, yhat) {
  if ((ncol(y) != ncol(yhat)) || (nrow(y) != nrow(yhat))) {
    stop("y and yhat must have equal dimensions")
  }

  sum(sapply(1:ncol(y), function(i) {
    -1*sum(y[,i]*log(yhat[,i] + 1e-10))
  }))
}


#' Function for evaluating linear model performance
#'
#' @param mdl.fit Linear model fit on training data
#' @param x.test Test design matrix
#' @param y.test Test response
#' @param groups Perturbation dictionary (in list format) or named vector
#' @param eval.metric Model evaluation metric
#' @param do.scale Whether to scale response before evaluation (only applies to pearson & spearman metrics)
#' @param min.cells Minimum cells per group to include in evaluation
#'
#' @export
#'
EvaluateLinearModel <- function(mdl.fit,
                                x.test,
                                y.test,
                                groups,
                                eval.metric,
                                do.scale = T,
                                min.cells = 3) {

  if (!is.list(groups)) {
    groups <- UnflattenGroups(groups)
  }
  groups <- lapply(groups, function(cells) intersect(cells, rownames(x.test)))

  ncells <- sapply(groups, length)
  groups.keep <- names(ncells[ncells >= min.cells])
  groups <- groups[groups.keep]

  if ((length(dim(y.test)) == 2) && (!eval.metric %in% c("pearson", "spearman"))) {
    stop("eval.metric must be either pearson, spearman for multi-gaussian response")
  }

  if ((length(dim(y.test)) == 1) && (eval.metric %in% c("pearson", "spearman"))) {
    stop("response must be multi-gaussian (2D) to use pearson or spearman metrics")
  }

  if (is.factor(y.test) && (eval.metric != "cross-entropy")) {
    stop("eval.metric must be cross entropy for multi-class response")
  }

  y.pred <- predict(mdl.fit, newx = x.test, type = "response")

  if (eval.metric %in% c("pearson", "spearman")) {
    y.pred <- y.pred[,,1]
    cell.names <- rownames(y.test)
    if (do.scale) {
      y.pred <- apply(y.pred, 2, scale)
      y.pred[is.na(y.pred) | is.nan(y.pred)] <- 0
      y.test <- apply(y.test, 2, scale)
      y.test[is.na(y.test) | is.nan(y.test)] <- 0
      rownames(y.pred) <- rownames(y.test) <- cell.names
    }

    pred.avg <- sapply(groups, function(cells) Matrix::colMeans(y.pred[cells,]))
    true.avg <- sapply(groups, function(cells) Matrix::colMeans(y.test[cells,]))

    groups.eval <- sapply(names(groups), function(i) {
      cor(pred.avg[,i], true.avg[,i], method = eval.metric)
    })
    names(groups.eval) <- names(groups)

  } else if (eval.metric == "cross-entropy") {
    y.pred <- y.pred[,,1]
    y_1hot <- model.matrix(~0 + y.test)
    colnames(y_1hot) <- gsub("y.test", "", colnames(y_1hot))
    rownames(y_1hot) <- names(y.test)

    y_1hot <- y_1hot[,colnames(y.pred)]

    groups.eval <- sapply(groups, function(cells) {
      cross_entropy(y.pred[cells,], y_1hot[cells,])
    })
  } else if (eval.metric == "euclidean") {
    y.pred <- y.pred[,1]
    groups.eval <- sapply(groups, function(cells) {
      sum((y.pred[cells] - y.test[cells])^2)
    })
  }

  df <- data.frame(group = names(groups.eval),
                   ncells = ncells[names(groups.eval)])
  df[[eval.metric]] <- groups.eval

  return(df)
}


#' Wrapper function for running a ElasticNet model on a single cell dataset
#'
#' @param obj Seurat object
#' @param pert.col Metadata column containing perturbation info (required)
#' @param batch.col Metadata column containing batch info (default: NULL)
#' @param size.factor.col Metadata column containing library size info (default: NULL)
#' @param mito.col Metadata column containing mitochondrial fraction info (default: NULL)
#' @param response.col Metadata column containing the desired response.
#'                     If NULL, the normalized data in the Seurat object will be used as the response.
#' @param features.use If response.col is NULL, a subset of features to use for the response
#' @param test.prop Proportion of cells held out for model evaluation
#' @param family GLM family to use for elasticnet (default: mgaussian)
#' @param alpha.seq Sequence of alpha values to test
#' @param plot.cv Plot cross validation results (default: True)
#' @param n.cores Number of cores to use (default: 4)
#' @param nlambda Number of lambda values to test (default: 10)
#' @param lambda.min.ratio Sets the minimum lambda value to test (default: 0.01)
#' @param nfolds Number of folds to cross validate over (default: 5)
#' @param seed Random seed for fold reproducibility
#' @param seq.lambda.pred Predict expression at each lambda value sequentially to save memory (default: F)
#' @param eval.metric Model evaluation metric
#'
#' @return Returns a list with the design matrix, response, train/test split
#'         cross validation results, fitted model, model coefficients, and evaluation metrics
#' @export
#'
RunLinearModel <- function(obj,
                           pert.col,
                           batch.col = NULL,
                           size.factor.col = NULL,
                           mito.col = NULL,
                           response.col = NULL,
                           features.use = NULL,
                           test.prop = 0.2,
                           family = NULL,
                           alpha.seq = c(0.1, 0.4, 0.7, 0.9, 0.95, 1),
                           plot.cv = T,
                           n.cores = 4,
                           nlambda = 10,
                           lambda.min.ratio = 0.01,
                           nfolds = 4,
                           seed = NULL,
                           seq.lambda.pred = F,
                           eval.metric = NULL) {
  ## Setup model
  mdl.list <- SetupLinearModel(obj,
                               pert.col = pert.col,
                               response.col = response.col,
                               batch.col = batch.col,
                               size.factor.col = size.factor.col,
                               mito.col = mito.col,
                               features.use = features.use,
                               test.prop = test.prop)

  train.ix <- mdl.list$train.cells
  test.ix <- mdl.list$test.cells
  perts <- as.factor(mdl.list$perts)

  if (is.null(family)) {
    family <- mdl.list$family
  }

  ## Optimize hyperparameters on training set
  if (mdl.list$family == "mgaussian") {
    y.train <- mdl.list$y[train.ix,]
    y.test <- mdl.list$y[test.ix,]
  } else {
    y.train <- mdl.list$y[train.ix]
    y.test <- mdl.list$y[test.ix]
  }


  cv <- RunCrossValidation(mdl.list$x[train.ix,],
                           y.train,
                           groups = perts[train.ix],
                           family = family,
                           nfolds = nfolds,
                           nlambda = nlambda,
                           n.cores = n.cores,
                           seed = seed,
                           plot = plot.cv,
                           lambda.min.ratio = lambda.min.ratio,
                           seq.lambda.pred = seq.lambda.pred)
  mdl.list$cv <- cv

  ## Run full model on training set
  fit <- glmnet(mdl.list$x[train.ix,],
                y.train,
                family = family,
                alpha = cv$alpha,
                lambda = cv$lambda,
                standardize = F)
  mdl.list$model <- fit
  mdl.list$coefs <- GetCoefMatrix(fit)

  ## Run covariates only model on training set
  cov.vars <- colnames(mdl.list$x)[!colnames(mdl.list$x) %in% levels(perts)]
  cov.fit <- glmnet(mdl.list$x[train.ix, cov.vars],
                    y.train,
                    family = family,
                    alpha = cv$alpha,
                    lambda = cv$lambda,
                    standardize = F)

  ## Evaluate model
  if (is.null(eval.metric)) {
    if (family == "mgaussian") {
      eval.metric <- "pearson"
    } else if (family == "multinomial") {
      eval.metric <- "cross-entropy"
    } else {
      eval.metric <- "euclidean"
    }
  }

  mdl.eval <- EvaluateLinearModel(fit,
                                  x.test = mdl.list$x[test.ix,],
                                  y.test = y.test,
                                  groups = perts,
                                  eval.metric = eval.metric)
  cov.eval <- EvaluateLinearModel(cov.fit,
                                  x.test = mdl.list$x[test.ix, cov.vars],
                                  y.test = y.test,
                                  groups = perts,
                                  eval.metric = eval.metric)

  mdl.eval[[paste0(eval.metric, "_reduced")]] <- cov.eval[[eval.metric]]
  if (family == "mgaussian") {
    mdl.eval[["rel_performance"]] <- log2((mdl.eval[[eval.metric]] + 1e-10)/(cov.eval[[eval.metric]] + 1e-10))
  } else {
    mdl.eval[["rel_performance"]] <- log2((cov.eval[[eval.metric]] + 1e-10)/(mdl.eval[[eval.metric]] + 1e-10))
  }
  mdl.eval <- mdl.eval[order(mdl.eval$rel_performance, decreasing = T),]
  mdl.list$evaluation <- mdl.eval

  return(mdl.list)
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



