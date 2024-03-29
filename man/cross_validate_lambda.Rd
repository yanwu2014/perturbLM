% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glmnet.R
\name{cross_validate_lambda}
\alias{cross_validate_lambda}
\title{Helper function for running cross validation on a sequence of lambda values for the mgaussian family}
\usage{
cross_validate_lambda(x, y, groups, alpha, lambda, folds, seq.lambda.pred)
}
\arguments{
\item{x}{Design matrix + covariates matrix}

\item{y}{Expression response}

\item{groups}{Perturbation dictionary (in list format) or named vector}

\item{alpha}{Alpha value}

\item{lambda}{lambda values to test}

\item{folds}{List of train/test splits}

\item{seq.lambda.pred}{Predict expression at each lambda value sequentially to save memory}

\item{metric}{Metric to use for evaluation (pearson, spearman)}

\item{family}{GLM family to use for elasticnet (default: mgaussian)}
}
\value{
Matrix of cross validation correlations per fold
}
\description{
Helper function for running cross validation on a sequence of lambda values for the mgaussian family
}
