% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_handling.R
\name{CleanControls}
\alias{CleanControls}
\title{Remove all cells in the control perturbation(s) that also received a non-control perturbation}
\usage{
CleanControls(group.list, ctrl.groups)
}
\arguments{
\item{group.list}{Perturbation dictionary (in list format)}

\item{ctrl.groups}{Control perturbation(s)}
}
\value{
Filtered perturbation dictionary
}
\description{
Remove all cells in the control perturbation(s) that also received a non-control perturbation
}
