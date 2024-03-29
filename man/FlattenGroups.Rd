% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_handling.R
\name{FlattenGroups}
\alias{FlattenGroups}
\title{Convert list of sample groups into a flat vector.
Removes all samples belonging to multiple groups}
\usage{
FlattenGroups(groups.list)
}
\arguments{
\item{groups.list}{List of sample groups}
}
\value{
A character vector of the group each sample belongs to. The vector names are the sample names.
}
\description{
Convert list of sample groups into a flat vector.
Removes all samples belonging to multiple groups
}
