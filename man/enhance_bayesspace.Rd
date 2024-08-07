% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enhance.R
\name{enhance_bayesspace}
\alias{enhance_bayesspace}
\title{enhance_bayesspace}
\usage{
enhance_bayesspace(x, test_num_clusters = seq(2, 12), num_clusters = NULL, ...)
}
\arguments{
\item{x}{Object, converted by \code{\link{as_SingleCellExperiment}}}

\item{test_num_clusters}{Range of number of clusters to be tested (\code{qs}).}

\item{num_clusters}{Optional picked number}

\item{...}{Additional paramters passed to \code{BayesSpace::spatialEnhance}}
}
\description{
Enhance BayesSpace data
}
\examples{
#enhance_bayesspace()
}
