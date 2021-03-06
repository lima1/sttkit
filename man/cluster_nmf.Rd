% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster.R
\name{cluster_nmf}
\alias{cluster_nmf}
\title{cluster_nmf}
\usage{
cluster_nmf(
  obj,
  rank,
  randomize = FALSE,
  variable_features = TRUE,
  max_features = NULL,
  ...
)
}
\arguments{
\item{obj}{Object, clustered by \code{\link{cluster_spatial}}.}

\item{rank}{Number of clusters}

\item{randomize}{Randomize data, useful for diagnostics and picking rank}

\item{variable_features}{If \code{TRUE}, only use variable features}

\item{max_features}{Reduce runtime by only using the top 
\code{max_features} features}

\item{...}{Additional parameters passed to the \code{nmf} function.}
}
\description{
Performes NMF clustering
}
\examples{
cluster_nmf
}
