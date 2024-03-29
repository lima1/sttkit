% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/giotto.R
\name{find_giotto_dwls_matrix}
\alias{find_giotto_dwls_matrix}
\title{find_giotto_dwls_matrix}
\usage{
find_giotto_dwls_matrix(
  singlecell_giotto,
  refdata,
  method = "scran",
  expression_values = "normalized",
  num_markers = 100
)
}
\arguments{
\item{singlecell_giotto}{List of Giotto objects}

\item{refdata}{Column with cell type information}

\item{method}{Name of the method used to find markers}

\item{expression_values}{Giotto expression values slot}

\item{num_markers}{Number of significant markers per cell type}
}
\description{
find_giotto_dwls_matrix
}
\examples{
#
}
