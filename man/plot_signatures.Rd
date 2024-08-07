% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_signatures}
\alias{plot_signatures}
\title{plot_signatures}
\usage{
plot_signatures(
  obj_spatial,
  file,
  gmt,
  nbin = 24,
  ctrl = 30,
  method = c("seurat", "mean"),
  width = 10,
  pdf = FALSE,
  png = FALSE,
  cells = NULL,
  zero_cutoff = NULL,
  assay = "Spatial",
  ...
)
}
\arguments{
\item{obj_spatial}{Seurat object with SpatialTranscriptomics data}

\item{file}{Output PDF}

\item{gmt}{GMT file with gene signatures or gene signature read by
\code{\link{read_signatures}}}

\item{nbin}{Argument of \code{Seurat::AddModuleScore}}

\item{ctrl}{Argument of \code{Seurat::AddModuleScore}}

\item{method}{Either use \code{Seurat::AddModuleScore} or simple mean}

\item{width}{Plot width}

\item{pdf}{Create PDF image files}

\item{png}{Create PNG image files}

\item{cells}{Plot only specified cells}

\item{zero_cutoff}{Cutoff defining zero. Defaults to half the number of genes in the signature.}

\item{assay}{Name of the assay corresponding to the initial input data.}

\item{...}{Arguments passed to \code{Seurat::SpatialPlot}.}
}
\description{
Plots signatures defined in a GMT file and applied to a Seurat object on HE slides
}
\examples{
plot_signatures()
}
