% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integrate.R
\name{integrate_spatial}
\alias{integrate_spatial}
\title{integrate_spatial}
\usage{
integrate_spatial(
  obj_spatial = NULL,
  references,
  features = 2000,
  min_features = 400,
  min_spots = 200,
  scale = TRUE,
  max_percent_mito = 33,
  min_max_counts = 3,
  min_detected = 2,
  reference_technology = "single_cell",
  skip_alternative_batch_corrections = c(),
  force,
  serialize = TRUE,
  prefix,
  verbose = FALSE
)
}
\arguments{
\item{obj_spatial}{Seurat object containing SpatialTranscriptomics data.
Can be \code{NULL}; in that case only references will be integrated.}

\item{references}{List of Seurat reference datasets}

\item{features}{Features to be used for integration}

\item{min_features}{Remove cells or spots with small number of detected genes}

\item{min_spots}{Remove or merge small samples with fewer spots}

\item{scale}{Scale features}

\item{max_percent_mito}{Remove cells or spots with high mitochondria expression.}

\item{min_max_counts}{Remove cells or spots with low maximum counts.}

\item{min_detected}{Remove features detected in fewer than \code{min_detected}
\code{references} (including \code{obj_spatial}).}

\item{reference_technology}{Used to fill the \code{meta.data@technology} field}

\item{skip_alternative_batch_corrections}{By default, harmony and fastmnn will
be run when installed. Provide methods to skip here as array.}

\item{force}{Recalculate, even when serialized objects are available}

\item{serialize}{Serialize output objects}

\item{prefix}{Prefix of output files}

\item{verbose}{Verbose Seurat output}
}
\description{
Integrate SpatialTranscriptomics data with (matched) scRNA-seq data
}
\examples{
integrate_spatial()
}
