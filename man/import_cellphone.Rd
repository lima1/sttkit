% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellphonedb.R
\name{import_cellphone}
\alias{import_cellphone}
\title{Import cell-cell interaction analysis using CellPhoneDB}
\usage{
import_cellphone(cellphone_outpath, orgdb, prefix, slot = "data", assay = NULL)
}
\arguments{
\item{orgdb}{orgdb Object}

\item{prefix}{Prefix of output files}

\item{slot}{Specific information to pull
(i.e. counts, data, scale.data,...)}

\item{assay}{Name of assay to pull data from}

\item{obj}{Seurat Object}
}
\value{
interaction pairs in GMT format
}
\description{
This function allows you to read CellPhoneDB output and rerank it based on local interactions.
}
\examples{

}