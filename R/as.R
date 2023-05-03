
#' as_SingleCellExperiment
#'
#' Convert Seurat Visium object for BayesSpace. Currently only works for single
#' sample objects.
#' @param object Seurat Visium object
#' @param ... Arguments passed to \code{spatialPreprocess}
#' @export as_SingleCellExperiment
#' @examples
#' # as_SingleCellExperiment(object)
as_SingleCellExperiment <- function(object, ...) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop("Install SingleCellExperiment.")
    }
    if (length(Images(object)) > 1) {
        stop("Currently only single sample objects supported.")
    }
    x <- as.SingleCellExperiment(object)
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("Install SummarizedExperiment.")
    }
    SummarizedExperiment::colData(x) <- cbind(
        SummarizedExperiment::colData(x),
        object@images[[1]]@coordinates[rownames(SummarizedExperiment::colData(x)), ]
    )
    if (requireNamespace("BayesSpace", quietly = TRUE)) {
        x <- BayesSpace::spatialPreprocess(x, ...)
    }
    return(x)
}

#' as_SpatialRNA
#'
#' Convert Seurat Visium object for spacexr

#' @param object Seurat Visium object
#' @param assay Seurat assay in \code{object}
#' @param slot Seurat slot in \code{object} and \code{assay}
#' @param ... Arguments passed to \code{spacexr::SpatialRNA}
#' @export as_SpatialRNA
#' @examples
#' # as_SpatialRNA(object)
as_SpatialRNA <- function(object, assay = "Spatial", slot = "counts",  ...) {
    if (!requireNamespace("spacexr", quietly = TRUE)) {
        stop("Install spacexr.")
    }
    if (length(Images(object)) > 1) {
        stop("Currently only single sample objects supported.")
    }
    counts <- GetAssayData(object, assay = assay, slot = slot)
    coords <- GetTissueCoordinates(object)

    return(spacexr::SpatialRNA(coords, counts, ...))
}

#' as_Reference
#'
#' Convert Seurat single cell object for spacexr

#' @param object Seurat object
#' @param refdata Column with cell type annotation in \code{object}
#' @param assay Seurat assay in \code{object}
#' @param slot Seurat slot in \code{object} and \code{assay}
#' @param ... Arguments passed to \code{spacexr::Reference}
#' @export as_Reference
#' @examples
#' # as_Reference(object)
as_Reference <- function(object, refdata, assay = "RNA", slot = "counts",  ...) {
    if (!requireNamespace("spacexr", quietly = TRUE)) {
        stop("Install spacexr.")
    }
    counts <- GetAssayData(object, assay = assay, slot = slot)
    cell_types_df <- FetchData(object, vars = refdata)
    cell_types <- as.factor(cell_types_df[, 1])
    names(cell_types) <- rownames(cell_types_df)

    return(spacexr::Reference(counts = counts, cell_types = cell_types, ...))
}
