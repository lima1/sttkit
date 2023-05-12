
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

#' as_AssayObject
#'
#' Convert deconvolution output object to a TransferPrediction object

#' @param object Currently supported is RCDT and Giotto (spatEnrObj) output
#' @importFrom data.table data.table dcast rbindlist
#' @importFrom methods as
#' @export as_AssayObject
#' @examples
#' # as_AssayObject(object)
as_AssayObject <- function(object) {
    if (is(object, "RCTD")) {
        if (!requireNamespace("spacexr", quietly = TRUE)) {
            stop("Install spacexr.")
        }
        .as_AssayObject_rcdt(object)
    } else if (is(object, "spatEnrObj")) {
        if (!requireNamespace("spacexr", quietly = TRUE)) {
            stop("Install spacexr.")
        }
        .as_AssayObject_giotto(object)
    } else {
        stop("Only spacexr::RCTD and Giotto::spatEnrObj objects supported")
    }
}

.as_AssayObject_rcdt <- function(object) {
    r <- object@results
    if (length(r) > 1) {
        if (!is.null(r[[1]]$sub_weights)) {
            sw <- rbindlist(lapply(seq_along(r), function(i)
                    data.table(
                        barcode = colnames(object@spatialRNA@counts)[i],
                        cell_type = names(r[[i]]$sub_weights),
                        weight = r[[i]]$sub_weights
                  )), fill = TRUE)
            sw$cell_type[is.na(sw$cell_type)] <- "unassigned"
            swd <- data.table::dcast(sw, barcode ~ cell_type, value.var = "weight", fill = 0)
            swm <- as.matrix(swd[, -1])
            rownames(swm) <- swd$barcode
            swm <- t(spacexr::normalize_weights(swm))
            swm <- rbind(swm, max = apply(swm[!rownames(swm) %in% "unassigned", ], 2, max))
            swm <- as(swm, "sparseMatrix")
            return(CreateAssayObject(data = swm))
        }
    } else if (length(r) == 1) {
        m <- t(spacexr::normalize_weights(as.matrix(r$weights)))
        m <- rbind(m, max = apply(m, 2, max))
        return(CreateAssayObject(data = m))
    }
}

.as_AssayObject_giotto <- function(object) {
    m <- as.matrix(object@enrichDT[,-1])
    rownames(m) <- object@enrichDT$cell_ID
    m <- t(m)
    m <- rbind(m, max = apply(m, 2, max))
    return(CreateAssayObject(data = m))
}

#' as_GiottoObject
#'
#' Convert Seurat Visium object to a GiottoObject

#' @param object Seurat Object
#' @param assay Seurat assay in \code{object}
#' @param slot Seurat slot in \code{object} and \code{assay}
#' @param ... Arguments passed to \code{Giotto::createGiottoObject}
#' @export as_GiottoObject
#' @examples
#' # as_GiottoObject(object)
as_GiottoObject <- function(object, assay = "Spatial", slot = "counts", ...) {
    if (!is(object, "Seurat")) {
        stop("Only Seurat objects supported")
    }
    if (!requireNamespace("Giotto", quietly = TRUE)) {
        stop("Install Giotto.")
    }
    if (!requireNamespace("magick", quietly = TRUE)) {
        stop("Install magick.")
    }
    if (length(Images(object)) > 1) {
        stop("Currently only single sample objects supported.")
    }
    raw_matrix <- GetAssayData(object, assay = assay, slot = slot)
    spatial_locs <- GetTissueCoordinates(object)
    spatial_locs <- spatial_locs / object@images[[1]]@scale.factors$lowres
    spatial_locs <-  spatial_locs[, c(2, 1)]
    colnames(spatial_locs) <- c("sdimx", "sdimy")
    spatial_locs$sdimy <- spatial_locs$sdimy * -1
    cell_metadata <-  data.table(
        cell_ID = rownames(object@images[[1]]@coordinates),
        object@images[[1]]@coordinates[, c("tissue", "row", "col")])

    colnames(cell_metadata) <- c("cell_ID", "in_tissue", "array_row", "array_col")
    mg_object <- Giotto::createGiottoImage(spatial_locs = spatial_locs,
        mg_object = magick::image_read(object@images[[1]]@image),
        name = Images(object)[1],
        scale_factor=object@images[[1]]@scale.factors$lowres)
    images <- list(mg_object)
    names(images) <- Images(object)[1]

    Giotto::createGiottoObject(
        expression = raw_matrix,
        expression_feat = "rna",
        spatial_locs = spatial_locs,
        cell_metadata = list(cell = list(rna = cell_metadata)),
        images = images,
        ...)
}    
