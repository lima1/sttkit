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

.as_AssayObject_giotto <- function(object) {
    m <- as.matrix(object@enrichDT[,-1])
    rownames(m) <- object@enrichDT$cell_ID
    m <- t(m)
    m <- rbind(m, max = apply(m, 2, max))
    return(CreateAssayObject(data = m))
}


#' as_spatEnrObj
#'
#' Convert \code{Seurat::AssayObject} to a \code{Giotto::spatEnrObj}

#' @param object Seurat \code{AssayObject}
#' @param slot Seurat slot in \code{object} 
#' @param ignore Features to be ignored
#' @param method Name of the method used to generate \code{object}
#' @param name Name of this enrichment object. 
#' @param spat_unit Giotto spatial unit.
#' @param feat_type Giotto feature type.
#' @param ... Additional arguments passed to \code{Giotto::new("spatEnrObj")}
#' @export as_spatEnrObj
#' @examples
#' # giotto_enrichment <- as_spatEnrObj(x$predictions)
#' # giotto_object <- as_GiottoObject(x)
#' # giotto_object <- set_spatial_enrichment(gobject = giotto_object,
#' #                                         spatenrichment = giotto_enrichment)
#' # spatDeconvPlot(giotto_object, radius = 150)
#' # 
as_spatEnrObj <- function(object, slot = "data", ignore = c("max", "unassigned"),
    method = "DWLS", name = method, spat_unit = "cell", feat_type = "rna", ...) {
    if (!is(object, "Assay")) {
        stop("Only Assay objects supported")
    }
    if (!requireNamespace("Giotto", quietly = TRUE)) {
        stop("Install Giotto.")
    }
    spot_proportion <- GetAssayData(object, slot = slot)
    spot_proportion <- spot_proportion[!rownames(spot_proportion) %in% ignore, ]
    deconvolutionDT <- data.table::data.table(cell_ID = colnames(spot_proportion))
    deconvolutionDT <- cbind(deconvolutionDT, as.data.table(t(spot_proportion)))
    new("spatEnrObj", name = name, method = method, enrichDT = deconvolutionDT,
        spat_unit = spat_unit, feat_type = feat_type, ...)
}
