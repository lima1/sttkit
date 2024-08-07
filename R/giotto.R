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
        stop("Install Giotto >= 4.0.")
    }
    if (!requireNamespace("GiottoClass", quietly = TRUE)) {
        stop("Install GiottoClass")
    }
    if (!requireNamespace("magick", quietly = TRUE)) {
        stop("Install magick.")
    }
    if (length(Images(object)) > 1) {
        stop("Currently only single sample objects supported.")
    }
    raw_matrix <- GetAssayData(object, assay = assay, slot = slot)
    spatial_locs <- GetTissueCoordinates(object)[, 1:2]
    spatial_locs <- spatial_locs / object@images[[1]]@scale.factors$lowres
    spatial_locs <-  spatial_locs[, c(2, 1)]
    colnames(spatial_locs) <- c("sdimx", "sdimy")
    spatial_locs$sdimy <- spatial_locs$sdimy * -1
    cell_metadata <- data.table(
        cell_ID = rownames(object@meta.data),
        object@meta.data)
    mg_object <- GiottoClass::createGiottoImage(spatial_locs = spatial_locs,
        mg_object = magick::image_read(object@images[[1]]@image),
        name = Images(object)[1],
        scale_factor=object@images[[1]]@scale.factors$lowres)
    images <- list(mg_object)
    names(images) <- Images(object)[1]

    gobject <- GiottoClass::createGiottoObject(
        expression = raw_matrix,
        expression_feat = "rna",
        spatial_locs = spatial_locs,
        cell_metadata = list(cell = list(rna = cell_metadata)),
        images = images,
        ...)

    if ("SCT" %in% Assays(object)) {
        norm_exp <- Seurat::GetAssayData(object = object, 
                        slot = "data", assay = "SCT")
        expr_obj <- new("exprObj", name = "normalized", exprMat = norm_exp,
                spat_unit = "cell", feat_type = "rna", provenance = "cell")
        gobject <- GiottoClass::set_expression_values(gobject = gobject, values = expr_obj, 
            set_defaults = FALSE)
    }
    if ("predictions" %in% Assays(object)) {
        enrObj <- as_spatEnrObj(object$predictions)
        gobject <- GiottoClass::set_spatial_enrichment(gobject, enrObj)
    }    
    return(gobject) 
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
#'
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
        stop("Install Giotto >= 4.0.")
    }
    spot_proportion <- GetAssayData(object, slot = slot)
    spot_proportion <- spot_proportion[!rownames(spot_proportion) %in% ignore, ]
    deconvolutionDT <- data.table::data.table(cell_ID = colnames(spot_proportion))
    deconvolutionDT <- cbind(deconvolutionDT, as.data.table(t(spot_proportion)))
    new("spatEnrObj", name = name, method = method, enrichDT = deconvolutionDT,
        spat_unit = spat_unit, feat_type = feat_type, ...)
}

#' find_giotto_dwls_matrix
#'
#'
#' @param singlecell_giotto List of Giotto objects 
#' @param refdata Column with cell type information
#' @param method Name of the method used to find markers
#' @param expression_values Giotto expression values slot
#' @param num_markers Number of significant markers per cell type
#' @export find_giotto_dwls_matrix
#' @examples
#' #
find_giotto_dwls_matrix <- function(singlecell_giotto, refdata, method = "scran", expression_values = "normalized",
                                    num_markers = 100) {

    if (!requireNamespace("Giotto", quietly = TRUE)) {
        stop("Install Giotto >= 4.0.")
    }
    if (!requireNamespace("GiottoClass", quietly = TRUE)) {
        stop("Install GiottoClass")
    }
    scran_markers_subclusters <- lapply(singlecell_giotto, Giotto::findMarkers_one_vs_all,
           method = method,
           expression_values = expression_values,
           cluster_column = refdata)
    
    sign_matrix <- lapply(seq_along(singlecell_giotto), function(i) {
        sig_scran <- unique(scran_markers_subclusters[[i]]$feats[which(scran_markers_subclusters[[i]]$ranking <= num_markers)])
        gene_metadata <- GiottoClass::fDataDT(singlecell_giotto[[i]])
        feats <- rownames(infile)
        feats <- unique(c(sig_scran, feats[feats %in% GiottoClass::fDataDT(singlecell_giotto[[i]])$feat_ID]))
        id <- GiottoClass::pDataDT(singlecell_giotto[[i]])[[opt$refdata]]
        exp <- Giotto::makeSignMatrixDWLS(singlecell_giotto[[i]], cell_type_vector = id, sign_gene = feats)
        list(matrix = exp, sig_feats = sig_scran)
    })
    return(sign_matrix)
}

#' calculate_giotto_spatial_correlation_cell_type
#'
#' Calculates spatial correlation of feature and cell type
#'
#' @param gobject giotto object
#' @param cell_types Cell type of interest
#' @param features Features of interest. If \code{NULL}, use all.
#' @export calculate_giotto_spatial_correlation_cell_type
#' @examples
#' #
calculate_giotto_spatial_correlation_cell_type <- function(gobject, cell_types, features = NULL) {
    if (!requireNamespace("GiottoClass", quietly = TRUE)) {
        stop("Install GiottoClass")
    }
    dt_net <- gobject@spatial_network$cell$Delaunay_network@networkDT
    dt_enr <- gobject@spatial_enrichment$cell$rna$DWLS@enrichDT
    expr_m <- GiottoClass::getExpression(gobject, "normalized", output = "matrix")
    expr_m <- expr_m[apply(expr_m, 1, function(x) sum(x > 0)) >= 3, ]
    if (is.null(features)) features <- rownames(expr_m)

    features <- features[features %in% rownames(expr_m)]
    .get_correlations <- function(expr, ct) {
        dt_net$from_expr <- expr[dt_net$from]
        dt_net$from_cell <- ct[dt_net$from]
        dt_net$to_expr <- expr[dt_net$to]
        dt_net$to_cell <- ct[dt_net$to]
        list(
            external = mean(cor(dt_net$from_cell, dt_net$to_expr),
                            cor(dt_net$to_cell, dt_net$from_expr)),
            internal = mean(cor(dt_net$from_cell, dt_net$from_expr),
                            cor(dt_net$to_cell, dt_net$to_expr))
        )
    }
    r <- rbindlist(lapply(cell_types, function(cell_type) 
             rbindlist(lapply(features, function(feature) {
                 ct <- dt_enr[[cell_type]]
                 names(ct) <- dt_enr$cell_ID
                 expr <- expr_m[feature, ]
                 cors <- .get_correlations(expr, ct)
                 data.table(feature = feature, cell_type = cell_type, external = cors$external, internal = cors$internal)
             }))
         ))
    
    return(r)
}    
