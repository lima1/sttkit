.as_AssayObject_scvi <- function(object) {
    m <- as.matrix(py_to_r(object$obsm$get("proportions")))
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
        stop("Install Giotto.")
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
    scran_markers_subclusters <- lapply(singlecell_giotto, findMarkers_one_vs_all,
           method = method,
           expression_values = expression_values,
           cluster_column = refdata)
    
    sign_matrix <- lapply(seq_along(singlecell_giotto), function(i) {
        sig_scran <- unique(scran_markers_subclusters[[i]]$feats[which(scran_markers_subclusters[[i]]$ranking <= num_markers)])
        gene_metadata <- fDataDT(singlecell_giotto[[i]])
        feats <- rownames(infile)
        feats <- unique(c(sig_scran, feats[feats %in% fDataDT(singlecell_giotto[[1]])$feat_ID]))
        id <- pDataDT(singlecell_giotto[[i]])[[opt$refdata]]
        exp <- makeSignMatrixDWLS(singlecell_giotto[[1]], cell_type_vector = id, sign_gene = feats)
        list(matrix = exp, sig_feats = sig_scran)
    })
    return(sign_matrix)
}
