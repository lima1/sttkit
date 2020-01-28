
#' normalize_spatial
#'
#' Normalize spatial transcriptomics Seurat object
#' @param obj Object, read by \code{\link{read_spatial}}.
#' @param nfeatures Number of variable features
#' @param scale Use \code{Seurat::ScaleData} to scale data.
#' @param center Use \code{Seurat::ScaleData} to center data.
#' @param correct_umi Currently only supported by \code{sctransform}
#' @param method Use \code{sctransform} to normalize and scale data,
#' or use \code{scran} or Seurat 2 style normalization and scaling.
#' @param regressout Regressout these features.
#' @param assay Name of the assay corresponding to the initial input data.
#' @param serialize Automatically serialize object
#' @param prefix Prefix of output files
#' @param ... Additional parameters passed to the normalization 
#' function.

#' @export normalize_spatial
#' @examples
#' normalize_spatial()

normalize_spatial <- function(obj, nfeatures = 2500, scale = TRUE, center = TRUE,
                              correct_umi = TRUE,
                              method = c("sctransform", "seurat2", "scran"), 
                              regressout = NULL, assay = "Spatial",
                              serialize = TRUE, prefix, ...) {
    regressout <- .check_regressout(obj, regressout)

    method <- match.arg(method)
    if (method == "sctransform" && requireNamespace("sctransform")) {
        flog.info("Using sctransform with %i features...", nfeatures)
        if (!is.null(regressout)) {
            flog.info("Regressing out %s.", paste(regressout, collapse = ", "))
        }
        min_cells <- min(Matrix::rowSums(GetAssayData(obj, "counts"))) + 1
        obj <- SCTransform(obj, variable.features.n = nfeatures, assay = assay,
            vars.to.regress = regressout, do.correct.umi = correct_umi,
            return.only.var.genes = FALSE, min_cells = min_cells, ...)
        if (serialize) .serialize(obj, prefix, "_scaled.rds")
        return(obj)    
    } else if (method == "scran" && requireNamespace("scran")) {
        flog.info("Using scran normalization...")
        sce <- SingleCellExperiment::SingleCellExperiment(assays = 
            list(counts = as.matrix(GetAssayData(obj, "counts")))) # read data from Seurat
        clusters <- NULL
        sizes <- seq(21, 101, 5)
        if (ncol(obj) < 200) {
            flog.warn("Not enough cells (%i) for clustered scran normalization.",
                ncol(obj))
            sizes <-seq(min(21, ncol(obj)), max(101, ncol(obj)))
        } else {    
            clusters <- scran::quickCluster(sce, min.size = 100)
        } 
        obj@meta.data$scran.cluster <- if (is.null(clusters)) "0" else clusters
        sce <- scran::computeSumFactors(sce, sizes = sizes, clusters = clusters)
        sce <- scater::normalize(sce, return_log = FALSE) # without(!) log transform
        obj <- NormalizeData(object = obj)
        obj@misc[["seurat_norm_data"]] = as.matrix(x = GetAssayData(obj)) # backup Seurat's norm data
        SetAssayData(obj, slot = "data", new.data = 
            log(x = SummarizedExperiment::assay(sce, "normcounts") + 1))
    } else {
        flog.info("Using standard log-normalization...")
        obj <- NormalizeData(object = obj, ...)
    }    
    flog.info("Finding %i variable features...", nfeatures)
    obj <- FindVariableFeatures(object = obj, selection.method = "vst",
        nfeatures = nfeatures, verbose = FALSE)
    if (serialize) .serialize(obj, prefix, "_unscaled.rds")
    if (!is.null(regressout)) {
        flog.info("Regressing out %s.", paste(regressout, collapse = ", "))
    }
    obj <- ScaleData(object = obj, vars.to.regress = regressout,
        do.scale = scale, do.center = center)
    if (serialize) {
        flog.info("Writing R data structure to %s...", paste0(prefix, "_scaled.rds"))
        .serialize(obj, prefix, "_scaled.rds")
    }    
    obj     
}

.check_regressout <- function(obj, regressout) {
    if (!is.null(regressout)) {
        regressout <- regressout[sapply(regressout, function(x) x %in% colnames(obj@meta.data))]
        if (!length(regressout)) regressout <- NULL
    }
    regressout
}        
