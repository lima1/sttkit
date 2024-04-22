
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
#' @param backend_method Method used by the backend. Currently only used by 
#' \code{sctransform} and \code{sctransform2} since the default assay is kept unchanged.
#' @param feature_filter If feature meta data contains 'included' column, only use 
#  these.
#' @param regressout Regressout these features.
#' @param cell_cycle_score Use \code{CellCycleScoring} to add S/G1/G2M scores
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
                              method = c("sctransform", "sctransform2", "seurat2", "scran"),
                              backend_method = "poisson",
                              feature_filter = TRUE,
                              cell_cycle_score = TRUE,
                              regressout = NULL, assay = "Spatial",
                              serialize = TRUE, prefix, ...) {
    regressout <- .check_regressout(obj, regressout)


    flog.info("Using standard log-normalization...")
    obj <- NormalizeData(object = obj, ...)
    scale_alt_assay <- NULL

    method <- match.arg(method)
    vst.flavor <- NULL
    if (method == "sctransform2") {
        method <- "sctransform"
        vst.flavor <- "v2"
    }
    feature_filter <- feature_filter && is(obj[[assay]][[]]$included, "logical")

    if (feature_filter && any(!obj[[assay]][[]]$included)) {
        flog.info("Removing %i flagged genes from assay %s",
            length(which(!obj[[assay]][[]]$included)), assay)

        obj <- subset(obj,
            features = rownames(obj[[assay]][[]][which(obj[[assay]][[]]$included), , drop = FALSE]))
    }
    if (method == "sctransform" && requireNamespace("sctransform")) {
        scale_alt_assay <- DefaultAssay(obj)

        flog.info("Using sctransform with %i features...", nfeatures)
        if (!is.null(regressout)) {
            flog.info("Regressing out %s.", paste(regressout, collapse = ", "))
        }
        min_cells <- min(Matrix::rowSums(GetAssayData(obj, slot = "counts"))) + 1
        conserve.memory <- ncol(obj) > 5000
        if ("vst.flavor" %in% names(formals(sctransform::vst))) {
            obj <- SCTransform(obj, variable.features.n = nfeatures, assay = assay,
                vars.to.regress = regressout, do.correct.umi = correct_umi,
                min_cells = min_cells, conserve.memory = conserve.memory,
                vst.flavor = vst.flavor, method = backend_method, ...)
        } else {
            obj <- SCTransform(obj, variable.features.n = nfeatures, assay = assay,
                vars.to.regress = regressout, do.correct.umi = correct_umi,
                min_cells = min_cells, conserve.memory = conserve.memory,
                method = backend_method, ...)
        }    

        if (cell_cycle_score) obj <- .add_cc_score(obj)
        if (serialize) .serialize(obj, prefix, "_scaled.rds")
        flog.info("Scaling alternative assay %s...", scale_alt_assay)
        obj <- ScaleData(object = obj, vars.to.regress = regressout,
            do.scale = scale, do.center = center, assay = scale_alt_assay)
        flog.info("Default assay is set to %s.", DefaultAssay(obj))
        return(obj)    
    } else if (method == "scran" && requireNamespace("scran")) {
        scale_alt_assay <- DefaultAssay(obj)
        flog.info("Using scran normalization...")
        sce <- SingleCellExperiment::SingleCellExperiment(assays =
            list(counts = as.matrix(GetAssayData(obj, slot = "counts")))) # read data from Seurat
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
        obj@misc[["seurat_norm_data"]] <- as.matrix(x = GetAssayData(obj)) # backup Seurat's norm data
        SetAssayData(obj, slot = "data", new.data =
            log(x = SummarizedExperiment::assay(sce, "normcounts") + 1))
    }
    flog.info("Finding %i variable features...", nfeatures)
    obj <- FindVariableFeatures(object = obj, selection.method = "vst",
        nfeatures = nfeatures, verbose = FALSE)
    if (cell_cycle_score) obj <- .add_cc_score(obj)
    if (serialize) .serialize(obj, prefix, "_unscaled.rds")
    if (!is.null(regressout)) {
        flog.info("Regressing out %s.", paste(regressout, collapse = ", "))
    }
    obj <- ScaleData(object = obj, vars.to.regress = regressout,
        do.scale = scale, do.center = center)
    if (!is.null(scale_alt_assay)) {
        flog.info("Scaling alternative assay %s...", scale_alt_assay)
        obj <- ScaleData(object = obj, vars.to.regress = regressout,
            do.scale = scale, do.center = center, assay = scale_alt_assay)
    }
    flog.info("Default assay is set to %s.", DefaultAssay(obj))
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
