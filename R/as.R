
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
        object@images[[1]]@coordinates[rownames(SummarizedExperiment::colData(x)),]
    )
    if (requireNamespace("BayesSpace", quietly = TRUE)) {
        x <- BayesSpace::spatialPreprocess(x, ...)
    }    
    return(x)
}    
