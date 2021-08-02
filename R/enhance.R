#' enhance_bayesspace
#'
#' Enhance BayesSpace data
#' @param x Object, converted by \code{\link{as_SingleCellExperiment}}
#' @param test_num_clusters Range of number of clusters to be tested (\code{qs}). 
#' @param num_clusters Optional picked number
#' @param force Recalculate, even when serialized objects are available
#' @param serialize Serialize output objects
#' @param prefix Prefix of output files
#' @param ... Additional paramters passed to \code{BayesSpace::spatialCluster}
#' @export enhance_bayesspace
#' @examples
#' #enhance_bayesspace()

enhance_bayesspace <- function(x, test_num_clusters = seq(2, 12),
                               num_clusters = NULL, ...) {
    flog.info("Enhancing %i clusters. This will take a while...", opt$num_clusters)
    ndata_enhanced <- spatialEnhance(ndata,
        q = opt$num_clusters,
        nrep = opt$num_iter,
        verbose = TRUE)
    flog.info("Writing R data structure to %s...", filename)
    saveRDS(ndata_enhanced, file = filename)
}    
